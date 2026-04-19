# =============================================================================
# AVONET imputation benchmark for BACE
# =============================================================================
#
# Purpose
# -------
# Evaluate BACE's imputation accuracy on a real comparative dataset -
# the AVONET bird trait database (Tobias et al. 2022, Ecol Lett 25:581)
# - under a known-masking protocol. Complements the controlled
# simulation benchmark (`dev/simulation_imputation_quality.R`).
#
#   Simulated  known truth, dialled phylogenetic signal, dialled
#              missingness mechanism. Answers "is BACE well-calibrated
#              under specific conditions I can control?"
#
#   AVONET     real bird trait data with a known mask, real phylogeny
#              (Hackett et al. 2008, Science 320:1763). Answers "does
#              BACE handle a representative empirical dataset the way
#              a comparative biologist actually uses it?"
#
# This script is designed to be re-run as BACE evolves (e.g. after
# changes to `.predict_bace`, priors, MCMC settings, or the default
# `ovr_categorical` behaviour). Each run writes a tagged, versioned
# output bundle that can be diffed against previous runs.
#
# Data
# ----
#   dev/testing_data/avonet_2000_masked.csv    2000 species x 8 traits;
#                                              ~10% of cells set to NA
#   dev/testing_data/avonet_2000_truth.csv     long form (species, trait,
#                                              true_value) for masked cells
#   dev/testing_data/Hackett_tree_2000.tre     matching 2000-tip phylogeny
#
# Traits (all continuous, log-transformed or appropriately scaled):
#   Mass, Wing.Length, Beak.Length_Culmen, Tarsus.Length, Tail.Length,
#   Range.Size, Centroid.Latitude, Centroid.Longitude.
#
# Metric glossary (consistent with simulation benchmark)
# ------------------------------------------------------
# All metrics are computed on HIDDEN CELLS ONLY and averaged across the
# n_final imputations. Coverage is computed from the full ensemble.
#
#   NRMSE       = RMSE(imp, true) / sd(complete column).
#                 Scale-free; 0 = perfect; 1 = as bad as marginal-mean.
#                 Normalised by marginal sd (Stekhoven & Buhlmann 2012;
#                 van Buuren 2018 FIMD §5.1), not hidden-subset sd.
#   MAE         = mean(|imp - true|) on the raw scale.
#   correlation = Pearson r(imp, true).
#   coverage95  = proportion of hidden cells whose true value lies
#                 inside the 95% PI [q_0.025, q_0.975] computed across
#                 the n_final imputations (Rubin 1987; van Buuren 2018
#                 FIMD §2.5). Well-calibrated -> ~ 0.95.
#
# Re-running across BACE versions
# -------------------------------
# 1. Set BACE_LOAD below to either:
#      "load_all"  -> uses the CURRENT source tree (this branch)
#      "library"   -> uses the globally-installed BACE
# 2. Configure MCMC / subset size / output dir if needed.
# 3. Run the script. Output lands in:
#      dev/benchmark_results/avonet/run_<tag>_<date>/
#    with the git commit hash or "uninstalled" tag so diffs are
#    traceable.
# 4. Use `dev/compare_avonet_runs.R` (TODO) to generate side-by-side
#    comparison tables across tagged runs.
# =============================================================================

# ---- Package loading --------------------------------------------------------
# Set to "load_all" to test the current (possibly uncommitted) source, or
# "library" to test the globally-installed BACE. The run output is tagged
# accordingly for later comparison.
BACE_LOAD <- "load_all"    # "load_all" or "library"

if (BACE_LOAD == "load_all") {
  devtools::load_all(quiet = TRUE)
} else {
  library(BACE)
}
library(ape)
library(MASS)

set.seed(2026)

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

# Subset (set to NA for the full 2000 species). Subsampling is useful for
# quick iteration; the full benchmark is the production configuration.
SUBSET_N <- NA

# MCMC budget (see simulation benchmark header for the Hadfield 2010 /
# Graham et al. 2007 / van Buuren 2018 rationale).
NITT    <- 50000
THIN    <- 25
BURNIN  <- 10000
RUNS    <- 10
N_FINAL <- 20
MAX_ATTEMPTS <- 2
N_CORES <- 4L

# Traits to impute (all 8 continuous columns of AVONET).
TRAIT_COLS <- c("Mass", "Wing.Length", "Beak.Length_Culmen",
                "Tarsus.Length", "Tail.Length", "Range.Size",
                "Centroid.Latitude", "Centroid.Longitude")

# =============================================================================
# 2. RUN TAGGING (traceability across reruns)
# =============================================================================

#' Build a run tag combining BACE source info + date + scope, so outputs
#' from successive re-runs are immediately distinguishable.
.make_run_tag <- function() {
  date_str <- format(Sys.Date(), "%Y%m%d")

  bace_tag <- if (BACE_LOAD == "load_all") {
    # Try to read the current git HEAD so the tag traces to a commit.
    head_sha <- tryCatch(
      substr(system("git rev-parse HEAD", intern = TRUE,
                     ignore.stderr = TRUE), 1, 8),
      error = function(e) NA_character_, warning = function(w) NA_character_
    )
    dirty <- tryCatch(
      length(system("git status --porcelain", intern = TRUE,
                     ignore.stderr = TRUE)) > 0,
      error = function(e) FALSE
    )
    if (is.na(head_sha) || !nzchar(head_sha)) "uncommitted"
    else paste0("src-", head_sha, if (dirty) "-dirty" else "")
  } else {
    paste0("pkg-", as.character(utils::packageVersion("BACE")))
  }

  scope_tag <- if (!is.na(SUBSET_N)) paste0("n", SUBSET_N) else "full"
  paste(bace_tag, scope_tag, date_str, sep = "_")
}

RUN_TAG <- .make_run_tag()
OUT_DIR <- file.path("dev", "benchmark_results", "avonet",
                     paste0("run_", RUN_TAG))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# =============================================================================
# 3. DATA + TREE PREPARATION
# =============================================================================

masked <- read.csv("dev/testing_data/avonet_2000_masked.csv",
                   stringsAsFactors = FALSE)
truth  <- read.csv("dev/testing_data/avonet_2000_truth.csv",
                   stringsAsFactors = FALSE)
tree   <- ape::read.tree("dev/testing_data/Hackett_tree_2000.tre")

# bace() expects the random-effect grouping column to be "Species".
names(masked)[names(masked) == "species"] <- "Species"
stopifnot(all(masked$Species %in% tree$tip.label))

#' Prep a phylogeny for MCMCglmm:
#'   (a) Force ultrametric if the tree is only ultrametric modulo
#'       floating-point drift (the published Hackett tree is).
#'   (b) Replace zero-length edges with a small positive value — these
#'       can appear after `ape::keep.tip` subsampling.
#' MCMCglmm::inverseA refuses non-ultrametric trees or zero edges.
.prep_tree <- function(tree) {
  if (!ape::is.ultrametric(tree)) {
    if (requireNamespace("phytools", quietly = TRUE)) {
      tree <- phytools::force.ultrametric(tree, method = "extend")
    } else {
      depths    <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
      extend    <- max(depths) - depths
      tip_edges <- match(seq_len(ape::Ntip(tree)), tree$edge[, 2])
      tree$edge.length[tip_edges] <- tree$edge.length[tip_edges] + extend
    }
  }
  if (any(tree$edge.length <= 0)) {
    pos <- tree$edge.length[tree$edge.length > 0]
    eps <- if (length(pos)) min(pos) / 1e3 else 1e-6
    tree$edge.length[tree$edge.length <= 0] <- eps
  }
  tree
}

# Optionally subsample BEFORE ultrametricizing — pruning can introduce
# zero-length edges that .prep_tree then cleans up.
if (!is.na(SUBSET_N) && SUBSET_N < nrow(masked)) {
  keep_species <- sample(masked$Species, SUBSET_N)
  masked <- masked[masked$Species %in% keep_species, ]
  truth  <- truth[truth$species_tip %in% keep_species, ]
  tree   <- ape::keep.tip(tree, keep_species)
}
tree <- .prep_tree(tree)

# =============================================================================
# 4. RUN BANNER
# =============================================================================

cat("===========================================================\n")
cat("  AVONET imputation benchmark\n")
cat("  Tag     :", RUN_TAG, "\n")
cat("  Output  :", OUT_DIR, "\n")
cat("  BACE    :", BACE_LOAD,
    if (BACE_LOAD == "library")
      paste0(" v", utils::packageVersion("BACE")) else "",
    "\n")
cat("  Species :", nrow(masked),
    if (is.na(SUBSET_N)) "(full)" else paste0("(subset of 2000)"),
    "\n")
cat("  Traits  :", length(TRAIT_COLS), "\n")
cat("  Masked  :", nrow(truth), "cells total (avg ~",
    round(nrow(truth) / length(TRAIT_COLS)), "per trait)\n")
cat("  MCMC    : nitt=", NITT, " thin=", THIN, " burnin=", BURNIN,
    " runs=", RUNS, " n_final=", N_FINAL, "\n", sep = "")
cat("===========================================================\n\n")

# =============================================================================
# 5. RUN BACE
# =============================================================================

t0 <- Sys.time()
fixformulas <- lapply(TRAIT_COLS, function(v) {
  others <- setdiff(TRAIT_COLS, v)
  paste(v, "~", paste(others, collapse = " + "))
})

res <- bace(
  fixformula     = fixformulas,
  ran_phylo_form = "~1|Species",
  phylo          = tree,
  data           = masked,
  nitt           = NITT, thin = THIN, burnin = BURNIN,
  runs           = RUNS, n_final = N_FINAL,
  species        = FALSE, verbose = TRUE,
  skip_conv      = FALSE, max_attempts = MAX_ATTEMPTS,
  n_cores        = N_CORES
)

runtime_min <- as.numeric(round(difftime(Sys.time(), t0, units = "mins"), 1))
cat("\nbace() finished in", runtime_min, "min\n")
cat("Converged:", res$converged, "  attempts:", res$n_attempts, "\n")

# =============================================================================
# 6. METRICS
# =============================================================================

#' Per-trait accuracy + calibration metrics. Point-estimate metrics
#' (NRMSE, MAE, cor) are computed per imputation then averaged.
#' coverage95 uses the full n_final ensemble (central 95% interval).
compute_avonet_metrics <- function(res, masked, truth, trait_cols) {
  rows <- lapply(trait_cols, function(v) {
    truth_v <- truth[truth$trait == v, ]
    idx     <- match(truth_v$species_tip, masked$Species)
    tv      <- truth_v$true_value
    sd_full <- sd(masked[[v]], na.rm = TRUE)

    imp_mat <- vapply(res$imputed_datasets,
                      function(d) as.numeric(d[[v]][idx]),
                      FUN.VALUE = numeric(length(tv)))
    if (!is.matrix(imp_mat)) imp_mat <- matrix(imp_mat, nrow = length(tv))

    per_imp_nrmse <- apply(imp_mat, 2,
                            function(iv) sqrt(mean((iv - tv)^2)) / sd_full)
    per_imp_mae   <- apply(imp_mat, 2, function(iv) mean(abs(iv - tv)))
    per_imp_cor   <- apply(imp_mat, 2, function(iv)
      if (length(unique(tv)) > 1 && length(unique(iv)) > 1)
        suppressWarnings(cor(iv, tv)) else NA_real_)

    lo <- apply(imp_mat, 1, stats::quantile, probs = 0.025,
                na.rm = TRUE, names = FALSE)
    hi <- apply(imp_mat, 1, stats::quantile, probs = 0.975,
                na.rm = TRUE, names = FALSE)

    data.frame(
      trait       = v,
      n_hidden    = length(tv),
      nrmse       = mean(per_imp_nrmse, na.rm = TRUE),
      mae         = mean(per_imp_mae),
      correlation = mean(per_imp_cor, na.rm = TRUE),
      coverage95  = mean(tv >= lo & tv <= hi),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summary_df <- compute_avonet_metrics(res, masked, truth, TRAIT_COLS)

cat("\n===========================================================\n")
cat("  AVONET results by trait\n")
cat("===========================================================\n")
print(summary_df, digits = 3, row.names = FALSE)

cat("\nAverages across traits:\n")
avg_row <- data.frame(
  trait = "(mean)",
  n_hidden = sum(summary_df$n_hidden),
  nrmse = mean(summary_df$nrmse, na.rm = TRUE),
  mae = NA_real_,  # MAE not comparable across traits with different scales
  correlation = mean(summary_df$correlation, na.rm = TRUE),
  coverage95 = mean(summary_df$coverage95, na.rm = TRUE)
)
print(avg_row, digits = 3, row.names = FALSE)

# =============================================================================
# 7. METADATA + OUTPUT
# =============================================================================

metadata <- list(
  run_tag         = RUN_TAG,
  date            = Sys.time(),
  bace_source     = BACE_LOAD,
  bace_version    = as.character(utils::packageVersion("BACE")),
  git_head        = tryCatch(
    system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE),
    error = function(e) NA_character_),
  git_dirty       = tryCatch(
    length(system("git status --porcelain", intern = TRUE,
                   ignore.stderr = TRUE)) > 0,
    error = function(e) NA),
  n_species       = nrow(masked),
  subset_n        = if (is.na(SUBSET_N)) "full" else SUBSET_N,
  traits          = TRAIT_COLS,
  mcmc            = list(nitt = NITT, thin = THIN, burnin = BURNIN,
                         runs = RUNS, n_final = N_FINAL,
                         max_attempts = MAX_ATTEMPTS),
  runtime_min     = runtime_min,
  converged       = res$converged,
  n_attempts      = res$n_attempts,
  seed            = 2026
)

saveRDS(list(metadata = metadata,
             summary  = summary_df,
             imputed_datasets = res$imputed_datasets),
        file.path(OUT_DIR, "results.rds"))
write.csv(summary_df, file.path(OUT_DIR, "summary.csv"),
          row.names = FALSE)
writeLines(jsonlite::toJSON(metadata, pretty = TRUE, auto_unbox = TRUE,
                            POSIXt = "ISO8601"),
           file.path(OUT_DIR, "metadata.json"))

# =============================================================================
# 8. DIAGNOSTIC PLOTS
# =============================================================================

plot_file <- file.path(OUT_DIR, "plots.pdf")
pdf(plot_file, width = 11, height = 8)

# True vs imputed (posterior mean) per trait, with 95% PI bars.
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
for (v in TRAIT_COLS) {
  truth_v <- truth[truth$trait == v, ]
  idx     <- match(truth_v$species_tip, masked$Species)
  tv      <- truth_v$true_value
  imp_mat <- vapply(res$imputed_datasets,
                    function(d) as.numeric(d[[v]][idx]),
                    FUN.VALUE = numeric(length(tv)))
  if (!is.matrix(imp_mat)) imp_mat <- matrix(imp_mat, nrow = length(tv))
  mid <- rowMeans(imp_mat)
  lo  <- apply(imp_mat, 1, stats::quantile, 0.025, na.rm = TRUE, names = FALSE)
  hi  <- apply(imp_mat, 1, stats::quantile, 0.975, na.rm = TRUE, names = FALSE)

  xl <- range(c(tv, lo, hi), na.rm = TRUE)
  plot(tv, mid, pch = 19, col = adjustcolor("steelblue", 0.4),
       xlim = xl, ylim = xl,
       xlab = "True value", ylab = "Posterior mean (imputed)",
       main = v, cex = 0.5)
  segments(tv, lo, tv, hi, col = adjustcolor("steelblue", 0.15))
  abline(0, 1, col = "red", lty = 2)
  mt <- summary_df[summary_df$trait == v, ]
  legend("topleft", bty = "n", cex = 0.8,
         legend = c(sprintf("NRMSE=%.2f r=%.2f", mt$nrmse, mt$correlation),
                    sprintf("cov95=%.0f%%", mt$coverage95 * 100)))
}
mtext(paste("AVONET imputation diagnostics  -", RUN_TAG),
      outer = TRUE, cex = 1.1)

# Coverage scorecard: bar chart of coverage95 with target line.
par(mfrow = c(1, 1), mar = c(6, 4, 3, 1))
bp <- barplot(summary_df$coverage95, names.arg = summary_df$trait,
              ylim = c(0, 1), col = "#6a5acd", border = NA,
              main = "95% PI coverage per trait",
              ylab = "Coverage", las = 2)
abline(h = 0.95, lty = 2, col = "red")
text(bp, summary_df$coverage95,
     sprintf("%.2f", summary_df$coverage95), pos = 3, cex = 0.8, xpd = NA)

dev.off()

cat("\nOutputs written to:", OUT_DIR, "\n")
cat("  - results.rds    full output + imputed datasets\n")
cat("  - summary.csv    per-trait metrics\n")
cat("  - metadata.json  run metadata (BACE version, MCMC settings, etc.)\n")
cat("  - plots.pdf      diagnostic plots\n\n")
cat("Done.\n")
