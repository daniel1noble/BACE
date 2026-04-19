# =============================================================================
# AVONET regression benchmark for the cov-fix branch
# =============================================================================
#
# Purpose
# -------
# The cov-fix branch changes how `.predict_bace(sample = TRUE)` produces
# posterior predictive draws for gaussian, poisson, threshold, and
# categorical (incl. OVR) response families. We want to confirm that:
#
#   (a) Point-estimate performance on the AVONET subset is within noise
#       of what BACE produced on main. This protects against silent
#       accuracy regressions introduced by the sampling changes.
#
#   (b) Calibration (95% PI coverage) is higher under cov-fix, because
#       the draws now include observation-level variance.
#
# Data
# ----
#   - `dev/testing_data/avonet_2000_masked.csv`: 2000 bird species x 8
#     continuous traits (log-transformed already); each trait has ~10%
#     of cells set to NA under an externally-defined masking scheme.
#   - `dev/testing_data/avonet_2000_truth.csv`: long-form (species_tip,
#     trait, true_value) giving the hidden true values.
#   - `dev/testing_data/Hackett_tree_2000.tre`: matching phylogeny.
#
# Metrics (on hidden cells only, averaged over n_final imputations)
# -----------------------------------------------------------------
#   NRMSE       RMSE / sd(complete column)
#   MAE         mean(|imp - true|)
#   correlation Pearson r (rank / pattern recovery)
#   coverage95  proportion of true values inside the 95% PI from the
#               n_final draws
#
# Run cost: ~20-40 min per trait at production MCMC with n_final = 20
#           on 2000 rows (8 traits * separate formulas = ~160 min total).
#           Reduce n_final or nitt to shorten.
# =============================================================================

# Load the cov-fix SOURCE (NOT the installed BACE) so this benchmark
# actually exercises the changes we care about.
devtools::load_all(quiet = TRUE)
library(ape)

set.seed(2026)

# ---- Config -----------------------------------------------------------------
# SUBSET_N: if not NA, subsample to this many species to keep runtime
# tractable during development. At 2000 species the ginverse is
# 2000x2000 and fitting is expensive. Set to NA for the full benchmark.
SUBSET_N <- 500

NITT    <- 30000
THIN    <- 20
BURNIN  <- 5000
RUNS    <- 5
N_FINAL <- 10
N_CORES <- 4L

OUT_DIR <- file.path("dev", "simulation_results", "avonet_regression")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---- Load data --------------------------------------------------------------
masked <- read.csv("dev/testing_data/avonet_2000_masked.csv",
                   stringsAsFactors = FALSE)
truth  <- read.csv("dev/testing_data/avonet_2000_truth.csv",
                   stringsAsFactors = FALSE)
tree   <- ape::read.tree("dev/testing_data/Hackett_tree_2000.tre")

#' Prepare a tree for MCMCglmm: ensure ultrametric (floating-point drift)
#' and no zero-length edges (subsampling can collapse short edges).
.prep_tree <- function(tree) {
  if (!ape::is.ultrametric(tree)) {
    if (requireNamespace("phytools", quietly = TRUE)) {
      tree <- phytools::force.ultrametric(tree, method = "extend")
    } else {
      depths   <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
      extend   <- max(depths) - depths
      tip_edges <- match(seq_len(ape::Ntip(tree)), tree$edge[, 2])
      tree$edge.length[tip_edges] <- tree$edge.length[tip_edges] + extend
    }
  }
  if (any(tree$edge.length <= 0)) {
    # Replace zero/negative edges with a small positive value proportional
    # to the median positive edge length. Preserves the tree topology
    # and ultrametric structure at the scale MCMCglmm needs.
    pos <- tree$edge.length[tree$edge.length > 0]
    eps <- if (length(pos)) min(pos) / 1e3 else 1e-6
    tree$edge.length[tree$edge.length <= 0] <- eps
  }
  tree
}

# Rename species -> Species to match the bace() formula convention.
names(masked)[names(masked) == "species"] <- "Species"
stopifnot(all(masked$Species %in% tree$tip.label))

# Subsample species if requested. Prune the tree to match.
if (!is.na(SUBSET_N) && SUBSET_N < nrow(masked)) {
  keep_species <- sample(masked$Species, SUBSET_N)
  masked <- masked[masked$Species %in% keep_species, ]
  truth  <- truth[truth$species_tip %in% keep_species, ]
  tree   <- ape::keep.tip(tree, keep_species)
  cat("Subsampled to", nrow(masked), "species.\n")
}

# Prep tree AFTER any subsampling — subtree pruning can introduce
# zero-length edges that MCMCglmm::inverseA() rejects.
tree <- .prep_tree(tree)

trait_cols <- c("Mass", "Wing.Length", "Beak.Length_Culmen", "Tarsus.Length",
                "Tail.Length", "Range.Size",
                "Centroid.Latitude", "Centroid.Longitude")

cat("===========================================================\n")
cat("  AVONET regression benchmark (cov-fix branch)\n")
cat("  Species:", nrow(masked), "   Traits:", length(trait_cols), "\n")
cat("  MCMC   : nitt=", NITT, " thin=", THIN, " burnin=", BURNIN,
    " runs=", RUNS, " n_final=", N_FINAL, "\n", sep = "")
cat("===========================================================\n\n")

# ---- Run bace() -------------------------------------------------------------
t0 <- Sys.time()
fixformulas <- lapply(trait_cols, function(v) {
  others <- setdiff(trait_cols, v)
  paste(v, "~", paste(others, collapse = " + "))
})

res <- bace(
  fixformula     = fixformulas,
  ran_phylo_form = "~1|Species",
  phylo          = tree,
  data           = masked,
  nitt           = NITT,
  thin           = THIN,
  burnin         = BURNIN,
  runs           = RUNS,
  n_final        = N_FINAL,
  species        = FALSE,
  verbose        = TRUE,
  skip_conv      = FALSE,
  max_attempts   = 2,
  n_cores        = N_CORES
)
cat("\nbace() finished in",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("Converged:", res$converged, "  attempts:", res$n_attempts, "\n")

# ---- Metrics ----------------------------------------------------------------
# For each trait: align the n_final imputations with the truth table,
# compute NRMSE, MAE, correlation, and 95% PI coverage.
results <- list()
for (v in trait_cols) {
  truth_v   <- truth[truth$trait == v, ]
  idx       <- match(truth_v$species_tip, masked$Species)
  tv        <- truth_v$true_value
  sd_full   <- sd(masked[[v]], na.rm = TRUE)

  imp_mat <- vapply(res$imputed_datasets,
                    function(d) as.numeric(d[[v]][idx]),
                    FUN.VALUE = numeric(length(tv)))
  if (!is.matrix(imp_mat)) imp_mat <- matrix(imp_mat, nrow = length(tv))

  mid <- rowMeans(imp_mat)
  lo  <- apply(imp_mat, 1, stats::quantile, probs = 0.025,
               na.rm = TRUE, names = FALSE)
  hi  <- apply(imp_mat, 1, stats::quantile, probs = 0.975,
               na.rm = TRUE, names = FALSE)

  per_imp_nrmse <- apply(imp_mat, 2,
                          function(iv) sqrt(mean((iv - tv)^2)) / sd_full)
  per_imp_mae   <- apply(imp_mat, 2, function(iv) mean(abs(iv - tv)))
  per_imp_cor   <- apply(imp_mat, 2,
                          function(iv) if (length(unique(tv)) > 1)
                                         cor(iv, tv) else NA_real_)

  results[[v]] <- data.frame(
    trait = v,
    n_hidden = length(tv),
    nrmse = mean(per_imp_nrmse, na.rm = TRUE),
    mae = mean(per_imp_mae),
    correlation = mean(per_imp_cor, na.rm = TRUE),
    coverage95 = mean(tv >= lo & tv <= hi)
  )
}
summary_df <- do.call(rbind, results)
rownames(summary_df) <- NULL

cat("\n===========================================================\n")
cat("  AVONET regression results (cov-fix)\n")
cat("===========================================================\n")
print(summary_df, digits = 3, row.names = FALSE)

saveRDS(list(bace = res, summary = summary_df),
        file.path(OUT_DIR, "avonet_regression_cov-fix.rds"))
write.csv(summary_df,
          file.path(OUT_DIR, "avonet_regression_cov-fix.csv"),
          row.names = FALSE)

cat("\nSaved raw + summary to:", OUT_DIR, "\n")
cat("Done.\n")
