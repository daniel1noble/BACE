# =============================================================================
# 00_benchmark_AVONET.R - Empirical AVONET benchmark (routine-run)
# =============================================================================
#
# Role in the benchmarking suite
# ------------------------------
#   00_benchmark_AVONET.R       THIS FILE. Empirical benchmark on the
#                                2000-species AVONET subset (routine).
#   01_benchmark_simulated.R    small-scale simulation, one replicate
#                                end-to-end (routine).
#   02_benchmark_simulated_full.R  full crossed-design simulation
#                                (NOT routine; periodic/pre-release).
#
# Purpose
# -------
# Evaluate BACE's imputation accuracy on the AVONET bird trait database
# (Tobias et al. 2022, Ecol Lett 25:581) under a known-masking protocol.
# Complements the controlled-simulation benchmarks (01, 02):
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
#   dev/testing_data/data/avonet_traits.rda    9993 species x 16 cols
#                                              (8 continuous + 2 ordinal +
#                                               3 categorical + 3 tax);
#                                              built by data-raw/make_avonet.R
#                                              from AVONET.csv
#   dev/testing_data/data/avonet_tree.rda      9993-tip Hackett phylogeny
#                                              (Hackett et al. 2008,
#                                              Science 320:1763)
#
# The benchmark subsamples SUBSET_N species at runtime (set.seed-controlled),
# then masks CONT_MISS_RATE of observed continuous cells and CAT_MISS_RATE of
# observed categorical cells. Truth values for the masked cells are stored
# in `truth` (continuous) / `cat_truth` (categorical) for scoring. This
# replaces the previously bundled avonet_2000_masked.csv / _truth.csv /
# Hackett_tree_2000.tre fixtures.
#
# Traits
# ------
# Continuous (imputed as gaussian):
#   Mass, Wing.Length, Beak.Length_Culmen, Tarsus.Length, Tail.Length,
#   Range.Size                            -> LOG-transformed before imputation
#   Centroid.Latitude, Centroid.Longitude -> imputed on raw scale
#
# Categorical (imputed as multinomial factor; BACE ovr_categorical
# default routes these through J binary threshold models):
#   Trophic.Level       4 levels (Scavenger, Omnivore, Herbivore,
#                       Carnivore). Unordered here (though arguably
#                       an ecological gradient).
#   Primary.Lifestyle   5 levels (Aerial, Aquatic, Generalist,
#                       Insessorial, Terrestrial).
#
# Why log-transform morphometric + Range.Size? These variables span
# 5+ orders of magnitude (Mass: 2 g to 33 kg; Range.Size: 3 to
# 118M km^2) and are right-skewed. Phylogenetic imputation assumes
# approximately Gaussian residuals (Brownian motion on the tree); that
# assumption is violently violated on the raw scale and a handful of
# large values dominate any fit. Log-transformation is the standard
# convention in comparative biology (Felsenstein 1985 Am Nat 125:1;
# Garland et al. 1993 Syst Biol 42:265; Tobias et al. 2022 Ecol Lett
# 25:581 specifically for AVONET). Latitudes / longitudes are
# angular and not log-transformed.
#
# Imputation runs in log space for the transformed traits; the imputed
# values are back-transformed with expm1() for comparison against the
# raw-scale truth. All metrics are also reported on the log scale
# (where Gaussian assumptions hold) for interpretability.
#
# Metric glossary (consistent with simulation benchmark)
# ------------------------------------------------------
# All metrics are computed on HIDDEN CELLS ONLY and averaged across the
# n_final imputations. Coverage / Brier use the full ensemble.
#
# Continuous:
#   NRMSE       = RMSE(imp, true) / sd(complete column).
#                 Scale-free; 0 = perfect; 1 = as bad as marginal-mean.
#                 Normalised by marginal sd (Stekhoven & Buhlmann 2012;
#                 van Buuren 2018 FIMD §5.1), not hidden-subset sd.
#   MAE_fit     = mean(|imp - true|) on the fit scale (log for
#                 LOG_TRAITS, raw otherwise).
#   MAE_raw     = back-transformed to raw scale via expm1() when
#                 applicable; interpretable in original units.
#   correlation = Pearson r(imp, true).
#   coverage95  = proportion of hidden cells whose true value lies
#                 inside the 95% PI [q_0.025, q_0.975] computed across
#                 the n_final imputations (Rubin 1987; van Buuren 2018
#                 FIMD §2.5). Well-calibrated -> ~ 0.95.
#
# Categorical:
#   accuracy          = proportion of hidden cells imputed with the
#                       correct class (majority vote across n_final).
#   balanced_accuracy = mean of per-class recall. Robust to class
#                       imbalance (e.g. Trophic.Level: Carnivore 56%,
#                       Scavenger 0.2% in AVONET).
#   brier             = multiclass Brier score averaged across hidden
#                       cells (Brier 1950 MWR 78:1; Gneiting & Raftery
#                       2007 JASA 102:359). Per-cell score =
#                       sum_k (p_k - 1{true=k})^2 where p_k is the
#                       empirical frequency of class k across the
#                       n_final imputations. Lower is better; 0 is
#                       perfect, 0.5 is uninformative for binary.
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
SUBSET_N <- 300

# MCMC budget (see simulation benchmark header for the Hadfield 2010 /
# Graham et al. 2007 / van Buuren 2018 rationale).
# Production settings:
#   NITT=50000, THIN=25, BURNIN=10000, RUNS=10, N_FINAL=20
NITT    <- 20000
THIN    <- 15
BURNIN  <- 4000
RUNS    <- 5
N_FINAL <- 10
MAX_ATTEMPTS <- 2
N_CORES <- 4L

# Continuous traits to impute.
CONT_TRAITS <- c("Mass", "Wing.Length", "Beak.Length_Culmen",
                 "Tarsus.Length", "Tail.Length", "Range.Size",
                 "Centroid.Latitude", "Centroid.Longitude")

# Continuous traits that must be log-transformed before imputation.
# See header for rationale / citations. `log1p` is used so that zeros
# (rare, but Range.Size can be tiny) don't produce -Inf.
LOG_TRAITS <- c("Mass", "Wing.Length", "Beak.Length_Culmen",
                "Tarsus.Length", "Tail.Length", "Range.Size")

# Categorical traits to impute. Pulled directly from avonet_traits.
# BACE will dispatch these via ovr_categorical = TRUE (default), i.e.
# J binary threshold models.
CAT_TRAITS <- c("Trophic.Level", "Primary.Lifestyle")

# Rates at which we hide trait cells. The continuous mask used to be
# pre-baked into avonet_2000_masked.csv (~10% per col); we now generate
# it on the fly from the bundled .rda so the scope is set by SUBSET_N
# rather than by the fixture file.
CONT_MISS_RATE <- 0.10
CAT_MISS_RATE  <- 0.10

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

# Load bundled trait data + tree (built by data-raw/make_avonet.R).
load("dev/testing_data/data/avonet_traits.rda")
load("dev/testing_data/data/avonet_tree.rda")

# Translate the rda's snake_case column names back to the AVONET-native
# names used downstream (Mass, Wing.Length, Trophic.Level, ...). This
# keeps the rest of the script readable in the conventions of the
# AVONET / phylo-comp literature without altering benchmark logic.
.rda_to_native <- c(
  mass_g                = "Mass",
  wing_length_mm        = "Wing.Length",
  beak_length_culmen_mm = "Beak.Length_Culmen",
  tarsus_length_mm      = "Tarsus.Length",
  tail_length_mm        = "Tail.Length",
  range_size_km2        = "Range.Size",
  centroid_lat          = "Centroid.Latitude",
  centroid_lon          = "Centroid.Longitude",
  trophic_level         = "Trophic.Level",
  primary_lifestyle     = "Primary.Lifestyle"
)
keep_cols <- intersect(names(.rda_to_native), colnames(avonet_traits))
masked <- avonet_traits[, keep_cols, drop = FALSE]
colnames(masked) <- unname(.rda_to_native[keep_cols])
masked$Species <- rownames(masked)
masked <- masked[, c("Species", CONT_TRAITS, CAT_TRAITS), drop = FALSE]
tree   <- avonet_tree
stopifnot(all(masked$Species %in% tree$tip.label))

#' Phylogenetic signal for a continuous trait — Pagel's lambda and
#' Blomberg's K on the tips that have an observed value.
#'
#' @param tree   phylo object with tip.label %in% names(values)
#' @param values named numeric vector (names == tip labels)
#' @return named list with lambda, K, both NA if computation fails
#' References: Pagel 1999, Nature 401:877; Blomberg et al. 2003,
#' Evolution 57:717; Münkemüller et al. 2012, MEE 3:743 for the
#' rationale for reporting both side-by-side.
.phylo_signal_cont <- function(tree, values) {
  if (!requireNamespace("phytools", quietly = TRUE))
    return(list(lambda = NA_real_, K = NA_real_))
  x <- values[!is.na(values)]
  if (length(x) < 10L || var(x) == 0)
    return(list(lambda = NA_real_, K = NA_real_))
  tr <- ape::keep.tip(tree, names(x))
  # phylosig() prints a lot of noise for non-ultrametric drift; suppress.
  lambda <- tryCatch(
    suppressWarnings(phytools::phylosig(tr, x, method = "lambda")$lambda),
    error = function(e) NA_real_)
  K <- tryCatch(
    suppressWarnings(phytools::phylosig(tr, x, method = "K")),
    error = function(e) NA_real_)
  K_val <- if (is.list(K)) K$K else K
  list(lambda = lambda, K = as.numeric(K_val))
}

#' Fritz & Purvis D statistic for a BINARY trait. D ≈ 0 under Brownian-
#' motion evolution; D ≈ 1 under random (phylogeny-free); D < 0 can
#' indicate more conserved than BM. Reference: Fritz & Purvis 2010,
#' Conservation Biology 24:1042.
#'
#' For multinomial traits we binarise OVR and average D across levels.
#' @return numeric: the D statistic (scalar for binary, mean(D) for
#'         multinomial), NA if caper unavailable or computation fails.
.phylo_signal_cat <- function(tree, values, species_names) {
  if (!requireNamespace("caper", quietly = TRUE))
    return(NA_real_)

  chr <- as.character(values)
  keep <- !is.na(chr) & nzchar(chr)
  chr <- chr[keep]
  sp  <- species_names[keep]
  if (length(chr) < 10L) return(NA_real_)
  classes <- unique(chr)
  if (length(classes) < 2L) return(NA_real_)

  # For each class, build a binary 0/1 trait and compute D. Average
  # across classes (mean across OVR binarisations; Münkemüller et al.
  # 2012 discuss the pragmatics of this for categorical data).
  Ds <- vapply(classes, function(cl) {
    binary <- as.integer(chr == cl)
    if (sum(binary) < 2 || sum(binary) > length(binary) - 2)
      return(NA_real_)       # caper::phylo.d fails on ~all-same
    df <- data.frame(Species = sp, y = binary)
    tr_cl <- ape::keep.tip(tree, sp)
    cd <- tryCatch(
      caper::comparative.data(tr_cl, df, names.col = "Species",
                              na.omit = FALSE, warn.dropped = FALSE),
      error = function(e) NULL)
    if (is.null(cd)) return(NA_real_)
    tryCatch(
      suppressMessages(caper::phylo.d(cd, binvar = y,
                                       permut = 100)$DEstimate),
      error = function(e) NA_real_)
  }, numeric(1))
  mean(Ds, na.rm = TRUE)
}

#' Prep a phylogeny for MCMCglmm:
#'   (a) Replace zero-length edges with a small positive value — these
#'       can appear after `ape::keep.tip` subsampling.
#'   (b) Re-force ultrametric: the zero-edge fix in (a) bumps individual
#'       tip edges up by epsilon and breaks ultrametricity, so we have
#'       to re-ultrametricize unconditionally afterwards.
#' MCMCglmm::inverseA refuses non-ultrametric trees or zero edges.
.prep_tree <- function(tree) {
  if (any(tree$edge.length <= 0)) {
    pos <- tree$edge.length[tree$edge.length > 0]
    eps <- if (length(pos)) min(pos) / 1e3 else 1e-6
    tree$edge.length[tree$edge.length <= 0] <- eps
  }
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
  tree
}

# Optionally subsample BEFORE ultrametricizing — pruning can introduce
# zero-length edges that .prep_tree then cleans up.
if (!is.na(SUBSET_N) && SUBSET_N < nrow(masked)) {
  keep_species <- sample(masked$Species, SUBSET_N)
  masked <- masked[masked$Species %in% keep_species, ]
  tree   <- ape::keep.tip(tree, keep_species)
}
tree <- .prep_tree(tree)

# Apply a reproducible CONT_MISS_RATE mask to the continuous traits and
# stash the true values in `truth` for downstream scoring. (Replaces
# the old pre-baked avonet_2000_masked.csv / avonet_2000_truth.csv
# fixtures.)
truth <- list()
for (v in CONT_TRAITS) {
  vals_full <- masked[[v]]
  known_idx <- which(!is.na(vals_full))
  n_mask    <- floor(length(known_idx) * CONT_MISS_RATE)
  mask_idx  <- sample(known_idx, n_mask)
  truth[[v]] <- data.frame(
    species_tip = masked$Species[mask_idx],
    trait       = v,
    true_value  = vals_full[mask_idx],
    stringsAsFactors = FALSE
  )
  masked[[v]][mask_idx] <- NA
}
truth <- do.call(rbind, truth)

# Log-transform the skewed continuous traits so bace() fits on a scale
# where Brownian / Gaussian residual assumptions are reasonable (see
# header for citations). truth$true_value stays on the raw scale; we
# handle the comparison on both scales below. log1p is used so exact
# zeros (rare but possible for Range.Size) don't produce -Inf.
for (v in LOG_TRAITS) {
  masked[[v]] <- log1p(masked[[v]])
}

# Apply a reproducible CAT_MISS_RATE mask to the categorical traits,
# preserving the true values in `cat_truth` for later scoring.
cat_truth <- list()
for (v in CAT_TRAITS) {
  vals_full <- masked[[v]]
  # If the source has its own NAs for this trait, keep them missing
  # (we can't score what we don't know).
  known_idx <- which(!is.na(vals_full))
  n_mask    <- floor(length(known_idx) * CAT_MISS_RATE)
  mask_idx  <- sample(known_idx, n_mask)
  cat_truth[[v]] <- data.frame(
    species_tip = masked$Species[mask_idx],
    trait       = v,
    true_value  = as.character(vals_full[mask_idx]),
    stringsAsFactors = FALSE
  )
  masked[[v]][mask_idx] <- NA
  # Force unordered factor for BACE's type detection.
  masked[[v]] <- factor(as.character(masked[[v]]))
}
cat_truth <- do.call(rbind, cat_truth)

# =============================================================================
# 4. RUN BANNER
# =============================================================================

ALL_TRAITS <- c(CONT_TRAITS, CAT_TRAITS)

# ---- Phylogenetic signal of each trait (pre-imputation) --------------------
# Compute Pagel's lambda + Blomberg's K on the observed (non-NA) values
# of each continuous trait, and Fritz & Purvis D for each categorical
# trait. These are REFERENCE VALUES for the dataset — they tell us how
# phylogenetically structured each trait is before we do any imputation.
# Stored in phylo_signal_df and also annotated on the diagnostic plots.
cat("Computing phylogenetic signal per trait (may take a minute)...\n")
phylo_signal_df <- do.call(rbind, c(
  lapply(CONT_TRAITS, function(v) {
    vals <- setNames(masked[[v]], masked$Species)
    s <- .phylo_signal_cont(tree, vals)
    data.frame(trait = v, lambda = s$lambda, K = s$K, D = NA_real_,
               stringsAsFactors = FALSE)
  }),
  lapply(CAT_TRAITS, function(v) {
    d <- .phylo_signal_cat(tree, masked[[v]], masked$Species)
    data.frame(trait = v, lambda = NA_real_, K = NA_real_, D = d,
               stringsAsFactors = FALSE)
  })
))
cat("\nPhylogenetic signal per trait:\n")
print(phylo_signal_df, digits = 3, row.names = FALSE)
cat("\n  lambda, K: Pagel 1999 / Blomberg et al. 2003 (continuous)\n")
cat("  D        : Fritz & Purvis 2010 (binary / multinomial OVR mean)\n")
cat("  Interpretation: lambda~1 / K~1 = BM-like; lambda~0 = no phylo\n")
cat("                  signal; D~0 = BM, D~1 = random, D<0 = conserved\n\n")

cat("===========================================================\n")
cat("  AVONET imputation benchmark\n")
cat("  Tag     :", RUN_TAG, "\n")
cat("  Output  :", OUT_DIR, "\n")
cat("  BACE    :", BACE_LOAD,
    if (BACE_LOAD == "library")
      paste0(" v", utils::packageVersion("BACE")) else "",
    "\n")
cat("  Species :", nrow(masked),
    if (is.na(SUBSET_N)) "(full)" else "(subset of 2000)", "\n")
cat("  Traits  :", length(ALL_TRAITS),
    "  (", length(CONT_TRAITS), "continuous +",
    length(CAT_TRAITS), "categorical)\n")
cat("  Masked  : continuous =", nrow(truth),
    "  categorical =", nrow(cat_truth), "\n")
cat("  MCMC    : nitt=", NITT, " thin=", THIN, " burnin=", BURNIN,
    " runs=", RUNS, " n_final=", N_FINAL, "\n", sep = "")
cat("===========================================================\n\n")

# =============================================================================
# 5. RUN BACE
# =============================================================================

t0 <- Sys.time()
fixformulas <- lapply(ALL_TRAITS, function(v) {
  others <- setdiff(ALL_TRAITS, v)
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

#' Per-trait accuracy + calibration metrics.
#'
#' For traits in LOG_TRAITS the imputed values (imp_mat) are already
#' in log space and the truth is raw; we compare in LOG space where the
#' model was fit (where Gaussian assumptions hold) and also report
#' raw-scale MAE using the back-transformed posterior mean, for
#' interpretability. For non-log traits we just compare in raw space.
#'
#' Point-estimate metrics (NRMSE, MAE, cor) are computed per imputation
#' then averaged. coverage95 uses the full n_final ensemble.
compute_avonet_metrics <- function(res, masked, truth, trait_cols,
                                    log_traits = character(0)) {
  rows <- lapply(trait_cols, function(v) {
    is_log   <- v %in% log_traits
    truth_v  <- truth[truth$trait == v, ]
    idx      <- match(truth_v$species_tip, masked$Species)
    tv_raw   <- truth_v$true_value

    # Build the imputation matrix in whatever scale bace() returned
    # (log space for LOG_TRAITS, raw for others — matches masked).
    imp_mat <- vapply(res$imputed_datasets,
                      function(d) as.numeric(d[[v]][idx]),
                      FUN.VALUE = numeric(length(tv_raw)))
    if (!is.matrix(imp_mat)) imp_mat <- matrix(imp_mat, nrow = length(tv_raw))

    # Comparison-scale truth (log1p(tv_raw) for LOG_TRAITS, else tv_raw).
    tv_cmp  <- if (is_log) log1p(tv_raw) else tv_raw
    sd_full <- sd(masked[[v]], na.rm = TRUE)  # marginal sd on fit scale

    per_imp_nrmse <- apply(imp_mat, 2,
                            function(iv) sqrt(mean((iv - tv_cmp)^2)) / sd_full)
    per_imp_mae_cmp <- apply(imp_mat, 2, function(iv) mean(abs(iv - tv_cmp)))
    per_imp_cor <- apply(imp_mat, 2, function(iv)
      if (length(unique(tv_cmp)) > 1 && length(unique(iv)) > 1)
        suppressWarnings(cor(iv, tv_cmp)) else NA_real_)

    # Raw-scale MAE: back-transform the posterior mean then compare to
    # raw truth. Gives an interpretable "grams" / "km^2" error for
    # log-transformed traits. Non-log traits: same as the fit-scale MAE.
    post_mean  <- rowMeans(imp_mat)
    post_mean_raw <- if (is_log) expm1(post_mean) else post_mean
    mae_raw    <- mean(abs(post_mean_raw - tv_raw))

    # Coverage on fit scale (the scale MCMCglmm residuals live on).
    lo <- apply(imp_mat, 1, stats::quantile, probs = 0.025,
                na.rm = TRUE, names = FALSE)
    hi <- apply(imp_mat, 1, stats::quantile, probs = 0.975,
                na.rm = TRUE, names = FALSE)

    data.frame(
      trait       = v,
      scale       = if (is_log) "log" else "raw",
      n_hidden    = length(tv_raw),
      nrmse       = mean(per_imp_nrmse, na.rm = TRUE),
      mae_fit     = mean(per_imp_mae_cmp),        # on fit scale (log or raw)
      mae_raw     = mae_raw,                      # always on raw scale
      correlation = mean(per_imp_cor, na.rm = TRUE),
      coverage95  = mean(tv_cmp >= lo & tv_cmp <= hi),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Per-trait categorical metrics (accuracy, balanced accuracy, Brier).
#' Consumes majority-vote class assignments across n_final imputations
#' as the point estimate; Brier uses the empirical per-class frequency
#' across the n_final ensemble.
compute_avonet_cat_metrics <- function(res, masked, cat_truth, cat_traits) {
  rows <- lapply(cat_traits, function(v) {
    tv_df  <- cat_truth[cat_truth$trait == v, ]
    idx    <- match(tv_df$species_tip, masked$Species)
    tv_chr <- as.character(tv_df$true_value)

    im <- vapply(res$imputed_datasets,
                 function(d) as.character(d[[v]][idx]),
                 FUN.VALUE = character(length(tv_chr)))
    if (!is.matrix(im)) im <- matrix(im, nrow = length(tv_chr))

    # Majority-vote point estimate per cell.
    vote <- apply(im, 1, function(row) {
      row <- row[!is.na(row)]
      if (length(row) == 0) return(NA_character_)
      tab <- table(row)
      top <- names(tab)[tab == max(tab)]
      if (length(top) == 1L) top else sample(top, 1L)
    })

    accuracy <- mean(vote == tv_chr, na.rm = TRUE)

    classes <- unique(tv_chr)
    recalls <- vapply(classes, function(cl) {
      is_cl <- tv_chr == cl
      if (!any(is_cl)) NA_real_ else mean(vote[is_cl] == cl, na.rm = TRUE)
    }, numeric(1))
    bal_acc <- mean(recalls, na.rm = TRUE)

    # Multiclass Brier from empirical n_final class frequencies.
    lvls <- if (is.factor(masked[[v]])) levels(masked[[v]]) else sort(unique(tv_chr))
    prob_mat <- sapply(lvls, function(cl) rowMeans(im == cl, na.rm = TRUE))
    y_indic  <- sapply(lvls, function(cl) as.integer(tv_chr == cl))
    brier    <- mean(rowSums((prob_mat - y_indic)^2))

    data.frame(
      trait             = v,
      scale             = "categorical",
      n_hidden          = length(tv_chr),
      accuracy          = accuracy,
      balanced_accuracy = bal_acc,
      brier             = brier,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summary_cont <- compute_avonet_metrics(res, masked, truth, CONT_TRAITS,
                                        log_traits = LOG_TRAITS)
summary_cat  <- compute_avonet_cat_metrics(res, masked, cat_truth, CAT_TRAITS)

cat("\n===========================================================\n")
cat("  AVONET results - continuous traits\n")
cat("===========================================================\n")
print(summary_cont, digits = 3, row.names = FALSE)

cat("\n===========================================================\n")
cat("  AVONET results - categorical traits\n")
cat("===========================================================\n")
print(summary_cat, digits = 3, row.names = FALSE)

cat("\nAverages across continuous traits (NRMSE / cor / coverage):\n")
avg_cont <- data.frame(
  trait = "(cont. mean)",
  n_hidden = sum(summary_cont$n_hidden),
  nrmse = mean(summary_cont$nrmse, na.rm = TRUE),
  correlation = mean(summary_cont$correlation, na.rm = TRUE),
  coverage95 = mean(summary_cont$coverage95, na.rm = TRUE)
)
print(avg_cont, digits = 3, row.names = FALSE)

cat("\nAverages across categorical traits (accuracy / Brier):\n")
avg_cat <- data.frame(
  trait = "(cat. mean)",
  n_hidden = sum(summary_cat$n_hidden),
  accuracy = mean(summary_cat$accuracy, na.rm = TRUE),
  balanced_accuracy = mean(summary_cat$balanced_accuracy, na.rm = TRUE),
  brier = mean(summary_cat$brier, na.rm = TRUE)
)
print(avg_cat, digits = 3, row.names = FALSE)

# Bind continuous + categorical into a unified summary for saving.
# Columns are a superset so summary_cat gets NA for continuous-only
# columns and vice versa.
all_cols <- union(names(summary_cont), names(summary_cat))
pad <- function(df, cols) {
  for (c in setdiff(cols, names(df))) df[[c]] <- NA
  df[, cols, drop = FALSE]
}
summary_df <- rbind(pad(summary_cont, all_cols), pad(summary_cat, all_cols))

# Join phylogenetic-signal columns so each trait row carries its
# lambda / K / D alongside the imputation-accuracy metrics.
summary_df <- merge(summary_df, phylo_signal_df,
                    by = "trait", all.x = TRUE, sort = FALSE)

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
  cont_traits     = CONT_TRAITS,
  log_traits      = LOG_TRAITS,
  cat_traits      = CAT_TRAITS,
  cat_miss_rate   = CAT_MISS_RATE,
  mcmc            = list(nitt = NITT, thin = THIN, burnin = BURNIN,
                         runs = RUNS, n_final = N_FINAL,
                         max_attempts = MAX_ATTEMPTS),
  runtime_min     = runtime_min,
  converged       = res$converged,
  n_attempts      = res$n_attempts,
  seed            = 2026
)

saveRDS(list(metadata        = metadata,
             summary         = summary_df,
             phylo_signal    = phylo_signal_df,
             imputed_datasets = res$imputed_datasets),
        file.path(OUT_DIR, "results.rds"))
write.csv(phylo_signal_df, file.path(OUT_DIR, "phylo_signal.csv"),
          row.names = FALSE)
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

# ---- Page 1: continuous traits -- true vs imputed, with 95% PI bars.
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
for (v in CONT_TRAITS) {
  is_log  <- v %in% LOG_TRAITS
  truth_v <- truth[truth$trait == v, ]
  idx     <- match(truth_v$species_tip, masked$Species)
  tv_raw  <- truth_v$true_value
  tv_cmp  <- if (is_log) log1p(tv_raw) else tv_raw

  imp_mat <- vapply(res$imputed_datasets,
                    function(d) as.numeric(d[[v]][idx]),
                    FUN.VALUE = numeric(length(tv_raw)))
  if (!is.matrix(imp_mat)) imp_mat <- matrix(imp_mat, nrow = length(tv_raw))
  mid <- rowMeans(imp_mat)
  lo  <- apply(imp_mat, 1, stats::quantile, 0.025, na.rm = TRUE, names = FALSE)
  hi  <- apply(imp_mat, 1, stats::quantile, 0.975, na.rm = TRUE, names = FALSE)

  xl <- range(c(tv_cmp, lo, hi), na.rm = TRUE)
  plot(tv_cmp, mid, pch = 19, col = adjustcolor("steelblue", 0.4),
       xlim = xl, ylim = xl,
       xlab = if (is_log) paste0("log1p(", v, ") true") else "True value",
       ylab = "Imputed (posterior mean)",
       main = paste0(v, if (is_log) " (log)" else ""), cex = 0.5)
  segments(tv_cmp, lo, tv_cmp, hi, col = adjustcolor("steelblue", 0.15))
  abline(0, 1, col = "red", lty = 2)
  mt <- summary_cont[summary_cont$trait == v, ]
  ps <- phylo_signal_df[phylo_signal_df$trait == v, ]
  legend("topleft", bty = "n", cex = 0.8,
         legend = c(sprintf("NRMSE=%.2f r=%.2f", mt$nrmse, mt$correlation),
                    sprintf("cov95=%.0f%%", mt$coverage95 * 100),
                    sprintf("lambda=%.2f K=%.2f",
                             ps$lambda, ps$K)))
}
mtext(paste("AVONET continuous traits  -", RUN_TAG,
            "  (lambda, K shown per trait)"),
      outer = TRUE, cex = 1.1)

# ---- Page 2: categorical traits -- confusion heatmaps.
if (length(CAT_TRAITS) > 0) {
  par(mfrow = c(1, length(CAT_TRAITS)), mar = c(5, 6, 3, 1),
      oma = c(0, 0, 2, 0))
  for (v in CAT_TRAITS) {
    tv_df <- cat_truth[cat_truth$trait == v, ]
    idx   <- match(tv_df$species_tip, masked$Species)
    tv_chr <- as.character(tv_df$true_value)
    im <- vapply(res$imputed_datasets,
                 function(d) as.character(d[[v]][idx]),
                 FUN.VALUE = character(length(tv_chr)))
    if (!is.matrix(im)) im <- matrix(im, nrow = length(tv_chr))
    vote <- apply(im, 1, function(row) {
      row <- row[!is.na(row)]
      if (length(row) == 0) NA_character_ else
        names(which.max(table(row)))
    })
    lvls <- if (is.factor(masked[[v]])) levels(masked[[v]])
            else sort(unique(tv_chr))
    tab     <- table(factor(tv_chr, lvls), factor(vote, lvls))
    rowprop <- sweep(tab, 1, pmax(rowSums(tab), 1), FUN = "/")

    image(x = seq_along(lvls), y = seq_along(lvls),
          z = t(rowprop[nrow(rowprop):1, , drop = FALSE]),
          col = grDevices::colorRampPalette(c("white", "#1b7837"))(64),
          zlim = c(0, 1), axes = FALSE,
          xlab = "Imputed class", ylab = "True class",
          main = sprintf("%s\n(accuracy=%.0f%%, D=%.2f)", v,
                          summary_cat$accuracy[summary_cat$trait == v] * 100,
                          phylo_signal_df$D[phylo_signal_df$trait == v]))
    axis(1, seq_along(lvls), lvls, las = 2, cex.axis = 0.8)
    axis(2, seq_along(lvls), rev(lvls), las = 1, cex.axis = 0.8)
    box()
    for (i in seq_along(lvls)) for (j in seq_along(lvls)) {
      count <- tab[i, j]; pct <- rowprop[i, j] * 100
      text(j, length(lvls) - i + 1,
           sprintf("%d\n(%.0f%%)", count, pct),
           col = if (pct > 50) "white" else "black", cex = 0.7)
    }
  }
  mtext("AVONET categorical traits - confusion (row %)",
        outer = TRUE, cex = 1.1)
}

# ---- Page 3: per-trait scorecard. Coverage (cont.) + accuracy (cat.)
par(mfrow = c(1, 2), mar = c(8, 4, 3, 1), oma = c(0, 0, 0, 0))
bp <- barplot(summary_cont$coverage95, names.arg = summary_cont$trait,
              ylim = c(0, 1), col = "#6a5acd", border = NA,
              main = "Continuous: 95% PI coverage",
              ylab = "Coverage", las = 2)
abline(h = 0.95, lty = 2, col = "red")
text(bp, summary_cont$coverage95,
     sprintf("%.2f", summary_cont$coverage95), pos = 3, cex = 0.8, xpd = NA)

if (length(CAT_TRAITS) > 0) {
  bp <- barplot(summary_cat$accuracy, names.arg = summary_cat$trait,
                ylim = c(0, 1), col = "#1b7837", border = NA,
                main = "Categorical: accuracy",
                ylab = "Accuracy", las = 2)
  text(bp, summary_cat$accuracy,
       sprintf("%.2f", summary_cat$accuracy), pos = 3, cex = 0.8, xpd = NA)
}

# ---- Page 4: phylogenetic signal vs imputation accuracy --------------------
# A biologically interesting diagnostic: does higher phylogenetic
# signal translate to better imputation? Continuous traits plot
# lambda vs correlation; categorical traits plot D vs accuracy.
# Included to answer the question the Münkemüller et al. 2012 and
# Nakagawa et al. 2022 reviews discuss — signal strength is
# conceptually linked to what phylogenetic imputation can recover.
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
cont_join <- merge(summary_cont, phylo_signal_df, by = "trait")
if (any(!is.na(cont_join$correlation) & !is.na(cont_join$lambda))) {
  plot(cont_join$lambda, cont_join$correlation, pch = 19,
       col = "steelblue", xlim = c(0, 1.05), ylim = c(-0.05, 1.05),
       xlab = "Pagel's lambda (pre-imputation)",
       ylab = "Correlation(imputed, true)",
       main = "Continuous: signal vs imputation quality")
  text(cont_join$lambda, cont_join$correlation,
       labels = cont_join$trait, pos = 3, cex = 0.7)
  abline(0, 1, col = "grey60", lty = 3)
}
if (length(CAT_TRAITS) > 0) {
  cat_join <- merge(summary_cat, phylo_signal_df, by = "trait")
  plot(cat_join$D, cat_join$accuracy, pch = 19, col = "#1b7837",
       xlim = c(-0.5, 1.5), ylim = c(0, 1.05),
       xlab = "Fritz & Purvis D (0=BM, 1=random, <0=conserved)",
       ylab = "Imputation accuracy",
       main = "Categorical: signal vs imputation quality")
  text(cat_join$D, cat_join$accuracy,
       labels = cat_join$trait, pos = 3, cex = 0.8)
  abline(v = c(0, 1), col = "grey60", lty = 3)
}

dev.off()

cat("\nOutputs written to:", OUT_DIR, "\n")
cat("  - results.rds    full output + imputed datasets\n")
cat("  - summary.csv    per-trait metrics\n")
cat("  - metadata.json  run metadata (BACE version, MCMC settings, etc.)\n")
cat("  - plots.pdf      diagnostic plots\n\n")
cat("Done.\n")
