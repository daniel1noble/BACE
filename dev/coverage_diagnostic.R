# =============================================================================
# Coverage diagnostic: theoretical PPD variance vs empirical imputation variance
# =============================================================================
#
# The question: when BACE's 95% PI across n_final imputations under-covers
# the true value, is it because:
#   (A) the n_final imputations' empirical variance is smaller than the
#       posterior predictive variance MCMCglmm estimates, i.e. the
#       sampling path is losing something (median-of-3, bad iteration
#       selection, etc.), OR
#   (B) the posterior predictive distribution itself is too narrow, i.e.
#       the fit is under-estimating residual + parameter variance?
#
# This script answers that by:
#   1. Loading the most recent 01_benchmark_simulated.R run's artefacts
#      (which include the fitted MCMCglmm models and the n_final imputations).
#   2. For each hidden cell of y (gaussian trait) :
#        a. "theoretical PPD": draw many samples of (eta + eps) directly
#           from the MCMCglmm posterior chain (Sol[i,] + VCV[i,"units"]).
#           Compute the sd across those samples.
#        b. "empirical imputation": compute the sd across the n_final
#           values in bace_result$imputed_datasets.
#   3. Compare per cell and on average.
#
# If (empirical sd << theoretical sd), the sampling path is the bottleneck
# and fixing it recovers coverage. If they match, the posterior itself
# is the issue and we need to widen priors or get more data.
# =============================================================================

devtools::load_all(quiet = TRUE)
library(ape); library(MASS)

art <- readRDS("dev/simulation_results/demo/demo_artefacts.rds")
bace_result   <- art$bace_result
complete_data <- art$sim$complete_data
names(complete_data)[names(complete_data) == "species"] <- "Species"

# Reproduce the mask used in STEP 3 of 01_benchmark_simulated.R. The
# demo does set.seed(42) at the top, then STEP 2 calls
# inject_missingness() three times (phylo_MAR, trait_MAR, trait_MNAR)
# for the mechanism diagnostics before STEP 3's single inject for the
# bace() run. Reproducing that exact sequence is brittle; easier to
# just reuse the mask stored in art$masks$trait_MAR which was
# generated with the SAME seed but on a different (low-signal) sim.
#
# Cleanest path: rerun inject_missingness with the same seed + sim
# and capture the mask. Deterministic replay.
set.seed(42)
SIGNAL <- list(high = 0.90, moderate = 0.60, low = 0.20)
DEP_STRENGTH <- 1.5
# sim_bace uses a lot of the RNG stream when constructing the two sims
# in STEPS 1-2; easiest to re-run the same code path exactly.
source_funcs <- function(path, fns) {
  src <- readLines(path)
  for (fn in fns) {
    s <- grep(paste0("^", fn, " <- function"), src)[1]
    d <- 0; e <- NA
    for (i in s:length(src)) {
      d <- d + sum(gregexpr("\\{", src[i])[[1]] > 0) -
               sum(gregexpr("\\}", src[i])[[1]] > 0)
      if (i > s && d == 0) { e <- i; break }
    }
    eval(parse(text = src[s:e]), envir = globalenv())
  }
}
source_funcs("dev/02_benchmark_simulated_full.R",
             c(".trait_to_numeric", ".calibrate_intercept",
               "inject_missingness", "make_phylo_signal"))

# Exactly mirror STEP 1 + STEP 2 of 01_benchmark_simulated.R up to the
# trait_MAR call in STEP 3.
phylo_signal <- make_phylo_signal("all_high")
sim_hi <- sim_bace(response_type = "gaussian",
                   predictor_types = c("binary","multinomial3","poisson","threshold3"),
                   var_names = c("y","x1","x2","x3","x4"),
                   phylo_signal = phylo_signal,
                   n_cases = 150, n_species = 60, beta_sparsity = 0.3,
                   missingness = c(0,0,0,0,0))
cd <- sim_hi$complete_data
names(cd)[names(cd) == "species"] <- "Species"

# STEP 2: low-signal sim used for mechanism diagnostic (same seed sequence).
sim_lo <- sim_bace(response_type = "gaussian",
                   predictor_types = c("binary","multinomial3","poisson","threshold3"),
                   var_names = c("y","x1","x2","x3","x4"),
                   phylo_signal = make_phylo_signal("all_low"),
                   n_cases = 150, n_species = 60, beta_sparsity = 0.3,
                   missingness = c(0,0,0,0,0))
cd_lo <- sim_lo$complete_data
names(cd_lo)[names(cd_lo) == "species"] <- "Species"
# Three mechanism masks in the order done in STEP 2.
.x1 <- inject_missingness(cd_lo, sim_lo$tree, "phylo_MAR",  0.35)
.x2 <- inject_missingness(cd_lo, sim_lo$tree, "trait_MAR",  0.35)
.x3 <- inject_missingness(cd_lo, sim_lo$tree, "trait_MNAR", 0.35)

# STEP 3: actual bace() mask.
miss <- inject_missingness(cd, sim_hi$tree, "trait_MAR", 0.35)
miss_mask <- miss$miss_mask
complete_data <- cd   # use this one consistently

# The fitted final-imputation models are on bace_result$final_results$all_models.
# Each entry is one of n_final final-run models per variable. We'll use run 1
# for illustration; the comparison holds regardless of which run.
final_models <- bace_result$final_results$all_models
y_model_idx  <- grep("^y\\b", names(final_models[[1]]))
if (length(y_model_idx) == 0) y_model_idx <- 1   # fall back
y_models <- lapply(final_models, function(run) run[[y_model_idx]])

# For the first run, extract posterior chains.
m <- y_models[[1]]
Sol <- as.matrix(m$Sol)
VCV <- as.matrix(m$VCV)
X   <- as.matrix(m$X)
Z   <- if (!is.null(m$Z)) as.matrix(m$Z) else NULL
W   <- if (!is.null(Z)) cbind(X, Z) else X
common <- intersect(colnames(W), colnames(Sol))

# Hidden cells of y.
idx_hidden <- which(miss_mask$y)
n_hidden   <- length(idx_hidden)

# ------------------------------------------------------------------
# (A) theoretical PPD: n_theoretical draws per hidden cell, from the
#     MCMCglmm posterior samples of (Sol, VCV).
# ------------------------------------------------------------------
n_theoretical <- min(500, nrow(Sol))
i_samps <- sample.int(nrow(Sol), n_theoretical)

# Theoretical draws on z-score scale (that's the scale MCMCglmm fits).
theoretical_z <- matrix(NA_real_, nrow = n_hidden, ncol = n_theoretical)
for (k in seq_len(n_theoretical)) {
  i <- i_samps[k]
  eta <- as.numeric(Sol[i, common, drop = FALSE] %*% t(W[idx_hidden, common, drop = FALSE]))
  eps <- rnorm(n_hidden, 0, sqrt(VCV[i, "units"]))
  theoretical_z[, k] <- eta + eps
}

# Back-transform to raw scale using the complete y's empirical mean/sd
# (that's what .data_prep used internally).
y_mean <- mean(complete_data$y)
y_sd   <- sd(complete_data$y)
theoretical_raw <- theoretical_z * y_sd + y_mean

sd_theoretical_z   <- apply(theoretical_z,   1, sd)
sd_theoretical_raw <- apply(theoretical_raw, 1, sd)

# Empirical draws across n_final imputations (raw scale already).
empirical_raw <- vapply(bace_result$imputed_datasets,
                        function(d) d$y[idx_hidden],
                        FUN.VALUE = numeric(n_hidden))
if (!is.matrix(empirical_raw)) empirical_raw <- matrix(empirical_raw, nrow = n_hidden)
sd_empirical <- apply(empirical_raw, 1, sd)

cat("sd(complete y, raw) =", round(y_sd, 3), "\n")
cat("sd(theoretical PPD, z-scale)   =", round(mean(sd_theoretical_z), 4), "\n")
cat("sd(theoretical PPD, raw-scale) =", round(mean(sd_theoretical_raw), 3),
    " (marginal y sd x theoretical z-sd)\n")
cat("sd(empirical imputations, raw) =", round(mean(sd_empirical), 3), "\n\n")

sd_theoretical <- sd_theoretical_raw
theoretical_draws_raw <- theoretical_raw

# ------------------------------------------------------------------
# Compare.
# ------------------------------------------------------------------
cat("=============================================================\n")
cat("  PPD variance diagnostic -- y (gaussian) hidden cells\n")
cat("=============================================================\n\n")
cat("n_hidden = ", n_hidden, "\n", sep = "")
cat("n_theoretical = ", n_theoretical, " samples per cell from model$Sol / VCV\n", sep = "")
cat("n_final       = ", ncol(empirical_draws),
    " imputations per cell (from imputed_datasets)\n\n", sep = "")

summary_df <- data.frame(
  quantity      = c("sd_theoretical (per-cell)",
                     "sd_empirical (per-cell)",
                     "ratio empirical / theoretical"),
  mean          = c(mean(sd_theoretical, na.rm = TRUE),
                     mean(sd_empirical, na.rm = TRUE),
                     mean(sd_empirical / sd_theoretical, na.rm = TRUE)),
  median        = c(median(sd_theoretical, na.rm = TRUE),
                     median(sd_empirical, na.rm = TRUE),
                     median(sd_empirical / sd_theoretical, na.rm = TRUE))
)
print(summary_df, digits = 3, row.names = FALSE)

cat("\nInterpretation:\n")
cat("  ratio ~ 1.0  -> empirical sd matches theoretical PPD; sampling\n")
cat("                  path is faithful. Coverage gap must come from\n")
cat("                  other sources (posterior too narrow in the fit).\n")
cat("  ratio < 1.0  -> empirical sd is SMALLER than the theoretical PPD.\n")
cat("                  Sampling path (e.g. median-of-3) is eating variance.\n")
cat("                  Expected ~0.72 under median-of-3; ~0.52 under\n")
cat("                  median-of-7; ~1.0 under single-draw.\n")

# Also check how often the TRUTH falls inside the theoretical 95% PI
# vs the empirical 95% PI.
true_y <- complete_data$y[idx_hidden]
th_lo <- apply(theoretical_draws_raw, 1, quantile, 0.025)
th_hi <- apply(theoretical_draws_raw, 1, quantile, 0.975)
emp_lo <- apply(empirical_draws, 1, quantile, 0.025)
emp_hi <- apply(empirical_draws, 1, quantile, 0.975)
cat("\nCoverage of true y (cells hidden by trait_MAR):\n")
cat("  theoretical 95% PI coverage = ",
    round(mean(true_y >= th_lo & true_y <= th_hi), 3), "\n", sep = "")
cat("  empirical 95% PI coverage   = ",
    round(mean(true_y >= emp_lo & true_y <= emp_hi), 3), "\n", sep = "")
