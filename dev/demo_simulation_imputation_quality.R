# =============================================================================
# Demonstration: What `simulation_imputation_quality.R` Actually Does
# =============================================================================
#
# Purpose
# -------
# The main script runs ~1500 BACE imputations headlessly (4 phylo scenarios
# x 3 missingness mechanisms x 125 replicates, with n_final = 20 inner
# imputations each). This file runs the same machinery on a SINGLE
# replicate, printing every intermediate artefact so a reader can
# confirm:
#
#   1. sim_bace() produces a complete dataset with known ground truth
#      (including an ordered categorical trait).
#   2. inject_missingness() generates masks whose empirical pattern
#      matches the intended mechanism (phylo-clustered, trait-driven,
#      value-driven).
#   3. bace() returns n_final imputed datasets with every NA filled in,
#      and reports convergence.
#   4. The accuracy metrics (NRMSE, MAE, log-MAE, correlation, accuracy,
#      balanced accuracy, ordinal MAE) match a hand calculation.
#   5. The calibration metrics (95% PI coverage for continuous/count;
#      Brier score for categorical/ordered) behave sensibly.
#
# Design references
#   Rubin 1976, Biometrika 63:581 - MCAR/MAR/MNAR framework.
#   Penone et al. 2014, MEE 5:961 - phylogenetic clustering of
#     missingness in trait databases.
#   van Buuren 2018, Flexible Imputation of Data 2e - pooled MI
#     calibration (§2.5), number of imputations (§2.8).
#   Gneiting & Raftery 2007, JASA 102:359 - Brier score as a proper
#     scoring rule for categorical forecasts.
# =============================================================================

library(BACE)
library(ape)
library(MASS)

set.seed(42)

DEMO_DIR <- file.path("dev", "simulation_results", "demo")
if (!dir.exists(DEMO_DIR)) dir.create(DEMO_DIR, recursive = TRUE)

# Source just the helpers from the main script so the demo reuses the
# same code paths end-to-end. This isolates the functions we need
# (inject_missingness, evaluate_imputation_ensemble, make_phylo_signal)
# without triggering the main loop.
source_funcs <- function(path, fn_names) {
  src <- readLines(path)
  for (fn in fn_names) {
    start_i <- grep(paste0("^", fn, " <- function"), src)[1]
    if (is.na(start_i)) stop("Function not found in ", path, ": ", fn)
    depth <- 0; end_i <- NA
    for (i in start_i:length(src)) {
      depth <- depth + sum(gregexpr("\\{", src[i])[[1]] > 0) -
                       sum(gregexpr("\\}", src[i])[[1]] > 0)
      if (i > start_i && depth == 0) { end_i <- i; break }
    }
    eval(parse(text = src[start_i:end_i]), envir = globalenv())
  }
}
MAIN <- file.path("dev", "simulation_imputation_quality.R")

# We also need the constants the helpers reference. Keep in sync with
# the main script - see inject_missingness() for the rationale.
SIGNAL       <- list(high = 0.90, moderate = 0.60, low = 0.20)
DEP_STRENGTH <- 1.5
source_funcs(MAIN, c(".trait_to_numeric", ".calibrate_intercept",
                     "inject_missingness", "evaluate_imputation_ensemble",
                     "make_phylo_signal"))

# -----------------------------------------------------------------------------
# STEP 1. Simulate complete data (all-high signal for a clean demo)
# -----------------------------------------------------------------------------
# Five variables, one of each type. All-high phylo signal keeps the demo
# focused on the MECHANISM axis rather than on signal effects.
phylo_signal <- make_phylo_signal("all_high")

sim <- sim_bace(
  response_type   = "gaussian",
  predictor_types = c("binary", "multinomial3", "poisson", "threshold3"),
  var_names       = c("y", "x1", "x2", "x3", "x4"),
  phylo_signal    = phylo_signal,
  n_cases         = 150,
  n_species       = 60,
  beta_sparsity   = 0.3,
  missingness     = c(0, 0, 0, 0, 0)
)
complete_data <- sim$complete_data
tree          <- sim$tree
names(complete_data)[names(complete_data) == "species"] <- "Species"

cat("\n--- STEP 1: Simulated complete data (all-high signal) ---\n")
cat("Rows x cols :", nrow(complete_data), "x", ncol(complete_data), "\n")
cat("Tree tips   :", ape::Ntip(tree), "\n")
cat("Column types:\n"); print(sapply(complete_data, function(x) class(x)[1]))
stopifnot(sum(is.na(complete_data)) == 0)

# -----------------------------------------------------------------------------
# STEP 2. Demonstrate each missingness MECHANISM using a LOW-signal sim
# -----------------------------------------------------------------------------
# Why a separate low-signal sim here? With all-high phylo signal, every
# trait is itself strongly clade-clustered, so trait_MAR and trait_MNAR
# will ALSO produce phylo-clustered masks (the clustering is inherited
# from the driver trait, not from the phylogeny directly). To cleanly
# separate "mechanism is phylogenetic" from "mechanism is trait-based",
# we use an all-low-signal sim here, then resume the all-high sim for
# the bace pipeline.
sim_lo <- sim_bace(
  response_type   = "gaussian",
  predictor_types = c("binary", "multinomial3", "poisson", "threshold3"),
  var_names       = c("y", "x1", "x2", "x3", "x4"),
  phylo_signal    = make_phylo_signal("all_low"),
  n_cases         = 150, n_species = 60, beta_sparsity = 0.3,
  missingness     = c(0, 0, 0, 0, 0)
)
cd_lo <- sim_lo$complete_data
names(cd_lo)[names(cd_lo) == "species"] <- "Species"

cat("\n--- STEP 2: Missingness masks under three mechanisms",
    "(low-signal sim) ---\n")

masks <- list(
  phylo_MAR  = inject_missingness(cd_lo, sim_lo$tree, "phylo_MAR",  0.35),
  trait_MAR  = inject_missingness(cd_lo, sim_lo$tree, "trait_MAR",  0.35),
  trait_MNAR = inject_missingness(cd_lo, sim_lo$tree, "trait_MNAR", 0.35)
)

cat("\nRealised per-variable missingness rates (target = 35%):\n")
for (m in names(masks)) {
  cat("  ", sprintf("%-11s", m), ":",
      paste0(round(colMeans(masks[[m]]$miss_mask) * 100), "%",
             collapse = " "), "\n")
}

# Diagnostic 1: species-level ICC of the mask. With LOW trait signal,
# traits are not themselves clade-clustered, so only phylo_MAR should
# give a high ICC. Interpret ICC loosely: 0 = no clustering, 1 = every
# species is either all-missing or all-observed.
species_icc <- function(mask, species) {
  vapply(mask, function(m) {
    agg <- aggregate(m ~ species, FUN = mean)$m
    if (var(m) == 0) NA_real_ else var(agg) / var(m)
  }, numeric(1))
}
cat("\nSpecies-level ICC of the mask (expect HIGH only for phylo_MAR):\n")
for (m in names(masks)) {
  cat("  ", sprintf("%-11s", m), ":",
      paste(sprintf("%s=%.2f", names(masks[[m]]$miss_mask),
                    species_icc(masks[[m]]$miss_mask, cd_lo$Species)),
            collapse = " "), "\n")
}

# Diagnostic 2: how does y's miss-probability relate to (a) y's value
# itself and (b) y's designed driver x3?
#   trait_MAR (driver = x3): cor(p, x3) should be strongly NEGATIVE
#     (low x3 -> high p); cor(p, y) only as strong as cor(y, x3).
#   trait_MNAR: cor(p, y) should be strongly NEGATIVE; cor(p, x3) only
#     inherited via cor(y, x3).
#   phylo_MAR: both should be near zero (p is driven by the tree, not
#     by the trait values).
cor_yx3 <- cor(cd_lo$y, cd_lo$x3)
cat("\nBackground cor(y, x3) in this sim:", round(cor_yx3, 2), "\n\n")
cat("Diagnostic correlations for p_miss(y):\n")
cat("  mechanism    cor(p, y)   cor(p, x3)\n")
for (m in names(masks)) {
  p  <- masks[[m]]$miss_probs$y
  cy <- cor(p, cd_lo$y)
  cx <- cor(p, cd_lo$x3)
  cat("  ", sprintf("%-11s  %+0.2f        %+0.2f", m, cy, cx), "\n")
}
cat("Interpretation: trait_MAR drives p(y) via x3, trait_MNAR via y,",
    "\nphylo_MAR via neither.\n")

# -----------------------------------------------------------------------------
# STEP 3. Pick one mechanism (trait_MAR) and run bace() end-to-end
# -----------------------------------------------------------------------------
# We use trait_MAR here rather than phylo_MAR for the pipeline walk-
# through. Reason: phylo_MAR's clade-clustered masks occasionally produce
# Species-factor configurations that hit a ginverse edge case inside
# MCMCglmm. run_one_sim() in the main script catches such errors and
# discards the replicate; a single-sim demo is not robust to them, so
# we pick the mechanism least likely to trigger it.
#
# MCMC budget
# -----------
# The settings below are deliberately at the "production" end for a
# single-replicate demo (i.e. several minutes of runtime) so that
# coverage / calibration results are NOT confounded by:
#
#   - Too few imputations: Rubin (1987) and Graham, Olchowski & Gilreath
#     (2007, Prev Sci 8:206) recommend m >= 20 at our fraction-of-
#     missing-information (~35%). With m = 6 the sample 2.5/97.5
#     quantiles are basically min/max and cannot achieve 95% coverage
#     by construction.
#
#   - Too-short chains: Hadfield (2010, JSS 33(2)) and Gelman & Rubin
#     (1992, Stat Sci 7:457) recommend an effective sample size of
#     >= 1000 for reliable posterior summaries. nitt = 50000, burnin
#     = 10000, thin = 25 gives 1600 retained draws per formula per
#     chain, which is safely above that threshold for well-mixed
#     gaussian/poisson models.
#
# If coverage is still far from the n_final-imposed ceiling (~90% with
# m = 20) after these settings, the remaining gap points at the
# .predict_bace sampling path (posterior-mean vs posterior-predictive
# draws) rather than at budget, and THAT becomes a BACE fix.
cat("\n--- STEP 3: bace() under trait_MAR (production MCMC budget) ---\n")
miss <- inject_missingness(complete_data, tree, "trait_MAR", 0.35)
miss_data <- miss$miss_data
miss_mask <- miss$miss_mask

for (v in c("y","x1","x2","x3","x4")) {
  stopifnot(all(is.na(miss_data[[v]]) == miss_mask[[v]]))
}
cat("Invariant OK: NA positions in miss_data exactly match mask.\n")

t0 <- Sys.time()
bace_result <- bace(
  fixformula     = list("y  ~ x1 + x2 + x3 + x4",
                        "x1 ~ y  + x2 + x3 + x4",
                        "x2 ~ y  + x1 + x3 + x4",
                        "x3 ~ y  + x1 + x2 + x4",
                        "x4 ~ y  + x1 + x2 + x3"),
  ran_phylo_form = "~1|Species",
  phylo          = tree,
  data           = miss_data,
  # Production-grade MCMC (see block above for rationale):
  # nitt/thin gives 1600 retained draws per formula per chain.
  nitt           = 50000, thin = 25, burnin = 10000,
  runs           = 10,  n_final = 20,
  species        = FALSE, verbose = FALSE,
  skip_conv      = FALSE, max_attempts = 2, n_cores = 1L
)
cat("bace() finished in",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
cat("Converged   :", bace_result$converged,
    "  (attempts =", bace_result$n_attempts, ")\n")
cat("n_final     :", length(bace_result$imputed_datasets), "\n")

# -----------------------------------------------------------------------------
# STEP 4. Side-by-side: true vs. imputed values (first imputation)
# -----------------------------------------------------------------------------
cat("\n--- STEP 4: True vs. imputed (imputation #1) ---\n")
imp1 <- bace_result$imputed_datasets[[1]]
for (v in c("y", "x1", "x2", "x3", "x4")) {
  idx <- miss_mask[[v]]
  tv  <- complete_data[[v]][idx]
  iv  <- imp1[[v]][idx]
  cat("\nVariable:", v, " (", class(complete_data[[v]])[1], ",",
      sum(idx), "hidden cells)\n")
  print(head(data.frame(true = tv, imputed = iv), 8), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# STEP 5. Hand-compute metrics on imputation #1 (point-estimate family)
# -----------------------------------------------------------------------------
cat("\n--- STEP 5: Hand-computed point-estimate metrics (imp #1) ---\n\n")

# y (gaussian): NRMSE (marginal-sd), MAE, correlation.
# NOTE: NRMSE divides by sd of the COMPLETE data, not the hidden subset.
# Under trait_MAR the hidden subset is a biased slice of y, so dividing
# by its (narrower) sd would inflate NRMSE spuriously. Normalising by
# the marginal sd is the modern convention (Stekhoven & Buhlmann 2012;
# van Buuren 2018 FIMD §5.1).
tv <- as.numeric(complete_data$y[miss_mask$y])
iv <- as.numeric(imp1$y[miss_mask$y])
sd_y_full <- sd(as.numeric(complete_data$y))
cat("y (gaussian):\n")
cat("  NRMSE =", round(sqrt(mean((iv - tv)^2)) / sd_y_full, 4),
    "  (RMSE / sd(complete y); 0 = perfect, 1 = as bad as marginal mean)\n")
cat("  MAE   =", round(mean(abs(iv - tv)), 4), "  (raw scale)\n")
cat("  cor   =", round(cor(iv, tv), 4), "\n")

# x3 (poisson): NRMSE (marginal-sd), MAE, log_mae, correlation
tv <- as.numeric(complete_data$x3[miss_mask$x3])
iv <- as.numeric(imp1$x3[miss_mask$x3])
sd_x3_full <- sd(as.numeric(complete_data$x3))
cat("\nx3 (poisson):\n")
cat("  NRMSE   =", round(sqrt(mean((iv - tv)^2)) / sd_x3_full, 4), "\n")
cat("  MAE     =", round(mean(abs(iv - tv)), 4), "\n")
cat("  log_mae =",
    round(mean(abs(log1p(pmax(iv, 0)) - log1p(pmax(tv, 0)))), 4), "\n")
cat("  cor     =", round(cor(iv, tv), 4), "\n")

# x1 (binary): accuracy + confusion
tv <- as.character(complete_data$x1[miss_mask$x1])
iv <- as.character(imp1$x1[miss_mask$x1])
cat("\nx1 (binary):\n")
cat("  accuracy =", round(mean(iv == tv), 4), "\n")
cat("  confusion (rows=true, cols=imp):\n"); print(table(true = tv, imp = iv))

# x4 (threshold3): ordinal MAE
lvls <- levels(complete_data$x4)
tv_chr <- as.character(complete_data$x4[miss_mask$x4])
iv_chr <- as.character(imp1$x4[miss_mask$x4])
cat("\nx4 (threshold3) levels:", paste(lvls, collapse = " < "), "\n")
cat("  accuracy    =", round(mean(iv_chr == tv_chr), 4), "\n")
cat("  ordinal MAE =", round(mean(abs(match(iv_chr, lvls) -
                                        match(tv_chr, lvls))), 4), "\n")

# -----------------------------------------------------------------------------
# STEP 6. Hand-compute CALIBRATION metrics from the n_final ensemble
# -----------------------------------------------------------------------------
# Coverage (continuous/count): for each hidden cell, form the central
# 95% interval of the n_final imputed values. Count the proportion of
# cells where the TRUE value falls inside.
#
# Brier (categorical): estimate class probabilities from the frequency
# across imputations; score against the 0/1 truth indicator.
cat("\n--- STEP 6: Hand-computed calibration metrics (ensemble) ---\n\n")

n_imp <- length(bace_result$imputed_datasets)

# y coverage
idx <- miss_mask$y
tv  <- as.numeric(complete_data$y[idx])
iv_mat <- sapply(bace_result$imputed_datasets, function(d)
                  as.numeric(d$y[idx]))
lo <- apply(iv_mat, 1, quantile, 0.025)
hi <- apply(iv_mat, 1, quantile, 0.975)
cat("y coverage95 =", round(mean(tv >= lo & tv <= hi), 4),
    "  (well calibrated ~ 0.95)\n")

# x3 coverage
idx <- miss_mask$x3
tv  <- as.numeric(complete_data$x3[idx])
iv_mat <- sapply(bace_result$imputed_datasets, function(d)
                  as.numeric(d$x3[idx]))
lo <- apply(iv_mat, 1, quantile, 0.025)
hi <- apply(iv_mat, 1, quantile, 0.975)
cat("x3 coverage95 =", round(mean(tv >= lo & tv <= hi), 4), "\n")

# x1 Brier
idx <- miss_mask$x1
tv_chr <- as.character(complete_data$x1[idx])
im_mat <- sapply(bace_result$imputed_datasets, function(d)
                  as.character(d$x1[idx]))
classes <- levels(complete_data$x1)
prob_mat <- sapply(classes, function(cl) rowMeans(im_mat == cl))
y_indic  <- sapply(classes, function(cl) as.integer(tv_chr == cl))
brier <- mean(rowSums((prob_mat - y_indic)^2))
cat("x1 brier     =", round(brier, 4),
    "  (0 = perfect; 0.5 = uninformative for binary)\n")

# -----------------------------------------------------------------------------
# STEP 7. Cross-check against the pipeline evaluator
# -----------------------------------------------------------------------------
cat("\n--- STEP 7: evaluate_imputation_ensemble() output ---\n")
var_types <- c(y = "gaussian", x1 = "categorical", x2 = "categorical",
               x3 = "count",    x4 = "ordered")
pipeline_scores <- evaluate_imputation_ensemble(
  complete_data, bace_result$imputed_datasets, miss_mask, var_types)
print(pipeline_scores, row.names = FALSE)

# -----------------------------------------------------------------------------
# STEP 8. Diagnostic plots (three pages, one question each)
# -----------------------------------------------------------------------------
# The plots are designed to answer three questions in order:
#   PAGE 1. "Are the three missingness mechanisms actually different?"
#   PAGE 2. "How close are imputed values to truth, per variable?"
#   PAGE 3. "Is the imputation uncertainty honest (well calibrated)?"
plot_file <- file.path(DEMO_DIR, "demo_true_vs_imputed.pdf")
pdf(plot_file, width = 11, height = 9)

# =============================================================================
# PAGE 1. Mechanism fingerprint
# =============================================================================
# For each mechanism, plot the per-species missingness rate for y against
# the tree-tip order (i.e. species arranged in phylogenetic order). Under
# phylo_MAR neighbouring tips share similar rates (the LOESS smooth
# visibly wiggles); under trait_MAR / trait_MNAR the rate is driven by
# trait value, not tree position, so phylo-order rates look flat-noisy.
par(mfrow = c(3, 1), mar = c(3.5, 4, 3, 1), oma = c(0, 0, 2.5, 0))
tree_lo    <- sim_lo$tree
tip_order  <- tree_lo$tip.label

for (m in names(masks)) {
  rates <- tapply(masks[[m]]$miss_mask$y, cd_lo$Species, mean)[tip_order]
  rates[is.na(rates)] <- 0   # species absent from any row (rare)
  plot(seq_along(rates), rates, pch = 19,
       col = adjustcolor("steelblue", 0.6),
       xlab = "Species (ordered along tree tips)",
       ylab = "Prop. missing for y",
       main = paste0(m, "   (target rate = 35%)"),
       ylim = c(0, 1))
  # LOESS smooth - wiggly under phylo_MAR (clade-correlated), flat otherwise.
  lw <- lowess(seq_along(rates), rates, f = 0.25)
  lines(lw, col = "red", lwd = 2)
  abline(h = 0.35, lty = 2, col = "grey50")
  legend("topright",
         legend = c("per-species rate", "LOESS smooth", "target rate"),
         col = c("steelblue", "red", "grey50"),
         pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 1),
         bty = "n", cex = 0.85)
}
mtext(paste("Page 1: mechanism fingerprint - per-species missingness for y",
            "(low-signal sim so clade clustering is visible)"),
      outer = TRUE, cex = 1.05)

# =============================================================================
# PAGE 2. True vs. imputed, per variable
# =============================================================================
# Continuous / count: scatter of true vs. posterior-mean imputed, with
# the 95% posterior predictive interval drawn as a vertical bar.
# Categorical / ordered: confusion-matrix heatmap with row percentages.
# For every panel we print the headline metric(s) in the title so the
# plot can stand alone.

var_types_demo <- c(y = "gaussian", x3 = "count",
                    x1 = "binary", x2 = "multinomial3",
                    x4 = "threshold3")

# Precompute per-variable imputation matrices + metrics once, reuse below.
imp_mats <- lapply(names(var_types_demo), function(v) {
  idx <- miss_mask[[v]]
  sapply(bace_result$imputed_datasets,
         function(d) if (var_types_demo[v] %in% c("gaussian", "count"))
                       as.numeric(d[[v]][idx])
                     else as.character(d[[v]][idx]))
})
names(imp_mats) <- names(var_types_demo)

layout(matrix(c(1, 2, 0,
                3, 4, 5), nrow = 2, byrow = TRUE))
par(mar = c(4, 4, 3, 1), oma = c(0, 0, 2.5, 0))

# --- Continuous / count scatters ---------------------------------------------
for (v in c("y", "x3")) {
  idx <- miss_mask[[v]]
  tv  <- as.numeric(complete_data[[v]][idx])
  im  <- imp_mats[[v]]
  mid <- rowMeans(im)
  lo  <- apply(im, 1, quantile, 0.025)
  hi  <- apply(im, 1, quantile, 0.975)

  sd_full <- sd(as.numeric(complete_data[[v]]), na.rm = TRUE)
  nrmse   <- sqrt(mean((mid - tv)^2)) / sd_full
  cov95   <- mean(tv >= lo & tv <= hi)
  cor_v   <- cor(mid, tv)

  plot(tv, mid, pch = 19, col = adjustcolor("steelblue", 0.55),
       xlim = range(c(tv, lo, hi)), ylim = range(c(tv, lo, hi)),
       xlab = "True value", ylab = "Posterior mean (imputed)",
       main = sprintf("%s (%s)\nNRMSE=%.2f  r=%.2f  cov95=%.0f%%",
                      v, var_types_demo[v], nrmse, cor_v, cov95 * 100))
  segments(tv, lo, tv, hi, col = adjustcolor("steelblue", 0.25))
  abline(0, 1, col = "red", lty = 2)
}

# --- Categorical / ordered confusion-matrix heatmaps -------------------------
plot_confusion <- function(v, label) {
  idx <- miss_mask[[v]]
  tv  <- as.character(complete_data[[v]][idx])
  im  <- imp_mats[[v]]
  # Majority-vote imputation per cell (most common value across n_final).
  vote <- apply(im, 1, function(row) names(sort(table(row),
                                               decreasing = TRUE))[1])

  lvls <- if (is.factor(complete_data[[v]]))
    levels(complete_data[[v]]) else sort(unique(tv))
  tab      <- table(factor(tv,  levels = lvls),
                    factor(vote, levels = lvls))
  rowprop  <- sweep(tab, 1, pmax(rowSums(tab), 1), FUN = "/")

  acc <- mean(tv == vote)

  image(x = seq_along(lvls), y = seq_along(lvls), z = t(rowprop[nrow(rowprop):1, , drop = FALSE]),
        col = grDevices::colorRampPalette(c("white", "#1b7837"))(64),
        zlim = c(0, 1), axes = FALSE,
        xlab = "Imputed class", ylab = "True class",
        main = sprintf("%s (%s)\naccuracy = %.0f%%", v, label, acc * 100))
  axis(1, at = seq_along(lvls), labels = lvls)
  axis(2, at = seq_along(lvls), labels = rev(lvls), las = 1)
  box()
  # Annotate each cell with count and row-percentage.
  for (i in seq_along(lvls)) for (j in seq_along(lvls)) {
    count <- tab[i, j]; pct <- rowprop[i, j] * 100
    text(j, length(lvls) - i + 1,
         sprintf("%d\n(%.0f%%)", count, pct),
         col = if (pct > 50) "white" else "black", cex = 0.8)
  }
}

plot_confusion("x1", "binary")
plot_confusion("x2", "multinomial3")
plot_confusion("x4", "threshold3 ordered")

mtext(paste("Page 2: imputed vs true per variable.",
            "Top: posterior mean +/- 95% PI.",
            "Bottom: confusion matrix (row %)."),
      outer = TRUE, cex = 1.05)

# =============================================================================
# PAGE 3. Calibration scorecard
# =============================================================================
# One bar per variable for each of the key metrics, with target
# reference lines where they exist. This is the single "summary view"
# a reader would share to communicate overall imputation quality.
par(mfrow = c(2, 2), mar = c(4, 5, 3, 1), oma = c(0, 0, 2.5, 0))

scores <- pipeline_scores   # from STEP 7
bar_by_metric <- function(metric, main, ylab, ref = NA, ylim = NULL) {
  sub <- scores[scores$metric == metric, , drop = FALSE]
  if (!nrow(sub)) return(invisible())
  sub <- sub[order(sub$variable), ]
  bp <- barplot(sub$value, names.arg = sub$variable, col = "#6a5acd",
                border = NA, main = main, ylab = ylab,
                ylim = if (is.null(ylim)) c(0, max(sub$value) * 1.15) else ylim,
                las = 1)
  text(bp, sub$value, sprintf("%.2f", sub$value),
       pos = 3, cex = 0.85, xpd = NA)
  if (!is.na(ref)) abline(h = ref, lty = 2, col = "red")
}

bar_by_metric("accuracy",
              "% imputed cells matching truth (categorical/ordered)",
              "Accuracy", ref = 1, ylim = c(0, 1))
bar_by_metric("nrmse",
              "NRMSE (continuous/count) - lower is better",
              "NRMSE", ref = 1)
bar_by_metric("coverage95",
              "95% posterior PI coverage (continuous/count)",
              "Coverage", ref = 0.95, ylim = c(0, 1))
bar_by_metric("brier",
              "Brier score (categorical/ordered) - lower is better",
              "Brier")

mtext("Page 3: calibration & accuracy scorecard for this replicate",
      outer = TRUE, cex = 1.05)

dev.off()
cat("\nDiagnostic plots saved to:", plot_file, "\n")

saveRDS(list(sim = sim, masks = masks, bace_result = bace_result,
             pipeline_scores = pipeline_scores),
        file.path(DEMO_DIR, "demo_artefacts.rds"))
cat("Artefacts saved to:", file.path(DEMO_DIR, "demo_artefacts.rds"), "\n")
cat("\nDemo complete.\n")
