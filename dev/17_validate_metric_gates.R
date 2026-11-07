# dev/17_validate_metric_gates.R
#
# Metric-level sanity gates for bace_validate() at a realistic scale
# (the unit tests check the ARITHMETIC on hand-made inputs; this checks the
# MAGNITUDES on real fits, where the borrowed pass bounds are stable:
# >= 60 species and prop >= 0.3 -- see tests/testthat/test-imputation_accuracy.R
# for the provenance of the bounds).
#
# Gates (averaged over REPS replicates):
#   G1 gaussian, strong signal: zRMSE < 0.8 * mean-filling baseline,
#      correlation with truth r > 0.30, coverage within 0.90-1.00.
#   G2 categorical (3 classes), strong signal: accuracy >= its
#      commonest-class baseline (imputation must not be WORSE than guessing).
#
# Run:  Rscript dev/17_validate_metric_gates.R [REPS]   (default 5)

suppressMessages({ devtools::load_all(quiet = TRUE); library(ape) })
reps <- {a <- commandArgs(trailingOnly = TRUE)
         if (length(a) >= 1) as.integer(a[1]) else 5L}

run_one <- function(seed, response_type) {
  set.seed(seed)
  sim <- sim_bace(response_type = response_type,
                  predictor_types = c("gaussian", "gaussian"),
                  var_names = c("y", "x1", "x2"),
                  phylo_signal = c(0.9, 0.9, 0.9),
                  n_cases = 60, n_species = 60,
                  missingness = c(0, 0, 0))
  v <- suppressWarnings(bace_validate(
    "y ~ x1 + x2", ran_phylo_form = "~ 1 | Species",
    phylo = sim$tree, data = sim$complete_data,
    prop = 0.3, runs = 3, nitt = 2000, burnin = 400, thin = 4,
    n_final = 20, skip_conv = TRUE, verbose = FALSE, seed = seed))
  v$table[v$table$variable == "y", ]
}

cat("gaussian replicates:\n")
g <- do.call(rbind, lapply(seq_len(reps), run_one, response_type = "gaussian"))
print(g[, c("zrmse", "zrmse_baseline", "r", "coverage")], digits = 2)
g1a <- mean(g$zrmse) < 0.8 * mean(g$zrmse_baseline)
g1b <- mean(g$r) > 0.30
g1c <- mean(g$coverage) >= 0.90 && mean(g$coverage) <= 1.00
cat(sprintf("G1 gaussian: zRMSE %.2f vs baseline %.2f [%s], r %.2f [%s], coverage %.2f [%s]\n",
            mean(g$zrmse), mean(g$zrmse_baseline),
            if (g1a) "PASS" else "FAIL",
            mean(g$r), if (g1b) "PASS" else "FAIL",
            mean(g$coverage), if (g1c) "PASS" else "FAIL"))

cat("\ncategorical replicates:\n")
k <- do.call(rbind, lapply(seq_len(reps) + 100L, run_one,
                           response_type = "multinomial3"))
print(k[, c("accuracy", "accuracy_baseline", "balanced_accuracy", "brier")],
      digits = 2)
g2 <- mean(k$accuracy) >= mean(k$accuracy_baseline)
cat(sprintf("G2 categorical: accuracy %.2f vs baseline %.2f [%s]; Brier %.2f\n",
            mean(k$accuracy), mean(k$accuracy_baseline),
            if (g2) "PASS" else "FAIL", mean(k$brier, na.rm = TRUE)))

cat(sprintf("\nVALIDATE METRIC GATES: %s\n",
            if (g1a && g1b && g1c && g2) "ALL PASS" else "FAIL"))
