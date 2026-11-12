# =================================================================================================
# reproducibility_probe.R -- run on EVERY machine that will run the campaign, then diff the output.
# Takes seconds. Uses no MCMC.
#
#   cd dev/simulation_2026-11 && Rscript reproducibility_probe.R
#
# WHY. Seed 1 of the same scenario gave different answers on two machines (2026-11: ordinal
# accuracy 0.0625 vs 0.28). A campaign whose results depend on the machine cannot be reproduced by
# a reader, and -- worse -- we cannot tell a harmless numerical difference from a real bug.
#
# This splits the question in two:
#   SECTION A  the generated DATA. R's RNG + ape only. No MCMC.
#   SECTION B  the environment: R, platform, maths library, package versions.
#
# READING IT
#   A differs            -> the SIMULATOR is not reproducing. The two machines are not even
#                           analysing the same datasets. Check the R version first (sample()
#                           changed in R 3.6.0, and the campaign does not fix RNGkind), then ape.
#                           This must be fixed before the campaign runs.
#   A matches, results
#   still differ         -> divergence is downstream in MCMCglmm. Compare its version before
#                           blaming the numerical libraries.
# =================================================================================================

suppressPackageStartupMessages({library(ape); library(phytools)})
source("R/01_scenarios.R")
source("R/02_simulate.R")

fmt <- function(x, d = 10) formatC(x, format = "f", digits = d)

cat("=============== SECTION A: generated data (no MCMC) ===============\n")
sc <- make_scenarios()
for (id in c("n100_ps0.7_typical_mcar", "n100_ps0.3_typical_clade")) {
  s <- sc[sc$id == id, ]
  d <- simulate_traits(s, seed = 1)
  cat("\nscenario: ", id, "   seed 1\n", sep = "")
  cat("  tree total branch length : ", fmt(sum(d$tree$edge.length)), "\n", sep = "")
  cat("  first 3 tip labels       : ", paste(head(d$tree$tip.label, 3), collapse = " "), "\n", sep = "")
  cat("  x1[1:3]                  : ", paste(fmt(head(d$truth$x1, 3)), collapse = " "), "\n", sep = "")
  cat("  mean(x1), sd(x1)         : ", fmt(mean(d$truth$x1)), " ", fmt(sd(d$truth$x1)), "\n", sep = "")
  cat("  table(ord)               : ",
      paste(names(table(d$truth$ord)), table(d$truth$ord), sep = "=", collapse = " "), "\n", sep = "")
  cat("  table(bin)               : ",
      paste(names(table(d$truth$bin)), table(d$truth$bin), sep = "=", collapse = " "), "\n", sep = "")
  cat("  realised lambda_x1       : ", fmt(d$realised$lambda_x1), "\n", sep = "")
  cat("  realised lambda_bin      : ", fmt(d$realised$lambda_bin), "\n", sep = "")
  m <- apply_missingness(d$truth, d$tree, s$mechanism, s$miss_frac, seed = 1)
  obs <- if (is.list(m) && !is.data.frame(m)) m[[1]] else m
  cat("  which cells hidden in x1 : ", paste(head(which(is.na(obs$x1)), 8), collapse = " "), "\n", sep = "")
}

cat("\n  Two bare RNG checks, independent of every package:\n")
set.seed(1); cat("    runif(3)      : ", paste(fmt(runif(3)), collapse = " "), "\n", sep = "")
set.seed(1); cat("    sample(100, 5): ", paste(sample(100, 5), collapse = " "), "\n", sep = "")

cat("\n=============== SECTION B: environment ===============\n")
si <- sessionInfo()
cat("  R           : ", si$R.version$version.string, "\n", sep = "")
cat("  platform    : ", si$platform, "\n", sep = "")
cat("  BLAS        : ", basename(si$BLAS %||% "unknown"), "\n", sep = "")
cat("  LAPACK      : ", basename(si$LAPACK %||% "unknown"), "\n", sep = "")
cat("  RNGkind     : ", paste(RNGkind(), collapse = " / "), "\n", sep = "")
for (p in c("BACE", "MCMCglmm", "ape", "phytools", "Matrix", "nlme", "coda"))
  cat(sprintf("  %-12s: %s\n", p,
      tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")))
