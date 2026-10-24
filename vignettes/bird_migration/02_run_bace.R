## =============================================================================
## Bird migration example, step 2: run BACE
## =============================================================================
##
## Runs one analysis (full, cv1-cv5 or mixed; see 00_setup.R) in three stages.
## Final imputations are independent given the chained phase, so each one can
## run in its own process. Run from the BACE repository root:
##
##   Rscript vignettes/bird_migration/02_run_bace.R full chain       # chained phase
##   Rscript vignettes/bird_migration/02_run_bace.R full final 1     # final imputation 1
##   ...                                                             # (in parallel)
##   Rscript vignettes/bird_migration/02_run_bace.R full final 50    # final imputation 50
##
## On a shared Linux server, pin the BLAS/OpenMP threads before R starts:
##   export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
## Each final imputation of the full data takes about 50 min and 2 GB of RAM.
## Results go to vignettes/bird_migration/results/tmp/ (not tracked); step 3
## (03_summarise.R) combines them.
## =============================================================================

args     <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]
mode     <- args[2]

suppressMessages({
  pkgload::load_all(".", quiet = TRUE)   # the BACE source in this repository
  library(ape)
})
source("vignettes/bird_migration/00_setup.R")
stopifnot(analysis %in% analyses, mode %in% c("chain", "final"))

tmp_dir <- "vignettes/bird_migration/results/tmp"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
a <- make_analysis_data(analysis)
analysis_id <- match(analysis, analyses)

if (mode == "chain") {
  runs  <- if (analysis == "mixed") chain_settings$runs_mixed else chain_settings$runs_ordinal
  chain <- bace_chain_step(a, runs = runs, seed = 1000 * analysis_id)
  saveRDS(chain, file.path(tmp_dir, sprintf("%s_chain.rds", analysis)))
  cat(sprintf("%s: chained phase done in %.1f min\n", analysis, chain$minutes))
}

if (mode == "final") {
  i <- as.integer(args[3])
  stopifnot(i >= 1, i <= n_final[[analysis]])
  chain <- readRDS(file.path(tmp_dir, sprintf("%s_chain.rds", analysis)))
  final <- bace_final_one(a, chain, seed = 1000 * analysis_id + 100 + i)
  saveRDS(final, file.path(tmp_dir, sprintf("%s_final%02d.rds", analysis, i)))
  cat(sprintf("%s: final imputation %d done in %.1f min\n", analysis, i, final$minutes))
}
