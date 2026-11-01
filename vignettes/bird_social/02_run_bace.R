## =============================================================================
## Bird social behaviour example, step 2: run BACE
## =============================================================================
##
## Runs one analysis (full, validation or mixed; see 00_setup.R) in two stages. Final
## imputations are independent given the chained phase, so each one can run in
## its own process. Run from the BACE repository root:
##
##   Rscript vignettes/bird_social/02_run_bace.R full chain       # chained phase
##   Rscript vignettes/bird_social/02_run_bace.R full final 1     # final imputation 1
##   ...                                                          # (in parallel)
##   Rscript vignettes/bird_social/02_run_bace.R full final 50    # final imputation 50
##
## On a shared Linux server, pin the BLAS/OpenMP threads before R starts:
##   export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
## The chained phase took 22 min (full; one model), 63 min (validation; a single
## formula string also fits models for the complete predictors) or 2.5 h
## (mixed); each final imputation about 55 min (full, validation) or 75 min (mixed).
## Results go to vignettes/bird_social/results/tmp/ (not tracked); step 3
## (03_summarise.R) combines them.
## =============================================================================

args     <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]
mode     <- args[2]

suppressMessages({
  pkgload::load_all(".", quiet = TRUE)   # the BACE source in this repository
                                         # (results in results/ were produced at git dc79be1)
  library(ape)
})
source("vignettes/bird_social/00_setup.R")
stopifnot(analysis %in% analyses, mode %in% c("chain", "final"))

tmp_dir <- "vignettes/bird_social/results/tmp"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
a  <- make_analysis_data(analysis)
cs <- chain_settings[[analysis]]
fs <- final_settings[[analysis]]

## -----------------------------------------------------------------------------
## Chained phase
## -----------------------------------------------------------------------------
if (mode == "chain") {
  set.seed(chain_seed[[analysis]])
  t0  <- Sys.time()
  fit <- bace_imp(fixformula = a$formulas, ran_phylo_form = "~ 1 | species",
                  phylo = a$tree, data = a$data, runs = cs$runs,
                  nitt = cs$nitt, burnin = cs$burnin, thin = cs$thin,
                  species = FALSE, verbose = TRUE)
  conv <- assess_convergence(fit, method = "summary")
  # Drop the heavy fitted models of the last chained iteration
  fit$models_last_run <- NULL
  fit$pred_list_last_run <- NULL
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(list(fit = fit, convergence = conv, minutes = minutes),
          file.path(tmp_dir, sprintf("%s_chain.rds", analysis)))
  print(conv)
  cat(sprintf("%s: chained phase done in %.1f min\n", analysis, minutes))
}

## -----------------------------------------------------------------------------
## One final (posterior-predictive) imputation, reduced to what we keep
## -----------------------------------------------------------------------------
if (mode == "final") {
  i <- as.integer(args[3])
  stopifnot(i >= 1, i <= n_final[[analysis]])
  chain <- readRDS(file.path(tmp_dir, sprintf("%s_chain.rds", analysis)))
  set.seed(final_seed[[analysis]] + i)
  t0  <- Sys.time()
  fin <- bace_final_imp(chain$fit, fixformula = a$formulas, ran_phylo_form = "~ 1 | species",
                        phylo = a$tree, nitt = fs$nitt, burnin = fs$burnin, thin = fs$thin,
                        n_final = 1, species = FALSE, verbose = FALSE, n_cores = 1,
                        seed = final_seed[[analysis]] + i)   # sweeps: package default
  completed  <- fin$all_datasets[[1]]
  fit_social <- fin$all_models[[1]][["social_max"]]
  out <- list(
    imputed_social = completed$social_max[a$social_rows],
    # probability of each category for the imputed rows, averaged over the
    # posterior draws of this imputation's model
    probs = as.matrix(fin$all_prob_preds[[1]][["social_max"]]),
    sol = fit_social$Sol[, seq_len(fit_social$Fixed$nfl), drop = FALSE],
    vcv = fit_social$VCV, cp = fit_social$CP,
    diagnostics = fin$diagnostics)
  if (analysis == "mixed") {
    fit_hwi <- fin$all_models[[1]][["hwi"]]
    out$imputed_hwi <- completed$hwi[a$hwi_rows]
    out$hwi_sol <- fit_hwi$Sol[, seq_len(fit_hwi$Fixed$nfl), drop = FALSE]
    out$hwi_vcv <- fit_hwi$VCV
    # The (scaled) hand-wing index that the social model of this imputation
    # conditioned on, at the rows where it was imputed
    out$hwi_used_by_social <- as.matrix(fit_social$X)[a$hwi_rows, "hwi"]
    out$sweeps <- fin$sweeps
  }
  out$minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(out, file.path(tmp_dir, sprintf("%s_final%02d.rds", analysis, i)))
  cat(sprintf("%s: final imputation %d done in %.1f min\n", analysis, i, out$minutes))
}
