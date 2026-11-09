# archive-2026-09/rubin_bace_chain.R
#
# RETIRED 2026-11-09. This was a WORKAROUND for a bug that no longer exists.
#
# In 2026-09 the shipped bace_final_imp() started all n_final runs from the same converged
# dataset, so the M imputed datasets shared the same point fills: between-imputation variance came
# out near zero and the intervals were too narrow. This study found that. mi_bace_chain() was the
# proposed fix, written outside the package so it could be tested without editing BACE -- run m
# starts from run m - 1, making a proper Gibbs chain.
#
# BACE fixed it internally on 2026-10-26 with `sweeps` (default 3), so bace() now does properly
# what this worked around. It also reached into BACE:::bace_final_imp(), a non-exported function,
# which is the kind of call that breaks silently on a refactor. The campaign calls bace() and
# nothing else.
#
# Requires fit_bace_mi() and mi_bace_shipped() from rubin_bace.R.

#' Chained BACE (the arm Meng's review B2 proposed; Shinichi approved it 2026-09-24). BACE's own final
#' step starts all n_final runs from the same converged dataset, so predictors later in the formula order
#' are identical across the M imputations. Here final run m starts from final dataset m - 1 instead
#' (run 1 from BACE's converged dataset), using BACE's own bace_final_imp() with n_final = 1 and the same
#' MCMC settings; BACE's code is called, never edited. It reuses the shipped fit's convergence phase
#' (fb$outb$initial_results), so the shipped and chained arms are paired on one fit. Cost: M extra
#' sweeps, run sequentially.
#'
#' @param fb a fit_bace_mi() result
#' @param df_miss the cell's missing-data data.frame
#' @param M number of chained datasets
#' @return list(datasets = M data.frames aligned like mi_bace_shipped(), models = per-run model lists,
#'   wall_s = chain wall time)
mi_bace_chain <- function(fb, df_miss, M = 20L) {
  obj <- fb$outb$initial_results
  k <- length(obj$data)
  raw <- vector("list", M); models <- vector("list", M)
  t0 <- Sys.time()
  for (m in seq_len(M)) {
    r <- BACE:::bace_final_imp(obj, fixformula = fb$fixformula, ran_phylo_form = "~1|Species",
                               phylo = fb$tree_b, nitt = fb$mcmc[["nitt"]], thin = fb$mcmc[["thin"]],
                               burnin = fb$mcmc[["burnin"]], n_final = 1L, species = FALSE,
                               verbose = FALSE, n_cores = 1L, ovr_categorical = fb$ovr,
                               sweeps = 1L)  # pin one sweep per link: the 2026-09 chain algorithm
                                             # (package default is now 3; see file header).
                                             # ovr follows the parent fit so chain links use the
                                             # same categorical construction as the shipped arm.
    raw[[m]] <- r$all_datasets[[1]]; models[[m]] <- r$all_models[[1]]
    obj$data[[k]] <- raw[[m]]   # the next run starts from this draw
  }
  wall_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(datasets = mi_bace_shipped(list(imputed_datasets = raw), df_miss), models = models, wall_s = wall_s)
}

