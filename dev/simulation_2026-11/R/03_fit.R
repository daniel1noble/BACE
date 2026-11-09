# =================================================================================================
# STEP 3 of 4 -- run the methods on one dataset.
#
# Three arms. Each one answers a different question:
#
#   bace           impute the gaps with BACE::bace()        -- the method the paper is about
#   complete_case  drop every species with any gap          -- what people do now
#   oracle         the truth, before we made any holes      -- the best anyone could possibly do
#
# BACE is called through bace(). We do not reach into bace_imp(), bace_final_imp() or any internal:
# the point is to test the function users actually call.
#
# NOTE ON SIZE: a fitted bace object is about 35 MB. A campaign of 54,000 of them cannot keep any
# of them, so this file returns the imputed datasets and nothing else, and the fit goes out of
# scope. Never saveRDS() a fit.
# =================================================================================================


#' Run BACE on one dataset.
#'
#' @param dataset   a simulate_dataset() result
#' @param n_final   number of completed datasets (BACE's n_final; the M of multiple imputation)
#' @param nitt,burnin,thin,runs  MCMC settings, passed straight through to bace()
#' @return list(imputed = list of completed data frames, converged = BACE's verdict, wall_s)
fit_bace <- function(dataset, n_final = 20L, nitt = 50000L, burnin = 10000L,
                     thin = 25L, runs = 5L) {

  ## STEP 1. One imputation model per trait: that trait on the left, every OTHER trait on the
  ##         right, plus `driver` and (where present) `aux`. The two helpers appear only on the
  ##         right -- they are fully observed, so there is nothing in them to impute, and giving
  ##         them a model of their own would be wasted compute.
  d        <- dataset$data
  helpers  <- intersect(c("driver", "aux"), names(d))
  traits   <- setdiff(names(d), c("Species", helpers))
  formulas <- lapply(traits, function(v)
    paste(v, "~", paste(c(setdiff(traits, v), helpers), collapse = " + ")))

  ## STEP 2. Impute. skip_conv = TRUE because a campaign cannot afford unbounded retries; the
  ##         convergence verdict is recorded per dataset instead and reported as a rate.
  t0  <- Sys.time()
  fit <- BACE::bace(fixformula    = formulas,
                    ran_phylo_form = "~1|Species",
                    phylo          = dataset$tree,
                    data           = d,
                    nitt = nitt, burnin = burnin, thin = thin, runs = runs,
                    n_final = n_final, n_cores = 1L,
                    verbose = FALSE, skip_conv = TRUE)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ## STEP 3. Take the completed datasets and let the fit go. Two fixes BACE's output needs:
  ##         the Species column comes back LAST (so never index by position), and count columns
  ##         come back as doubles, which would make a second pass treat them as gaussian.
  imputed <- BACE::get_imputed_data(fit, format = "list")
  imputed <- lapply(imputed, function(x) {
    x <- x[, names(d), drop = FALSE]
    for (v in names(d)) if (is.integer(d[[v]])) x[[v]] <- as.integer(round(x[[v]]))
    x
  })

  ## The per-cell 95% intervals (Rubin's t across the M imputations) are pulled out HERE, while
  ## the fit is still in scope, because step 4 must never be handed a 35 MB object.
  intervals <- tryCatch(BACE::get_imputation_intervals(fit), error = function(e) NULL)

  list(imputed = imputed, intervals = intervals,
       converged = isTRUE(fit$converged), wall_s = wall)
}


#' The complete-case arm: drop every species with a gap.
#'
#' It is allowed to decline. Under MCAR with six traits at 30% each, only about 0.7^6 = 12% of
#' species survive, and sometimes too few to fit anything. That is a result about complete-case
#' analysis, not a failure of the campaign, so it returns a status and a survivor count rather
#' than a number that means nothing.
#'
#' @param dataset     a simulate_dataset() result
#' @param min_species refuse below this many survivors
#' @return list(data = the surviving rows or NULL, n_complete, frac_kept, status)
fit_complete_case <- function(dataset, min_species = 10L) {
  keep <- stats::complete.cases(dataset$data)
  n_cc <- sum(keep)
  out  <- list(n_complete = n_cc,
               frac_kept  = n_cc / nrow(dataset$data),
               status     = if (n_cc < min_species) "too_few" else "scored")
  out$data <- if (out$status == "scored") dataset$data[keep, , drop = FALSE] else NULL
  out
}
