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


`%||%` <- function(x, y) if (is.null(x)) y else x

#' Run an expression and count BACE's poisson rate-cap warnings, per variable.
#'
#' BACE caps an imputed poisson rate at 3 * max(observed) + 5 to stop runaway draws, and warns
#' each time it does. The cap biases counts DOWN exactly when the hidden counts are the large
#' ones -- driven and clade missingness hide those, and at strong signal the counts span orders
#' of magnitude -- so how often it bit is a result, not noise. It only lived in log files until
#' 2026-11-14; now it travels with every dataset. Every other warning passes through untouched.
#'
#' @return list(value = what expr returned, rate_cap = named list, one entry per capped variable:
#'   n_warnings, n_rates (rates truncated, summed over warnings), cap (the last cap value))
with_rate_cap_count <- function(expr) {
  events <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    m <- regmatches(msg, regexec(
      "rate guard bound on ([0-9]+) predicted rate\\(s\\) for '([^']+)' \\(cap = ([^ ]+) =", msg))[[1]]
    if (length(m) == 4L) {
      v <- m[3]
      prev <- events[[v]] %||% list(n_warnings = 0L, n_rates = 0L, cap = NA_real_)
      events[[v]] <<- list(n_warnings = prev$n_warnings + 1L,
                           n_rates    = prev$n_rates + as.integer(m[2]),
                           cap        = as.numeric(m[4]))
      invokeRestart("muffleWarning")
    }
  })
  list(value = value, rate_cap = events)
}

#' Run BACE on one dataset.
#'
#' RESUMING. `checkpoint_dir` makes bace() save the convergence chain when it finishes and each
#' final run as it completes, so an interrupted dataset picks up instead of starting over.
#' Measured on a small fit: the resumed run took 5% of the time and gave an IDENTICAL result.
#' It is worth having on a three-week campaign, where wall-time limits and node failures are
#' ordinary -- at n = 500 a lost dataset is 2.2 hours of work.
#'
#' Two things to be clear about:
#'   - It does NOT help with the two known BACE crashes. Those are deterministic: the same
#'     scenario and seed fails the same way every time, so a resume just re-hits the bug.
#'   - The checkpoint is BIG -- 51.6 MB even for a tiny fit, so much more at production settings.
#'     54,000 of them cannot accumulate, so this function DELETES the directory on success. Only
#'     interrupted datasets leave one behind, which is exactly what should be resumed.
#'
#' `seed` is passed whether or not checkpointing is on: resume requires it, and it makes every
#' dataset exactly reproducible, which is worth having regardless.
#'
#' @param dataset   a simulate_dataset() result
#' @param n_final   number of completed datasets (BACE's n_final; the M of multiple imputation)
#' @param nitt,burnin,thin,runs  MCMC settings, passed straight through to bace()
#' @param seed      makes the run reproducible, and is required for resuming
#' @param checkpoint_dir  where to save progress, or NULL for none
#' @return list(imputed = list of completed data frames, intervals, converged = BACE's verdict,
#'   conv_detail, wall_s, rate_cap = how often BACE's poisson rate cap bit, per variable)
fit_bace <- function(dataset, n_final = 20L, nitt = 50000L, burnin = 10000L,
                     thin = 25L, runs = 7L, seed = 1L, checkpoint_dir = NULL,
                     final_nitt = NULL, final_burnin = NULL, final_thin = NULL) {

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
  run <- with_rate_cap_count(BACE::bace(fixformula    = formulas,
                    ran_phylo_form = "~1|Species",
                    phylo          = dataset$tree,
                    data           = d,
                    nitt = nitt, burnin = burnin, thin = thin, runs = runs,
                    n_final = n_final, n_cores = 1L,
                    verbose = FALSE, skip_conv = TRUE,
                    seed = seed, checkpoint_dir = checkpoint_dir,
                    ## The two phases do different jobs, so they can have different chains.
                    ## The chained phase only needs the imputed values to SETTLE; the final runs
                    ## supply the draws that are pooled for inference. NULL reuses the chained
                    ## settings, which is what BACE did before 1.0.0.
                    final_nitt = final_nitt, final_burnin = final_burnin,
                    final_thin = final_thin))
  fit  <- run$value
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

  ## The fit succeeded, so the checkpoint has done its job and must go -- at a few hundred MB
  ## each, keeping them would fill the disk long before the campaign finished.
  if (!is.null(checkpoint_dir) && dir.exists(checkpoint_dir))
    unlink(checkpoint_dir, recursive = TRUE)

  ## Keep WHY as well as whether. "converged = FALSE" on its own is unactionable: it could mean
  ## the chained loop genuinely has not settled, or that one variable of eight is unsettled, or
  ## that the assessment could not run at all. The per-variable verdicts and the attempt count
  ## are small, so they travel with the result.
  conv_detail <- tryCatch({
    cv <- fit$convergence
    ## BACE 1.0.1's stability check names what it found per variable: "ok", or one of the
    ## pathologies "drift" (curable with more runs), "oscillation", "stuck", "unassessable".
    ## notes are advisory only (e.g. mild alternation) and do not affect the verdict.
    verdicts <- cv$verdicts %||% NULL
    v <- cv$method_results$summary$var_converged
    list(n_attempts = fit$n_attempts %||% NA_integer_,
         n_vars = length(v), n_vars_converged = sum(v %in% TRUE, na.rm = TRUE),
         vars = v,
         verdicts = verdicts,
         notes = Filter(nzchar, cv$notes %||% character(0)))
  }, error = function(e) NULL)

  list(imputed = imputed, intervals = intervals,
       converged = isTRUE(fit$converged), conv_detail = conv_detail, wall_s = wall,
       rate_cap = run$rate_cap)
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
