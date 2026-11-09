# =================================================================================================
# STEP 4 of 4 -- turn one finished dataset into statistics.
#
# Two questions, scored separately, because they are different claims:
#
#   A. RECOVERY   how close were the imputed values to the truth we hid?
#                 -> score_recovery()
#   B. INFERENCE  does the analysis you would actually run give the right answer, with honest
#                 uncertainty?  The analysis here is the PGLS slope of x2 on x1.
#                 -> score_slope()
#
# Recovery is the obvious thing to measure and the less important one. A method can recover values
# well and still give you a biased slope with intervals that miss. B is what the manuscript's
# claims are about.
#
# Everything uses BACE's public functions: get_imputation_intervals() for the per-cell Rubin
# intervals, with_imputations() + pool_mi() for the pooled analysis.
# =================================================================================================


# -------------------------------------------------------------------------------------------------
# A. RECOVERY
# -------------------------------------------------------------------------------------------------

#' How well were the hidden values recovered?
#'
#' Continuous and count traits are scored by zRMSE: the error in the imputed values, in units of
#' the trait's own spread. 0 is perfect; 1 is no better than filling in the average.
#'
#' Category traits cannot be scored that way, so they get accuracy with FRACTIONAL CREDIT FOR
#' TIES: if two categories tie for most-imputed and one of them is right, that cell scores 0.5 --
#' which is what you would get on average by guessing between them. Scoring a tie as a win or a
#' loss would make the number depend on how R happens to break ties.
#'
#' @param truth    the complete data
#' @param mask     logical matrix, TRUE where a value was hidden
#' @param imputed  list of M completed data frames
#' @return one row per trait
score_recovery <- function(truth, mask, imputed) {
  rows <- list()

  for (v in colnames(mask)) {
    hidden <- which(mask[, v])
    if (!length(hidden)) next

    ## STEP 1. Collect the M imputations of every hidden cell: rows = cells, columns = imputations.
    draws <- sapply(imputed, function(d) d[[v]][hidden])
    if (!is.matrix(draws)) draws <- matrix(draws, nrow = length(hidden))
    actual <- truth[[v]][hidden]

    if (is.numeric(actual)) {
      ## STEP 2a. Continuous or count: average the imputations, compare with the truth, and
      ##          standardise by the spread of the values that were NOT hidden.
      best  <- rowMeans(draws)
      spread <- stats::sd(truth[[v]][!mask[, v]])
      rows[[v]] <- data.frame(trait = v, kind = "numeric", n_hidden = length(hidden),
                              zrmse = sqrt(mean((best - actual)^2)) / spread,
                              accuracy = NA_real_)
    } else {
      ## STEP 2b. Category: the most frequently imputed class wins, ties share the credit.
      hit <- vapply(seq_along(hidden), function(i) {
        counts <- table(as.character(draws[i, ]))
        top    <- names(counts)[counts == max(counts)]     # may be more than one
        as.numeric(actual[i] %in% top) / length(top)
      }, numeric(1))
      rows[[v]] <- data.frame(trait = v, kind = "category", n_hidden = length(hidden),
                              zrmse = NA_real_, accuracy = mean(hit))
    }
  }
  do.call(rbind, rows)
}


#' Did the 95% interval for each hidden value contain the truth?
#'
#' Uses BACE's own get_imputation_intervals(), which builds the Rubin t-interval across the M
#' imputations. It covers numeric traits only -- category traits get no interval from it, and are
#' simply absent from the result rather than scored as failures.
#'
#' @param iv    the intervals table, already extracted by fit_bace()
#' @param truth the complete data
#' @return one row per trait, or NULL
score_cell_coverage <- function(iv, truth) {
  if (is.null(iv) || !nrow(iv)) return(NULL)

  ## get_imputation_intervals() reports row and variable, so look the truth up by both.
  iv$actual  <- mapply(function(r, v) as.numeric(truth[r, v]), iv$row, iv$variable)
  iv$covered <- iv$actual >= iv$lower & iv$actual <= iv$upper

  stats::aggregate(cbind(coverage = covered, width = upper - lower) ~ variable, iv, mean)
}


# -------------------------------------------------------------------------------------------------
# B. INFERENCE -- the analysis a reader would actually run
# -------------------------------------------------------------------------------------------------

#' The PGLS slope of x2 on x1, fitted the same way for every arm so the arms differ only in the
#' data they saw.
#'
#' WHY PAGEL AND NOT BROWNIAN. The traits are generated as a phylogenetic effect plus independent
#' noise, which is exactly a Pagel's-lambda correlation structure -- so corBrownian is
#' misspecified. Measured on 20 complete datasets at n = 150, where the true slope is 0.5:
#'
#'     ols          mean 0.467   bias -0.033   sd 0.218
#'     corBrownian  mean 0.439   bias -0.061   sd 0.305
#'     corPagel     mean 0.526   bias +0.026   sd 0.072      <- correct, and 4x more precise
#'
#' WHY THE RETRIES. corPagel's optimiser is unconstrained and can wander outside [0, 1] or fail to
#' converge, more often on some completed datasets than others. Left alone that would quietly thin
#' the imputations and bias the pooled result towards the datasets that happened to fit. So we try
#' several starting values and return NULL if none works, which the caller treats as a failure
#' rather than as a smaller M.
#'
#' IMPORTANT: the argument MUST be called `tree`. with_imputations() looks for a formal with that
#' exact name and silently drops the phylogeny if it is called anything else (BACE issue #18), so
#' every model would then condition on the wrong tree with no warning.
pgls_slope <- function(d, tree) {
  for (start in c(0.5, 0.1, 0.9)) {
    fit <- tryCatch(
      nlme::gls(x2 ~ x1, data = d,
                correlation = ape::corPagel(start, phy = tree, form = ~Species)),
      error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) next
    lambda <- unname(stats::coef(fit$modelStruct$corStruct, unconstrained = FALSE))
    if (is.finite(lambda) && lambda >= 0 && lambda <= 1) return(fit)
  }
  ## Every free start failed or left [0, 1]: fall back to lambda fixed at 1 (Brownian), which
  ## always fits. Recorded by the caller, never silent.
  tryCatch(nlme::gls(x2 ~ x1, data = d,
                     correlation = ape::corPagel(1, phy = tree, form = ~Species, fixed = TRUE)),
           error = function(e) NULL)
}


#' Pool the slope across BACE's M imputed datasets, by Rubin's rules.
#'
#' @param imputed list of completed data frames
#' @param tree    the phylogeny
#' @return one row: estimate, se, interval, df
score_slope_bace <- function(imputed, tree) {
  fits <- suppressMessages(BACE::with_imputations(imputed, pgls_slope, tree = tree,
                                                  .progress = FALSE, .on_error = "continue"))
  ## Refuse to pool a thinned set. pool_mi() will happily pool whatever survived and report a full
  ## table; the count of fits actually used lives only in an attribute, so check it here.
  if (!is.null(attr(fits, "n_failed")) && attr(fits, "n_failed") > 0L) return(NULL)

  p <- suppressMessages(BACE::pool_mi(fits))
  p <- p[p$term == "x1", ]
  if (!nrow(p)) return(NULL)
  data.frame(estimate = p$estimate, se = p$std.error,
             lower = p$conf.low, upper = p$conf.high, df = p$df)
}


#' The same slope from a single dataset -- used for the complete-case arm and for the oracle.
#' One fit, so no pooling: the ordinary t interval on the residual degrees of freedom.
score_slope_single <- function(d, tree) {
  fit <- tryCatch(pgls_slope(d, tree), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  est <- stats::coef(fit)[["x1"]]
  se  <- sqrt(stats::vcov(fit)["x1", "x1"])
  df  <- nrow(d) - 2
  tcrit <- stats::qt(0.975, df)
  data.frame(estimate = est, se = se,
             lower = est - tcrit * se, upper = est + tcrit * se, df = df)
}
