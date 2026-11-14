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
      ## For a COUNT the same error on the log(1 + count) scale. At strong signal the counts
      ## span orders of magnitude and a handful of giants dominate the raw-scale error, so the
      ## raw number says how the tail came out and the log number how the typical species did.
      ## Both are reported; neither replaces the other. Continuous traits get NA here.
      zrmse_log <- if (is.integer(truth[[v]])) {
        sqrt(mean((log1p(best) - log1p(actual))^2)) / stats::sd(log1p(truth[[v]][!mask[, v]]))
      } else NA_real_
      rows[[v]] <- data.frame(trait = v, kind = "numeric", n_hidden = length(hidden),
                              zrmse = sqrt(mean((best - actual)^2)) / spread,
                              zrmse_log = zrmse_log,
                              accuracy = NA_real_)
    } else {
      ## STEP 2b. Category: the most frequently imputed class wins, ties share the credit.
      hit <- vapply(seq_along(hidden), function(i) {
        counts <- table(as.character(draws[i, ]))
        top    <- names(counts)[counts == max(counts)]     # may be more than one
        as.numeric(actual[i] %in% top) / length(top)
      }, numeric(1))
      rows[[v]] <- data.frame(trait = v, kind = "category", n_hidden = length(hidden),
                              zrmse = NA_real_, zrmse_log = NA_real_, accuracy = mean(hit))
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
#' several starting values before falling back.
#'
#' THE FALLBACK MUST GO TO THE NEAREST EDGE, NOT ALWAYS TO 1. An unconstrained lambda-hat below 0
#' says the data carry LESS phylogenetic structure than lambda = 0; the constrained estimate on
#' [0, 1] is then lambda = 0, because the profile likelihood is unimodal and still falling as it
#' leaves the interval. Until 2026-11-10 every out-of-range estimate was refitted with lambda
#' FIXED AT 1 -- the maximum-structure extreme, the opposite end of the dial from what the data
#' said. Since the free optimiser only ever strayed BELOW 0 in the runs we measured, that turned
#' "no phylogenetic signal" into "full Brownian motion", corrupting both the point estimate and
#' the interval, and it fired in roughly two of every three BACE replicates at phylo signal 0.3
#' (the BACE arm fits M imputations, so it gets M chances to trip).
#'
#' When no free start yields a finite lambda at all there is nothing to clamp towards, so the two
#' boundaries are compared by likelihood rather than one of them being assumed. Both fits are
#' REML with identical fixed effects and identical (zero) free correlation parameters, so their
#' log-likelihoods are comparable.
#'
#' The returned fit carries attributes `lambda_used` and `lambda_path` ("free", "clamped_lo",
#' "clamped_hi" or "profiled") so the caller can record how the slope was obtained.
#'
#' IMPORTANT: the argument MUST be called `tree`. with_imputations() looks for a formal with that
#' exact name. Up to BACE 1.0.0 a differently named argument meant the phylogeny was dropped in
#' silence and every model conditioned on the wrong tree (BACE issue #18); from 1.0.1 that is an
#' error instead, so this name is still required -- it just fails loudly now rather than quietly.
pgls_slope <- function(d, tree) {

  tag <- function(fit, lam, path) {
    if (is.null(fit)) return(NULL)
    attr(fit, "lambda_used") <- lam
    attr(fit, "lambda_path") <- path
    fit
  }

  fixed_fit <- function(lam) {
    tryCatch(nlme::gls(x2 ~ x1, data = d,
                       correlation = ape::corPagel(lam, phy = tree, form = ~Species,
                                                   fixed = TRUE)),
             error = function(e) NULL, warning = function(w) NULL)
  }

  ## 1. Try the free (unconstrained) optimiser from several starts.
  stray <- NA_real_          # last finite lambda that landed outside [0, 1]
  for (start in c(0.5, 0.1, 0.9)) {
    fit <- tryCatch(
      nlme::gls(x2 ~ x1, data = d,
                correlation = ape::corPagel(start, phy = tree, form = ~Species)),
      error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) next
    lambda <- unname(stats::coef(fit$modelStruct$corStruct, unconstrained = FALSE))
    if (is.finite(lambda) && lambda >= 0 && lambda <= 1) return(tag(fit, lambda, "free"))
    if (is.finite(lambda)) stray <- lambda
  }

  ## 2. It strayed out of [0, 1]: refit at the NEAREST edge, not always at 1.
  if (is.finite(stray)) {
    lam  <- if (stray < 0) 0 else 1
    path <- if (stray < 0) "clamped_lo" else "clamped_hi"
    f <- fixed_fit(lam)
    if (!is.null(f)) return(tag(f, lam, path))
  }

  ## 3. No finite free estimate at all: let the likelihood choose the edge.
  f0 <- fixed_fit(0)
  f1 <- fixed_fit(1)
  if (is.null(f0) && is.null(f1)) return(NULL)
  if (is.null(f1)) return(tag(f0, 0, "profiled"))
  if (is.null(f0)) return(tag(f1, 1, "profiled"))
  ll0 <- tryCatch(as.numeric(stats::logLik(f0)), error = function(e) -Inf)
  ll1 <- tryCatch(as.numeric(stats::logLik(f1)), error = function(e) -Inf)
  if (ll0 >= ll1) tag(f0, 0, "profiled") else tag(f1, 1, "profiled")
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
