#' Pool downstream model fits across imputations with Rubin's rules
#'
#' Combine coefficient estimates from \code{M} model fits -- one per imputed
#' dataset -- into a single table using Rubin's rules (Rubin 1987). The pooled
#' standard errors account for \emph{both} within-imputation sampling variance
#' and between-imputation variance, so downstream inference propagates the
#' uncertainty introduced by imputation.
#'
#' @details
#' Let \eqn{\hat\theta_i} and \eqn{U_i} be the coefficient vector and its
#' variance for fit \eqn{i}. With \eqn{\bar\theta = M^{-1}\sum_i\hat\theta_i},
#' within-variance \eqn{W = M^{-1}\sum_i U_i}, and between-variance
#' \eqn{B = (M-1)^{-1}\sum_i(\hat\theta_i-\bar\theta)^2}, the total variance is
#' \eqn{T = W + (1 + 1/M)B} and the pooled SE is \eqn{\sqrt T}. The relative
#' increase in variance is \eqn{r = (1+1/M)B/W}, and the fraction of missing
#' information is \eqn{\mathrm{fmi} = (r + 2/(\nu+3))/(r+1)}. Degrees of freedom
#' use the classical Rubin (1987) formula \eqn{\nu_{old}=(M-1)(1+1/r)^2}, or the
#' Barnard & Rubin (1999) small-sample correction when \code{df_fun} supplies
#' the complete-data residual df (recommended for the modest species counts of
#' comparative data). This matches \code{mice::pool}.
#'
#' \strong{MCMCglmm fits are accepted} (unlike some MI tools). For MCMCglmm the
#' default extractors take the posterior mean and posterior covariance of the
#' fixed effects as \eqn{(\hat\theta_i, U_i)} -- a normal approximation to the
#' exact stacked posterior. This is offered mainly for the \code{fmi}
#' diagnostic and standard MI reporting; for coefficient inference the
#' principled Bayesian combiner is \code{\link{pool_posteriors}} (stacking),
#' which needs no normality assumption. \strong{Never pool variance components /
#' phylogenetic signal / heritability with \code{pool_mi()}} -- their posteriors
#' are skewed and bounded; use \code{\link{pool_posteriors}} for those.
#'
#' @param fits A list of fits of length \code{M >= 2} (e.g. the output of
#'   \code{\link{with_imputations}}). Any class implementing \code{coef()} and
#'   \code{vcov()} works; \code{MCMCglmm} fits are handled specially.
#' @param conf.level Confidence level for the pooled interval (default 0.95).
#' @param coef_fun Optional extractor of a named coefficient vector from one
#'   fit. \code{NULL} (default) uses \code{stats::coef}, or a posterior-mean
#'   extractor for MCMCglmm fits.
#' @param vcov_fun Optional extractor of the coefficient covariance matrix.
#'   \code{NULL} (default) uses \code{stats::vcov}, or the posterior covariance
#'   for MCMCglmm fits.
#' @param df_fun Optional function returning the complete-data residual df from
#'   one fit; enables the Barnard & Rubin (1999) df. \code{NULL} uses Rubin
#'   (1987) df.
#'
#' @return A data.frame of class \code{"bace_pooled_mi"} with one row per
#'   coefficient: \code{term}, \code{estimate}, \code{std.error}, \code{df},
#'   \code{statistic}, \code{p.value}, \code{conf.low}, \code{conf.high},
#'   \code{fmi}, \code{riv}.
#'
#' @references
#' Rubin DB (1987). \emph{Multiple Imputation for Nonresponse in Surveys}. Wiley.
#'
#' Barnard J, Rubin DB (1999). Small-sample degrees of freedom with multiple
#' imputation. \emph{Biometrika} 86(4): 948-955.
#'
#' van Buuren S (2018). \emph{Flexible Imputation of Data}, 2nd ed., Sec 2.3-2.4.
#'
#' @seealso \code{\link{with_imputations}}, \code{\link{pool_posteriors}}
#' @examples \dontrun{
#' fits <- with_imputations(res, function(d) lm(y ~ x1 + x2, data = d))
#' pool_mi(fits)
#' }
#' @export
pool_mi <- function(fits, conf.level = 0.95,
                    coef_fun = NULL, vcov_fun = NULL, df_fun = NULL) {

  if (!is.list(fits)) stop("`fits` must be a list of model fits.", call. = FALSE)

  is_err <- vapply(fits, function(x)
    inherits(x, "try-error") || inherits(x, "bace_mi_error"), logical(1))
  if (any(is_err)) {
    warning(sprintf("Dropping %d fit(s) that failed in with_imputations().",
                    sum(is_err)), call. = FALSE)
    fits <- fits[!is_err]
  }

  M <- length(fits)
  if (M < 2L) stop("Need at least 2 fits to pool; got ", M, ".", call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      conf.level <= 0 || conf.level >= 1) {
    stop("`conf.level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  # Choose extractors. MCMCglmm -> posterior summaries + a guidance NOTE.
  if (inherits(fits[[1]], "MCMCglmm")) {
    if (is.null(coef_fun)) coef_fun <- .mcmcglmm_fixef_mean
    if (is.null(vcov_fun)) vcov_fun <- .mcmcglmm_fixef_vcov
    message(
      "pool_mi(): MCMCglmm fits detected. Rubin's rules on posterior summaries ",
      "is a normal approximation to the stacked posterior; it is provided for ",
      "the fmi diagnostic and standard MI reporting. For coefficient inference ",
      "the exact Bayesian combiner is pool_posteriors(). Do NOT use pool_mi() ",
      "for variance components / phylogenetic signal.")
  } else {
    if (is.null(coef_fun)) coef_fun <- stats::coef
    if (is.null(vcov_fun)) vcov_fun <- stats::vcov
  }

  # ---- Extract coefficients ----
  coefs <- lapply(fits, coef_fun)
  ok <- vapply(coefs, function(x) is.numeric(x) && !is.null(names(x)), logical(1))
  if (!all(ok)) {
    stop("`coef_fun()` must return a named numeric vector for every fit.",
         call. = FALSE)
  }
  nm_ref <- names(coefs[[1]])
  nm_ok  <- vapply(coefs, function(x) identical(names(x), nm_ref), logical(1))
  if (!all(nm_ok)) {
    off <- which(!nm_ok)[1]
    stop("Coefficient names differ across fits (Rubin's rules require a common ",
         "term set). First offending fit: index ", off, ".\n  reference: ",
         paste(nm_ref, collapse = ", "), "\n  offender:  ",
         paste(names(coefs[[off]]), collapse = ", "), call. = FALSE)
  }
  coef_mat <- do.call(rbind, lapply(coefs, function(x) unname(x[nm_ref])))

  # ---- Within-imputation variances (diagonal of vcov) ----
  vars_mat <- matrix(NA_real_, nrow = M, ncol = length(nm_ref))
  for (i in seq_len(M)) {
    V <- vcov_fun(fits[[i]])
    if (is.null(V) || !is.matrix(V) || any(dim(V) < length(nm_ref))) {
      stop("`vcov_fun()` did not return a valid matrix for fit ", i, ".",
           call. = FALSE)
    }
    dV <- diag(V)
    if (!is.null(names(dV)) && all(nm_ref %in% names(dV))) {
      vars_mat[i, ] <- dV[nm_ref]
    } else if (length(dV) == length(nm_ref)) {
      vars_mat[i, ] <- dV
    } else {
      stop("`vcov_fun()` returned a matrix of the wrong size for fit ", i, ".",
           call. = FALSE)
    }
  }

  # ---- Rubin's rules ----
  theta_bar <- colMeans(coef_mat)
  W         <- colMeans(vars_mat)
  B         <- apply(coef_mat, 2, stats::var)            # (M-1) denominator
  total_var <- W + (1 + 1 / M) * B
  se_pool   <- sqrt(total_var)

  r      <- ifelse(W > 0, (1 + 1 / M) * B / W, Inf)
  lambda <- ifelse(total_var > 0, (1 + 1 / M) * B / total_var, 1)

  v_old <- ifelse(is.finite(r) & r > 0, (M - 1) * (1 + 1 / r)^2, M - 1)

  if (!is.null(df_fun)) {
    nu_com <- tryCatch(df_fun(fits[[1]]), error = function(e) NA_real_)
    if (is.numeric(nu_com) && length(nu_com) == 1L && is.finite(nu_com) &&
        nu_com > 0) {
      v_obs <- ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - lambda)
      v_bar <- 1 / (1 / v_old + 1 / v_obs)
    } else {
      v_bar <- v_old
    }
  } else {
    v_bar <- v_old
  }

  fmi <- ifelse(is.finite(r), (r + 2 / (v_bar + 3)) / (r + 1), 1)

  t_stat <- theta_bar / se_pool
  p_val  <- ifelse(se_pool > 0, 2 * stats::pt(-abs(t_stat), df = v_bar), NA_real_)
  alpha  <- 1 - conf.level
  q      <- stats::qt(1 - alpha / 2, df = v_bar)

  out <- data.frame(
    term      = nm_ref,
    estimate  = theta_bar,
    std.error = se_pool,
    df        = v_bar,
    statistic = t_stat,
    p.value   = p_val,
    conf.low  = theta_bar - q * se_pool,
    conf.high = theta_bar + q * se_pool,
    fmi       = fmi,
    riv       = r,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "m")          <- M
  attr(out, "conf.level") <- conf.level
  class(out) <- c("bace_pooled_mi", "data.frame")
  out
}

# ---- MCMCglmm fixed-effect extractors --------------------------------------

#' Column indices of the fixed effects in an MCMCglmm Sol matrix
#' @keywords internal
.mcmcglmm_fixef_idx <- function(model) {
  sol <- as.matrix(model$Sol)
  nfl <- tryCatch(model$Fixed$nfl, error = function(e) NULL)
  if (!is.null(nfl) && is.finite(nfl) && nfl > 0L && nfl <= ncol(sol)) {
    return(seq_len(nfl))            # fixed effects come first in Sol
  }
  # Fallback: drop random-effect BLUP columns (named "<term>.<level>").
  which(!grepl("\\.", colnames(sol), perl = TRUE))
}

#' @keywords internal
.mcmcglmm_fixef_mean <- function(model) {
  sol <- as.matrix(model$Sol)
  idx <- .mcmcglmm_fixef_idx(model)
  stats::setNames(colMeans(sol[, idx, drop = FALSE]), colnames(sol)[idx])
}

#' @keywords internal
.mcmcglmm_fixef_vcov <- function(model) {
  sol <- as.matrix(model$Sol)
  idx <- .mcmcglmm_fixef_idx(model)
  V <- stats::cov(sol[, idx, drop = FALSE])
  dimnames(V) <- list(colnames(sol)[idx], colnames(sol)[idx])
  V
}

#' @export
print.bace_pooled_mi <- function(x, digits = 4, ...) {
  cat(sprintf("Pooled estimates from %d imputations (Rubin's rules)\n",
              attr(x, "m")))
  cat(sprintf("Confidence level: %.0f%%\n\n", 100 * attr(x, "conf.level")))
  df_show <- x
  num <- vapply(df_show, is.numeric, logical(1))
  df_show[num] <- lapply(df_show[num],
                         function(v) formatC(v, digits = digits, format = "fg",
                                             flag = "#"))
  class(df_show) <- "data.frame"
  print(df_show, row.names = FALSE)
  max_fmi <- suppressWarnings(max(x$fmi, na.rm = TRUE))
  if (is.finite(max_fmi) && max_fmi > 0.5) {
    cat(sprintf(
      "\nNote: max fmi = %.2f. Consider more imputations for stable SEs.\n",
      max_fmi))
  }
  invisible(x)
}
