# code/rubin_covariate.R
#
# THE AUXILIARY COVARIATE -- one well-measured variable that helps impute the others.
#
# WHAT IT IS: think body size. You have it for every species, you never have to impute it, and it
# tells you something about the trait that IS missing. That is the practical advice the paper
# wants to give: when phylogenetic signal is weak, a good auxiliary variable can carry the
# imputation instead.
#
# WHERE IT COMES FROM: build_covs() in campaign_gnn_off_lib.R already builds exactly this. For one
# covariate targeting one trait it is
#
#     cov  =  rho_cov * (the target trait's latent)  +  sqrt(1 - rho_cov^2) * e
#
# where e is noise that itself evolves on the tree, so the covariate looks like a real comparative
# trait and not a random column. The construction is correlated exactly rho_cov with its target,
# has unit variance, and is a valid (positive-definite) joint distribution for any rho_cov in
# (-1, 1).
#
# WHY ONE COVARIATE AND NOT A CORRELATION APPLIED TO EVERYTHING. The campaign's --rho flag fills a
# compound-symmetry matrix: it sets EVERY trait-trait pair at once, all 21 of them at K = 7. That
# is a different experiment, and it cannot be pushed to 0.8 anyway. For one variable correlated r
# with K traits that are themselves mutually correlated rho, the joint matrix is positive definite
# only while
#
#     r^2  <  (1 + (K - 1) * rho) / K
#
# At K = 7, rho = 0.5 that is r < 0.756, so "0.8 with everything" is not a design we could run
# even if we wanted it. 0.8 with ONE trait is both runnable and the realistic case.
#
# WHAT IT IS NOT: it is not the MAR missingness driver. The driver `d1` is a separate, fully
# observed column, hard-wired at correlation 0.35, whose job is to CAUSE the gaps. The auxiliary's
# job is to help FILL them. Two different variables doing two different jobs.
#
# THE FAILURE MODE THIS FILE EXISTS TO CATCH. A covariate that is built but never reaches the
# imputation model costs nothing, crashes nothing and changes no output -- the run completes and
# silently reports a covariate-free result. That is not hypothetical: it is what the campaign has
# been doing since it was written, because rubin_cell.R called make_cell() without n_cov. So every
# run now records the covariate's REALISED correlation with its target, and the test suite asserts
# it comes back where we asked for it.


#' Did the covariate actually turn out the way we asked?
#'
#' Measures the realised correlation between each covariate and the trait it targets, on the
#' observed data scale -- which is the scale a reader cares about and is NOT the same as rho_cov
#' once a latent has been turned into a count or a category.
#'
#' @param covs the covariate data.frame (make_cell()$covs), or NULL if none was built
#' @param truth the complete trait data.frame (make_cell()$truth)
#' @param target_names trait name(s) each covariate targets, recycled to length(covs)
#' @param rho_cov the correlation we asked for, recorded alongside for comparison
#' @return a one-row-per-covariate data.frame, or NULL when there is no covariate
covariate_check <- function(covs, truth, target_names = "c1", rho_cov = NA_real_) {

  ## STEP 1. No covariate is a perfectly legitimate state (two of the three regimes have none),
  ##         so say so quietly rather than erroring.
  if (is.null(covs) || !ncol(covs)) return(NULL)

  ## STEP 2. Line each covariate up with the trait it is supposed to be a proxy for.
  tg <- rep_len(target_names, ncol(covs))

  ## STEP 3. For each one, correlate it with its target on the OBSERVED scale. A factor has to be
  ##         made numeric first: an ordered factor by its rank, an unordered one is not meaningfully
  ##         correlated at all and is reported as NA rather than as a misleading number.
  as_numeric_trait <- function(x) {
    if (is.numeric(x)) return(as.numeric(x))
    if (is.ordered(x)) return(as.numeric(x))
    if (is.factor(x) && nlevels(x) == 2L) return(as.numeric(x) - 1)
    rep(NA_real_, length(x))
  }

  do.call(rbind, lapply(seq_len(ncol(covs)), function(j) {
    y <- as_numeric_trait(truth[[tg[j]]])
    x <- as.numeric(covs[[j]][match(rownames(truth), rownames(covs))])
    data.frame(cov = names(covs)[j], target = tg[j], rho_cov_requested = rho_cov,
               rho_cov_realised = if (all(is.na(y))) NA_real_ else stats::cor(x, y, use = "complete.obs"))
  }))
}


#' Attach the covariates to a dataset so BACE can use them as predictors.
#'
#' This is the one step that was missing. It does two things and nothing else:
#'   - puts the covariate columns beside the traits, so the model can see them;
#'   - returns their names, so the caller can add them to the RIGHT-HAND side of every trait's
#'     formula WITHOUT ever building a formula that has a covariate as the response.
#'
#' That second point is the whole distinction between a predictor and a target. A covariate is
#' fully observed, so there is nothing to impute; giving it its own imputation model would be
#' wasted compute at best and would let it drift at worst.
#'
#' @param df the model data.frame (traits, plus a Species column)
#' @param covs the covariate data.frame, or NULL
#' @param row_names the row order to align the covariates to
#' @return list(df = df with covariate columns appended, cov_names = character vector)
attach_covs <- function(df, covs, row_names) {

  ## STEP 1. Nothing to attach: hand back the data unchanged and an empty name list, so every
  ##         caller can use one code path whether or not a covariate exists.
  if (is.null(covs) || !ncol(covs)) return(list(df = df, cov_names = character(0)))

  ## STEP 2. Align by species name, never by position. The covariates were built in tree-tip order
  ##         and the model data.frame may not be.
  aligned <- covs[match(row_names, rownames(covs)), , drop = FALSE]
  rownames(aligned) <- row_names

  ## STEP 3. Append, and report the names back.
  list(df = cbind(df, aligned), cov_names = names(covs))
}
