# code/rubin_discrete.R
#
# SCORING for the discrete traits: binary, ordinal and categorical.
#
# The continuous traits are scored on how far the imputed number is from the truth. A category
# cannot be scored that way, so these traits get their own measures:
#
#   accuracy        how often the most frequently imputed category was the right one, with
#                   FRACTIONAL credit for ties -- if two categories tie and one is correct, that
#                   cell scores 0.5, which is what you would get on average by guessing between
#                   them. Without that, a tie would be scored as a win or a loss at random.
#   set coverage    the discrete version of a confidence interval: the smallest SET of categories
#                   holding 95% of the imputed draws, and whether the truth was in it.
#   set size        how big that set had to be. A set of "all three categories" covers the truth
#                   every time and says nothing, so coverage is only meaningful beside this.
#   Brier, ECE      how well calibrated the imputed category probabilities are.
#
# One edge case, handled the same way for every method: a trait with only ONE observed class has
# nothing to predict from. Every arm then imputes that single class (fill_degenerate), and the
# filled cells are counted rather than quietly scored as successes.
#
# Assumes campaign_gnn_off_lib.R and rubin_lib.R are sourced; rubin_cell.R does this.
#
# (The castor Mk machinery that used to live here filled the discrete traits for the frequentist
# comparison arms. Those arms are not part of the 2026-11 campaign -- it tests BACE -- so it moved
# to archive-2026-09/rubin_discrete_castor.R, intact and still callable if the competitor
# benchmark is revived.)

disc_traits_of <- function(truth) names(truth)[vapply(truth, is.factor, logical(1))]
#' A trait with a single observed class: impute that class in every missing cell, for every arm alike.
#' @return list(sets, n_filled = number of missing cells filled per trait, per dataset)
fill_degenerate <- function(sets, df_miss, traits) {
  n_filled <- setNames(integer(length(traits)), traits)
  for (v in traits) {
    cls <- unique(as.character(df_miss[[v]][!is.na(df_miss[[v]])]))
    if (length(cls) != 1L) next
    n_filled[v] <- sum(is.na(df_miss[[v]]))
    for (m in seq_along(sets)) {
      x <- as.character(sets[[m]][[v]]); x[is.na(x)] <- cls
      sets[[m]][[v]] <- factor(x, levels = levels(df_miss[[v]]), ordered = is.ordered(df_miss[[v]]))
    }
  }
  list(sets = sets, n_filled = n_filled)
}

#' Per-value scores of one arm on the discrete traits, from its M draws only (every arm on the same basis).
#' accuracy: of the modal class; brier: multi-class Brier score of the draw frequencies (sum over classes, mean over
#' cells); ece: top-label expected calibration error, 10 equal-width bins of the modal frequency;
#' set_coverage / set_size: the smallest set of classes holding >= 95% of the draws; mae_class: mean absolute class
#' distance of the modal class (ordinal only); frac_unanimous: all M draws agree. Ties (for the modal class and at the
#' edge of the 95% set) are counted fractionally: the expected value under a random tie-break.
#' A trait with one realised class in the complete data has nothing to predict and is skipped (its absence is counted
#' at aggregation). NA draws are dropped cell by cell; cells whose draws are all NA are dropped and counted.
score_discrete <- function(arm, sets, truth, mask, traits = disc_traits_of(truth)) {
  tol <- 1e-12; rows <- list()
  for (v in traits) {
    idx <- which(mask[, v]); if (!length(idx)) next
    if (nlevels(droplevels(truth[[v]])) < 2L) next
    lev <- levels(truth[[v]]); K <- length(lev)
    D <- vapply(sets, function(s) as.integer(factor(as.character(s[[v]][idx]), levels = lev)), integer(length(idx)))
    D <- matrix(D, nrow = length(idx))                 # cells x M
    n_na <- sum(is.na(D))
    keep <- rowSums(!is.na(D)) > 0L; n_all_na <- sum(!keep)
    if (!any(keep)) next
    idx <- idx[keep]; D <- D[keep, , drop = FALSE]
    Pm <- t(apply(D, 1, function(d) tabulate(d[!is.na(d)], K) / sum(!is.na(d))))
    tr <- as.integer(truth[[v]][idx]); Y <- diag(K)[tr, , drop = FALSE]
    top <- apply(Pm, 1, max); is_top <- abs(Pm - top) <= tol
    acc <- rowSums(is_top & Y == 1) / rowSums(is_top)
    in_set <- vapply(seq_along(idx), function(i) {
      p <- Pm[i, ]; ps <- sort(p, decreasing = TRUE)
      k <- which(cumsum(ps) >= 0.95 - 1e-9)[1]; thr <- ps[k]
      n_gt <- sum(p > thr + tol); n_eq <- sum(abs(p - thr) <= tol)
      pt <- p[tr[i]]
      c(if (pt > thr + tol) 1 else if (abs(pt - thr) <= tol) (k - n_gt) / n_eq else 0, k)
    }, numeric(2))
    bins <- cut(top, seq(0, 1, 0.1), include.lowest = TRUE)
    ece <- sum(tapply(seq_along(idx), bins, function(j) length(j) * abs(mean(acc[j]) - mean(top[j]))), na.rm = TRUE) /
      length(idx)
    mae <- if (is.ordered(truth[[v]])) mean(vapply(seq_along(idx), function(i) mean(abs(which(is_top[i, ]) - tr[i])), numeric(1))) else NA_real_
    rows[[v]] <- data.frame(arm = arm, trait = v, n_cells = length(idx), accuracy = mean(acc),
                            brier = mean(rowSums((Pm - Y)^2)), ece = ece,
                            set_coverage = mean(in_set[1, ]), set_size = mean(in_set[2, ]),
                            mae_class = mae, frac_unanimous = mean(top == 1), n_na = n_na, n_cells_all_na = n_all_na)
  }
  do.call(rbind, rows)
}

#' PGLS slope of y on x (Pagel's lambda by REML, the eigenbasis estimator of est_pgls_slope_fast) for any two
#' columns; x may be a two-level factor, coded 0/1 (second level = 1).
est_pgls_xy <- function(df, tree, y, x, eig = NULL) {
  d2 <- data.frame(row.names = rownames(df), c1 = df[[x]], c2 = df[[y]])
  if (is.factor(d2$c1)) d2$c1 <- as.numeric(d2$c1 == levels(d2$c1)[2L])
  est_pgls_slope_fast(d2, tree, eig = eig)
}

#' GLS slope of y on x at a fixed Pagel lambda (eigenbasis of pagel_eigen()).
.gls_slope_at <- function(y, x, eig, lam) {
  X <- cbind(1, x); Xs <- crossprod(eig$U, X); ys <- as.numeric(crossprod(eig$U, y))
  w <- 1 / (lam * eig$d + (1 - lam))
  solve(crossprod(Xs, w * Xs), crossprod(Xs, w * ys))[2, 1]
}

#' Rubin-pooled downstream estimand on the discrete trait: PGLS slope of c1 on bin (bin = "yes" coded 1).
#' Two targets are stored for coverage. The population target (filled at aggregation from the Monte Carlo of the
#' complete-data estimate) is not what any per-dataset interval estimates at lambda = 1, rho = 0.5, where the
#' dataset-specific value varies more than the standard error (review of 2026-09-26: complete-data coverage 0.31 to
#' 0.56). The per-dataset target, target_cond = rho times the GLS slope (at the true lambda) of the bin liability L5
#' on bin, is E[slope | tree, L5] with lambda known; complete-data coverage of it is about 0.96 in every cell.
#' Undefined (NULL) when bin has one class in the complete data, or is constant in every imputed dataset (one
#' observed class, filled by fill_degenerate()).
#' @param L5 the bin liability (make_cell()$L[, 5]) in rownames(truth) order; @param lambda,rho the cell's values
score_discrete_estimand <- function(arm, sets, truth, tree, eig, n, L5 = NULL, lambda = NA_real_, rho = NA_real_) {
  if (nlevels(droplevels(truth$bin)) < 2L) return(NULL)
  if (all(vapply(sets, function(s) nlevels(droplevels(s$bin)) < 2L && !anyNA(s$bin), logical(1)))) return(NULL)
  if (is.null(eig)) eig <- pagel_eigen(tree, rownames(truth))
  sl <- lapply(sets, est_pgls_xy, tree = tree, y = "c1", x = "bin", eig = eig)
  ok <- vapply(sl, function(s) is.finite(s$estimate) && is.finite(s$variance), logical(1))
  if (sum(ok) < 2L) stop(sprintf("%s: fewer than 2 analysable imputations for c1 ~ bin", arm))
  ps <- rubin_pool(vapply(sl[ok], `[[`, numeric(1), "estimate"), vapply(sl[ok], `[[`, numeric(1), "variance"),
                   df_com = n - 2)
  ref <- est_pgls_xy(truth, tree, y = "c1", x = "bin", eig = eig)
  # the complete-data interval (t with n - 2 df), the coverage baseline, as the continuous "complete" row
  ref_se <- sqrt(ref$variance); tq <- stats::qt(0.975, n - 2)
  target_cond <- if (!is.null(L5) && is.finite(lambda) && is.finite(rho))
    rho * .gls_slope_at(L5, as.numeric(truth$bin == levels(truth$bin)[2L]), eig, lambda) else NA_real_
  data.frame(arm = arm, estimand = "slope_c1_bin", estimate = ps$estimate, se = ps$se, lower = ps$lower,
             upper = ps$upper, df = ps$df, fmi = ps$fmi, m_ok = sum(ok), complete_data = ref$estimate,
             complete_se = ref_se, complete_lower = ref$estimate - tq * ref_se,
             complete_upper = ref$estimate + tq * ref_se, target_cond = target_cond)
}
