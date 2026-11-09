# archive-2026-09/rubin_discrete_castor.R
#
# RETIRED 2026-11-09. Castor Mk models that filled the discrete traits for the FREQUENTIST
# comparison arms (freqA / freqB). The 2026-11 campaign tests BACE, so these arms are off and this
# code is unused. Moved out of code/rubin_discrete.R intact rather than deleted: if the competitor
# benchmark is revived, source this file beside rubin_discrete.R and it works as before.
#
# Original header follows.
#
# Frequentist discrete arms mirror the continuous freq A / freq B:
#   freq B analogue: one castor Mk fit (castor::hsp_mk_model), then M joint conditional draws of
#     the missing tips at the fitted rates.
#   freq A analogue: for each m, a parametric bootstrap of the rates, then one joint conditional
#     draw at the refitted Q.
# Rate models: the default CASTOR_FLEX (all rates different, ARD, for binary and categorical;
# stepwise all-rates-different, SRD, for ordinal). Its stationary distribution can match the
# observed class frequencies, so at low phylogenetic signal it falls back to them, as the
# continuous arm's Pagel lambda falls back to the sample mean. Equal rates (ER) cannot: at low
# lambda its rate runs to the bound and it imputes uniform classes. Simulation v1's set
# (ER, ER, SUEDE) is kept as CASTOR_V1 for a comparison arm.
# mk_joint_draw is forward filtering, backward sampling on the tree: exact draws from the joint
# distribution of all missing tips given the observed tips.
#
# Requires disc_traits_of() from rubin_discrete.R.

CASTOR_FLEX <- c(binary = "ARD", categorical = "ARD", ordinal = "SRD")
CASTOR_V1   <- c(binary = "ER",  categorical = "ER",  ordinal = "SUEDE")
.trait_kind <- function(f) if (is.ordered(f)) "ordinal" else if (nlevels(f) == 2L) "binary" else "categorical"

.edge_P <- function(Q, t) {
  lapply(t, function(s) { P <- as.matrix(Matrix::expm(Q * s)); P[P < 0] <- 0; P / rowSums(P) })
}

#' Joint conditional draws of the missing tip states of one Mk trait.
#' @param tree ape phylo; @param tip integer states 1..K in tree$tip.label order, NA = missing
#' @param Q K x K transition-rate matrix; @param root_prior length-K probabilities; @param M number of draws
#' @return M x Ntip integer matrix of complete tip states (observed tips kept). Stops if the observed tips have zero
#'   likelihood under Q (possible for a bootstrap refit with rates at exactly 0), instead of returning NA draws.
mk_joint_draw <- function(tree, tip, Q, root_prior, M) {
  K <- nrow(Q); nt <- length(tree$tip.label); nn <- tree$Nnode
  tr <- ape::reorder.phylo(tree, "postorder")
  P <- .edge_P(Q, tr$edge.length)
  lik <- matrix(1, nt + nn, K)
  obs <- !is.na(tip); lik[which(obs), ] <- 0; lik[cbind(which(obs), tip[obs])] <- 1
  for (i in seq_len(nrow(tr$edge))) {
    p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
    lik[p, ] <- lik[p, ] * as.numeric(P[[i]] %*% lik[c, ])
    mx <- max(lik[p, ])
    if (!is.finite(mx) || mx <= 0) stop("the observed tips have zero likelihood under Q")
    lik[p, ] <- lik[p, ] / mx              # rescale against underflow; sampling ignores the scale
  }
  root <- nt + 1L
  draw_cat <- function(pr) {                 # pr: M x K unnormalised probabilities, one row per draw
    rs <- rowSums(pr)
    if (any(!is.finite(rs) | rs <= 0)) stop("the observed tips have zero likelihood under Q")
    cp <- t(apply(pr / rs, 1, cumsum))
    1L + rowSums(cp < stats::runif(nrow(pr)))
  }
  st <- matrix(NA_integer_, M, nt + nn)
  st[, root] <- draw_cat(matrix(root_prior * lik[root, ], M, K, byrow = TRUE))
  for (i in rev(seq_len(nrow(tr$edge)))) {   # reverse postorder visits every parent before its children
    p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
    st[, c] <- draw_cat(P[[i]][st[, p], , drop = FALSE] * matrix(lik[c, ], M, K, byrow = TRUE))
  }
  st[, seq_len(nt), drop = FALSE]
}

.empirical_prior <- function(tip, K) { f <- tabulate(tip[!is.na(tip)], K); f / sum(f) }

.castor_fit <- function(tree, tip, K, rate_model) {
  h <- castor::hsp_mk_model(tree, tip_states = tip, Nstates = K, rate_model = rate_model, Ntrials = 3,
                            root_prior = "empirical", include_likelihoods = FALSE, Nthreads = 1)
  if (!isTRUE(h$success)) stop("castor hsp_mk_model failed: ", h$error %||% "unknown")
  h$transition_matrix
}
`%||%` <- function(x, y) if (is.null(x)) y else x

#' castor arms for every discrete trait of a cell.
#' @param cell a make_cell() result; @param M draws; @param proper TRUE = freq A analogue (bootstrap rates)
#' @param base_sets optional list of M data.frames (e.g. the freq A continuous draws) whose discrete columns are
#'   filled; default copies of df_miss
#' @param models named rate models for binary, categorical and ordinal traits (CASTOR_FLEX or CASTOR_V1)
#' @return list(datasets, diag = per-trait data.frame: rate model, degenerate flag, largest exit rate of the fit and
#'   median over the refits used, bootstrap re-simulations (< 2 observed classes), failed refits, rejected refits
#'   (observed tips impossible under Q*), draws that fell back to the fitted Q, error)
mi_castor <- function(cell, M, proper, base_sets = NULL, models = CASTOR_FLEX) {
  tree <- cell$tree; df <- cell$df_miss
  sets <- base_sets %||% rep(list(df), M)
  dg <- list()
  for (v in disc_traits_of(cell$truth)) {
    f <- df[[v]]; lev <- levels(f); K <- length(lev)
    tip <- as.integer(f)[match(tree$tip.label, rownames(df))]
    rm <- models[[.trait_kind(f)]]
    miss_sp <- tree$tip.label[is.na(tip)]
    row <- data.frame(trait = v, rate_model = rm, degenerate = FALSE, q_hat = NA_real_, q_star_med = NA_real_,
                      n_resim = 0L, n_refit_fail = 0L, n_reject = 0L, n_fallback = 0L, error = NA_character_)
    if (length(unique(tip[!is.na(tip)])) < 2L) {   # castor cannot fit; leave NA for fill_degenerate()
      row$degenerate <- TRUE; dg[[v]] <- row; next
    }
    res <- tryCatch({
      Q <- .castor_fit(tree, tip, K, rm); prior <- .empirical_prior(tip, K)
      n_resim <- 0L; n_fail <- 0L; n_reject <- 0L; n_fallback <- 0L
      if (!proper) {
        draws <- mk_joint_draw(tree, tip, Q, prior, M); q_star <- rep(max(-diag(Q)), M)
      } else {
        draws <- matrix(NA_integer_, M, length(tip)); q_star <- numeric(M)
        obs <- !is.na(tip)
        for (m in seq_len(M)) {
          dm <- NULL; tries <- 0L
          # A bootstrap Q* is kept only if the observed tips are possible under it (a refit can put rates at exactly
          # 0); otherwise it is redrawn. After 20 failed tries the fitted Q is used. All three events are counted.
          while (is.null(dm) && tries < 20L) {
            tries <- tries + 1L
            repeat {                           # a bootstrap sample needs >= 2 observed classes to refit
              sim <- castor::simulate_mk_model(tree, Q, root_probabilities = prior, include_nodes = FALSE)$tip_states
              tip_b <- ifelse(obs, sim, NA_integer_)
              if (length(unique(tip_b[obs])) >= 2L || n_resim > 50L * M) break
              n_resim <- n_resim + 1L
            }
            Qm <- tryCatch(.castor_fit(tree, tip_b, K, rm), error = function(e) NULL)
            if (is.null(Qm)) { n_fail <- n_fail + 1L; next }
            dm <- tryCatch(mk_joint_draw(tree, tip, Qm, prior, 1L), error = function(e) NULL)
            if (is.null(dm)) n_reject <- n_reject + 1L
          }
          if (is.null(dm)) { Qm <- Q; dm <- mk_joint_draw(tree, tip, Q, prior, 1L); n_fallback <- n_fallback + 1L }
          draws[m, ] <- dm; q_star[m] <- max(-diag(Qm))
        }
      }
      list(draws = draws, q_hat = max(-diag(Q)), q_star = q_star, n_resim = n_resim, n_fail = n_fail,
           n_reject = n_reject, n_fallback = n_fallback)
    }, error = function(e) e)
    if (inherits(res, "error")) { row$error <- conditionMessage(res); dg[[v]] <- row; next }
    draws <- res$draws; colnames(draws) <- tree$tip.label
    for (m in seq_len(M)) {
      x <- as.character(sets[[m]][[v]]); x[match(miss_sp, rownames(df))] <- lev[draws[m, miss_sp]]
      sets[[m]][[v]] <- factor(x, levels = lev, ordered = is.ordered(f))
    }
    row$q_hat <- res$q_hat; row$q_star_med <- stats::median(res$q_star)
    row$n_resim <- res$n_resim; row$n_refit_fail <- res$n_fail; row$n_reject <- res$n_reject
    row$n_fallback <- res$n_fallback
    dg[[v]] <- row
  }
  list(datasets = sets, diag = do.call(rbind, dg))
}

