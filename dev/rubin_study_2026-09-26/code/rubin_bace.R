# script/rubin_bace.R
#
# BACE arms for the Rubin's-rules MI-arm comparison lane. Assumes campaign_gnn_off_lib.R has
# already been sourced (make_cell(), run_bace()'s data-prep conventions, bace_diagnostics());
# rubin_cell.R sources it first. BACE itself is installed as a package and never edited here.
#
# 2026-10 campaign: re-run on the CORRECTED package. The original (2026-09) campaign's
# "BACE as shipped" arm measured the pre-fix package, in which every final imputation started
# from one converged dataset whose missing entries held point fills, so the M datasets shared
# those fills (improper MI; found by this study, independently reported in PR #15, and fixed in
# BACE on 2026-10-26). The installed package now draws each final dataset from an independent run
# with `sweeps` internal Gibbs passes (default 3) -- the design this study's own discussion
# recommended ("independent chains, or several sweeps between retained imputations").
#
# Arms on top of one BACE::bace() fit:
#   mi_bace_shipped() -- BACE's own n_final imputed datasets, unmodified. Against the corrected
#                        package this IS the released algorithm; display label "Current BACE"
#                        (the internal arm code stays "bace" so tags/aggregation are unchanged).
#   mi_bace_chain()   -- the 2026-09 chained variant, kept for comparison: M sequential
#                        bace_final_imp(n_final = 1, sweeps = 1) calls, run m starting from
#                        dataset m - 1 (consecutive states of a single Gibbs chain, one sweep
#                        apart). sweeps = 1L is pinned so this arm remains exactly the algorithm
#                        studied in 2026-09 now that the package default is 3.
#
# The 2026-09 "bace_resid" negative-control arm (a second post-hoc residual draw, probing the
# old defect) is removed with the defect; its code and the premise-correction notes live in git
# history (import commit eaca7e3).

#' Fit BACE with n_final = M, reusing run_bace()'s data-prep conventions (script/campaign_gnn_off_lib.R
#' L601-653: Species column from rownames, one formula per trait regressed on all the others,
#' n_cores = 1L, skip_conv = TRUE, ovr_categorical = TRUE).
#'
#' @param cell a make_cell() result (uses cell$tree, cell$df_miss)
#' @param M number of BACE final imputation datasets (BACE's n_final)
#' @param nitt,burnin,thin,runs passed straight to BACE::bace()
#' @return list(outb = the "bace_complete" object, wall_s = fit wall time in seconds,
#'   converged = BACE's own convergence verdict, ess_med = median effective sample size over the
#'   final fits' fixed effects, diag = the full bace_diagnostics() list -- same extraction G9b/v1 use)
fit_bace_mi <- function(cell, M = 20L, nitt, burnin, thin, runs) {
  tree_b <- cell$tree
  if (any(tree_b$edge.length == 0)) tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
  # Input cleaning (2026-09-24, campaign): at lambda = 1 with fixed thresholds most datasets leave a discrete
  # level unobserved, and BACE passes the unused level to MCMCglmm, which stops at once with "Mixed model
  # equations singular". A careful user drops unused levels, so do that; a discrete trait with fewer than two
  # observed classes cannot be modelled and is left out of BACE's model (discrete traits are not scored).
  # A no-op on datasets without empty levels. Every change is recorded in diag$input_fix.
  df_b <- cell$df_miss; input_fix <- character(0)
  for (v in names(df_b)) if (is.factor(df_b[[v]])) {
    obs_lv <- unique(as.character(df_b[[v]][!is.na(df_b[[v]])]))
    if (length(obs_lv) < 2L) {
      df_b[[v]] <- NULL; input_fix <- c(input_fix, sprintf("%s: left out (%d observed class)", v, length(obs_lv)))
    } else if (nlevels(droplevels(df_b[[v]])) < nlevels(df_b[[v]])) {
      input_fix <- c(input_fix, sprintf("%s: dropped %d empty level(s)", v, nlevels(df_b[[v]]) - nlevels(droplevels(df_b[[v]]))))
      df_b[[v]] <- droplevels(df_b[[v]])
    }
  }
  df_b$Species <- rownames(cell$df_miss)
  all_traits <- setdiff(names(df_b), "Species")
  fixformula <- lapply(all_traits, function(v)
    paste0(v, " ~ ", paste(setdiff(all_traits, v), collapse = " + ")))

  t0 <- Sys.time()
  outb <- BACE::bace(fixformula = fixformula, ran_phylo_form = "~1|Species", phylo = tree_b,
                     data = df_b, nitt = nitt, burnin = burnin, thin = thin, runs = runs,
                     n_final = M, n_cores = 1L, verbose = FALSE, skip_conv = TRUE,
                     ovr_categorical = TRUE)
  wall_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  d <- bace_diagnostics(outb)   # same convergence/ESS extraction v1/G9b use (campaign_gnn_off_lib.R)
  d$input_fix <- input_fix
  list(outb = outb, wall_s = wall_s, converged = d$converged, ess_med = d$ess_med, diag = d,
       fixformula = fixformula, tree_b = tree_b, mcmc = c(nitt = nitt, burnin = burnin, thin = thin))
}

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
                               verbose = FALSE, n_cores = 1L, ovr_categorical = TRUE,
                               sweeps = 1L)  # pin one sweep per link: the 2026-09 chain algorithm
                                             # (package default is now 3; see file header)
    raw[[m]] <- r$all_datasets[[1]]; models[[m]] <- r$all_models[[1]]
    obj$data[[k]] <- raw[[m]]   # the next run starts from this draw
  }
  wall_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(datasets = mi_bace_shipped(list(imputed_datasets = raw), df_miss), models = models, wall_s = wall_s)
}

#' BACE's own n_final imputed datasets, unchanged, aligned to df_miss's row order/column
#' names/classes (BACE's internal Species column dropped; storage class only, never values, is
#' coerced back to df_miss's -- e.g. BACE returns a rounded double for count traits, df_miss has
#' them as integer).
#'
#' @param outb a "bace_complete" object (fit_bace_mi()$outb)
#' @param df_miss the cell's missing-data data.frame (make_cell()$df_miss)
#' @return list of length M, each a data.frame shaped like df_miss
mi_bace_shipped <- function(outb, df_miss) {
  lapply(outb$imputed_datasets, function(d) {
    for (v in setdiff(names(df_miss), names(d))) d[[v]] <- df_miss[[v]][match(rownames(d), rownames(df_miss))]  # left-out trait
    d <- d[rownames(df_miss), names(df_miss), drop = FALSE]
    for (v in names(df_miss)) {
      if (is.integer(df_miss[[v]]) && !is.integer(d[[v]])) {
        d[[v]] <- as.integer(round(d[[v]]))
      } else if (is.factor(df_miss[[v]]) && !identical(levels(d[[v]]), levels(df_miss[[v]]))) {
        d[[v]] <- factor(as.character(d[[v]]), levels = levels(df_miss[[v]]),
                         ordered = is.ordered(df_miss[[v]]))
      }
    }
    d
  })
}
