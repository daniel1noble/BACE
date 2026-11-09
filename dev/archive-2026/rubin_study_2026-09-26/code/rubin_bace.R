# script/rubin_bace.R
#
# BACE arms for the Rubin's-rules MI-arm comparison lane.
#
# REQUIRES, sourced before this file:
#   campaign_gnn_off_lib.R  make_cell(), run_bace()'s data-prep conventions, bace_diagnostics()
#   rubin_covariate.R       attach_covs(), which puts the auxiliary covariate into the model data
# rubin_cell.R sources both first. BACE itself is installed as a package and never edited here.
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
#
# The 2026-09 "bace_resid" negative-control arm (a second post-hoc residual draw, probing the
# old defect) is removed with the defect; its code and the premise-correction notes live in git
# history (import commit eaca7e3).

#' Fit BACE with n_final = M, reusing run_bace()'s data-prep conventions (script/campaign_gnn_off_lib.R
#' L601-653: Species column from rownames, one formula per trait regressed on all the others,
#' n_cores = 1L, skip_conv = TRUE, ovr_categorical = FALSE).
#'
#' @param cell a make_cell() result (uses cell$tree, cell$df_miss)
#' @param M number of BACE final imputation datasets (BACE's n_final)
#' @param nitt,burnin,thin,runs passed straight to BACE::bace()
#' @return list(outb = the "bace_complete" object, wall_s = fit wall time in seconds,
#'   converged = BACE's own convergence verdict, ess_med = median effective sample size over the
#'   final fits' fixed effects, diag = the full bace_diagnostics() list -- same extraction G9b/v1 use)
fit_bace_mi <- function(cell, M = 20L, nitt, burnin, thin, runs, ovr_categorical = FALSE,
                        covs = cell$covs) {
  ## STEP 1. The tree. MCMCglmm cannot take a zero-length branch, so nudge any to a tiny positive
  ##         number. A no-op on almost every tree.
  tree_b <- cell$tree
  if (any(tree_b$edge.length == 0)) tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
  ## STEP 2. Input cleaning, before the model ever sees the data.
  # (2026-09-24, campaign): at lambda = 1 with fixed thresholds most datasets leave a discrete
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
  ## STEP 3. BACE wants the species as a column, not as row names.
  df_b$Species <- rownames(cell$df_miss)

  ## STEP 4. Decide what is a TRAIT (gets imputed) and what is a PREDICTOR (helps, never imputed).
  ##         Everything left in df_b is a trait. The auxiliary covariate, if there is one, is
  ##         attached here but deliberately kept OUT of all_traits -- see rubin_covariate.R.
  all_traits <- setdiff(names(df_b), "Species")
  att <- attach_covs(df_b, covs, rownames(cell$df_miss))
  df_b <- att$df; cov_names <- att$cov_names

  ## STEP 5. One imputation model per trait: that trait on the left, every OTHER trait plus every
  ##         covariate on the right. Because cov_names never appears on the left, no covariate is
  ##         ever given a model of its own.
  fixformula <- lapply(all_traits, function(v)
    paste0(v, " ~ ", paste(c(setdiff(all_traits, v), cov_names), collapse = " + ")))

  t0 <- Sys.time()
  outb <- BACE::bace(fixformula = fixformula, ran_phylo_form = "~1|Species", phylo = tree_b,
                     data = df_b, nitt = nitt, burnin = burnin, thin = thin, runs = runs,
                     n_final = M, n_cores = 1L, verbose = FALSE, skip_conv = TRUE,
                     ovr_categorical = ovr_categorical)   # default FALSE = multinomial (matches the package)
                     # since 2026-11-06: the paired coverage experiment showed one-vs-rest (TRUE) costs cat3
                     # ~0.021 set coverage and ~0.02 accuracy; bin/ord unchanged. The 2026-10 campaign results
                     # were produced with TRUE -- never resume an old cache dir across this change.
                     # coverage_experiment.R passes the argument explicitly to vary it.
  wall_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  d <- bace_diagnostics(outb)   # same convergence/ESS extraction v1/G9b use (campaign_gnn_off_lib.R)
  d$input_fix <- input_fix
  list(outb = outb, wall_s = wall_s, converged = d$converged, ess_med = d$ess_med, diag = d,
       fixformula = fixformula, tree_b = tree_b, mcmc = c(nitt = nitt, burnin = burnin, thin = thin),
       cov_names = cov_names, df_b = df_b,   # the exact data BACE saw, so a caller can check it
       ovr = ovr_categorical)
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
