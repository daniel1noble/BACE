# Numerical correctness tests.
# These go beyond "function returns the right shape" to verify that
# BACE actually recovers known structure under controlled simulation:
#   - imputed values track ground truth on continuous traits
#   - BACE outperforms a column-mean baseline when phylogenetic
#     signal is real
#   - phylo_signal_summary recovers a known phylo signal within
#     tolerance
#   - posterior pooling preserves dimensionality and produces
#     sensible variance
#
# Every test is deterministic (set.seed at the top). Tolerances are
# tuned to catch regressions but loose enough that the small-sample
# noise inherent in MCMC and tree simulation doesn't produce flakes.
# All tests skip on CRAN to keep CRAN check fast.

# ---- Shared simulation helpers --------------------------------------------

# Generate a small phylogeny + a continuous response with known phylo
# signal, optionally with a single continuous predictor. Returns the
# ground-truth (complete) data, plus a copy with a fraction of the
# response masked at random.
make_continuous_fixture <- function(n_species = 40,
                                     phylo_signal = 0.8,
                                     n_predictors = 2,  # 2+ -- generate
                                                        # _default_beta_matrix
                                                        # requires it
                                     miss_frac = 0.30,
                                     seed = 2026) {
  set.seed(seed)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "gaussian",
    predictor_types = rep("gaussian", n_predictors),
    phylo_signal    = rep(phylo_signal, n_predictors + 1),
    n_cases         = n_species,
    n_species       = n_species
  )))
  truth  <- out$complete_data
  masked <- out$data
  resp_col <- out$params$var_names[1]
  # sim_bace produces a 'Species' column (post-fix, matches BACE
  # convention used elsewhere in the package).

  # Apply a clean MCAR mask on the response only -- sim_bace's
  # missingness=NULL leaves the response complete, which is what we
  # want for a clean truth-vs-imputed comparison.
  set.seed(seed + 1L)
  obs_idx  <- which(!is.na(masked[[resp_col]]))
  mask_idx <- sample(obs_idx, floor(length(obs_idx) * miss_frac))
  masked[[resp_col]][mask_idx] <- NA

  list(
    truth      = truth,
    masked     = masked,
    tree       = out$tree,
    resp_col   = resp_col,
    pred_cols  = out$params$var_names[-1],
    mask_idx   = mask_idx
  )
}

# ---- Continuous imputation: recovers truth + beats baseline ---------------

test_that("BACE imputed values correlate with truth (high phylo signal)", {
  testthat::skip_on_cran()
  # n_species = 120, not 40. At 40 species a 30 per cent mask leaves only 12
  # hidden cells, and a correlation on 12 points is far too noisy to assert
  # anything: measured over 20 seeds it ran from -0.29 to +0.93 (mean 0.487,
  # sd 0.410) and fell below this test's own threshold in 7 of 20. The test
  # was passing on a lucky seed, not on a property of the code -- it failed
  # the moment the simulator's draws changed for an unrelated reason, while
  # the old and new simulators recover equally well (mean r 0.510 vs 0.528
  # over 30 seeds, p = 0.85).
  #
  # With 120 species the mask leaves 36 cells and the same measurement is
  # stable: mean 0.724, sd 0.153, minimum 0.406 over 20 seeds, none below the
  # threshold. The extra species cost a few seconds and buy a test that means
  # something.
  fix <- make_continuous_fixture(n_species = 120,
                                 phylo_signal = 0.8,
                                 n_predictors = 1,
                                 miss_frac = 0.30,
                                 seed = 2026)

  set.seed(2026)
  res <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = paste(fix$resp_col, "~", fix$pred_cols[1]),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$tree,
    data           = fix$masked,
    runs           = 3,
    nitt = 2000, thin = 5, burnin = 500,
    verbose = FALSE
  )))

  # Take the last imputation as the point estimate.
  imp_last <- res$data[[length(res$data)]]
  truth_v  <- fix$truth[[fix$resp_col]][fix$mask_idx]
  imp_v    <- imp_last[[fix$resp_col]][fix$mask_idx]

  expect_true(all(!is.na(imp_v)))
  # With phylo_signal = 0.8 and a real predictor, BACE correlates with truth
  # on the hidden cells at about r = 0.72 (20-seed mean, sd 0.15). The 0.30
  # threshold sits roughly 2.8 sd below that, so it flags a real regression
  # without flaking on MCMC noise.
  r <- suppressWarnings(stats::cor(imp_v, truth_v))
  expect_gt(r, 0.30)
})

test_that("BACE beats column-mean baseline on continuous trait at high signal", {
  testthat::skip_on_cran()
  # Use a larger fixture (80 spp) + longer chain so MCMC noise
  # doesn't dominate the lift-vs-baseline comparison. At 40 spp the
  # signal/noise ratio is too thin to robustly detect a sub-5% lift.
  fix <- make_continuous_fixture(n_species = 80,
                                 phylo_signal = 0.8,
                                 n_predictors = 2,
                                 miss_frac = 0.30,
                                 seed = 31415)

  set.seed(31415)
  res <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = paste(fix$resp_col, "~",
                           paste(fix$pred_cols, collapse = " + ")),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$tree,
    data           = fix$masked,
    runs           = 4,
    nitt = 4000, thin = 5, burnin = 1000,
    verbose = FALSE
  )))

  # Average over the last 3 imputations to reduce per-iteration MCMC
  # noise -- still much fewer than n_final but stable enough for the
  # lift assertion.
  imp_runs <- res$data[(length(res$data) - 2):length(res$data)]
  imp_v <- rowMeans(vapply(imp_runs,
                            function(d) d[[fix$resp_col]][fix$mask_idx],
                            numeric(length(fix$mask_idx))))
  truth_v <- fix$truth[[fix$resp_col]][fix$mask_idx]

  bace_rmse     <- sqrt(mean((imp_v - truth_v)^2))
  baseline_pred <- mean(fix$masked[[fix$resp_col]], na.rm = TRUE)
  baseline_rmse <- sqrt(mean((baseline_pred - truth_v)^2))

  # At 80 spp + phylo signal 0.8 + 2 predictors + 4-run chain,
  # BACE typically beats baseline by 15-50%. We assert >= 5% to
  # make this robust to MCMC noise but still catch a substantive
  # regression.
  expect_lt(bace_rmse, baseline_rmse * 0.95)
})

# ---- Phylogenetic signal recovery -----------------------------------------

test_that("phylo_signal_summary recovers high simulated lambda", {
  testthat::skip_on_cran()

  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace_gaussian(
    n_predictors = 1, n_cases = 60, n_species = 60,
    phylo_signal = 0.9, beta_sparsity = 1.0  # no covariate effect
  )))
  resp <- out$params$var_names[1]
  d <- out$complete_data
  d$Species <- d$Species  # already the species column; explicit

  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = d,
    tree        = out$tree,
    species_col = "Species",
    variables   = resp,
    quick       = TRUE,
    verbose     = FALSE
  )))

  expect_s3_class(res, "phylo_signal")
  expect_true("table" %in% names(res))
  expect_true("lambda" %in% colnames(res$table))
  lam <- res$table$lambda[res$table$variable == resp]
  expect_false(is.na(lam))
  # Phylo signal 0.9 should give a moderate-to-high lambda. Gibbs
  # / MCMC noise at n=60 routinely produces lambda in (0.3, 1.0);
  # set a loose floor.
  expect_gt(as.numeric(lam), 0.30)
})

# ---- Posterior pooling: pooled variance is reasonable ---------------------

test_that("pool_posteriors output has expected dimensions and finite stats", {
  testthat::skip_on_cran()

  fix <- make_continuous_fixture(n_species = 30, phylo_signal = 0.5,
                                 n_predictors = 1, miss_frac = 0.30,
                                 seed = 99)

  set.seed(99)
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = paste(fix$resp_col, "~", fix$pred_cols[1]),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$tree,
    data           = fix$masked,
    runs           = 2, nitt = 800, thin = 5, burnin = 200,
    verbose = FALSE
  )))
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = paste(fix$resp_col, "~", fix$pred_cols[1]),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$tree,
    nitt = 800, thin = 5, burnin = 200,
    n_final = 3, verbose = FALSE
  )))
  pooled <- suppressWarnings(suppressMessages(pool_posteriors(step2)))

  # Structural checks
  expect_s3_class(pooled, "bace_pooled")
  expect_true(!is.null(pooled$models))
  expect_true(length(pooled$models) >= 1L)

  # Numerical checks on the pooled model for the response variable.
  resp_model <- pooled$models[[fix$resp_col]]
  expect_s3_class(resp_model$Sol, "mcmc")
  sol_mat <- as.matrix(resp_model$Sol)
  # Pooled chain should have rows = sum(per-imputation iterations),
  # which is positive and finite.
  expect_gt(nrow(sol_mat), 0L)
  expect_true(all(is.finite(sol_mat)))
  # Posterior variance of intercept should be > 0 (otherwise pooling
  # is degenerate).
  intercept_col <- which(colnames(sol_mat) == "(Intercept)")
  if (length(intercept_col)) {
    expect_gt(stats::var(sol_mat[, intercept_col]), 0)
  }
})

# ---- Convergence check gives real verdicts on real fits --------------------

test_that("the stability check passes a healthy real fit and reports verdicts", {
  testthat::skip_on_cran()
  fix <- make_continuous_fixture(n_species = 25, phylo_signal = 0.5,
                                 n_predictors = 1, miss_frac = 0.30,
                                 seed = 7)

  set.seed(7)
  fit <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = paste(fix$resp_col, "~", fix$pred_cols[1]),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$tree,
    data           = fix$masked,
    runs           = 6, nitt = 4000, thin = 5, burnin = 1000,
    verbose        = FALSE
  )))
  conv <- suppressWarnings(assess_convergence(fit))

  # The pre-1.1 test here ("flips OK as chains lengthen") asserted only the
  # S3 class — it could not fail. This one asserts the verdict: a healthy
  # gaussian fit with a decent chain must PASS the repaired check (the
  # pre-1.1 criteria failed settled chains like this one ~97% of the time,
  # which is the bug this file now guards against). Deterministic by seed.
  expect_s3_class(conv, "bace_convergence")
  expect_true(conv$converged)
  expect_true(all(conv$verdicts == "ok"))
  # and the monitored series are actually populated for every variable
  expect_gte(ncol(conv$summary_stats), 2)
  expect_identical(nrow(conv$summary_stats), 6L)
})
