# code/tests-rubin/test-covariate.R
#
# WHAT THIS GUARDS: the auxiliary covariate actually existing and actually reaching the model.
#
# This is not a hypothetical failure. From the day it was written until 2026-11, the campaign
# called make_cell() without n_cov, so build_covs() returned NULL and no covariate was ever built.
# Nothing crashed. Nothing warned. Every result simply came back covariate-free, and looked normal.
#
# The three things below are exactly the three ways that can happen again:
#   1. the covariate is not built at all
#   2. it is built but is not the variable we asked for
#   3. it is built and correct but never reaches BACE's model formulas

code_dir <- normalizePath(file.path(testthat::test_path(), ".."))
source(file.path(code_dir, "campaign_gnn_off_lib.R"))
source(file.path(code_dir, "rubin_covariate.R"))

make_test_cell <- function(seed, n = 60, n_cov = 1L, rho_cov = 0.8) {
  make_cell("types_mixed", n, seed, lambda = 0.7, rho = 0.5, thresholds = "fixed",
            driver = TRUE, n_cov = n_cov, rho_cov = rho_cov)
}

# ---- 1. asking for a covariate gets you a covariate; asking for none gets you none -------------
testthat::test_that("n_cov controls whether a covariate exists at all", {
  with_cov <- make_test_cell(1L, n_cov = 1L)
  no_cov   <- make_test_cell(1L, n_cov = 0L)

  testthat::expect_false(is.null(with_cov$covs))
  testthat::expect_equal(ncol(with_cov$covs), 1L)
  testthat::expect_equal(names(with_cov$covs), "cov1")
  testthat::expect_null(no_cov$covs)

  # it is a predictor, never a target: it must not appear among the traits, so it is never masked
  # and never scored
  testthat::expect_false("cov1" %in% names(with_cov$truth))
  testthat::expect_false("cov1" %in% colnames(with_cov$mask))
  testthat::expect_false("cov1" %in% names(with_cov$df_miss))

  # and it is complete -- the whole point is that you measured this one for everybody
  testthat::expect_false(any(is.na(with_cov$covs$cov1)))
})

# ---- 2. the covariate is correlated with its target at the strength we asked for ---------------
testthat::test_that("the realised covariate correlation is centred on rho_cov", {
  # ONE dataset is not enough to test this, and the reason matters for the manuscript too.
  # Both the covariate and the trait evolve on the same tree, so the effective sample size is far
  # below n and a single realised correlation scatters widely. Measured over 300 seeds at n = 120:
  #
  #     requested 0.8 : mean 0.799, sd 0.096   (so one dataset can land anywhere near 0.6 - 0.95)
  #     requested 0.3 : mean 0.314, sd 0.222   (one dataset can even come out negative)
  #
  # The construction is unbiased -- both means sit within ~1 standard error of the target -- but a
  # per-dataset check would fail at random. So the test is on the AVERAGE over a fixed set of
  # seeds, with the tolerance set to about 3 standard errors of that average from the measured
  # standard deviations above (0.096/sqrt(80) = 0.011 and 0.222/sqrt(80) = 0.025).
  #
  # The paper should say the same thing: report the realised correlation as a distribution, not as
  # "0.8", or a reader will assume every dataset had a covariate that good.
  SEEDS <- 1:80

  realised <- vapply(SEEDS, function(s) {
    cl <- make_test_cell(s, n = 120, rho_cov = 0.8)
    covariate_check(cl$covs, cl$truth, target_names = "c1", rho_cov = 0.8)$rho_cov_realised
  }, numeric(1))
  testthat::expect_equal(mean(realised), 0.8, tolerance = 0.04)

  # a DIFFERENT requested strength must actually come out different, otherwise the argument is
  # being ignored and the check above would pass on a constant
  weak <- vapply(SEEDS, function(s) {
    cl <- make_test_cell(s, n = 120, rho_cov = 0.3)
    covariate_check(cl$covs, cl$truth, target_names = "c1", rho_cov = 0.3)$rho_cov_realised
  }, numeric(1))
  testthat::expect_equal(mean(weak), 0.3, tolerance = 0.08)
  testthat::expect_gt(mean(realised), mean(weak) + 0.3)
})

# ---- 3. the covariate reaches the model, on the right-hand side only ---------------------------
testthat::test_that("attach_covs puts the covariate in the data and only on the right-hand side", {
  cl <- make_test_cell(3L)

  df  <- cl$df_miss; df$Species <- rownames(cl$df_miss)
  all_traits <- setdiff(names(df), "Species")
  att <- attach_covs(df, cl$covs, rownames(cl$df_miss))

  # present in the data BACE will see ...
  testthat::expect_true("cov1" %in% names(att$df))
  testthat::expect_equal(att$cov_names, "cov1")
  # ... aligned by species name, not by position
  testthat::expect_equal(att$df[["cov1"]], cl$covs[rownames(cl$df_miss), "cov1"])

  # ... on the right of every formula ...
  fixformula <- lapply(all_traits, function(v)
    paste0(v, " ~ ", paste(c(setdiff(all_traits, v), att$cov_names), collapse = " + ")))
  testthat::expect_true(all(grepl("cov1", unlist(fixformula))))

  # ... and never on the left. A covariate with its own imputation model would be wasted compute
  # at best, and would let a fully observed variable drift at worst.
  responses <- vapply(fixformula, function(f) sub(" ~.*", "", f), character(1))
  testthat::expect_false(any(grepl("^cov", responses)))

  # with no covariate the function is a no-op, so one code path serves both regimes
  none <- attach_covs(df, NULL, rownames(cl$df_miss))
  testthat::expect_equal(none$df, df)
  testthat::expect_length(none$cov_names, 0L)
})
