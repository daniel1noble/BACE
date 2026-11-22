# A factor level that no observed row has must not reach MCMCglmm's threshold sampler.
#
# THE BUG (2026-11-22). Seven datasets of the 2026-11 simulation crashed R inside MCMCglmm
# ("bus error, invalid alignment"), all at strong signal with whole clades hidden. Each had an
# ordinal level present in the factor but absent from every observed row (seed 79 of
# n100_ps0.9_good_covariate_clade: observed 0/7/24/34). Through BACE's own fitting path the
# empty level kept = crash, dropped = fit OK, imputed 9/38/53. So BACE drops it before fitting.

test_that(".drop_unobserved_levels drops empty levels, warns once, refuses a one-level response", {
  d <- data.frame(y = factor(c("a", "b", "c", NA, "b", "c", NA), levels = c("a", "b", "c", "d"),
                             ordered = TRUE), x = rnorm(7))
  # "a" is observed once, "d" never; "d" goes, the order is kept
  expect_warning(out <- BACE:::.drop_unobserved_levels(d, "y", "threshold"), "no observed value \\(d\\)")
  expect_identical(levels(out$y), c("a", "b", "c"))
  expect_true(is.ordered(out$y))
  expect_identical(as.character(out$y), as.character(d$y))
  # silent when asked
  expect_silent(BACE:::.drop_unobserved_levels(d, "y", "threshold", warn = FALSE))
  # nothing to drop: unchanged, no warning
  full <- data.frame(y = factor(c("a", "b", "b"), levels = c("a", "b")), x = 1:3)
  expect_silent(same <- BACE:::.drop_unobserved_levels(full, "y", "threshold"))
  expect_identical(same, full)
  # not a factor, or not a categorical type: untouched
  num <- data.frame(y = c(1, 2, NA), x = 1:3)
  expect_identical(BACE:::.drop_unobserved_levels(num, "y", "poisson"), num)
  # observed values in one level only: refuse, naming the empty ones
  one <- data.frame(y = factor(c("a", "a", NA, NA), levels = c("a", "b", "c")), x = 1:4)
  expect_error(BACE:::.drop_unobserved_levels(one, "y", "threshold"), "only one of its levels")
})

# The crash itself (seed 79 of n100_ps0.9_good_covariate_clade, campaign settings) is too slow
# for a test; it was reproduced by hand on 2026-11-22. These tests pin the behaviour the guard
# gives: a warning, imputations that never take the empty level, and a column that keeps all
# its levels -- in the chained phase, in the final phase, and for a categorical response.

test_that("an ordinal with an unobserved level is left out of its model: warning, never imputed, levels kept", {
  skip_on_cran()
  set.seed(5)
  sim <- sim_bace(response_type = "gaussian", predictor_types = c("gaussian", "threshold4"),
                  var_names = c("y", "x1", "ord"), beta_matrix = matrix(c(0, 0, 0.5, 0), 2, 2, byrow = TRUE),
                  beta_resp = list(x1 = 0.5, ord = 0.3), phylo_signal = c(0.5, 0.5, 0.5),
                  missingness = c(0, 0.2, 0), n_cases = 60, n_species = 60)
  d <- sim$data
  # hide EVERY row of the lowest ordinal level: the level stays in the factor with no data
  lowest <- levels(d$ord)[1]
  d$ord[d$ord == lowest] <- NA
  expect_true(lowest %in% levels(d$ord)); expect_equal(sum(d$ord == lowest, na.rm = TRUE), 0)

  expect_warning(
    fit <- suppressMessages(bace_imp(fixformula = c("ord ~ y + x1", "x1 ~ y + ord"),
                                     ran_phylo_form = "~ 1 | Species", phylo = sim$tree, data = d,
                                     nitt = 600, burnin = 100, thin = 4, runs = 2, verbose = FALSE)),
    "no observed value")
  imputed <- fit$data[[length(fit$data)]]
  # the column keeps its full set of levels, the imputations never use the empty one
  expect_identical(levels(imputed$ord), levels(d$ord))
  expect_false(any(is.na(imputed$ord)))
  expect_equal(sum(imputed$ord == lowest), 0)

  # the FINAL phase has the guard too: no final dataset imputes the empty level
  full <- suppressWarnings(suppressMessages(
    bace(fixformula = c("ord ~ y + x1", "x1 ~ y + ord"), ran_phylo_form = "~ 1 | Species",
         phylo = sim$tree, data = d, runs = 2, nitt = 600, burnin = 100, thin = 4,
         n_final = 3, skip_conv = TRUE, verbose = FALSE)))
  finals <- get_imputed_data(full, format = "list")
  expect_length(finals, 3L)
  for (f in finals) {
    expect_identical(levels(f$ord), levels(d$ord))
    expect_equal(sum(f$ord == lowest), 0)
  }
})

test_that("a categorical response with an unobserved level runs on both categorical paths", {
  skip_on_cran()
  # The multinomial prior is sized from the number of levels; it must be the number the
  # model is fitted with (the empty level dropped), or MCMCglmm refuses the prior.
  set.seed(8)
  sim <- sim_bace(response_type = "multinomial4", predictor_types = "gaussian",
                  var_names = c("cat", "x1"), beta_resp = list(x1 = 0.5),
                  phylo_signal = c(0.5, 0.5), missingness = c(0, 0), n_cases = 60, n_species = 60)
  d <- sim$data
  gone <- levels(d$cat)[3]
  d$cat[d$cat == gone] <- NA                                   # that level is now unobserved
  d$cat[sample(which(!is.na(d$cat)), 8)] <- NA                 # plus a few ordinary gaps
  expect_true(gone %in% levels(d$cat)); expect_equal(sum(d$cat == gone, na.rm = TRUE), 0)
  for (ovr in c(FALSE, TRUE)) {
    expect_warning(
      fit <- suppressMessages(bace_imp(fixformula = "cat ~ x1", ran_phylo_form = "~ 1 | Species",
                                       phylo = sim$tree, data = d, nitt = 600, burnin = 100, thin = 4,
                                       runs = 2, verbose = FALSE, ovr_categorical = ovr)),
      "no observed value")
    imputed <- fit$data[[length(fit$data)]]
    expect_identical(levels(imputed$cat), levels(d$cat))
    expect_false(any(is.na(imputed$cat)))
    expect_equal(sum(imputed$cat == gone), 0)
  }
})
