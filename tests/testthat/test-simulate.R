# Tests for the simulation framework (R/simulate_simBACE.R + R/simulate_auxiliary.r).
# These exercise the public API at small scale -- enough to drive coverage
# through every documented response type and the helper utilities.

# ----------------------------------------------------------------------------
# sim_bace + convenience wrappers
# ----------------------------------------------------------------------------

test_that("sim_bace_gaussian returns a list with expected structure", {
  set.seed(2026)
  out <- suppressMessages(sim_bace_gaussian(
    n_predictors = 3, n_cases = 60, n_species = 20,
    phylo_signal = 0.5, beta_sparsity = 0.5
  ))
  expect_type(out, "list")
  expect_named(out, c("data", "complete_data", "tree", "params",
                       "random_effects"))
  expect_s3_class(out$tree, "phylo")
  expect_s3_class(out$data, "data.frame")
  expect_equal(nrow(out$data), 60L)
  # response + 3 predictors + Species column
  expect_true(ncol(out$data) >= 4L)
  expect_equal(out$params$response_type, "gaussian")
  expect_length(out$params$predictor_types, 3L)
})

test_that("sim_bace_poisson produces integer-valued response", {
  set.seed(2026)
  out <- suppressMessages(sim_bace_poisson(
    n_predictors = 2, n_cases = 50, n_species = 15
  ))
  resp_col <- out$params$var_names[1]
  resp <- out$complete_data[[resp_col]]
  expect_true(is.integer(resp) || all(resp == as.integer(resp), na.rm = TRUE))
  expect_equal(out$params$response_type, "poisson")
})

test_that("sim_bace_binary produces a 2-level factor response", {
  set.seed(2026)
  out <- suppressMessages(sim_bace_binary(
    n_predictors = 2, n_cases = 50, n_species = 15
  ))
  resp_col <- out$params$var_names[1]
  resp <- out$complete_data[[resp_col]]
  expect_s3_class(resp, "factor")
  expect_equal(nlevels(resp), 2L)
})

test_that("sim_bace handles mixed predictor types", {
  set.seed(2026)
  out <- suppressMessages(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "binary", "poisson"),
    phylo_signal    = c(0.3, 0.3, 0.3, 0.3),
    n_cases         = 40, n_species = 12
  ))
  expect_named(out, c("data", "complete_data", "tree", "params",
                       "random_effects"))
  # Predictor columns should match their declared types
  for (i in seq_along(out$params$predictor_types)) {
    col <- out$params$var_names[i + 1]
    pt  <- out$params$predictor_types[i]
    val <- out$complete_data[[col]]
    if (pt == "gaussian")      expect_true(is.numeric(val))
    else if (pt == "binary")   expect_s3_class(val, "factor")
    else if (pt == "poisson")  expect_true(is.integer(val) ||
                                            all(val == as.integer(val),
                                                na.rm = TRUE))
  }
})

test_that("sim_bace applies missingness when requested", {
  set.seed(2026)
  out <- suppressMessages(sim_bace(
    response_type = "gaussian",
    predictor_types = c("gaussian", "gaussian"),
    phylo_signal  = c(0.2, 0.2, 0.2),
    n_cases = 40, n_species = 12,
    missingness = c(0.0, 0.30, 0.30)
  ))
  # Predictor columns should have ~30% NA; response should have 0%.
  preds <- out$params$var_names[-1]
  for (p in preds) {
    expect_true(sum(is.na(out$data[[p]])) > 0L)
  }
  # complete_data has no NAs
  expect_equal(sum(is.na(out$complete_data)), 0L)
})

test_that("print_sim_bace_summary runs without error", {
  set.seed(2026)
  out <- suppressMessages(sim_bace_gaussian(n_predictors = 2,
                                            n_cases = 30, n_species = 10))
  expect_output(print_sim_bace_summary(out), "Simulation Summary")
})

# ----------------------------------------------------------------------------
# Auxiliary helpers (R/simulate_auxiliary.r)
# ----------------------------------------------------------------------------

test_that("sim_tree returns a list with phylo + case mapping", {
  set.seed(2026)
  out <- sim_tree(n_species = 20, birth = 0.6, death = 0.3, n_cases = 50)
  expect_type(out, "list")
  expect_named(out, c("tree", "cases", "species"))
  expect_s3_class(out$tree, "phylo")
  # n_cases > n_species so we expect <= n_species distinct species
  expect_true(length(out$species) <= 20L)
  expect_equal(ape::Ntip(out$tree), length(out$species))
})

test_that("generate_default_beta_matrix returns a square matrix with expected sparsity", {
  set.seed(2026)
  b <- generate_default_beta_matrix(n_predictors = 5, sparsity = 0.5)
  expect_true(is.matrix(b))
  expect_equal(dim(b), c(5L, 5L))
  # Lower-triangular structure: upper triangle (incl diag) is zero
  expect_true(all(b[upper.tri(b, diag = TRUE)] == 0))
})

test_that("mnom_liab2cat samples categories from softmax", {
  set.seed(2026)
  liability <- matrix(rnorm(100 * 2), nrow = 100, ncol = 2)
  cats <- c("A", "B", "C")
  out <- mnom_liab2cat(liability, cats)
  expect_length(out, 100L)
  expect_true(all(out %in% cats))
})

test_that("ordinal_liab2cat returns integer levels in 1..K", {
  set.seed(2026)
  liab <- rnorm(100)
  out <- ordinal_liab2cat(liab, n_cats = 4)
  expect_true(is.integer(out))
  expect_length(out, 100L)
  expect_true(all(out >= 1L & out <= 4L))
})
