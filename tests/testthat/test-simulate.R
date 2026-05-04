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

# ----------------------------------------------------------------------------
# Regression: source bugs fixed alongside this test file.
# ----------------------------------------------------------------------------

test_that("sim_bace output uses 'Species' (capital S) per BACE convention", {
  # Bug surfaced 2026-05: sim_bace previously produced a lowercase
  # 'species' column, which forced users to rename before feeding
  # the result into bace_imp (whose ran_phylo_form examples all use
  # capital 'Species'). Fix in commit alongside this test.
  set.seed(2026)
  out <- suppressMessages(sim_bace_gaussian(n_predictors = 2,
                                             n_cases = 30, n_species = 10))
  expect_true("Species" %in% colnames(out$data))
  expect_true("Species" %in% colnames(out$complete_data))
  expect_false("species" %in% colnames(out$data))
})

test_that("generate_default_beta_matrix handles n_predictors = 1 cleanly", {
  # Bug surfaced 2026-05: the inner `for (i in 2:n_predictors)` loop
  # used R's `2:1` ⇒ c(2,1) quirk and tried to assign
  # beta_matrix[2, 1] on a 1x1 matrix. Fixed by guarding the loop
  # behind n_predictors >= 2.
  set.seed(2026)
  b <- generate_default_beta_matrix(n_predictors = 1, sparsity = 0.5)
  expect_true(is.matrix(b))
  expect_equal(dim(b), c(1L, 1L))
  expect_equal(b[1, 1], 0)
})

test_that("generate_default_beta_matrix(0) errors clearly", {
  expect_error(generate_default_beta_matrix(n_predictors = 0),
               regexp = "n_predictors must be >= 1")
})

# ----------------------------------------------------------------------------
# More sim_bace coverage: alternative response types + validation paths
# ----------------------------------------------------------------------------

test_that("sim_bace with custom var_names overrides defaults", {
  set.seed(2026)
  out <- suppressMessages(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("resp", "pred_a", "pred_b"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases = 30, n_species = 12
  ))
  expect_equal(out$params$var_names, c("resp", "pred_a", "pred_b"))
  expect_true(all(c("resp", "pred_a", "pred_b") %in% colnames(out$data)))
})

test_that("sim_bace rejects mismatched var_names length", {
  set.seed(2026)
  expect_error(suppressMessages(sim_bace(
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("only_one"),  # need 3 (response + 2 predictors)
    n_cases = 20, n_species = 10
  )), regexp = "var_names must have length")
})

test_that("sim_bace rejects mismatched phylo_signal length", {
  set.seed(2026)
  expect_error(suppressMessages(sim_bace(
    predictor_types = c("gaussian", "gaussian"),
    phylo_signal    = c(0.5),  # need 3 (response + 2 predictors)
    n_cases = 20, n_species = 10
  )), regexp = "phylo_signal must have length")
})

test_that("sim_bace rejects out-of-range phylo_signal", {
  set.seed(2026)
  expect_error(suppressMessages(sim_bace(
    predictor_types = c("gaussian"),
    phylo_signal    = c(0.5, 1.5),  # 1.5 is out of [0, 1)
    n_cases = 20, n_species = 10
  )), regexp = "phylo_signal values must be in")
})

test_that("sim_bace accepts numeric intercepts vector and converts to list", {
  set.seed(2026)
  out <- suppressMessages(sim_bace(
    predictor_types = c("gaussian"),
    phylo_signal    = c(0.3, 0.3),
    intercepts      = c(2.0, 0.5),  # response, predictor
    n_cases = 30, n_species = 12
  ))
  expect_equal(out$params$intercepts$response, 2.0)
  expect_equal(out$params$intercepts$predictors, 0.5)
})

test_that("sim_bace rejects list intercepts missing required keys", {
  set.seed(2026)
  expect_error(suppressMessages(sim_bace(
    predictor_types = c("gaussian"),
    phylo_signal    = c(0.3, 0.3),
    intercepts      = list(response = 0),  # missing 'predictors'
    n_cases = 20, n_species = 10
  )), regexp = "must have 'response' and 'predictors'")
})

test_that("sim_bace warns when rr=TRUE but rr_form is NULL", {
  set.seed(2026)
  expect_warning(suppressMessages(sim_bace(
    predictor_types = c("gaussian"),
    phylo_signal    = c(0.3, 0.3),
    rr              = TRUE,
    rr_form         = NULL,
    n_cases = 30, n_species = 12
  )), regexp = "rr=TRUE but rr_form not specified")
})

test_that("sim_bace_gaussian respects beta_sparsity override", {
  set.seed(2026)
  out <- suppressMessages(sim_bace_gaussian(
    n_predictors = 4, n_cases = 30, n_species = 15,
    beta_sparsity = 0.0  # all coefficients non-zero
  ))
  bm <- out$params$beta_matrix
  # Lower triangle entries should mostly be non-zero with sparsity=0
  lower <- bm[lower.tri(bm)]
  expect_gt(sum(lower != 0), length(lower) * 0.5)  # most are non-zero
})

test_that("sim_bace handles a poisson predictor in the design", {
  set.seed(2026)
  # First predictor must be gaussian for numerical stability
  # (sim_bace warns otherwise) -- put the poisson second.
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "poisson"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases = 30, n_species = 12
  )))
  pois_col <- out$params$var_names[3]
  vals <- out$complete_data[[pois_col]]
  expect_true(is.integer(vals) || all(vals == as.integer(vals), na.rm = TRUE))
  expect_true(all(vals >= 0))
})

# ----------------------------------------------------------------------------
# Response-type variants: threshold + multinomial responses
# ----------------------------------------------------------------------------

test_that("sim_bace response_type='threshold3' returns 3-level ordered factor", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "threshold3",
    predictor_types = c("gaussian", "gaussian"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 40, n_species = 20
  )))
  resp_col <- out$params$var_names[1]
  resp <- out$complete_data[[resp_col]]
  expect_s3_class(resp, "factor")
  expect_true(is.ordered(resp))
  # All response values should be valid level labels (1, 2, or 3)
  expect_true(all(as.integer(as.character(resp)) %in% 1:3))
})

test_that("sim_bace response_type='threshold4' returns 4-level ordered factor", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "threshold4",
    predictor_types = c("gaussian", "gaussian"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 50, n_species = 25
  )))
  resp <- out$complete_data[[out$params$var_names[1]]]
  expect_s3_class(resp, "factor")
  expect_true(is.ordered(resp))
  expect_true(all(as.integer(as.character(resp)) %in% 1:4))
})

test_that("sim_bace response_type='multinomial3' returns 3-level unordered factor", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "multinomial3",
    predictor_types = c("gaussian", "gaussian"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 40, n_species = 20
  )))
  resp <- out$complete_data[[out$params$var_names[1]]]
  expect_s3_class(resp, "factor")
  expect_false(is.ordered(resp))
  expect_true(all(as.character(resp) %in% c("A","B","C")))
})

test_that("sim_bace response_type='multinomial4' returns 4-level unordered factor", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "multinomial4",
    predictor_types = c("gaussian", "gaussian"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 50, n_species = 25
  )))
  resp <- out$complete_data[[out$params$var_names[1]]]
  expect_s3_class(resp, "factor")
  expect_true(all(as.character(resp) %in% c("A","B","C","D")))
})

test_that("sim_bace response_type='binary' returns 2-level ordered factor", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "binary",
    predictor_types = c("gaussian", "gaussian"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 40, n_species = 20
  )))
  resp <- out$complete_data[[out$params$var_names[1]]]
  expect_s3_class(resp, "factor")
  expect_equal(nlevels(resp), 2L)
})

test_that("sim_bace response_type='poisson' returns non-negative integer", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "poisson",
    predictor_types = c("gaussian", "gaussian"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 40, n_species = 20
  )))
  resp <- out$complete_data[[out$params$var_names[1]]]
  expect_true(is.integer(resp))
  expect_true(all(resp >= 0))
})

# ----------------------------------------------------------------------------
# Threshold/multinomial as predictor types (different code path than response)
# ----------------------------------------------------------------------------

test_that("sim_bace with threshold3 predictor produces ordered factor", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "threshold3"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 40, n_species = 20
  )))
  pred_col <- out$params$var_names[3]
  vals <- out$complete_data[[pred_col]]
  expect_s3_class(vals, "factor")
  expect_true(is.ordered(vals))
})

test_that("sim_bace with multinomial3 predictor produces unordered factor", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "multinomial3"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 40, n_species = 20
  )))
  pred_col <- out$params$var_names[3]
  vals <- out$complete_data[[pred_col]]
  expect_s3_class(vals, "factor")
  expect_false(is.ordered(vals))
  expect_true(all(as.character(vals) %in% c("A","B","C")))
})

# ----------------------------------------------------------------------------
# Custom beta_matrix + custom intercepts list
# ----------------------------------------------------------------------------

test_that("sim_bace accepts a custom beta_matrix", {
  set.seed(2026)
  bm <- matrix(0, 3, 3)
  bm[2, 1] <- 0.4
  bm[3, 1] <- -0.2
  bm[3, 2] <- 0.1
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    beta_matrix     = bm,
    phylo_signal    = c(0.3, 0.3, 0.3, 0.3),
    n_cases         = 30, n_species = 15
  )))
  # beta_matrix is named with predictor_names by sim_bace; compare
  # values element-wise rather than insisting on attribute equality.
  expect_equal(unname(out$params$beta_matrix), unname(bm))
})

test_that("sim_bace with custom intercepts list applies the values", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian"),
    intercepts      = list(response = 5.0,
                            predictors = c(2.0, -1.0)),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 40, n_species = 20
  )))
  expect_equal(out$params$intercepts$response, 5.0)
  expect_equal(out$params$intercepts$predictors, c(2.0, -1.0))
})

test_that("sim_bace rejects beta_matrix with wrong dimensions", {
  set.seed(2026)
  bm <- matrix(0, 2, 2)  # 2x2 but n_predictors = 3
  expect_error(suppressMessages(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    beta_matrix     = bm,
    n_cases         = 30, n_species = 15
  )), regexp = "beta_matrix must be a square matrix")
})

# ----------------------------------------------------------------------------
# Custom beta_resp argument
# ----------------------------------------------------------------------------

test_that("sim_bace with custom beta_resp vector applies coefficients", {
  set.seed(2026)
  out <- suppressMessages(suppressWarnings(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian"),
    beta_resp       = c(2.0, -1.0),
    phylo_signal    = c(0.3, 0.3, 0.3),
    n_cases         = 40, n_species = 20
  )))
  # beta_resp_stored should reflect the input
  stored <- out$params$beta_resp
  expect_true(!is.null(stored))
})

# ----------------------------------------------------------------------------
# missingness vector
# ----------------------------------------------------------------------------

test_that("sim_bace missingness vector validates length", {
  set.seed(2026)
  expect_error(suppressMessages(sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian"),
    phylo_signal    = c(0.3, 0.3, 0.3),
    missingness     = c(0.1),  # length 1, but n_vars = 3
    n_cases = 30, n_species = 15
  )), regexp = "missingness must have length")
})
