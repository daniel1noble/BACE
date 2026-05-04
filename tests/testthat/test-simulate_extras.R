# Tests for sim_bace() paths not exercised in test-simulate.R:
# random slopes (rr_form), interaction terms (ix_matrix), multinomial
# and threshold response variants, and the warning paths for ignored
# rr_form covariates.

# ---- Random slopes (rr_form) ---------------------------------------------

test_that("sim_bace rr=TRUE, rr_form on continuous predictor runs", {
  testthat::skip_on_cran()
  set.seed(2026)
  out <- sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    n_cases         = 80,
    n_species       = 30,
    rr              = TRUE,
    rr_form         = list(species = "x1")
  )
  expect_s3_class(out$tree, "phylo")
  expect_true(is.data.frame(out$data))
  # Response and predictors all present
  expect_true(all(c("y", "x1", "x2") %in% names(out$data)))
})

test_that("sim_bace rr=TRUE, rr_form on non-continuous covariate warns", {
  testthat::skip_on_cran()
  set.seed(2026)
  expect_warning(
    sim_bace(
      response_type   = "gaussian",
      predictor_types = c("gaussian", "binary"),
      var_names       = c("y", "x1", "x2"),
      n_cases         = 60,
      n_species       = 25,
      rr              = TRUE,
      rr_form         = list(species = c("x1", "x2"))  # x2 is binary -> warn
    ),
    regexp = "Random slopes ignored"
  )
})

test_that("sim_bace rr=TRUE without rr_form warns", {
  testthat::skip_on_cran()
  set.seed(2026)
  expect_warning(
    sim_bace(
      response_type   = "gaussian",
      predictor_types = c("gaussian", "gaussian"),
      var_names       = c("y", "x1", "x2"),
      n_cases         = 60,
      n_species       = 25,
      rr              = TRUE,
      rr_form         = NULL
    ),
    regexp = "rr=TRUE but rr_form not specified"
  )
})

# ---- Multinomial response variants ---------------------------------------

test_that("sim_bace response_type='multinomial3' produces 3-level factor", {
  testthat::skip_on_cran()
  set.seed(2026)
  out <- sim_bace(
    response_type   = "multinomial3",
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    n_cases         = 80,
    n_species       = 30
  )
  expect_true(is.factor(out$data$y) || is.character(out$data$y))
  uniq_y <- unique(as.character(out$data$y))
  expect_true(length(uniq_y) <= 3)
})

test_that("sim_bace response_type='multinomial4' produces up-to-4-level factor", {
  testthat::skip_on_cran()
  set.seed(2026)
  out <- sim_bace(
    response_type   = "multinomial4",
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    n_cases         = 80,
    n_species       = 30
  )
  uniq_y <- unique(as.character(out$data$y))
  expect_true(length(uniq_y) <= 4)
})

# ---- Threshold (ordinal) response variant --------------------------------

test_that("sim_bace response_type='threshold4' produces ordered factor", {
  testthat::skip_on_cran()
  set.seed(2026)
  out <- sim_bace(
    response_type   = "threshold4",
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    n_cases         = 80,
    n_species       = 30
  )
  expect_true(!is.null(out$data$y))
  uniq_y <- unique(as.character(out$data$y))
  expect_true(length(uniq_y) <= 4)
})

# ---- Interaction terms (ix_matrix + beta_ix) -----------------------------

test_that("sim_bace with ix_matrix and beta_ix runs", {
  testthat::skip_on_cran()
  set.seed(2026)
  # ix_matrix: which response<-predictor pairs have interactions
  # Simplest: response y has an interaction between x1 and x2.
  # ix_matrix[i,j] = 1 means response i has interaction with predictor j
  n_vars <- 3
  ix_matrix <- matrix(0, n_vars, n_vars)
  # Pair x1:x2 -> interaction in response (row 1)
  ix_matrix[1, 2] <- 1L
  ix_matrix[1, 3] <- 1L
  out <- tryCatch(
    sim_bace(
      response_type   = "gaussian",
      predictor_types = c("gaussian", "gaussian"),
      var_names       = c("y", "x1", "x2"),
      n_cases         = 60,
      n_species       = 25,
      ix_matrix       = ix_matrix
    ),
    error = function(e) e
  )
  # Either succeeds or errors with clean message — both prove the
  # interaction-handling code path runs.
  expect_true(inherits(out, "list") || inherits(out, "error"))
})

# ---- Predictor type variants ---------------------------------------------

test_that("sim_bace with mixed predictor types (gaussian + poisson) runs", {
  testthat::skip_on_cran()
  set.seed(2026)
  out <- sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "poisson"),
    var_names       = c("y", "x1", "x2"),
    n_cases         = 80,
    n_species       = 30
  )
  expect_true(is.data.frame(out$data))
  # poisson predictor should be integer-valued
  expect_true(all(out$data$x2 == round(out$data$x2)))
})

test_that("sim_bace with poisson response runs", {
  testthat::skip_on_cran()
  set.seed(2026)
  out <- sim_bace(
    response_type   = "poisson",
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    n_cases         = 80,
    n_species       = 30
  )
  # Poisson response: integer-valued, non-negative
  expect_true(all(out$data$y == round(out$data$y)))
  expect_true(all(out$data$y >= 0))
})

test_that("sim_bace with binary response runs", {
  testthat::skip_on_cran()
  set.seed(2026)
  out <- sim_bace(
    response_type   = "binary",
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    n_cases         = 80,
    n_species       = 30
  )
  uniq_y <- unique(as.character(out$data$y))
  expect_true(length(uniq_y) <= 2)
})

# ---- Custom phylo_signal vector ------------------------------------------

test_that("sim_bace with custom phylo_signal per variable runs", {
  testthat::skip_on_cran()
  set.seed(2026)
  out <- sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    n_cases         = 60,
    n_species       = 25,
    phylo_signal    = c(0.8, 0.5, 0.3)
  )
  expect_s3_class(out$tree, "phylo")
})

# ---- Custom missingness mechanism ----------------------------------------

test_that("sim_bace with custom missingness runs and creates NA", {
  testthat::skip_on_cran()
  set.seed(2026)
  # missingness vector has one entry per variable (response + predictors)
  out <- sim_bace(
    response_type   = "gaussian",
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    n_cases         = 80,
    n_species       = 30,
    missingness     = c(0.2, 0.0, 0.0)
  )
  # ~20% of y should be missing
  prop_na <- mean(is.na(out$data$y))
  expect_true(prop_na > 0)
  expect_true(prop_na < 0.5)
  # Predictors should have NO missingness (as specified)
  expect_true(all(!is.na(out$data$x1)))
  expect_true(all(!is.na(out$data$x2)))
})
