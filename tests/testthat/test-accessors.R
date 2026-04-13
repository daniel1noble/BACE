# Tests for get_pooled_model() and get_imputed_data() accessor functions

# ---- Helper: reusable test fixtures ------------------------------------
# Uses the same setup pattern as test-bace_wrapper.R
setup_test_data <- function() {
  set.seed(12345)

  phylo <- ape::rtree(20)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length / max(ape::node.depth.edgelength(phylo))

  data <- data.frame(
    y       = rpois(20, lambda = 5),
    x1      = rnorm(20, 10, 2),
    x2      = factor(rep(c("A", "B"), 10)),
    x3      = rnorm(20, 5, 1),
    Species = phylo$tip.label
  )
  data$y[c(1, 5, 10, 15)]  <- NA
  data$x1[c(2, 7, 12)]     <- NA
  data$x3[c(3, 8, 13, 18)] <- NA

  list(phylo = phylo, data = data)
}

# Build a bace_final and bace_complete fixture (expensive; cached per file)
build_bace_objects <- function() {
  ts <- setup_test_data()

  initial <- bace_imp(
    fixformula    = list("y ~ x1 + x3", "x1 ~ x3", "x3 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo  = ts$phylo,
    data   = ts$data,
    runs   = 3,
    nitt   = 500,
    thin   = 2,
    burnin = 100,
    verbose = FALSE
  )

  final <- bace_final_imp(
    bace_object    = initial,
    fixformula     = list("y ~ x1 + x3", "x1 ~ x3", "x3 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo  = ts$phylo,
    n_final = 3,
    nitt    = 500,
    thin    = 2,
    burnin  = 100,
    verbose = FALSE
  )

  pooled <- pool_posteriors(final)

  # Assemble a minimal bace_complete object (same structure as bace())
  bace_complete <- list(
    pooled_models    = pooled,
    final_results    = final,
    imputed_datasets = final$all_datasets,
    initial_results  = initial,
    convergence      = NULL,
    converged        = TRUE,
    n_attempts       = 1L,
    call             = match.call()
  )
  class(bace_complete) <- "bace_complete"

  list(
    initial       = initial,
    final         = final,
    pooled        = pooled,
    bace_complete = bace_complete,
    test_data     = ts
  )
}

# ========================================================================
# Tests for get_pooled_model()
# ========================================================================

test_that("get_pooled_model extracts a single variable from bace_complete", {
  skip_on_cran()
  objs <- build_bace_objects()

  y_model <- get_pooled_model(objs$bace_complete, variable = "y")
  expect_s3_class(y_model, "MCMCglmm")
  expect_s3_class(y_model, "bace_pooled_MCMCglmm")
  expect_true("Sol"  %in% names(y_model))
  expect_true("VCV"  %in% names(y_model))
  expect_true("BACE_pooling" %in% names(y_model))
})

test_that("get_pooled_model extracts a single variable from bace_pooled", {
  skip_on_cran()
  objs <- build_bace_objects()

  x1_model <- get_pooled_model(objs$pooled, variable = "x1")
  expect_s3_class(x1_model, "MCMCglmm")
  expect_true("Sol" %in% names(x1_model))
})

test_that("get_pooled_model returns all models when variable is NULL", {
  skip_on_cran()
  objs <- build_bace_objects()

  all_models <- get_pooled_model(objs$bace_complete)
  expect_true(is.list(all_models))
  expect_true(all(c("y", "x1", "x3") %in% names(all_models)))
  expect_s3_class(all_models$y, "MCMCglmm")
})

test_that("get_pooled_model errors on invalid variable name", {
  skip_on_cran()
  objs <- build_bace_objects()

  expect_error(
    get_pooled_model(objs$bace_complete, variable = "nonexistent"),
    "not found"
  )
})

test_that("get_pooled_model errors on non-character variable", {
  skip_on_cran()
  objs <- build_bace_objects()

  expect_error(
    get_pooled_model(objs$bace_complete, variable = 1),
    "single character string"
  )
  expect_error(
    get_pooled_model(objs$bace_complete, variable = c("y", "x1")),
    "single character string"
  )
})

test_that("get_pooled_model errors on wrong object class", {
  expect_error(
    get_pooled_model(list(a = 1), variable = "y"),
    "bace_complete.*bace_pooled"
  )
  expect_error(
    get_pooled_model(data.frame(x = 1), variable = "y"),
    "bace_complete.*bace_pooled"
  )
})

# ========================================================================
# Tests for get_imputed_data()
# ========================================================================

test_that("get_imputed_data returns a list from bace_complete", {
  skip_on_cran()
  objs <- build_bace_objects()

  imp <- get_imputed_data(objs$bace_complete, format = "list")
  expect_true(is.list(imp))
  expect_equal(length(imp), 3)  # n_final = 3
  expect_true(is.data.frame(imp[[1]]))
})

test_that("get_imputed_data returns a list from bace_final", {
  skip_on_cran()
  objs <- build_bace_objects()

  imp <- get_imputed_data(objs$final, format = "list")
  expect_true(is.list(imp))
  expect_equal(length(imp), 3)
})

test_that("get_imputed_data default format is list", {
  skip_on_cran()
  objs <- build_bace_objects()

  imp <- get_imputed_data(objs$bace_complete)
  expect_true(is.list(imp))
  expect_false(is.data.frame(imp))
})

test_that("get_imputed_data data.frame format stacks with .imputation column", {
  skip_on_cran()
  objs <- build_bace_objects()

  imp_df <- get_imputed_data(objs$bace_complete, format = "data.frame")
  expect_true(is.data.frame(imp_df))
  expect_true(".imputation" %in% names(imp_df))

  # The .imputation column should have values 1..n_final
  expect_equal(sort(unique(imp_df$.imputation)), 1:3)

  # Total rows = n_final * rows_per_dataset
  n_rows_each <- nrow(objs$bace_complete$imputed_datasets[[1]])
  expect_equal(nrow(imp_df), 3 * n_rows_each)
})

test_that("get_imputed_data data.frame preserves column names", {
  skip_on_cran()
  objs <- build_bace_objects()

  original_cols <- names(objs$bace_complete$imputed_datasets[[1]])
  imp_df <- get_imputed_data(objs$bace_complete, format = "data.frame")
  expect_true(all(original_cols %in% names(imp_df)))
})

test_that("get_imputed_data errors on wrong object class", {
  expect_error(
    get_imputed_data(list(a = 1)),
    "bace_complete.*bace_final"
  )
  expect_error(
    get_imputed_data(data.frame(x = 1)),
    "bace_complete.*bace_final"
  )
})

test_that("get_imputed_data errors on invalid format argument", {
  skip_on_cran()
  objs <- build_bace_objects()

  expect_error(
    get_imputed_data(objs$bace_complete, format = "matrix"),
    "should be one of"
  )
})
