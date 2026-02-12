# Tests for bace wrapper function and related components
# Tests for bace_final_imp, pool_posteriors, and bace wrapper

# Setup test data ----
setup_test_data <- function() {
  set.seed(12345)
  
  # Create a simple phylogenetic tree with 20 tips
  phylo <- ape::rtree(20)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  
  # Make ultrametric by simple rescaling (sufficient for tests)
  # This avoids needing phytools dependency
  phylo$edge.length <- phylo$edge.length / max(ape::node.depth.edgelength(phylo))
  
  # Create test data with multiple variable types
  data <- data.frame(
    y = rpois(20, lambda = 5),  # Response variable (count)
    x1 = rnorm(20, 10, 2),       # Continuous predictor
    x2 = factor(rep(c("A", "B"), 10)),  # Categorical predictor
    x3 = rnorm(20, 5, 1),        # Another continuous
    Species = phylo$tip.label
  )
  
  # Introduce missing data
  data$y[c(1, 5, 10, 15)] <- NA
  data$x1[c(2, 7, 12)] <- NA
  data$x3[c(3, 8, 13, 18)] <- NA
  
  list(phylo = phylo, data = data)
}


# Test bace_final_imp ----
test_that("bace_final_imp runs successfully with valid bace object", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run initial bace_imp
  initial <- bace_imp(
    fixformula = list("y ~ x1 + x3", "x1 ~ x3", "x3 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,  # Small number for testing
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Run final imputation
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = list("y ~ x1 + x3", "x1 ~ x3", "x3 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 5,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Check output structure
  expect_s3_class(final, "bace_final")
  expect_true("all_models" %in% names(final))
  expect_true("all_datasets" %in% names(final))
  expect_equal(length(final$all_models), 5)
  expect_equal(length(final$all_datasets), 5)
  
  # Check that models were saved
  expect_true("y" %in% names(final$all_models[[1]]))
  expect_true("x1" %in% names(final$all_models[[1]]))
  expect_true("x3" %in% names(final$all_models[[1]]))
})


test_that("bace_final_imp throws error with non-bace object", {
  expect_error(
    bace_final_imp(
      bace_object = list(data = "not a bace object"),
      fixformula = "y ~ x",
      ran_phylo_form = "~ 1 | Species",
      phylo = ape::rtree(10)
    ),
    "must be an object of class 'bace'"
  )
})


test_that("bace_final_imp works with single formula", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Only y has missing data for this test
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  # Run initial bace_imp
  initial <- bace_imp(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Run final imputation
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 3,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  expect_s3_class(final, "bace_final")
  expect_equal(length(final$all_models), 3)
})


# Test pool_posteriors ----
test_that("pool_posteriors pools multiple imputation models correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run initial bace_imp
  initial <- bace_imp(
    fixformula = list("y ~ x1 + x3", "x1 ~ x3", "x3 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Run final imputation
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = list("y ~ x1 + x3", "x1 ~ x3", "x3 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 5,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Pool posteriors
  pooled <- pool_posteriors(final)
  
  # Check output structure
  expect_s3_class(pooled, "bace_pooled")
  expect_true("models" %in% names(pooled))
  expect_true("n_imputations" %in% names(pooled))
  expect_equal(pooled$n_imputations, 5)
  
  # Check pooled models are proper MCMCglmm objects
  expect_true("y" %in% names(pooled$models))
  expect_true("x1" %in% names(pooled$models))
  expect_true("x3" %in% names(pooled$models))
  
  # Check that pooled models have MCMCglmm class
  y_model <- pooled$models$y
  expect_s3_class(y_model, "MCMCglmm")
  expect_s3_class(y_model, "bace_pooled_MCMCglmm")
  expect_true("Sol" %in% names(y_model))
  expect_true("VCV" %in% names(y_model))
  expect_true("BACE_pooling" %in% names(y_model))
  
  # Check that total samples equals samples per imputation * n_imputations
  n_samples_per <- nrow(final$all_models[[1]]$y$Sol)
  expect_equal(nrow(y_model$Sol), n_samples_per * 5)
  
  # Check BACE_pooling metadata
  expect_equal(y_model$BACE_pooling$n_imputations, 5)
  expect_equal(y_model$BACE_pooling$n_samples_per_imputation, n_samples_per)
})


test_that("pool_posteriors works with specific variable selection", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run initial and final
  initial <- bace_imp(
    fixformula = list("y ~ x1 + x3", "x1 ~ x3", "x3 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = list("y ~ x1 + x3", "x1 ~ x3", "x3 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 3,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Pool only one variable
  pooled <- pool_posteriors(final, variable = "y")
  
  expect_s3_class(pooled, "bace_pooled")
  expect_equal(length(pooled$models), 1)
  expect_true("y" %in% names(pooled$models))
  expect_false("x1" %in% names(pooled$models))
  
  # Check MCMCglmm class
  expect_s3_class(pooled$models$y, "MCMCglmm")
})


test_that("pool_posteriors throws error with non-bace_final object", {
  expect_error(
    pool_posteriors(list(all_models = "not valid")),
    "must be an object of class 'bace_final'"
  )
})


test_that("pool_posteriors throws error with invalid variable name", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Only y has missing data
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  initial <- bace_imp(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  expect_error(
    pool_posteriors(final, variable = "nonexistent"),
    "not found in models"
  )
})


# Test bace wrapper function ----
test_that("bace wrapper runs complete workflow successfully", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run complete bace workflow
  # Use small values for speed
  result <- bace(
    fixformula = list("y ~ x1 + x3", "x1 ~ x3", "x3 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    n_final = 3,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE,
    plot = FALSE
  )
  
  # Check output structure
  expect_s3_class(result, "bace_complete")
  expect_true("pooled_models" %in% names(result))
  expect_true("final_results" %in% names(result))
  expect_true("initial_results" %in% names(result))
  expect_true("convergence" %in% names(result))
  expect_true("converged" %in% names(result))
  
  # Check that we have pooled models
  expect_s3_class(result$pooled_models, "bace_pooled")
  expect_s3_class(result$pooled_models$models[[1]], "MCMCglmm")
  expect_true(length(result$pooled_models$models) > 0)
  
  # Check convergence tracking
  expect_type(result$converged, "logical")
  expect_type(result$n_attempts, "double")
})


test_that("bace wrapper handles non-convergence appropriately", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run with very few iterations to likely not converge
  # Use minimum 3 runs for convergence assessment
  result <- bace(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,  # Minimum for convergence assessment
    n_final = 2,
    nitt = 300,
    thin = 2,
    burnin = 50,
    verbose = FALSE,
    plot = FALSE,
    max_attempts = 2  # Limit attempts
  )
  
  # Should still complete even without convergence
  expect_s3_class(result, "bace_complete")
  expect_true("pooled_models" %in% names(result))
})


test_that("bace wrapper works with single formula", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Only y has missing data
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  result <- bace(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    n_final = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE,
    plot = FALSE
  )
  
  expect_s3_class(result, "bace_complete")
  expect_true("y" %in% names(result$pooled_models$models))
  
  # Verify MCMCglmm class
  expect_s3_class(result$pooled_models$models$y, "MCMCglmm")
})


test_that("bace wrapper respects verbose and plot parameters", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Only y has missing data
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  # Test verbose = FALSE (should run silently)
  expect_silent({
    result <- bace(
      fixformula = "y ~ x1 + x3",
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 3,  # Minimum for convergence assessment
      n_final = 2,
      nitt = 500,
      thin = 2,
      burnin = 100,
      verbose = FALSE,
      plot = FALSE
    )
  })
  
  expect_s3_class(result, "bace_complete")
})


# Test print methods ----
test_that("print.bace_final works correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Only y has missing data
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  initial <- bace_imp(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 3,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Test that print doesn't error
  expect_output(print(final), "BACE Final Imputation Results")
})


test_that("print.bace_pooled works correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  initial <- bace_imp(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  pooled <- pool_posteriors(final)
  
  # Test that print doesn't error
  expect_output(print(pooled), "BACE Pooled Posterior Results")
})


test_that("print.bace_pooled_MCMCglmm works correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  initial <- bace_imp(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  pooled <- pool_posteriors(final)
  
  # Print should call MCMCglmm print method
  expect_output(print(pooled$models$y), "BACE Pooled MCMCglmm Model")
})


test_that("summary.bace_pooled_MCMCglmm works correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  initial <- bace_imp(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  pooled <- pool_posteriors(final)
  
  # Summary should print header and call MCMCglmm summary method
  output <- capture.output(summ <- summary(pooled$models$y))
  # Check for the header text (allowing for spaces and formatting)
  expect_true(any(grepl("Pooled MCMCglmm Summary", output)))
})


test_that("print.bace_complete works correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  result <- bace(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,  # Minimum for convergence assessment
    n_final = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE,
    plot = FALSE
  )
  
  expect_output(print(result), "BACE Complete Analysis Results")
})


test_that("summary.bace_pooled works correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$x1[is.na(data$x1)] <- mean(data$x1, na.rm = TRUE)
  data$x3[is.na(data$x3)] <- mean(data$x3, na.rm = TRUE)
  
  initial <- bace_imp(
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1 + x3",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 2,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  pooled <- pool_posteriors(final)
  
  # Get summary - returns standard MCMCglmm summary
  output <- capture.output(summ <- summary(pooled$models$y))
  
  # Should be MCMCglmm summary class
  expect_s3_class(summ, "summary.MCMCglmm")
  
  # Should contain standard MCMCglmm summary components
  expect_true("solutions" %in% names(summ))
  expect_true("Gcovariances" %in% names(summ))
})
