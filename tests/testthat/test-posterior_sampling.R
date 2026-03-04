# Tests for posterior sampling functionality in pool_posteriors and bace wrapper

# Setup test data ----
setup_sampling_test_data <- function() {
  set.seed(99999)
  
  # Create a simple phylogenetic tree with 20 tips
  phylo <- ape::rtree(20)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  
  # Make ultrametric by simple rescaling
  phylo$edge.length <- phylo$edge.length / max(ape::node.depth.edgelength(phylo))
  
  # Create test data
  data <- data.frame(
    y = rpois(20, lambda = 5),  # Count response
    x1 = rnorm(20, 10, 2),       # Continuous predictor
    x2 = rnorm(20, 5, 1),        # Another continuous
    Species = phylo$tip.label
  )
  
  # Introduce missing data
  data$y[c(1, 5, 10, 15)] <- NA
  data$x1[c(2, 7, 12)] <- NA
  data$x2[c(3, 8, 13)] <- NA
  
  list(phylo = phylo, data = data)
}


# Test pool_posteriors with sampling ----
test_that("pool_posteriors with sample_size reduces posterior samples correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_sampling_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run initial bace_imp
  initial <- bace_imp(
    fixformula = list("y ~ x1 + x2", "x1 ~ x2", "x2 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    nitt = 1000,
    thin = 2,
    burnin = 200,
    verbose = FALSE
  )
  
  # Run final imputation
  n_final <- 5
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = list("y ~ x1 + x2", "x1 ~ x2", "x2 ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = n_final,
    nitt = 1000,
    thin = 2,
    burnin = 200,
    verbose = FALSE
  )
  
  # Get original sample size per imputation
  original_n_samples <- nrow(final$all_models[[1]]$y$Sol)
  
  # Pool with sampling
  sample_size <- 50
  pooled <- pool_posteriors(final, sample_size = sample_size)
  
  # Check output structure
  expect_s3_class(pooled, "bace_pooled")
  expect_true("models" %in% names(pooled))
  
  # Check that pooled y model has correct number of samples
  y_model <- pooled$models$y
  expect_s3_class(y_model, "MCMCglmm")
  expect_s3_class(y_model, "bace_pooled_MCMCglmm")
  
  # Total samples should be sample_size * n_final
  expected_total <- sample_size * n_final
  expect_equal(nrow(y_model$Sol), expected_total)
  expect_equal(nrow(y_model$VCV), expected_total)
  
  # Check BACE_pooling metadata
  expect_true("BACE_pooling" %in% names(y_model))
  expect_equal(y_model$BACE_pooling$n_imputations, n_final)
  expect_equal(y_model$BACE_pooling$n_samples_per_imputation, sample_size)
  expect_equal(y_model$BACE_pooling$original_samples_per_imputation, original_n_samples)
  expect_equal(y_model$BACE_pooling$total_samples, expected_total)
  expect_true(y_model$BACE_pooling$sampled)
  expect_true(y_model$BACE_pooling$pooled)
})


test_that("pool_posteriors without sample_size uses all samples", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_sampling_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run workflow
  initial <- bace_imp(
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    nitt = 600,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  n_final <- 3
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = n_final,
    nitt = 600,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Get original sample size
  original_n_samples <- nrow(final$all_models[[1]]$y$Sol)
  
  # Pool without sampling (sample_size = NULL)
  pooled <- pool_posteriors(final, sample_size = NULL)
  
  # Check that all samples are used
  y_model <- pooled$models$y
  expect_equal(nrow(y_model$Sol), original_n_samples * n_final)
  
  # Check metadata shows no sampling was used
  expect_false(y_model$BACE_pooling$sampled)
  expect_equal(y_model$BACE_pooling$n_samples_per_imputation, original_n_samples)
  expect_equal(y_model$BACE_pooling$original_samples_per_imputation, original_n_samples)
})


test_that("pool_posteriors handles sample_size larger than available samples", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_sampling_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run with small nitt to get few samples
  initial <- bace_imp(
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    nitt = 300,
    thin = 5,
    burnin = 50,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 3,
    nitt = 300,
    thin = 5,
    burnin = 50,
    verbose = FALSE
  )
  
  # Get actual sample size (should be small due to high thin)
  original_n_samples <- nrow(final$all_models[[1]]$y$Sol)
  
  # Request more samples than available
  large_sample_size <- original_n_samples + 100
  pooled <- pool_posteriors(final, sample_size = large_sample_size)
  
  # Should use all available samples (no sampling applied)
  y_model <- pooled$models$y
  expect_equal(nrow(y_model$Sol), original_n_samples * 3)
  expect_false(y_model$BACE_pooling$sampled)
})


test_that("pool_posteriors sampling works for multiple variables", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_sampling_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run workflow with multiple variables
  initial <- bace_imp(
    fixformula = list("y ~ x1", "x1 ~ x2", "x2 ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    nitt = 800,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  n_final <- 4
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = list("y ~ x1", "x1 ~ x2", "x2 ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = n_final,
    nitt = 800,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Pool with sampling
  sample_size <- 75
  pooled <- pool_posteriors(final, sample_size = sample_size)
  
  # Check all variables are correctly sampled
  for (var in c("y", "x1", "x2")) {
    model <- pooled$models[[var]]
    expect_equal(nrow(model$Sol), sample_size * n_final,
                 info = paste("Variable:", var))
    expect_true(model$BACE_pooling$sampled,
                info = paste("Variable:", var))
    expect_equal(model$BACE_pooling$n_samples_per_imputation, sample_size,
                 info = paste("Variable:", var))
  }
})


test_that("pool_posteriors with sample_size=1 works correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_sampling_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run workflow
  initial <- bace_imp(
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  n_final <- 10
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = n_final,
    nitt = 500,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Pool with just 1 sample per imputation (edge case)
  pooled <- pool_posteriors(final, sample_size = 1)
  
  y_model <- pooled$models$y
  expect_equal(nrow(y_model$Sol), n_final)  # 1 sample * 10 imputations
  expect_true(y_model$BACE_pooling$sampled)
  expect_equal(y_model$BACE_pooling$n_samples_per_imputation, 1)
})


# Test bace wrapper with sample_size ----
test_that("bace wrapper passes sample_size correctly to pool_posteriors", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_sampling_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run bace with sample_size
  sample_size <- 60
  result <- bace(
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,  # Need at least 3 for convergence assessment
    n_final = 4,
    nitt = 800,
    thin = 2,
    burnin = 100,
    verbose = FALSE,
    plot = FALSE,
    sample_size = sample_size,
    skip_conv = TRUE  # Skip convergence retries for speed
  )
  
  # Check that result has pooled models
  expect_s3_class(result, "bace_complete")
  expect_true("pooled_models" %in% names(result))
  
  # Check that pooled model used sampling
  y_model <- result$pooled_models$models$y
  expect_equal(nrow(y_model$Sol), sample_size * 4)
  expect_true(y_model$BACE_pooling$sampled)
  expect_equal(y_model$BACE_pooling$n_samples_per_imputation, sample_size)
})


test_that("bace wrapper without sample_size uses all samples", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_sampling_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run bace without sample_size (should use all)
  result <- bace(
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,  # Need at least 3 for convergence assessment
    n_final = 3,
    nitt = 600,
    thin = 2,
    burnin = 100,
    verbose = FALSE,
    plot = FALSE,
    sample_size = NULL,
    skip_conv = TRUE
  )
  
  # Check that no sampling was used
  y_model <- result$pooled_models$models$y
  expect_false(y_model$BACE_pooling$sampled)
  
  # Total samples should be (nitt - burnin) / thin * n_final
  expected_samples_per <- (600 - 100) / 2
  expect_equal(nrow(y_model$Sol), expected_samples_per * 3)
})


# Test print methods with sampling ----
test_that("print methods display sampling information correctly", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_sampling_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Run workflow with sampling
  initial <- bace_imp(
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    nitt = 600,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 3,
    nitt = 600,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Pool with sampling
  pooled <- pool_posteriors(final, sample_size = 50)
  
  # Test print.bace_pooled - should not error and should mention sampling
  output <- capture.output(print(pooled))
  expect_true(any(grepl("50", output)))  # Sample size should appear
  expect_true(any(grepl("sampled", output, ignore.case = TRUE)))
  
  # Test print.bace_pooled_MCMCglmm - should not error
  y_model <- pooled$models$y
  output_model <- capture.output(print(y_model))
  expect_true(length(output_model) > 0)
  expect_true(any(grepl("sampled", output_model, ignore.case = TRUE)))
})


# Test validation and edge cases ----
test_that("pool_posteriors validates sample_size parameter", {
  skip_on_cran()
  
  # Setup
  test_setup <- setup_sampling_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Create a minimal final object
  initial <- bace_imp(
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    nitt = 400,
    thin = 2,
    burnin = 50,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 2,
    nitt = 400,
    thin = 2,
    burnin = 50,
    verbose = FALSE
  )
  
  # Test with sample_size = 0 (should warn and use all samples)
  expect_warning(
    result <- pool_posteriors(final, sample_size = 0),
    "sample_size must be >= 1"
  )
  expect_s3_class(result, "bace_pooled")
  expect_false(result$models$y$BACE_pooling$sampled)
  
  # Test with negative sample_size (should warn and use all samples)
  expect_warning(
    result <- pool_posteriors(final, sample_size = -5),
    "sample_size must be >= 1"
  )
  expect_s3_class(result, "bace_pooled")
  
  # Test with non-integer sample_size (should warn and round down)
  expect_warning(
    result <- pool_posteriors(final, sample_size = 50.7),
    "must be an integer.*Rounding down"
  )
  expect_s3_class(result, "bace_pooled")
  expect_equal(result$models$y$BACE_pooling$n_samples_per_imputation, 50)
  
  # Test with non-numeric sample_size (should error)
  expect_error(
    pool_posteriors(final, sample_size = "50"),
    "must be a single numeric value"
  )
  
  # Test with vector sample_size (should error)
  expect_error(
    pool_posteriors(final, sample_size = c(50, 100)),
    "must be a single numeric value"
  )
})


test_that("sampling maintains posterior structure for different model types", {
  skip_on_cran()
  
  # Setup with different response types
  set.seed(54321)
  phylo <- ape::rtree(20)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length / max(ape::node.depth.edgelength(phylo))
  
  data <- data.frame(
    y_count = rpois(20, lambda = 5),       # Count
    y_continuous = rnorm(20, 10, 2),       # Gaussian
    x1 = rnorm(20, 5, 1),
    Species = phylo$tip.label
  )
  
  # Add missingness
  data$y_count[c(1, 5, 10)] <- NA
  data$y_continuous[c(2, 7, 12)] <- NA
  
  # Run workflow
  initial <- bace_imp(
    fixformula = list("y_count ~ x1", "y_continuous ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 3,
    nitt = 600,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  final <- bace_final_imp(
    bace_object = initial,
    fixformula = list("y_count ~ x1", "y_continuous ~ x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    n_final = 3,
    nitt = 600,
    thin = 2,
    burnin = 100,
    verbose = FALSE
  )
  
  # Pool with sampling
  pooled <- pool_posteriors(final, sample_size = 40)
  
  # Check both model types maintained structure after sampling
  expect_s3_class(pooled$models$y_count, "MCMCglmm")
  expect_s3_class(pooled$models$y_continuous, "MCMCglmm")
  
  # Both should have correct number of samples
  expect_equal(nrow(pooled$models$y_count$Sol), 40 * 3)
  expect_equal(nrow(pooled$models$y_continuous$Sol), 40 * 3)
  
  # Both should have metadata indicating sampling
  expect_true(pooled$models$y_count$BACE_pooling$sampled)
  expect_true(pooled$models$y_continuous$BACE_pooling$sampled)
})
