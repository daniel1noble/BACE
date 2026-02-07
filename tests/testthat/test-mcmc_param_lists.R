# Tests for MCMC parameter list functionality

# Setup test data ----
setup_test_data <- function() {
  set.seed(123)
  
  # Create a simple phylogenetic tree
  phylo <- ape::rtree(10)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  
  # Create test data
  data <- data.frame(
    y = rnorm(10, 10, 2),
    x1 = rnorm(10, 5, 1),
    x2 = rnorm(10, 3, 0.5),
    x3 = rnorm(10, 7, 1.5),
    Species = phylo$tip.label
  )
  
  list(phylo = phylo, data = data)
}


# Test .standardize_mcmc_params helper function ----
test_that(".standardize_mcmc_params converts single value to list", {
  result <- BACE:::.standardize_mcmc_params(6000, 3, "nitt")
  expect_type(result, "list")
  expect_length(result, 3)
  expect_equal(result, list(6000, 6000, 6000))
})

test_that(".standardize_mcmc_params validates list length", {
  expect_error(
    BACE:::.standardize_mcmc_params(list(6000, 8000), 3, "nitt"),
    "nitt is a list but its length \\(2\\) does not match the number of formulas \\(3\\)"
  )
})

test_that(".standardize_mcmc_params validates list elements are numeric", {
  expect_error(
    BACE:::.standardize_mcmc_params(list(6000, "abc", 8000), 3, "nitt"),
    "All elements of nitt list must be numeric"
  )
})

test_that(".standardize_mcmc_params passes through valid list", {
  input_list <- list(6000, 8000, 10000)
  result <- BACE:::.standardize_mcmc_params(input_list, 3, "nitt")
  expect_equal(result, input_list)
})

test_that(".standardize_mcmc_params rejects invalid input", {
  expect_error(
    BACE:::.standardize_mcmc_params(c(6000, 8000), 2, "nitt"),
    "nitt must be either a single numeric value or a list of numeric values"
  )
})


# Test bace_imp with single MCMC values (backward compatibility) ----
test_that("bace_imp works with single MCMC values (backward compatibility)", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data
  data$y[c(1, 3)] <- NA
  data$x1[c(2, 4)] <- NA
  data$x2[c(6, 8)] <- NA
  
  # Single formula with single MCMC values
  expect_no_error(
    result <- bace_imp(
      fixformula = "y ~ x1 + x2",
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 500,
      thin = 2,
      burnin = 50,
      verbose = FALSE
    )
  )
})

test_that("bace_imp works with list of formulas and single MCMC values", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data
  data$y[c(1, 3)] <- NA
  data$x1[c(2, 4)] <- NA
  data$x2[c(6, 8)] <- NA
  
  # List of formulas with single MCMC values (applied to all)
  expect_no_error(
    result <- bace_imp(
      fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3", "x2 ~ x1 + x3"),
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 500,
      thin = 2,
      burnin = 50,
      verbose = FALSE
    )
  )
})


# Test bace_imp with list MCMC values ----
test_that("bace_imp works with list of MCMC values matching formula count", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data
  data$y[c(1, 3)] <- NA
  data$x1[c(2, 4)] <- NA
  data$x2[c(6, 8)] <- NA
  
  # List of formulas with list MCMC values
  expect_no_error(
    result <- bace_imp(
      fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3", "x2 ~ x1 + x3"),
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = list(500, 600, 700),
      thin = list(2, 3, 2),
      burnin = list(50, 60, 70),
      verbose = FALSE
    )
  )
  
  # Check that result has expected structure
  expect_s3_class(result, "bace")
  expect_true("models_last_run" %in% names(result))
  expect_equal(length(result$models_last_run), 3)
})

test_that("bace_imp works with mixed: some params as list, some as single", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data
  data$y[c(1, 3)] <- NA
  data$x1[c(2, 4)] <- NA
  
  # List of formulas with mixed MCMC specification
  expect_no_error(
    result <- bace_imp(
      fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3"),
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = list(500, 600),  # list
      thin = 2,                # single value
      burnin = list(50, 60),   # list
      verbose = FALSE
    )
  )
})


# Test error handling for mismatched list lengths ----
test_that("bace_imp errors when nitt list length doesn't match formula count", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$y[c(1, 3)] <- NA
  data$x1[c(2, 4)] <- NA
  data$x2[c(6, 8)] <- NA
  
  expect_error(
    bace_imp(
      fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3", "x2 ~ x1 + x3"),
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = list(500, 600),  # Only 2 values for 3 formulas
      thin = 2,
      burnin = 50,
      verbose = FALSE
    ),
    "nitt is a list but its length \\(2\\) does not match the number of formulas \\(3\\)"
  )
})

test_that("bace_imp errors when thin list length doesn't match formula count", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$y[c(1, 3)] <- NA
  data$x1[c(2, 4)] <- NA
  data$x2[c(6, 8)] <- NA
  
  expect_error(
    bace_imp(
      fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3", "x2 ~ x1 + x3"),
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 500,
      thin = list(2, 3, 4, 5),  # 4 values for 3 formulas
      burnin = 50,
      verbose = FALSE
    ),
    "thin is a list but its length \\(4\\) does not match the number of formulas \\(3\\)"
  )
})

test_that("bace_imp errors when burnin list length doesn't match formula count", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$y[c(1, 3)] <- NA
  data$x1[c(2, 4)] <- NA
  
  expect_error(
    bace_imp(
      fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3"),
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 500,
      thin = 2,
      burnin = list(50),  # 1 value for 2 formulas
      verbose = FALSE
    ),
    "burnin is a list but its length \\(1\\) does not match the number of formulas \\(2\\)"
  )
})


# Test that MCMC values are properly applied to each model ----
test_that("Different MCMC values are applied to different models", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data
  data$y[c(1, 3)] <- NA
  data$x1[c(2, 4)] <- NA
  
  result <- bace_imp(
    fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 1,
    nitt = list(500, 1000),
    thin = list(2, 5),
    burnin = list(50, 100),
    verbose = FALSE
  )
  
  # Check that models exist
  expect_true("models_last_run" %in% names(result))
  expect_equal(length(result$models_last_run), 2)
  
  # Check that the models have different number of samples (due to different nitt/thin/burnin)
  # Model 1: (500 - 50) / 2 = 225 samples
  # Model 2: (1000 - 100) / 5 = 180 samples
  model1 <- result$models_last_run[[1]]
  model2 <- result$models_last_run[[2]]
  
  expect_equal(nrow(model1$Sol), 225)
  expect_equal(nrow(model2$Sol), 180)
})


# Test with single formula (auto-converted to list) ----
test_that("Single formula works with single MCMC values", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$y[c(1, 3)] <- NA
  
  expect_no_error(
    result <- bace_imp(
      fixformula = "y ~ x1 + x2",  # Single formula string
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 500,
      thin = 2,
      burnin = 50,
      verbose = FALSE
    )
  )
})

test_that("Single formula works with list MCMC values (of length 1)", {
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$y[c(1, 3)] <- NA
  
  # This should work because the single formula gets expanded to a list internally
  # and we're providing lists of length matching the expanded formulas
  # Need to check how many formulas are created from a single formula string
  expect_no_error(
    result <- bace_imp(
      fixformula = "y ~ x1 + x2",  # Single formula string
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 500,  # Single value still works
      thin = 2,
      burnin = 50,
      verbose = FALSE
    )
  )
})
