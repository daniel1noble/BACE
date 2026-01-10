# Tests for bace_imp.R

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


# Test error handling for missing formulas ----
test_that("bace_imp throws error when variables with NA lack formulas in list", {
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data in multiple variables
  data$y[c(1, 3, 5)] <- NA
  data$x1[c(2, 4)] <- NA
  data$x2[c(6, 8)] <- NA
  data$x3[c(7, 9)] <- NA
  
  # Provide formulas for only some variables with missing data
  # y, x1, and x2 have formulas, but x3 does not
  fixformula_list <- list(
    "y ~ x1 + x2",
    "x1 ~ x2 + x3",
    "x2 ~ x1"
  )
  
  # Should throw an error because x3 has missing data but no formula
  expect_error(
    bace_imp(
      fixformula = fixformula_list,
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 100,
      verbose = FALSE
    )
  )
})

test_that("bace_imp succeeds when all variables with NA have formulas", {
  # Setup
  test_setup <- setup_test_data()
       phylo <- test_setup$phylo
        data <- test_setup$data
  
  # Introduce missing data in multiple variables
   data$y[c(1, 3)] <- NA
  data$x1[c(2, 4)] <- NA
  data$x2[c(6, 8)] <- NA
  
  # Provide formulas for all variables with missing data
  fixformula_list <- list(   "y ~ x1 + x2",
							"x1 ~ x2 + x3",
							"x2 ~ x1 + x3")
  
  # Should NOT throw an error because all variables with NA have formulas
  expect_no_error(
    bace_imp(
      fixformula = fixformula_list,
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 100,
      thin = 1,
      burnin = 10,
      verbose = FALSE
    )
  )
})

test_that("bace_imp handles single formula (not a list) correctly", {
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data only in y
  data$y[c(1, 3, 5)] <- NA
  
  # Single formula string (not a list)
  # This should work even if x1, x2, x3 don't have formulas
  # because they don't have missing data
  expect_no_error(
    bace_imp(
      fixformula = "y ~ x1 + x2",
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 100,
      thin = 1,
      burnin = 10,
      verbose = FALSE
    )
  )
})

test_that("bace_imp throws error when variable in formula not in data", {
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data
  data$y[c(1, 3)] <- NA
  
  # Formula references a variable that doesn't exist
  expect_error(
    bace_imp(
      fixformula = "y ~ x1 + x_nonexistent",
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      verbose = FALSE
    ),
    regexp = "The following variables are not in the dataframe.*x_nonexistent"
  )
})

test_that("bace_imp throws error when cluster variable has missing data", {
  # Setup
  test_setup <- setup_test_data()
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data in response and cluster variable
  data$y[c(1, 3)] <- NA
  data$Species[5] <- NA
  
  # Should throw an error because Species has missing data
  expect_error(
    bace_imp(
      fixformula = "y ~ x1 + x2",
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      verbose = FALSE
    ),
    regexp = "There are missing data in the random effect cluster variable.*Species"
  )
})
