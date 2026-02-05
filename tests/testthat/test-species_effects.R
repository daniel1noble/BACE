# Tests for species random effects decomposition functionality

# Setup test data with replicated species ----
setup_replicated_species_data <- function(n_species = 10, n_reps = 3) {
  set.seed(42)
  
  # Create a phylogenetic tree
  phylo <- ape::rtree(n_species)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  
  # Create replicated data for each species
  n_total <- n_species * n_reps
  
  data <- data.frame(
    y = rnorm(n_total, 10, 2),
    x1 = rnorm(n_total, 5, 1),
    x2 = rnorm(n_total, 3, 0.5),
    Species = rep(phylo$tip.label, each = n_reps)
  )
  
  list(phylo = phylo, data = data)
}

# Test 1: species=FALSE works as before (backward compatibility) ----
test_that("bace_imp with species=FALSE works as before (single phylo effect)", {
  test_setup <- setup_replicated_species_data(n_species = 5, n_reps = 2)
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data
  data$y[c(1, 3, 5)] <- NA
  data$x1[c(5, 8, 10)] <- NA
  
  # Run with species=FALSE (default)
  result <- bace_imp(
    fixformula = "y ~ x1 + x2",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 2,
    nitt = 100,
    thin = 1,
    burnin = 10,
    species = FALSE,
    verbose = FALSE
  )
  
  # Check that result is a bace object
  expect_s3_class(result, "bace")
  
  # Check that imputed values exist
  expect_true(!any(is.na(result$data$Iter_2$y)))
  
  # Check that phylo_ran structure is as expected (single cluster)
  expect_true("cluster" %in% names(result$phylo_ran))
})

# Test 2: species=TRUE with insufficient replication throws error ----
test_that("species=TRUE throws error with insufficient replicated species", {
  # Create data with no replication (each species appears once)
  set.seed(42)
  phylo <- ape::rtree(10)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  
  data <- data.frame(
    y = rnorm(10, 10, 2),
    x1 = rnorm(10, 5, 1),
    Species = phylo$tip.label
  )
  
  data$y[c(1, 3)] <- NA
  
  # Should throw error - no replicated species
  expect_error(
    bace_imp(
      fixformula = "y ~ x1",
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 100,
      species = TRUE,
      verbose = FALSE
    ),
    regexp = "fewer than 2 species having replicated observations"
  )
})

# Test 3: species=TRUE with low replication gives warning ----
test_that("species=TRUE gives warning with low proportion of replicated species", {
  # Create data where only 2 out of 10 species are replicated (20%)
  set.seed(42)
  phylo <- ape::rtree(10)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  
  # 8 species with 1 obs each, 2 species with 2 obs each = 12 total
  species_vec <- c(rep(phylo$tip.label[1:8], each = 1),
                   rep(phylo$tip.label[9:10], each = 2))
  
  data <- data.frame(
    y = rnorm(12, 10, 2),
    x1 = rnorm(12, 5, 1),
    Species = species_vec
  )
  
  data$y[c(1, 3)] <- NA
  
  # Should give warning about low replication
  expect_warning(
    bace_imp(
      fixformula = "y ~ x1",
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 100,
      thin = 1,
      burnin = 10,
      species = TRUE,
      verbose = FALSE
    ),
    regexp = "may be unreliable with limited replication"
  )
})

# Test 4: species=TRUE with good replication works ----
test_that("species=TRUE with sufficient replication runs successfully", {
  test_setup <- setup_replicated_species_data(n_species = 8, n_reps = 3)
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data
  data$y[c(1, 5, 10, 15)] <- NA
  data$x1[c(2, 6)] <- NA
  
  # Should run without error
  result <- expect_no_error(
    bace_imp(
      fixformula = "y ~ x1 + x2",
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 2,
      nitt = 100,
      thin = 1,
      burnin = 10,
      species = TRUE,
      verbose = FALSE
    )
  )
  
  # Check that result is a bace object
  expect_s3_class(result, "bace")
  
  # Check that imputed values exist
  expect_true(!any(is.na(result$data$Iter_2$y)))
  expect_true(!any(is.na(result$data$Iter_2$x1)))
})

# Test 5: .build_formula_string_random with species=FALSE ----
test_that(".build_formula_string_random with species=FALSE returns single formula", {
  result <- .build_formula_string_random("~ 1 | Species", species = FALSE)
  
  # Should be a single formula
  expect_true(inherits(result, "formula"))
  expect_false(is.list(result))
  
  # Should extract just the cluster variable
  expect_equal(as.character(result), c("~", "Species"))
})

# Test 6: .build_formula_string_random with species=TRUE returns list ----
test_that(".build_formula_string_random with species=TRUE returns list of two formulas", {
  result <- .build_formula_string_random("~ 1 | Species", species = TRUE)
  
  # Should be a list
  expect_true(is.list(result))
  expect_false(inherits(result, "formula"))
  
  # Should have phylo and species elements
  expect_true(all(c("phylo", "species") %in% names(result)))
  
  # Both should be formulas
  expect_true(inherits(result$phylo, "formula"))
  expect_true(inherits(result$species, "formula"))
  
  # For random intercept only, both should be simple
  expect_equal(as.character(result$phylo), c("~", "Species"))
  expect_equal(as.character(result$species), c("~", "Species"))
})

# Test 7: .build_formula_string_random with random slopes ----
test_that(".build_formula_string_random handles random slopes correctly", {
  # Random slope specification - note the pipe separator is required
  result <- .build_formula_string_random("~ us(1 + x1) | Species", species = TRUE)
  
  expect_true(is.list(result))
  expect_true(all(c("phylo", "species") %in% names(result)))
  
  # Phylo should have the slope structure
  phylo_char <- paste(as.character(result$phylo), collapse = " ")
  expect_true(grepl("us|:", phylo_char))
  
  # Species should be intercept only
  expect_equal(as.character(result$species), c("~", "Species"))
})

# Test 8: .make_prior handles n_rand=2 ----
test_that(".make_prior handles n_rand=2 for dual random effects", {
  # Test gaussian
  prior_1 <- .make_prior(n_rand = 1, type = "gaussian")
  prior_2 <- .make_prior(n_rand = 2, type = "gaussian")
  
  expect_equal(length(prior_1$G), 1)
  expect_equal(length(prior_2$G), 2)
  
  # Test poisson
  prior_pois_2 <- .make_prior(n_rand = 2, type = "poisson")
  expect_equal(length(prior_pois_2$G), 2)
  
  # Test ordinal
  prior_ord_2 <- .make_prior(n_rand = 2, type = "ordinal")
  expect_equal(length(prior_ord_2$G), 2)
})

# Test 9: .make_prior validates n_rand ----
test_that(".make_prior validates n_rand is 1 or 2", {
  expect_error(
    .make_prior(n_rand = 0, type = "gaussian"),
    regexp = "n_rand must be 1.*or 2"
  )
  
  expect_error(
    .make_prior(n_rand = 3, type = "gaussian"),
    regexp = "n_rand must be 1.*or 2"
  )
})

# Test 10: Integration test with multiple variables ----
test_that("species=TRUE works with multiple variables in formula list", {
  test_setup <- setup_replicated_species_data(n_species = 6, n_reps = 4)
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  # Introduce missing data in multiple variables
  data$y[c(1, 5, 9)] <- NA
  data$x1[c(2, 6, 10)] <- NA
  data$x2[c(3, 7, 11)] <- NA
  
  # Provide formulas for all variables with missing data
  fixformula_list <- list(
    "y ~ x1 + x2",
    "x1 ~ x2 + y",
    "x2 ~ x1 + y"
  )
  
  # Should run successfully
  result <- expect_no_error(
    bace_imp(
      fixformula = fixformula_list,
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 2,
      nitt = 100,
      thin = 1,
      burnin = 10,
      species = TRUE,
      verbose = FALSE
    )
  )
  
  # Check all variables are imputed
  expect_true(!any(is.na(result$data$Iter_2$y)))
  expect_true(!any(is.na(result$data$Iter_2$x1)))
  expect_true(!any(is.na(result$data$Iter_2$x2)))
})

# Test 11: Verify message about replication is shown when verbose=TRUE ----
test_that("species=TRUE with verbose=TRUE shows replication message", {
  test_setup <- setup_replicated_species_data(n_species = 5, n_reps = 3)
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$y[c(1, 3)] <- NA
  
  # Capture messages
  expect_message(
    bace_imp(
      fixformula = "y ~ x1",
      ran_phylo_form = "~ 1 | Species",
      phylo = phylo,
      data = data,
      runs = 1,
      nitt = 50,
      thin = 1,
      burnin = 5,
      species = TRUE,
      verbose = TRUE
    ),
    regexp = "Species effect decomposition enabled"
  )
})

# Test 12: Backward compatibility - old code still works ----
test_that("Existing code without species argument still works (backward compatibility)", {
  test_setup <- setup_replicated_species_data(n_species = 5, n_reps = 2)
  phylo <- test_setup$phylo
  data <- test_setup$data
  
  data$y[c(1, 3)] <- NA
  
  # Call without species argument - should default to FALSE
  result <- expect_no_error(
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
  
  expect_s3_class(result, "bace")
  expect_true(!any(is.na(result$data$Iter_1$y)))
})
