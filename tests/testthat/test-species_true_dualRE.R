# Tests for species=TRUE dual-RE code paths in model_functions.R
# (lines 25-192). Requires replicated species observations and
# triggers the dual phylo + non-phylo species random-effect setup.

build_replicated_fixture <- function(seed = 2026, n_species = 12, reps = 3) {
  set.seed(seed)
  phylo <- ape::rtree(n_species)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))

  # Replicated observations per species
  n <- n_species * reps
  data <- data.frame(
    y = rnorm(n, 0, 1),
    x = rnorm(n, 0, 1),
    Species = rep(phylo$tip.label, each = reps),
    stringsAsFactors = FALSE
  )
  data$y[c(2, 5, 8, 14, 20)] <- NA
  list(phylo = phylo, data = data)
}

test_that("bace_imp(species=TRUE) on replicated data runs", {
  testthat::skip_on_cran()
  fix <- build_replicated_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = 2,
    nitt = 400, thin = 5, burnin = 100,
    species        = TRUE,
    verbose        = FALSE
  )))
  expect_s3_class(step1, "bace")
})

test_that("bace_imp(species=TRUE) errors on insufficient replication", {
  testthat::skip_on_cran()
  set.seed(2026)
  n <- 12
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  # No replicates — every species has 1 observation
  data <- data.frame(
    y = rnorm(n, 0, 1),
    x = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y[c(2, 5)] <- NA
  expect_error(suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo,
    data           = data,
    runs           = 2,
    nitt = 400, thin = 5, burnin = 100,
    species        = TRUE,
    verbose        = FALSE
  ))), regexp = "Cannot decompose")
})

test_that("bace_final_imp(species=TRUE) on replicated data runs", {
  testthat::skip_on_cran()
  fix <- build_replicated_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = 2,
    nitt = 400, thin = 5, burnin = 100,
    species        = TRUE,
    verbose        = FALSE
  )))
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2,
    species        = TRUE,
    verbose        = FALSE
  )))
  expect_s3_class(step2, "bace_final")
})

test_that("bace(species=TRUE, skip_conv=TRUE) on replicated data runs", {
  testthat::skip_on_cran()
  fix <- build_replicated_fixture()
  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = 2, n_final = 2,
    nitt = 400, thin = 5, burnin = 100,
    species        = TRUE,
    skip_conv      = TRUE,
    max_attempts   = 1, n_cores = 1L,
    verbose        = FALSE
  )))
  expect_s3_class(res, "bace_complete")
})

test_that("phylo_signal_summary(species=TRUE) on replicated categorical runs", {
  testthat::skip_on_cran()
  set.seed(2026)
  n_species <- 15
  reps <- 3
  phylo <- ape::rtree(n_species)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  n <- n_species * reps
  data <- data.frame(
    y_cat = factor(sample(c("A", "B", "C"), n, replace = TRUE),
                   levels = c("A", "B", "C")),
    Species = rep(phylo$tip.label, each = reps),
    stringsAsFactors = FALSE
  )
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = data,
    tree        = phylo,
    species_col = "Species",
    variables   = "y_cat",
    species     = TRUE,
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
})

test_that("phylo_signal_summary(species=TRUE) on replicated threshold runs", {
  testthat::skip_on_cran()
  set.seed(2026)
  n_species <- 15
  reps <- 3
  phylo <- ape::rtree(n_species)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  n <- n_species * reps
  data <- data.frame(
    y_ord = factor(sample(c("1", "2", "3", "4"), n, replace = TRUE),
                   levels = c("1", "2", "3", "4"), ordered = TRUE),
    Species = rep(phylo$tip.label, each = reps),
    stringsAsFactors = FALSE
  )
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = data,
    tree        = phylo,
    species_col = "Species",
    variables   = "y_ord",
    species     = TRUE,
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
})
