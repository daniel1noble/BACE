# Additional phylo_signal_summary tests targeting the species=TRUE
# dual-RE path, ordinal/categorical paths, and quick mode parameter
# variants. Boosts coverage beyond what test-phylo_signal_summary.R
# already exercises.

# ---- Shared fixture ------------------------------------------------------

build_signal_fixture <- function(seed = 2026, n = 30) {
  set.seed(seed)
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    cont   = rnorm(n, 0, 1),
    bin    = factor(sample(c("no", "yes"), n, replace = TRUE),
                     levels = c("no", "yes")),
    ord    = factor(sample(c("1","2","3","4"), n, replace = TRUE),
                     levels = c("1","2","3","4"), ordered = TRUE),
    cat    = factor(sample(c("A","B","C"), n, replace = TRUE),
                     levels = c("A","B","C")),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  list(data = data, tree = phylo)
}

test_that("phylo_signal_summary on a continuous trait returns lambda + K", {
  testthat::skip_on_cran()
  fix <- build_signal_fixture(seed = 2026)
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "Species",
    variables   = "cont",
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
  expect_true("lambda" %in% colnames(res$table))
  expect_true("K" %in% colnames(res$table))
  lam <- res$table$lambda[res$table$variable == "cont"]
  expect_true(is.finite(lam))
})

test_that("phylo_signal_summary on a binary trait returns D statistic", {
  testthat::skip_on_cran()
  fix <- build_signal_fixture(seed = 2026, n = 40)
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "Species",
    variables   = "bin",
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
  expect_true("D" %in% colnames(res$table))
})

test_that("phylo_signal_summary on an ordinal trait runs without error", {
  testthat::skip_on_cran()
  fix <- build_signal_fixture(seed = 2026, n = 40)
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "Species",
    variables   = "ord",
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
  expect_true("ord" %in% res$table$variable)
})

test_that("phylo_signal_summary on a categorical trait runs without error", {
  testthat::skip_on_cran()
  fix <- build_signal_fixture(seed = 2026, n = 40)
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "Species",
    variables   = "cat",
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
  expect_true("cat" %in% res$table$variable)
})

test_that("phylo_signal_summary handles multiple variables in one call", {
  testthat::skip_on_cran()
  fix <- build_signal_fixture(seed = 2026, n = 30)
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "Species",
    variables   = c("cont", "bin"),
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
  expect_true(all(c("cont", "bin") %in% res$table$variable))
})

test_that("phylo_signal_summary species=TRUE adds species random effect", {
  testthat::skip_on_cran()
  set.seed(2026)
  n_sp <- 20
  phylo <- ape::rtree(n_sp)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  # Replicated observations so species=TRUE can decompose
  data <- data.frame(
    cont   = rnorm(2 * n_sp, 0, 1),
    Species = rep(phylo$tip.label, each = 2),
    stringsAsFactors = FALSE
  )
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = data,
    tree        = phylo,
    species_col = "Species",
    variables   = "cont",
    species     = TRUE,
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
  # With species=TRUE, R_species_mean should be present + finite-or-NA
  expect_true("R_species_mean" %in% colnames(res$table))
})

test_that("phylo_signal_summary stores models when keep_models=TRUE", {
  testthat::skip_on_cran()
  fix <- build_signal_fixture(seed = 2026, n = 25)
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "Species",
    variables   = "cont",
    keep_models = TRUE,
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_true(!is.null(res$models))
  expect_true(length(res$models) >= 1L)
})

test_that("phylo_signal_summary keep_models=FALSE returns NULL models", {
  testthat::skip_on_cran()
  fix <- build_signal_fixture(seed = 2026, n = 25)
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "Species",
    variables   = "cont",
    keep_models = FALSE,
    quick       = TRUE,
    verbose     = FALSE
  )))
  expect_true(is.null(res$models) || length(res$models) == 0L)
})

test_that("phylo_signal_summary errors when species_col is missing from data", {
  fix <- build_signal_fixture(seed = 2026, n = 20)
  expect_error(suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "NonexistentCol",
    variables   = "cont",
    quick       = TRUE,
    verbose     = FALSE
  ))))
})

test_that("phylo_signal_summary errors on all-NA variable (logical class)", {
  # All-NA cells become class 'logical' which doesn't fit any of the
  # imputable types. Function errors with a referenced-variable message.
  fix <- build_signal_fixture(seed = 2026, n = 20)
  fix$data$cont <- NA  # Make all values missing -> class becomes logical
  expect_error(suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "Species",
    variables   = "cont",
    quick       = TRUE,
    verbose     = FALSE
  ))), regexp = "cont|class|classify")
})

test_that("phylo_signal_summary produces print-able output", {
  testthat::skip_on_cran()
  fix <- build_signal_fixture(seed = 2026, n = 25)
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data        = fix$data,
    tree        = fix$tree,
    species_col = "Species",
    variables   = "cont",
    quick       = TRUE,
    verbose     = FALSE
  )))
  # print method should not error
  expect_no_error(capture.output(print(res)))
})
