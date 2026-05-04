# Tests for convergence_diagnostics paths not exercised elsewhere:
# categorical/threshold use_all_data=TRUE branch, plot helpers with
# missing inputs, plot_acf_convergence, plot_density_convergence,
# plot_imputation_convergence categorical path.

build_cat_fixture <- function(seed = 2026, n = 25) {
  set.seed(seed)
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y_cat = factor(sample(c("A", "B", "C"), n, replace = TRUE),
                   levels = c("A", "B", "C")),
    y_bin = factor(sample(c("no", "yes"), n, replace = TRUE),
                   levels = c("no", "yes")),
    x = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y_cat[c(2, 5, 8)] <- NA
  data$y_bin[c(3, 7, 11)] <- NA
  list(phylo = phylo, data = data)
}

# ---- categorical/threshold use_all_data branch ---------------------------

test_that("assess_convergence on binary trait with use_all_data=TRUE", {
  testthat::skip_on_cran()
  fix <- build_cat_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_bin ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 4, nitt = 400, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  out <- suppressWarnings(assess_convergence(step1, method = "summary",
                                              use_all_data = TRUE))
  expect_s3_class(out, "bace_convergence")
})

test_that("assess_convergence on categorical trait with use_all_data=FALSE", {
  testthat::skip_on_cran()
  fix <- build_cat_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_cat ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 4, nitt = 400, thin = 5, burnin = 100,
    ovr_categorical = FALSE,
    verbose        = FALSE
  )))
  out <- suppressWarnings(assess_convergence(step1, method = "summary",
                                              use_all_data = FALSE))
  expect_s3_class(out, "bace_convergence")
})

test_that("assess_convergence energy method on categorical runs", {
  testthat::skip_on_cran()
  fix <- build_cat_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_bin ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 4, nitt = 300, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  out <- suppressWarnings(assess_convergence(step1, method = "energy"))
  expect_s3_class(out, "bace_convergence")
})

test_that("assess_convergence wasserstein method on categorical runs", {
  testthat::skip_on_cran()
  fix <- build_cat_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_bin ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 4, nitt = 300, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  out <- suppressWarnings(assess_convergence(step1, method = "wasserstein"))
  expect_s3_class(out, "bace_convergence")
})

test_that("assess_convergence method='all' on categorical runs", {
  testthat::skip_on_cran()
  fix <- build_cat_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_bin ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 4, nitt = 300, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  out <- suppressWarnings(assess_convergence(step1, method = "all"))
  expect_s3_class(out, "bace_convergence")
  # method='all' aggregates votes across summary/energy/wasserstein
  expect_true(!is.null(out$diagnostics$votes))
})

# ---- plot helpers --------------------------------------------------------

test_that("plot_trace_convergence errors with no summary_stats", {
  fake_conv <- list(summary_stats = NULL,
                    method_results = list(),
                    types = list())
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_trace_convergence(fake_conv),
               regexp = "No summary statistics")
})

test_that("plot_pct_change_convergence errors with no summary_stats", {
  fake_conv <- list(summary_stats = NULL,
                    method_results = list(),
                    types = list())
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_pct_change_convergence(fake_conv),
               regexp = "No summary statistics")
})

test_that("plot_acf_convergence errors with no summary_stats", {
  fake_conv <- list(summary_stats = NULL,
                    method_results = list(),
                    types = list())
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_acf_convergence(fake_conv),
               regexp = "No summary statistics")
})

test_that("plot_energy_convergence errors with no energy results", {
  fake_conv <- list(method_results = list(energy = NULL),
                    types = list())
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_energy_convergence(fake_conv),
               regexp = "No energy distance results")
})

test_that("plot_density_convergence errors with no summary_stats", {
  fake_conv <- list(summary_stats = NULL,
                    method_results = list(),
                    types = list())
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_density_convergence(fake_conv),
               regexp = "No summary statistics")
})
