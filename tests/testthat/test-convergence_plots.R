# Tests for R/convergence_plots.R -- each plot function should produce
# graphics output without erroring on a representative bace_convergence
# object. We pipe everything through a null PDF device so the tests run
# headless.

# ---- Shared fixture ---------------------------------------------------------

build_conv_fixture <- function(seed = 2026) {
  set.seed(seed)
  phylo <- ape::rtree(12)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y  = rnorm(12, 0, 1),
    x1 = rnorm(12, 0, 1),
    x2 = rnorm(12, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y[c(1, 4, 7)] <- NA
  data$x1[c(2, 5)]   <- NA
  bo <- bace_imp(
    fixformula     = list("y ~ x1 + x2", "x1 ~ y + x2"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo,
    data           = data,
    runs           = 5,
    nitt = 600, thin = 5, burnin = 200,
    verbose = FALSE
  )
  suppressWarnings(assess_convergence(bo, method = "all"))
}

with_null_device <- function(expr) {
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

# ---- plot_* functions -------------------------------------------------------

test_that("plot_trace_convergence runs without error", {
  testthat::skip_on_cran()
  cv <- suppressMessages(suppressWarnings(build_conv_fixture()))
  expect_no_error(with_null_device(plot_trace_convergence(cv)))
})

test_that("plot_density_convergence runs without error", {
  testthat::skip_on_cran()
  cv <- suppressMessages(suppressWarnings(build_conv_fixture()))
  expect_no_error(with_null_device(plot_density_convergence(cv)))
})

test_that("plot_pct_change_convergence runs without error", {
  testthat::skip_on_cran()
  cv <- suppressMessages(suppressWarnings(build_conv_fixture()))
  expect_no_error(with_null_device(plot_pct_change_convergence(cv)))
})

test_that("plot_acf_convergence runs without error", {
  testthat::skip_on_cran()
  cv <- suppressMessages(suppressWarnings(build_conv_fixture()))
  expect_no_error(with_null_device(plot_acf_convergence(cv)))
})

test_that("plot_energy_convergence runs without error", {
  testthat::skip_on_cran()
  cv <- suppressMessages(suppressWarnings(build_conv_fixture()))
  # Only meaningful if energy method ran in fixture
  if (!is.null(cv$method_results$energy)) {
    expect_no_error(with_null_device(plot_energy_convergence(cv)))
  } else {
    succeed("energy method not present; skipping plot test")
  }
})

test_that("plot_wasserstein_convergence runs without error", {
  testthat::skip_on_cran()
  cv <- suppressMessages(suppressWarnings(build_conv_fixture()))
  if (!is.null(cv$method_results$wasserstein)) {
    expect_no_error(with_null_device(plot_wasserstein_convergence(cv)))
  } else {
    succeed("wasserstein method not present; skipping plot test")
  }
})

test_that("plot_convergence_summary runs without error", {
  testthat::skip_on_cran()
  cv <- suppressMessages(suppressWarnings(build_conv_fixture()))
  expect_no_error(with_null_device(plot_convergence_summary(cv)))
})

test_that("plot.bace_convergence dispatches by type", {
  testthat::skip_on_cran()
  cv <- suppressMessages(suppressWarnings(build_conv_fixture()))
  for (tp in c("trace", "density", "acf", "pct_change", "energy")) {
    if (tp == "energy" && is.null(cv$method_results$energy)) next
    expect_no_error(with_null_device(plot(cv, type = tp)))
  }
})

test_that("plot_trace_convergence errors when summary_stats missing", {
  fake <- list()
  class(fake) <- "bace_convergence"
  expect_error(with_null_device(plot_trace_convergence(fake)),
               regexp = "No summary statistics")
})

# ---- plot_bace_imputation_comparison ---------------------------------------

test_that("plot_bace_imputation_comparison runs without error", {
  testthat::skip_on_cran()
  set.seed(2026)
  phylo <- ape::rtree(12)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y  = rnorm(12, 0, 1),
    x1 = rnorm(12, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y[c(1, 4, 7)] <- NA
  bo <- bace_imp(
    fixformula     = "y ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo, data = data,
    runs = 3, nitt = 400, thin = 5, burnin = 100, verbose = FALSE
  )
  expect_no_error(with_null_device(
    plot_bace_imputation_comparison(bo, variable = "y")))
})
