# Tests for R/convergence_diagnostics.R
# Builds a tiny bace_imp result once and runs each diagnostic method
# on it. Exercises assess_convergence() through every documented
# `method` and the print.bace_convergence S3 method.

# ---- Shared fixture: a tiny bace_imp result with multiple iterations ----

build_tiny_bace_object <- function(seed = 2026) {
  set.seed(seed)
  phylo <- ape::rtree(12)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))

  data <- data.frame(
    y = rnorm(12, 0, 1),
    x1 = rnorm(12, 0, 1),
    x2 = rnorm(12, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  # Holes in y so there's something to impute / converge.
  data$y[c(1, 4, 7)]   <- NA
  data$x1[c(2, 5)]     <- NA

  bace_imp(
    fixformula     = list("y ~ x1 + x2", "x1 ~ y + x2"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo,
    data           = data,
    runs           = 5,
    nitt           = 600, thin = 5, burnin = 200,
    verbose        = FALSE
  )
}

# ---- assess_convergence: input validation ----

test_that("assess_convergence rejects non-bace input", {
  expect_error(assess_convergence(list()),
               regexp = "must be an object of class 'bace'")
})

test_that("assess_convergence warns when too few iterations", {
  # Build a fake bace object with only 2 iterations
  fake <- list(
    data     = list(NULL, data.frame(y = c(1,2)), data.frame(y = c(1,3))),
    types    = list(y = "gaussian"),
    miss_dat = data.frame(rowname = "1", colname = "y",
                           stringsAsFactors = FALSE)
  )
  class(fake) <- "bace"
  expect_warning(
    res <- assess_convergence(fake, method = "summary", min_iterations = 5),
    regexp = "Insufficient iterations|Only"
  )
  expect_true(is.na(res$converged))
})

# ---- assess_convergence: each method runs end-to-end ----

skip_if_slow <- function() {
  testthat::skip_on_cran()
}

test_that("assess_convergence(method='summary') returns expected structure", {
  skip_if_slow()
  bo <- suppressMessages(suppressWarnings(build_tiny_bace_object()))
  out <- suppressWarnings(assess_convergence(bo, method = "summary"))
  expect_type(out, "list")
  expect_s3_class(out, "bace_convergence")
  expect_true(!is.null(out$converged) || !is.null(out$method_results))
})

test_that("assess_convergence summary sub-criteria each run", {
  skip_if_slow()
  bo <- suppressMessages(suppressWarnings(build_tiny_bace_object()))
  for (m in c("summary.acf", "summary.percentage",
              "summary.trend", "summary.geweke")) {
    out <- suppressWarnings(assess_convergence(bo, method = m))
    expect_s3_class(out, "bace_convergence")
  }
})

test_that("assess_convergence(method='energy') runs", {
  skip_if_slow()
  bo <- suppressMessages(suppressWarnings(build_tiny_bace_object()))
  out <- suppressWarnings(assess_convergence(bo, method = "energy"))
  expect_s3_class(out, "bace_convergence")
})

test_that("assess_convergence(method='wasserstein') runs", {
  skip_if_slow()
  bo <- suppressMessages(suppressWarnings(build_tiny_bace_object()))
  out <- suppressWarnings(assess_convergence(bo, method = "wasserstein"))
  expect_s3_class(out, "bace_convergence")
})

test_that("assess_convergence(method='all') runs", {
  skip_if_slow()
  bo <- suppressMessages(suppressWarnings(build_tiny_bace_object()))
  out <- suppressWarnings(assess_convergence(bo, method = "all"))
  expect_s3_class(out, "bace_convergence")
  # 'all' should populate multiple method_results entries
  expect_true(length(out$method_results) >= 2L)
})

# ---- variable filtering + edge cases ----

test_that("assess_convergence filters to specified variables", {
  skip_if_slow()
  bo <- suppressMessages(suppressWarnings(build_tiny_bace_object()))
  out <- suppressWarnings(
    assess_convergence(bo, method = "summary", variables = "y"))
  expect_s3_class(out, "bace_convergence")
})

test_that("assess_convergence stops if no specified variable has missingness", {
  skip_if_slow()
  bo <- suppressMessages(suppressWarnings(build_tiny_bace_object()))
  expect_error(
    suppressWarnings(assess_convergence(bo, method = "summary",
                                         variables = "nonexistent_var")),
    regexp = "no missing data|None of"
  )
})

# ---- print.bace_convergence S3 method ----

test_that("print.bace_convergence runs without error", {
  skip_if_slow()
  bo <- suppressMessages(suppressWarnings(build_tiny_bace_object()))
  out <- suppressWarnings(assess_convergence(bo, method = "summary"))
  expect_output(print(out), ".+")  # any output
})
