# Final coverage push: parallel bace_final_imp, plot helper error
# paths, and plot_bace_imputation_comparison on a categorical fixture.

build_final_push_fixture <- function(seed = 2026, n = 12) {
  set.seed(seed)
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y = rnorm(n, 0, 1),
    x = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y[c(2, 5, 8)] <- NA
  list(phylo = phylo, data = data)
}

# ---- bace_final_imp parallel path (n_cores > 1) ---------------------------

test_that("bace_final_imp(n_cores=2) parallel path runs (skipped on Windows)", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")  # mclapply not available on Windows
  fix <- build_final_push_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 2, nitt = 400, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2,
    n_cores        = 2L,
    verbose        = FALSE
  )))
  expect_s3_class(step2, "bace_final")
})

test_that("bace_final_imp(n_cores=2, verbose=TRUE) prints completion messages", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  fix <- build_final_push_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 2, nitt = 400, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  out <- capture.output(suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2,
    n_cores        = 2L,
    verbose        = TRUE
  ))))
  combined <- paste(out, collapse = "\n")
  expect_true(nchar(combined) > 0L)
})

# ---- plot helpers: variables argument and invalid-variable error path ----

build_conv_obj <- function() {
  fix <- build_final_push_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 4, nitt = 400, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  suppressWarnings(assess_convergence(step1, method = "summary"))
}

test_that("plot_trace_convergence with valid variables arg runs", {
  testthat::skip_on_cran()
  conv <- build_conv_obj()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(suppressMessages(
    plot_trace_convergence(conv, variables = "y")
  ))
})

test_that("plot_trace_convergence with non-matching variables errors", {
  testthat::skip_on_cran()
  conv <- build_conv_obj()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(suppressMessages(
    plot_trace_convergence(conv, variables = "nonexistent_var")
  ), regexp = "No valid variables")
})

test_that("plot_pct_change_convergence with non-matching variables errors", {
  testthat::skip_on_cran()
  conv <- build_conv_obj()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(suppressMessages(
    plot_pct_change_convergence(conv, variables = "nonexistent_var")
  ), regexp = "No valid variables")
})

# ---- plot_bace_imputation_comparison input validation ------------------------

test_that("plot_bace_imputation_comparison errors on non-bace input", {
  expect_error(plot_bace_imputation_comparison(list(),
                                            variable = "y"),
               regexp = "class 'bace'")
})

test_that("plot_bace_imputation_comparison errors on missing variable", {
  testthat::skip_on_cran()
  fix <- build_final_push_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 3, nitt = 300, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  expect_error(plot_bace_imputation_comparison(step1, variable = "no_such_var"),
               regexp = "not found")
})

test_that("plot_bace_imputation_comparison errors on variable with no missing data", {
  testthat::skip_on_cran()
  fix <- build_final_push_fixture()
  # x has no NAs
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 3, nitt = 300, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  expect_error(plot_bace_imputation_comparison(step1, variable = "x"),
               regexp = "No missing data")
})

# ---- plot_bace_imputation_comparison categorical path ------------------------

test_that("plot_bace_imputation_comparison works on a binary trait", {
  testthat::skip_on_cran()
  set.seed(2026)
  n <- 20
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y_bin = factor(sample(c("no", "yes"), n, replace = TRUE),
                   levels = c("no", "yes")),
    x = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y_bin[c(2, 5, 8)] <- NA
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_bin ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo, data = data,
    runs           = 4, nitt = 300, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(suppressMessages(
    plot_bace_imputation_comparison(step1, variable = "y_bin")
  ))
})

# ---- plot_wasserstein_convergence error paths ----------------------------

test_that("plot_wasserstein_convergence errors with no wasserstein results", {
  fake_conv <- list(method_results = list(wasserstein = NULL),
                    types = list())
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_wasserstein_convergence(fake_conv),
               regexp = "No Wasserstein")
})
