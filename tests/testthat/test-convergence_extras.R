# Tests for convergence_diagnostics paths not exercised elsewhere:
# categorical/threshold series (proportion + binomial spread), the
# use_all_data=TRUE branch, defunct method errors on real fits, and the
# plot helpers' input validation.

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

# ---- categorical/threshold series through the stability check --------------

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
  out <- suppressWarnings(assess_convergence(step1, use_all_data = TRUE))
  expect_s3_class(out, "bace_convergence")
  expect_true(all(out$verdicts %in% c("ok", "drift", "oscillation", "stuck")))
  # the monitored series is a proportion: it must live in [0, 1]
  expect_true(all(out$summary_stats$y_bin >= 0 & out$summary_stats$y_bin <= 1))
  # and its spread series is the binomial sd, also bounded
  expect_true(all(out$sd_stats$y_bin >= 0 & out$sd_stats$y_bin <= 0.5 + 1e-12))
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
  out <- suppressWarnings(assess_convergence(step1, use_all_data = FALSE))
  expect_s3_class(out, "bace_convergence")
  expect_true(all(out$summary_stats$y_cat >= 0 & out$summary_stats$y_cat <= 1))
})

test_that("the defunct distribution methods error on a real fit too", {
  testthat::skip_on_cran()
  fix <- build_cat_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_bin ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo, data = fix$data,
    runs           = 4, nitt = 300, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  expect_error(assess_convergence(step1, method = "energy"), regexp = "defunct")
  expect_error(assess_convergence(step1, method = "wasserstein"), regexp = "defunct")
  expect_error(assess_convergence(step1, method = "all"), regexp = "defunct")
})

# ---- plot helpers ----------------------------------------------------------

test_that("plot_trace_convergence errors with no summary_stats", {
  fake_conv <- list(summary_stats = NULL,
                    method_results = list(),
                    types = list())
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_trace_convergence(fake_conv),
               regexp = "No summary statistics")
})

test_that("plot_steps_convergence errors with no summary_stats", {
  fake_conv <- list(summary_stats = NULL,
                    method_results = list(),
                    types = list())
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_steps_convergence(fake_conv),
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

test_that("plot_density_convergence errors with no summary_stats", {
  fake_conv <- list(summary_stats = NULL,
                    method_results = list(),
                    types = list())
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_density_convergence(fake_conv),
               regexp = "No summary statistics")
})

# ---- the defunct plot helpers point at their replacements ------------------

test_that("the pre-1.1 plot helpers are defunct with pointers", {
  fake_conv <- list(summary_stats = data.frame(iteration = 1:3, y = 1:3))
  class(fake_conv) <- c("bace_convergence", "list")
  expect_error(plot_pct_change_convergence(fake_conv), regexp = "defunct")
  expect_error(plot_energy_convergence(fake_conv), regexp = "defunct")
  expect_error(plot_wasserstein_convergence(fake_conv), regexp = "defunct")
  expect_error(plot(fake_conv, type = "pct_change"), regexp = "defunct")
  expect_error(plot(fake_conv, type = "energy"), regexp = "defunct")
})
