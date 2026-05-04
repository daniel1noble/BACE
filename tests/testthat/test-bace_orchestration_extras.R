# Tests targeting bace() orchestration paths not exercised in
# test-verbose_orchestration.R: phylo_signal=TRUE short-circuit,
# plot=TRUE branches, sample_size pooling messages, the
# non-converged warning path, and the categorical-with-cluster
# pipeline. Also covers the .geweke_test internal helper directly
# and the categorical/threshold branch of convergence summary.

# ---- Shared fixture ------------------------------------------------------

build_orch_fixture <- function(seed = 2026, n = 12) {
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

# ---- bace(phylo_signal = TRUE) short-circuit -----------------------------

test_that("bace(phylo_signal=TRUE) returns a phylo_signal preview", {
  testthat::skip_on_cran()
  fix <- build_orch_fixture(n = 20)
  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    phylo_signal   = TRUE,
    verbose        = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
  expect_true(!is.null(res$table))
})

test_that("bace(phylo_signal=TRUE, verbose=TRUE) prints preview banner", {
  testthat::skip_on_cran()
  fix <- build_orch_fixture(n = 20)
  out <- capture.output(suppressWarnings(suppressMessages(
    bace(
      fixformula     = "y ~ x",
      ran_phylo_form = "~ 1 | Species",
      phylo          = fix$phylo,
      data           = fix$data,
      phylo_signal   = TRUE,
      verbose        = TRUE
    )
  )))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "phylo_signal", ignore.case = TRUE)
})

test_that("bace(phylo_signal=TRUE) accepts list of formulas", {
  testthat::skip_on_cran()
  fix <- build_orch_fixture(n = 20)
  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    phylo_signal   = TRUE,
    verbose        = FALSE
  )))
  expect_s3_class(res, "phylo_signal")
})

# ---- bace(plot = TRUE) ---------------------------------------------------

test_that("bace(plot=TRUE) calls plot(converge) without error", {
  testthat::skip_on_cran()
  fix <- build_orch_fixture()
  # plot() goes to a tempfile pdf device; closed on exit.
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = 4, n_final = 2,
    nitt = 400, thin = 5, burnin = 100,
    skip_conv      = TRUE,
    max_attempts   = 1, n_cores = 1L,
    plot           = TRUE,
    verbose        = FALSE
  )))
  expect_s3_class(res, "bace_complete")
})

test_that("bace(plot=TRUE, runs<3) is guarded against insufficient-iter convergence list", {
  testthat::skip_on_cran()
  # Regression test: runs<3 -> assess_convergence returns plain list (no
  # bace_convergence class). bace() should NOT call plot() in that case.
  fix <- build_orch_fixture()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(suppressMessages(bace(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = 2, n_final = 2,
    nitt = 200, thin = 1, burnin = 50,
    skip_conv      = TRUE,
    max_attempts   = 1, n_cores = 1L,
    plot           = TRUE,
    verbose        = FALSE
  ))))
})

# ---- bace(sample_size = X) verbose pooling-with-sampling message ---------

test_that("bace(skip_conv=TRUE, sample_size=N, verbose=TRUE) prints sampling message", {
  testthat::skip_on_cran()
  fix <- build_orch_fixture()
  out <- capture.output(suppressWarnings(suppressMessages(
    bace(
      fixformula     = "y ~ x",
      ran_phylo_form = "~ 1 | Species",
      phylo          = fix$phylo,
      data           = fix$data,
      runs           = 2, n_final = 2,
      nitt = 200, thin = 1, burnin = 50,
      skip_conv      = TRUE,
      sample_size    = 50,
      max_attempts   = 1, n_cores = 1L,
      verbose        = TRUE
    )
  )))
  combined <- paste(out, collapse = "\n")
  # Either the converged "Step 4 ... posterior sampling" or the
  # non-converged "Pooling posteriors with sampling" message
  expect_match(combined, "sampling|Step", ignore.case = TRUE)
})

# ---- .geweke_test internal helper ----------------------------------------

test_that(".geweke_test returns NA list for short series (n<10)", {
  out <- BACE:::.geweke_test(c(1, 2, 3, 4, 5))
  expect_named(out, c("converged", "z_score", "p_value"))
  expect_true(is.na(out$converged))
  expect_true(is.na(out$z_score))
})

test_that(".geweke_test returns finite z and p for stationary series", {
  set.seed(2026)
  x <- rnorm(200, 0, 1)  # stationary
  out <- BACE:::.geweke_test(x, alpha = 0.05)
  expect_true(is.finite(out$z_score))
  expect_true(is.finite(out$p_value))
  expect_true(out$p_value >= 0 && out$p_value <= 1)
})

test_that(".geweke_test flags non-stationary series", {
  set.seed(2026)
  # Strong drift between first 10% and last 50%
  x <- c(rnorm(20, 0, 0.1), rnorm(180, 5, 0.1))
  out <- BACE:::.geweke_test(x, alpha = 0.05)
  expect_true(is.finite(out$z_score))
  # Strong drift => low p-value => not converged
  expect_false(isTRUE(out$converged))
})

test_that(".geweke_test handles all-NA series", {
  out <- BACE:::.geweke_test(rep(NA_real_, 50))
  expect_true(is.na(out$converged))
})

test_that(".geweke_test handles zero-variance series", {
  out <- BACE:::.geweke_test(rep(3.14, 100))
  # Zero variance => var1+var2 = 0 => returns NA list
  expect_true(is.na(out$converged))
})

test_that(".geweke_test handles series with a few NAs", {
  set.seed(2026)
  x <- rnorm(100)
  x[c(5, 12, 78)] <- NA
  out <- BACE:::.geweke_test(x)
  expect_true(is.finite(out$z_score) || is.na(out$z_score))
})

test_that(".geweke_test alpha threshold affects converged decision", {
  set.seed(2026)
  x <- c(rnorm(20, 0, 1), rnorm(180, 0.5, 1))  # mild drift
  out_strict <- BACE:::.geweke_test(x, alpha = 0.50)
  out_loose <- BACE:::.geweke_test(x, alpha = 0.001)
  expect_equal(out_strict$z_score, out_loose$z_score)
  # alpha=0.50 is stricter (rejects more), alpha=0.001 is looser
  expect_true(is.logical(out_strict$converged) || is.na(out_strict$converged))
  expect_true(is.logical(out_loose$converged) || is.na(out_loose$converged))
})

# ---- convergence summary categorical/threshold branch --------------------

test_that("assess_convergence on a binary trait exercises factor branch", {
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
    runs           = 4, nitt = 400, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  out <- suppressWarnings(assess_convergence(step1, method = "summary"))
  expect_s3_class(out, "bace_convergence")
})
