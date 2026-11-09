# Tests targeting bace() orchestration paths not exercised in
# test-verbose_orchestration.R: phylo_signal=TRUE short-circuit,
# plot=TRUE branches, sample_size pooling messages, the
# non-converged warning path, and the categorical-with-cluster
# pipeline. Also covers the .conv_assess_series internal helper directly
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

# ---- .conv_assess_series internal helper edge cases ----------------------

test_that(".conv_assess_series handles an all-equal series without error", {
  out <- BACE:::.conv_assess_series(rep(3.14, 10), rep(0, 10), n_cells = 20)
  expect_identical(out$verdict, "ok")
})

test_that(".conv_assess_series survives 3-point series (minimum assessable)", {
  out <- BACE:::.conv_assess_series(c(1, 1.1, 0.9), rep(1, 3), n_cells = 20)
  expect_true(out$verdict %in% c("ok", "drift", "oscillation", "stuck"))
})

test_that(".conv_assess_series reports the Monte-Carlo noise floor", {
  out <- BACE:::.conv_assess_series(rnorm(10), rep(1, 10), n_cells = 50)
  expect_equal(unname(out$stats["ref"]), sqrt(2 / 50))
})

test_that(".conv_assess_series flags a loud two-mode flip even with sd noise", {
  set.seed(1)
  m <- 2 * (-1)^(1:12) + rnorm(12, 0, 0.1)
  out <- BACE:::.conv_assess_series(m, rep(1, 12) * exp(rnorm(12, 0, 0.05)),
                                    n_cells = 30)
  expect_identical(out$verdict, "oscillation")
})

test_that(".conv_assess_series drift beats oscillation ordering sanity", {
  # a drifting chain whose steps also alternate slightly must come out
  # drift or oscillation, never ok
  set.seed(2)
  m <- 0.6 * (1:12) + 0.2 * (-1)^(1:12)
  out <- BACE:::.conv_assess_series(m, rep(1, 12), n_cells = 30)
  expect_true(out$verdict %in% c("drift", "oscillation"))
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
  out <- suppressWarnings(assess_convergence(step1))
  expect_s3_class(out, "bace_convergence")
})
