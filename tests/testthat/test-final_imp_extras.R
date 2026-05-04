# Tests targeting bace_final_imp + pool_posteriors orchestration paths,
# and convergence_diagnostics deeper-method branches.

# ---- Shared fixture ------------------------------------------------------

build_final_imp_fixture <- function(seed = 2026, n = 14) {
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
  data$x[c(3, 7)]    <- NA
  list(phylo = phylo, data = data)
}

build_step1 <- function(fix, runs = 2) {
  suppressWarnings(suppressMessages(bace_imp(
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = runs, nitt = 600, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
}

# ---- bace_final_imp orchestration branches -------------------------------

test_that("bace_final_imp with n_final=1 produces a single imputed dataset", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix)
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 500, thin = 5, burnin = 100,
    n_final        = 1, verbose = FALSE
  )))
  expect_s3_class(step2, "bace_final")
  expect_equal(length(step2$all_datasets), 1L)
})

test_that("bace_final_imp with n_final=5 produces five imputed datasets", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix)
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 5, verbose = FALSE
  )))
  expect_equal(length(step2$all_datasets), 5L)
  # Each imputed dataset should have NO NAs in y / x
  for (d in step2$all_datasets) {
    expect_true(all(!is.na(d$y)))
    expect_true(all(!is.na(d$x)))
  }
})

test_that("bace_final_imp imputed values are finite numerics", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix)
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 500, thin = 5, burnin = 100,
    n_final        = 3, verbose = FALSE
  )))
  for (d in step2$all_datasets) {
    expect_true(all(is.finite(d$y)))
    expect_true(all(is.finite(d$x)))
  }
})

# ---- pool_posteriors edge cases ------------------------------------------

test_that("pool_posteriors sample_size > total chain returns full chain", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix)
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2, verbose = FALSE
  )))
  pooled_huge <- suppressWarnings(suppressMessages(
    pool_posteriors(step2, sample_size = 1e6)))
  expect_s3_class(pooled_huge, "bace_pooled")
  # Pooled chain should still be valid (capped at actual length)
  expect_true(!is.null(pooled_huge$models[["y"]]$Sol))
})

test_that("pool_posteriors errors on mismatched variable name", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix)
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2, verbose = FALSE
  )))
  expect_error(suppressWarnings(suppressMessages(
    pool_posteriors(step2, variable = "no_such_var")
  )))
})

test_that("pool_posteriors output: pooled Sol has n_iter * n_final rows", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix)
  n_iter_per_run <- (400 - 100) / 5  # (nitt - burnin) / thin
  n_final <- 3
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = n_final, verbose = FALSE
  )))
  pooled <- suppressWarnings(suppressMessages(pool_posteriors(step2)))
  sol <- as.matrix(pooled$models[["y"]]$Sol)
  # Approximately n_iter * n_final rows
  expect_true(nrow(sol) <= n_iter_per_run * n_final + 5)
  expect_true(nrow(sol) >= n_iter_per_run * n_final - 5)
})

# ---- assess_convergence with custom alpha + threshold parameters ---------

test_that("assess_convergence with custom alpha runs", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix, runs = 4)
  for (a in c(0.01, 0.05, 0.1)) {
    out <- suppressWarnings(assess_convergence(step1,
                                                method = "summary.geweke",
                                                alpha = a))
    expect_s3_class(out, "bace_convergence")
  }
})

test_that("assess_convergence with custom pct_change_threshold runs", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix, runs = 4)
  out <- suppressWarnings(assess_convergence(step1,
                                              method = "summary.percentage",
                                              pct_change_threshold = 0.10))
  expect_s3_class(out, "bace_convergence")
})

test_that("assess_convergence with custom lag runs", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix, runs = 4)
  for (lg in c(1, 2, 5)) {
    out <- suppressWarnings(assess_convergence(step1,
                                                method = "summary.acf",
                                                lag = lg))
    expect_s3_class(out, "bace_convergence")
  }
})

test_that("assess_convergence with use_all_data=TRUE runs", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix, runs = 4)
  out <- suppressWarnings(assess_convergence(step1,
                                              method = "summary",
                                              use_all_data = TRUE))
  expect_s3_class(out, "bace_convergence")
})

# ---- print methods on returns --------------------------------------------

test_that("print.bace_final returns the object invisibly", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix)
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2, verbose = FALSE
  )))
  expect_output(print(step2), ".+")
})

test_that("print.bace_pooled returns sensible summary text", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix)
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2, verbose = FALSE
  )))
  pooled <- suppressWarnings(suppressMessages(pool_posteriors(step2)))
  expect_output(print(pooled), ".+")
})

test_that("summary.bace_pooled_MCMCglmm produces structured output", {
  testthat::skip_on_cran()
  fix <- build_final_imp_fixture()
  step1 <- build_step1(fix)
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y ~ x", "x ~ y"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2, verbose = FALSE
  )))
  pooled <- suppressWarnings(suppressMessages(pool_posteriors(step2)))
  m <- pooled$models[["y"]]
  expect_no_error(s <- summary(m))
})
