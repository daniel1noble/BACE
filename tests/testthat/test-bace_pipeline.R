# Phase C-extended: end-to-end bace() / bace_final_imp() / pool_posteriors()
# tests permuting non-default arguments to drive the still-uncovered
# branches in model_functions.R + bace.R.
#
# Each test runs the full pipeline on a tiny tree (10-15 species) with
# minimal MCMC. They are slow by unit-test standards (~10-30s each)
# so they are gated with skip_on_cran(). The GHA workflow runs with
# NOT_CRAN=true so they execute in CI.

# ---- Shared tiny-mixed-types fixture --------------------------------------

build_mixed_pipeline_fixture <- function(seed = 2026, n = 12) {
  set.seed(seed)
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y_cont = rnorm(n, 0, 1),
    y_bin  = factor(sample(c("no","yes"), n, replace = TRUE),
                    levels = c("no","yes")),
    x1     = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y_cont[c(2, 5)] <- NA
  data$y_bin[c(3, 7)]  <- NA
  list(phylo = phylo, data = data)
}

# ---- bace() with skip_conv = TRUE ------------------------------------------

test_that("bace() runs end-to-end with skip_conv=TRUE", {
  testthat::skip_on_cran()
  fix <- build_mixed_pipeline_fixture(seed = 2026, n = 12)

  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = list("y_cont ~ y_bin + x1",
                          "y_bin ~ y_cont + x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = 2, n_final = 2,
    nitt = 500, thin = 5, burnin = 100,
    skip_conv      = TRUE,
    max_attempts   = 1,
    n_cores        = 1L,
    verbose        = FALSE
  )))
  # Structure
  expect_s3_class(res, "bace_complete")
  expect_true(!is.null(res$pooled_models))
  expect_true(!is.null(res$imputed_datasets))

  # Numerical: every imputed dataset should have NO NAs in y_cont /
  # y_bin (the columns we asked to impute), and imputed continuous
  # values should be finite + within 5 sd of the column mean.
  for (imp in res$imputed_datasets) {
    expect_true(all(!is.na(imp$y_cont)))
    expect_true(all(!is.na(imp$y_bin)))
    expect_true(all(is.finite(imp$y_cont)))
    expect_true(all(as.character(imp$y_bin) %in% c("no", "yes")))
    col_mean <- mean(fix$data$y_cont, na.rm = TRUE)
    col_sd   <- stats::sd(fix$data$y_cont, na.rm = TRUE)
    # 5-sd envelope is intentionally generous; we want to catch
    # imputed values that explode (e.g. 1e10) but not flake on
    # legitimate posterior tails.
    expect_true(all(abs(imp$y_cont - col_mean) < 5 * col_sd + 5))
  }

  # Pooled posterior: Sol entries should be finite numerics.
  for (m in res$pooled_models$models) {
    sol <- as.matrix(m$Sol)
    expect_true(all(is.finite(sol)))
  }
})

# ---- bace() with species = TRUE (dual random effects) ---------------------

test_that("bace() runs with species=TRUE (dual phylo + species effects)", {
  testthat::skip_on_cran()
  # species=TRUE requires multiple observations per species so the
  # phylo-vs-species variance can be decomposed. Build a fixture with
  # 2 reps per species.
  set.seed(2026)
  n_sp <- 12
  phylo <- ape::rtree(n_sp)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y_cont = rnorm(2 * n_sp, 0, 1),
    x1     = rnorm(2 * n_sp, 0, 1),
    Species = rep(phylo$tip.label, each = 2),
    stringsAsFactors = FALSE
  )
  data$y_cont[c(2, 5, 9)] <- NA

  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = "y_cont ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo,
    data           = data,
    runs           = 2, n_final = 2,
    nitt = 500, thin = 5, burnin = 100,
    species        = TRUE,
    skip_conv      = TRUE,
    max_attempts   = 1, n_cores = 1L,
    verbose        = FALSE
  )))
  expect_s3_class(res, "bace_complete")
})

# ---- bace() with phylo_signal = TRUE (short-circuit) ----------------------

test_that("bace(phylo_signal=TRUE) returns the signal preview without imputing", {
  testthat::skip_on_cran()
  fix <- build_mixed_pipeline_fixture(seed = 2026, n = 14)

  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = "y_cont ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    phylo_signal   = TRUE,
    verbose        = FALSE
  )))
  # Should return a phylo_signal-like object, not a bace_complete.
  expect_false(inherits(res, "bace_complete"))
})

# ---- bace_final_imp + pool_posteriors with sample_size --------------------

test_that("pool_posteriors with sample_size subsamples MCMC draws", {
  testthat::skip_on_cran()
  fix <- build_mixed_pipeline_fixture(seed = 2026, n = 14)

  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_cont ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs = 2, nitt = 600, thin = 5, burnin = 100,
    verbose = FALSE
  )))
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = "y_cont ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 600, thin = 5, burnin = 100,
    n_final = 3, verbose = FALSE
  )))

  # sample_size limits each pooled chain to that many draws
  pooled <- suppressWarnings(suppressMessages(
    pool_posteriors(step2, sample_size = 30)))
  expect_s3_class(pooled, "bace_pooled")
  sol <- as.matrix(pooled$models[["y_cont"]]$Sol)
  # Pooled chain should have approximately sample_size * n_final rows
  # (or simply be capped at <= sample_size * n_final).
  expect_true(nrow(sol) <= 30 * 3)
})

# ---- pool_posteriors single-variable extraction ---------------------------

test_that("pool_posteriors(variable=...) returns just one model", {
  testthat::skip_on_cran()
  fix <- build_mixed_pipeline_fixture(seed = 2026, n = 14)

  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = list("y_cont ~ y_bin + x1",
                          "y_bin ~ y_cont + x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs = 2, nitt = 500, thin = 5, burnin = 100,
    verbose = FALSE
  )))
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = list("y_cont ~ y_bin + x1",
                          "y_bin ~ y_cont + x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    nitt = 500, thin = 5, burnin = 100,
    n_final = 2, verbose = FALSE
  )))

  pooled_one <- suppressWarnings(suppressMessages(
    pool_posteriors(step2, variable = "y_cont")))
  expect_s3_class(pooled_one, "bace_pooled")
  expect_true(length(pooled_one$models) == 1L)
  expect_true("y_cont" %in% names(pooled_one$models))
})

# ---- print methods on bace_complete / bace_final / bace_pooled ------------

test_that("print methods on bace pipeline objects run without error", {
  testthat::skip_on_cran()
  fix <- build_mixed_pipeline_fixture(seed = 2026, n = 12)

  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = "y_cont ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs = 2, n_final = 2,
    nitt = 500, thin = 5, burnin = 100,
    skip_conv      = TRUE,
    max_attempts   = 1, n_cores = 1L,
    verbose        = FALSE
  )))
  expect_output(print(res),       ".+")
  expect_output(print(res$final_results), ".+")
  expect_output(print(res$pooled_models), ".+")
})
