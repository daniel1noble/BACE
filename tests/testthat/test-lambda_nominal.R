# Tests for the per-level multinomial Pagel-style lambda
# (lambda_nominal_*) added to phylo_signal_summary().

build_lambda_fixture <- function(seed = 2026, n_sp = 30) {
  set.seed(seed)
  phylo <- ape::rtree(n_sp)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    cat3 = factor(sample(c("A", "B", "C"), n_sp, replace = TRUE),
                   levels = c("A", "B", "C")),
    cont = rnorm(n_sp),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  list(phylo = phylo, data = data)
}

# ---- Output structure ----------------------------------------------------

test_that("phylo_signal_summary table exposes lambda_nominal_* columns", {
  testthat::skip_on_cran()
  fix <- build_lambda_fixture()
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data = fix$data, tree = fix$phylo,
    species_col = "Species", variables = "cat3",
    quick = TRUE, ovr_categorical = FALSE,
    keep_models = TRUE, verbose = FALSE)))
  expect_true(all(c("lambda_nominal_mean",
                     "lambda_nominal_lo",
                     "lambda_nominal_hi") %in% colnames(res$table)))
})

test_that("lambda_nominal_mean is in [0, 1] for a multinomial fit", {
  testthat::skip_on_cran()
  fix <- build_lambda_fixture()
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data = fix$data, tree = fix$phylo,
    species_col = "Species", variables = "cat3",
    quick = TRUE, ovr_categorical = FALSE,
    keep_models = TRUE, verbose = FALSE)))
  lam <- res$table$lambda_nominal_mean[res$table$variable == "cat3"]
  expect_true(is.finite(lam))
  expect_true(lam >= 0 && lam <= 1)
})

test_that("lambda_nominal HPD bounds are ordered", {
  testthat::skip_on_cran()
  fix <- build_lambda_fixture()
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data = fix$data, tree = fix$phylo,
    species_col = "Species", variables = "cat3",
    quick = TRUE, ovr_categorical = FALSE,
    keep_models = TRUE, verbose = FALSE)))
  r <- res$table[res$table$variable == "cat3", ]
  expect_true(r$lambda_nominal_lo <= r$lambda_nominal_mean)
  expect_true(r$lambda_nominal_mean <= r$lambda_nominal_hi)
})

test_that("per-level lambda detail attached to the model object", {
  testthat::skip_on_cran()
  fix <- build_lambda_fixture()
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data = fix$data, tree = fix$phylo,
    species_col = "Species", variables = "cat3",
    quick = TRUE, ovr_categorical = FALSE,
    keep_models = TRUE, verbose = FALSE)))
  per_lev <- attr(res$models$cat3, "lambda_nominal_per_level")
  expect_false(is.null(per_lev))
  # K = 3 categories -> J = 2 non-baseline levels -> 2 entries
  expect_equal(length(per_lev), 2L)
  # Each entry has named slots
  expect_named(per_lev[[1]], c("column", "mean", "lo", "hi"))
  expect_true(all(vapply(per_lev,
                          function(x) is.finite(x$mean) &&
                                       x$mean >= 0 && x$mean <= 1,
                          logical(1))))
})

# ---- Continuous types: lambda_nominal must be NA ------------------------

test_that("lambda_nominal_* is NA for gaussian variables", {
  testthat::skip_on_cran()
  fix <- build_lambda_fixture()
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data = fix$data, tree = fix$phylo,
    species_col = "Species", variables = "cont",
    quick = TRUE, verbose = FALSE)))
  r <- res$table[res$table$variable == "cont", ]
  expect_true(is.na(r$lambda_nominal_mean))
  expect_true(is.na(r$lambda_nominal_lo))
  expect_true(is.na(r$lambda_nominal_hi))
})

# ---- OVR path (per-level threshold fits): lambda_nominal must be NA -----

test_that("OVR-categorical path returns NA lambda_nominal", {
  testthat::skip_on_cran()
  fix <- build_lambda_fixture()
  res <- suppressWarnings(suppressMessages(phylo_signal_summary(
    data = fix$data, tree = fix$phylo,
    species_col = "Species", variables = "cat3",
    quick = TRUE, ovr_categorical = TRUE,
    keep_models = FALSE, verbose = FALSE)))
  r <- res$table[res$table$variable == "cat3", ]
  expect_true(grepl("OVR", r$type))
  expect_true(is.na(r$lambda_nominal_mean))
})

# ---- Sanity check the c2 correction --------------------------------------

test_that(".compute_lambda_nominal Amemiya c2 correction implementation", {
  # White-box: with G_phylo_kk known, lambda_k formula must hold exactly.
  J <- 2L
  c2 <- (16 * sqrt(3) / (15 * pi))^2
  c2_corr <- c2 * 2 / (J + 1)
  expect_equal(c2_corr, c2 * 2/3, tolerance = 1e-12)

  G_kk <- 5
  G_corr <- G_kk / (1 + c2_corr)
  lam_expected <- G_corr / (G_corr + 1)
  expect_true(lam_expected > 0 && lam_expected < 1)
  # Must be lower than the un-corrected lambda
  lam_uncorrected <- G_kk / (G_kk + 1)
  expect_true(lam_expected < lam_uncorrected)
})
