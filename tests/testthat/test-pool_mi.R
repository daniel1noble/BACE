# Tests for the Rubin's-rules MI pathway: with_imputations() + pool_mi().

# ---- with_imputations: dispatch + structure --------------------------------

make_imputed_list <- function(M = 6, n = 60, seed = 1) {
  set.seed(seed)
  base <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  base$y <- 1 + 0.5 * base$x1 - 0.3 * base$x2 + rnorm(n, 0, 0.5)
  lapply(seq_len(M), function(m) {
    d <- base
    idx <- sample(n, floor(n / 4))
    d$y[idx] <- d$y[idx] + rnorm(length(idx), 0, 0.4)  # vary imputed cells
    d
  })
}

test_that("with_imputations accepts a list of data.frames and returns bace_mi_fits", {
  dl <- make_imputed_list()
  fits <- with_imputations(dl, function(d) lm(y ~ x1 + x2, data = d),
                           .progress = FALSE)
  expect_s3_class(fits, "bace_mi_fits")
  expect_length(fits, length(dl))
  expect_true(all(vapply(fits, inherits, logical(1), "lm")))
})

test_that("with_imputations captures per-fit errors; pool_mi drops them", {
  dl <- make_imputed_list(M = 5)
  counter <- 0L
  fits <- withCallingHandlers(
    with_imputations(dl, function(d) {
      counter <<- counter + 1L
      if (counter == 2L) stop("boom")     # second imputation fails
      lm(y ~ x1 + x2, data = d)
    }, .progress = FALSE),
    warning = function(w) invokeRestart("muffleWarning"))
  expect_s3_class(fits, "bace_mi_fits")
  expect_equal(attr(fits, "n_failed"), 1L)
  expect_true(inherits(fits[[2]], "bace_mi_error"))
  # pool_mi drops the failure (4 good fits remain) with a warning.
  expect_warning(p <- pool_mi(fits), "Dropping 1 fit")
  expect_equal(attr(p, "m"), 4L)
})

# ---- pool_mi: golden test vs mice::pool ------------------------------------

test_that("pool_mi matches mice::pool to machine precision (lm, Barnard-Rubin df)", {
  skip_if_not_installed("mice")
  dl   <- make_imputed_list(M = 8)
  fits <- with_imputations(dl, function(d) lm(y ~ x1 + x2, data = d),
                           .progress = FALSE)

  bace_tab <- as.data.frame(pool_mi(fits, df_fun = stats::df.residual))
  mice_obj <- mice::pool(fits)
  mice_tab <- summary(mice_obj)

  expect_equal(bace_tab$estimate,  mice_tab$estimate,  tolerance = 1e-10)
  expect_equal(bace_tab$std.error, mice_tab$std.error, tolerance = 1e-10)
  expect_equal(bace_tab$df,        mice_tab$df,        tolerance = 1e-8)
  expect_equal(bace_tab$fmi,       mice_obj$pooled$fmi, tolerance = 1e-10)
})

test_that("pool_mi returns the documented columns and sane values", {
  dl   <- make_imputed_list()
  fits <- with_imputations(dl, function(d) lm(y ~ x1 + x2, data = d),
                           .progress = FALSE)
  p <- pool_mi(fits)
  expect_s3_class(p, "bace_pooled_mi")
  expect_setequal(colnames(p),
    c("term", "estimate", "std.error", "df", "statistic", "p.value",
      "conf.low", "conf.high", "fmi", "riv"))
  expect_true(all(p$std.error > 0))
  expect_true(all(p$fmi >= 0 & p$fmi <= 1))
  expect_true(all(p$riv >= 0))
  # Pooled SE exceeds the naive single-fit SE (imputation uncertainty added).
  se1 <- sqrt(diag(vcov(fits[[1]])))
  expect_true(all(p$std.error >= se1 * 0.9))
})

# ---- pool_mi: edge cases ---------------------------------------------------

test_that("pool_mi errors on < 2 fits and on mismatched term sets", {
  dl   <- make_imputed_list(M = 3)
  fits <- with_imputations(dl, function(d) lm(y ~ x1 + x2, data = d),
                           .progress = FALSE)
  expect_error(pool_mi(fits[1]), "at least 2")

  mixed <- list(lm(y ~ x1, data = dl[[1]]), lm(y ~ x1 + x2, data = dl[[2]]))
  expect_error(pool_mi(mixed), "differ across fits")
})

# ---- pool_mi: MCMCglmm path ------------------------------------------------

test_that("pool_mi accepts MCMCglmm, extracts only fixed effects, warns", {
  testthat::skip_on_cran()
  set.seed(7)
  n  <- 30
  sp <- paste0("s", seq_len(n))
  mk <- function(shift) {
    d <- data.frame(x = rnorm(n), Species = sp)
    d$y <- 1 + 0.8 * d$x + rnorm(n, 0, 0.5) + shift
    d
  }
  fit_mm <- function(d) suppressWarnings(MCMCglmm::MCMCglmm(
    y ~ x, random = ~ Species, data = d, family = "gaussian",
    verbose = FALSE, pr = TRUE, nitt = 1200, burnin = 200, thin = 2,
    prior = list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002)))))

  fits <- lapply(c(0, 0.1, -0.1), function(s) fit_mm(mk(s)))

  # Extractor returns ONLY fixed effects (no Species.* BLUPs).
  cf <- BACE:::.mcmcglmm_fixef_mean(fits[[1]])
  expect_setequal(names(cf), c("(Intercept)", "x"))

  expect_message(p <- pool_mi(fits), "pool_posteriors")
  expect_s3_class(p, "bace_pooled_mi")
  expect_setequal(p$term, c("(Intercept)", "x"))
  expect_true(all(is.finite(p$estimate)))
  expect_true(all(p$std.error > 0))
})
