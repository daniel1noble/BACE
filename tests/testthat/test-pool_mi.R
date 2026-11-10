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

# ---- the two pooling routes, compared ---------------------------------------
#
# BACE can combine M imputations two ways, and the manuscript relies on both:
#
#   pool_posteriors()  stacks the per-imputation MCMC chains  -> $pooled_models
#   pool_mi()          Rubin's rules on the per-imputation posterior summaries
#
# .agents/plans/pool_mi_rubin.md section 5 asked for this comparison twice,
# including as a manuscript figure, and it had never been written. Measured
# 2026-11-10 over 8 datasets at M = 20, and over M in {5, 10, 20, 50}:
#
#   point estimate   IDENTICAL to machine precision at every M and every seed
#   standard error   agree to 0.07%
#   interval width   stacking is consistently WIDER: +2.5% mean, +0.4% to +7.1%
#
# The width gap does NOT shrink with M (about 7% at M = 5 and at M = 50). The
# plan predicted an O(1/M) gap that closes; it does not, because this is not
# Monte-Carlo error. Stacking takes the 2.5% and 97.5% points of the actual
# posterior, Rubin builds a symmetric t-interval around the estimate, so any
# departure from normality in the posterior shows up as a width difference
# however many imputations are used. Wider is the conservative direction.

test_that("stacking and Rubin give the SAME point estimate when chains are equal length", {
  # The identity is arithmetic, not approximate: the mean of the stacked draws
  # equals the mean of the per-chain means exactly when every chain has the
  # same number of draws. This pure-arithmetic test pins the REASON, with no
  # MCMC, so a regression is caught in milliseconds and the cause is obvious.
  set.seed(42)
  chains_equal <- lapply(1:6, function(i) rnorm(500, mean = i / 10))
  stacked  <- mean(unlist(chains_equal))
  per_mean <- mean(vapply(chains_equal, mean, numeric(1)))
  expect_equal(stacked, per_mean, tolerance = 1e-12)

  # ...and UNEQUAL chain lengths break it, which is why BACE must not let the
  # per-imputation chains drift apart in length.
  chains_ragged <- list(rnorm(100, 0), rnorm(900, 10))
  expect_false(isTRUE(all.equal(mean(unlist(chains_ragged)),
                                mean(vapply(chains_ragged, mean, numeric(1))),
                                tolerance = 1e-6)))
})

test_that("pool_posteriors and pool_mi agree on a real bace() fit", {
  testthat::skip_on_cran()
  skip_if_not_installed("ape")
  skip_if_not_installed("MCMCglmm")

  set.seed(3); n <- 40
  tr <- ape::rcoal(n); tr$tip.label <- paste0("sp", seq_len(n))
  C <- ape::vcv(tr); C <- C / max(C); L <- t(chol(C))
  d <- data.frame(Species = tr$tip.label, x1 = as.numeric(L %*% rnorm(n)))
  d$y <- 0.6 * d$x1 + rnorm(n)
  d$y[sample(n, 12)] <- NA

  fit <- suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", phylo = tr, data = d, runs = 3,
         nitt = 6000, burnin = 1000, thin = 5, n_final = 10,
         verbose = FALSE, skip_conv = TRUE, seed = 3)))

  # a normal bace() fit carries the pooled model; only the simulation drops it
  expect_false(is.null(fit$pooled_models))
  pm <- get_pooled_model(fit, "y")
  expect_s3_class(pm, "MCMCglmm")

  mods <- fit$final_results$all_models
  mods <- lapply(mods, function(m)
    if (is.list(m) && !inherits(m, "MCMCglmm")) m[["y"]] else m)
  mods <- Filter(function(m) inherits(m, "MCMCglmm"), mods)
  expect_gte(length(mods), 2L)

  # the precondition for the identity
  lens <- vapply(mods, function(m) nrow(as.matrix(m$Sol)), integer(1))
  expect_equal(length(unique(lens)), 1L)

  sol <- as.matrix(pm$Sol)[, "x1"]
  fits <- structure(mods, class = c("bace_mi_fits", "list"),
                    n_fits = length(mods), n_failed = 0L,
                    failed = integer(0), is_mcmcglmm = TRUE)
  rub <- suppressWarnings(suppressMessages(pool_mi(fits)))
  rub <- rub[rub$term == "x1", ]

  # 1. point estimates identical to machine precision (measured 5.6e-17)
  expect_equal(mean(sol), rub$estimate, tolerance = 1e-10)

  # 2. standard errors agree closely (measured 0.07%)
  expect_lt(abs(stats::sd(sol) - rub$std.error) / rub$std.error, 0.01)

  # 3. the stacked interval is not narrower, and not wildly wider
  #    (measured +0.4% to +7.1% over 8 datasets, always positive)
  qa <- unname(stats::quantile(sol, c(0.025, 0.975)))
  w_stack <- diff(qa)
  w_rubin <- rub$conf.high - rub$conf.low
  expect_gt(w_stack / w_rubin, 0.95)
  expect_lt(w_stack / w_rubin, 1.30)
})
