# Regression tests for the 2026-08 performance fixes.
#
# Root-cause background (investigation of the Study B failures; see
# .agents/roadmap.md):
#
# 1. bace_final_imp fit every final-phase model on the convergence chain's
#    point-imputed dataset with the fills treated as OBSERVED data (response
#    NAs never reinstated). That is improper multiple imputation (Rubin 1987
#    Ch. 4): the residual variance deflates (~30% of "observations" sit on the
#    fitted surface), the BLUPs anchor to the fills, and the final
#    "posterior-predictive draws" barely move — collapsing between-imputation
#    variance and producing severe interval under-coverage for gaussian
#    (0.76) and poisson (0.69).
#
# 2. The poisson point estimate was colMeans(exp(Liab)) — a mean of
#    lognormals dominated by single chain excursions (one exp(20) draw out of
#    800 adds ~1e6). Combined with (1) and the fixed pmin(rate, 1e6) clip,
#    this froze million-count imputations into every final dataset
#    (max standardised bias 217,838 in Study B).
#
# 3. The poisson priors were copied from gaussian: IW(V=1, nu=2) on the
#    log-link units variance forbids the typical overdispersion scale (~0.1),
#    and the PX working prior alpha.V = 1e4 permits log-scale random-effect
#    SDs of ~100.

# ---- 1. bace_final_imp reinstates response NAs (proper MI) -----------------

test_that("bace_final_imp fits each final model with response NAs reinstated", {
  testthat::skip_on_cran()
  set.seed(4090)
  n <- 30
  tr <- ape::rtree(n)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  C <- ape::vcv(tr, corr = TRUE)
  x <- as.numeric(MASS::mvrnorm(1, rep(0, n), C))
  u <- as.numeric(MASS::mvrnorm(1, rep(0, n), C))
  y <- 1 + 1.5 * x + u + rnorm(n, 0, 0.6)
  d <- data.frame(y = y, x = x, Species = tr$tip.label,
                  stringsAsFactors = FALSE)
  miss_rows <- sample(n, 9)
  d$y[miss_rows] <- NA

  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
    phylo = tr, data = d, runs = 2,
    nitt = 1200, thin = 2, burnin = 300, verbose = FALSE)))

  # The converged dataset is fully imputed; the load-bearing question is
  # whether the final phase re-blanks the originally missing response cells
  # before each fit (data augmentation) or fits on the point fills as if
  # observed (improper MI, the pre-fix behaviour). Spy on .model_fit.
  expect_false(anyNA(init$data[[length(init$data)]]$y))

  na_counts <- integer(0)
  real_fit <- BACE:::.model_fit
  testthat::local_mocked_bindings(
    .model_fit = function(data, ...) {
      na_counts <<- c(na_counts, sum(is.na(data$y)))
      real_fit(data = data, ...)
    },
    .package = "BACE"
  )

  fin <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = tr, n_final = 2,
    nitt = 1200, thin = 2, burnin = 300, verbose = FALSE)))

  # One fit per final run, each seeing exactly the originally-missing cells
  # as NA (the pre-fix behaviour recorded 0 NAs on every fit).
  expect_length(na_counts, 2L)
  expect_true(all(na_counts == length(miss_rows)))

  # And the imputations must still fill every originally missing cell.
  for (ds in fin$all_datasets) expect_false(anyNA(ds$y))
})

# ---- 2. Poisson point estimate + data-adaptive rate guard ------------------

mock_liab_model <- function(liab) list(Liab = coda::mcmc(liab))

test_that(".pred_count post_median is the exp of the median liability and is
           robust to single chain excursions", {
  set.seed(587)
  n_iter <- 100; n_obs <- 4
  liab <- matrix(rnorm(n_iter * n_obs, 1, 0.4), n_iter, n_obs)
  liab[1, 2] <- 20.4   # one excursion draw (exp = 7.3e8), as in Study B rep 587
  out <- BACE:::.pred_count(mock_liab_model(liab))

  expect_equal(out$post_median,
               as.numeric(exp(apply(liab, 2, stats::median))),
               tolerance = 1e-10)
  # The lognormal mean explodes from the single excursion; the median must not.
  expect_gt(out$post_mean[2], exp(20.4) / n_iter * 0.9)   # ~6.6e6
  expect_lt(out$post_median[2], 20)
})

test_that(".predict_bace poisson: excursion draws cannot produce catastrophic
           counts (data-adaptive rate cap replaces the fixed 1e6 clip)", {
  set.seed(587)
  n_iter <- 60; n_obs <- 6
  liab <- matrix(rnorm(n_iter * n_obs, 0.5, 0.3), n_iter, n_obs)
  liab[, 5] <- 20.4   # every draw for cell 5 is an excursion (worst case)
  model <- mock_liab_model(liab)

  # Observed counts small (max 2, as in Study B rep 587); two cells missing.
  y_obs <- c(0L, 1L, 2L, 1L, NA, NA)
  dat_prep <- list(data.frame(y = y_obs), NULL)
  cap <- 3 * 2 + 5   # 3 * max observed + 5 = 11

  # sample = FALSE: median-based point estimate, capped.
  pt <- BACE:::.predict_bace(model, dat_prep, response_var = "y",
                             type = "poisson", sample = FALSE)
  expect_true(all(pt$pred_values <= cap))

  # sample = TRUE: the rate is capped at 11, so a count anywhere near the old
  # rpois(1e6) regime is impossible (P(Pois(11) > 60) < 1e-24).
  for (i in 1:20) {
    dr <- BACE:::.predict_bace(model, dat_prep, response_var = "y",
                               type = "poisson", sample = TRUE)
    expect_true(all(dr$pred_values < 100))
  }
})

# ---- 3. Family-specific priors ---------------------------------------------

test_that("gaussian prior: weak Hadfield R-structure, PX G with alpha.V 1e4", {
  p <- BACE:::.make_prior(n_rand = 1, type = "gaussian")
  expect_equal(p$R$V, 1)
  expect_equal(p$R$nu, 0.002)
  expect_equal(p$G[[1]]$alpha.V, 1e4 * diag(1))
})

test_that("poisson prior: log-link-scaled R-structure and tight PX working
           prior (not the gaussian copy)", {
  p <- BACE:::.make_prior(n_rand = 1, type = "poisson")
  expect_equal(p$R$V, 0.1)
  expect_equal(p$R$nu, 1)
  expect_equal(p$G[[1]]$alpha.V, 1 * diag(1))
})
