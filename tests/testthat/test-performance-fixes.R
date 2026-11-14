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

  # sample = FALSE: median-based point estimate, capped. The cap binding is
  # the whole point of this test, so it must also SAY so -- truncating the
  # posterior predictive in silence is what we are guarding against.
  expect_warning(
    pt <- BACE:::.predict_bace(model, dat_prep, response_var = "y",
                               type = "poisson", sample = FALSE),
    regexp = "rate guard bound")
  expect_true(all(pt$pred_values <= cap))

  # sample = TRUE: the rate is capped at 11, so a count anywhere near the old
  # rpois(1e6) regime is impossible (P(Pois(11) > 60) < 1e-24).
  expect_warning(
    BACE:::.predict_bace(model, dat_prep, response_var = "y",
                         type = "poisson", sample = TRUE),
    regexp = "rate guard bound")
  for (i in 1:20) {
    dr <- suppressWarnings(
      BACE:::.predict_bace(model, dat_prep, response_var = "y",
                           type = "poisson", sample = TRUE))
    expect_true(all(dr$pred_values < 100))
  }

  # and a cell whose rates never approach the cap must NOT warn
  liab_ok <- matrix(rnorm(60 * 3, 0.5, 0.3), 60, 3)
  expect_silent(
    BACE:::.predict_bace(mock_liab_model(liab_ok),
                         list(data.frame(y = c(1L, 2L, NA)), NULL),
                         response_var = "y", type = "poisson", sample = TRUE))
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


# ---- 2b. The rate cap is an option (1.0.5): off for experiments, never silent ----

test_that("bace_options(poisson_rate_cap) scales the cap, Inf switches it off, bad values refuse", {
  old <- getOption("BACE.poisson_rate_cap"); on.exit(options(BACE.poisson_rate_cap = old), add = TRUE)
  set.seed(587)
  n_iter <- 60; n_obs <- 6
  liab <- matrix(rnorm(n_iter * n_obs, 0.5, 0.3), n_iter, n_obs)
  liab[, 5] <- 20.4                                  # every draw for cell 5 is an excursion
  model <- mock_liab_model(liab)
  dat_prep <- list(data.frame(y = c(0L, 1L, 2L, 1L, NA, NA)), NULL)

  # default multiplier 3: cap 11, and the message names the multiplier
  expect_equal(bace_options()$poisson_rate_cap, 3)
  expect_warning(pt3 <- BACE:::.predict_bace(model, dat_prep, response_var = "y",
                                             type = "poisson", sample = FALSE),
                 "cap = 11 = 3 \\* max\\(observed\\) \\+ 5")
  expect_true(all(pt3$pred_values <= 11))

  # multiplier 10: cap 25
  bace_options(poisson_rate_cap = 10)
  expect_warning(pt10 <- BACE:::.predict_bace(model, dat_prep, response_var = "y",
                                              type = "poisson", sample = FALSE),
                 "cap = 25 = 10 \\* max")
  expect_true(all(pt10$pred_values <= 25))
  expect_equal(max(pt10$pred_values), 25)

  # Inf: the cap is off -- the excursion cell comes through at exp(20.4), no warning,
  # on the point path and on the draw path (the draw stays finite)
  bace_options(poisson_rate_cap = Inf)
  expect_silent(ptI <- BACE:::.predict_bace(model, dat_prep, response_var = "y",
                                            type = "poisson", sample = FALSE))
  expect_equal(max(ptI$pred_values), round(exp(20.4)))
  expect_silent(dr <- BACE:::.predict_bace(model, dat_prep, response_var = "y",
                                           type = "poisson", sample = TRUE))
  expect_true(all(is.finite(dr$pred_values)))
  expect_gt(max(dr$pred_values), 1e8)
  expect_equal(BACE:::.poisson_rate_cap(c(0, 2), Inf), .Machine$integer.max)
  expect_equal(BACE:::.poisson_rate_cap(c(0, 2), 3), 11)
  expect_equal(BACE:::.poisson_rate_cap(numeric(0), 3), 1e6)

  # validation
  expect_error(bace_options(poisson_rate_cap = 0.5), ">= 1")
  expect_error(bace_options(poisson_rate_cap = "off"), ">= 1")
  expect_error(bace_options(poisson_rate_cap = c(3, 4)), ">= 1")

  # the checkpoint fingerprint carries it, so a resume cannot mix cap settings
  bace_options(poisson_rate_cap = 7)
  expect_equal(BACE:::.ckpt_common()$poisson_rate_cap, 7)
  bace_options(.reset = TRUE)
  expect_equal(BACE:::.ckpt_common()$poisson_rate_cap, 3)
})
