# get_imputation_intervals(): the Rubin t-interval and the model-based predictive interval.
#
# WHY THE SECOND ONE EXISTS (2026-11-22). In the 2026-11 simulation the t-interval covered hidden
# counts 0.933 overall and 0.756 of counts above 20: a symmetric interval on a whole-number,
# right-skewed quantity puts its lower bound below zero and its upper bound too low where the
# big counts are. Quantiles of the models' own posterior predictive draws covered 0.986 and 0.923
# (Shinichi, 24 datasets, 729 hidden counts). MCMCglmm keeps one latent column per data row, in
# data order, NA-response rows included, for gaussian, poisson and threshold families and with
# replicated species rows (checked directly), so the draws can be matched to cells by row.

make_fit <- function() {
  set.seed(11)
  sim <- sim_bace(response_type = "poisson", predictor_types = c("gaussian", "gaussian"),
                  var_names = c("cnt", "x1", "x2"), beta_resp = list(x1 = 0.8, x2 = 0.2),
                  beta_matrix = matrix(0, 2, 2), phylo_signal = c(0.5, 0.5, 0.5),
                  missingness = c(0.25, 0.2, 0), n_cases = 80, n_species = 80)
  fit <- suppressWarnings(suppressMessages(
    bace(fixformula = c("cnt ~ x1 + x2", "x1 ~ cnt + x2"), ran_phylo_form = "~ 1 | Species",
         phylo = sim$tree, data = sim$data, runs = 2, nitt = 800, burnin = 200, thin = 4,
         n_final = 8, skip_conv = TRUE, verbose = FALSE)))
  list(sim = sim, fit = fit)
}

test_that("the predictive interval has the right shape and matches the Rubin mean", {
  skip_on_cran()
  z <- make_fit(); fit <- z$fit
  rub <- get_imputation_intervals(fit)
  set.seed(1)
  pre <- get_imputation_intervals(fit, method = "predictive")

  # same cells, same order, labelled
  expect_identical(pre$row, rub$row)
  expect_identical(pre$variable, rub$variable)
  expect_identical(unique(pre$method), "predictive")
  expect_identical(unique(rub$method), "rubin")

  # counts: whole-number bounds, never negative, ordered around the mean
  cnt <- pre[pre$variable == "cnt", ]
  expect_true(all(cnt$lower == round(cnt$lower)))
  expect_true(all(cnt$upper == round(cnt$upper)))
  expect_true(all(cnt$lower >= 0))
  expect_true(all(cnt$lower <= cnt$mean & cnt$mean <= cnt$upper))

  # continuous: on the data scale, so the two methods estimate the same mean. The Rubin
  # mean is the average of only M = 8 random draws, so it scatters around the predictive
  # mean (which averages hundreds); agreement is to within that noise, not exact.
  x1_r <- rub[rub$variable == "x1", ]; x1_p <- pre[pre$variable == "x1", ]
  expect_gt(stats::cor(x1_r$mean, x1_p$mean), 0.6)
  expect_lt(max(abs(x1_r$mean - x1_p$mean)) / stats::sd(z$sim$complete_data$x1), 1.5)
  # and no systematic offset between the two
  expect_lt(abs(mean(x1_r$mean - x1_p$mean)) / stats::sd(z$sim$complete_data$x1), 0.3)
})

test_that("the predictive draws are matched to the right cells", {
  skip_on_cran()
  # cnt depends strongly on x1 (slope 0.8 on the log scale), so if the latent columns
  # were matched to the wrong rows the predictive mean would not track the truth.
  z <- make_fit(); fit <- z$fit; truth <- z$sim$complete_data
  set.seed(2)
  pre <- get_imputation_intervals(fit, method = "predictive", variable = "cnt")
  rub <- get_imputation_intervals(fit, variable = "cnt")
  # the truth correlation is a sanity floor (0.4-0.9 across seeds for correct code); the
  # agreement with the Rubin mean is the discriminating check -- nine wrong row mappings
  # never exceeded 0.41 on it (review, 2026-11-22)
  expect_gt(stats::cor(log1p(pre$mean), log1p(truth$cnt[pre$row])), 0.3)
  expect_gt(stats::cor(log1p(pre$mean), log1p(rub$mean)), 0.6)   # Rubin mean = 8 draws, noisy
  # and the interval covers the hidden truth at a sensible rate on this one dataset
  covered <- mean(truth$cnt[pre$row] >= pre$lower & truth$cnt[pre$row] <= pre$upper)
  expect_gt(covered, 0.8)
})

test_that("the predictive interval refuses what it cannot do, by name", {
  skip_on_cran()
  z <- make_fit(); fit <- z$fit
  broken <- fit
  broken$final_results$all_models[[1]]$cnt$Liab <- NULL
  expect_error(get_imputation_intervals(broken, method = "predictive", variable = "cnt"),
               "carries no latent values")
  no_models <- fit; no_models$final_results$all_models <- NULL
  expect_error(get_imputation_intervals(no_models, method = "predictive"), "final models")
  expect_error(get_imputation_intervals(fit, method = "nonsense"), "arg")
})

# ---- the arithmetic, pinned on a hand-built object ------------------------------------------
#
# A bace_final with known latent values: two fits, six rows, a count and a continuous variable.
# Everything below is checkable by hand, so a wrong back-transform, a wrong row, a
# deterministic draw in place of a Poisson one, or an ignored cap each fail a line.

make_synthetic <- function(count_latent_row5 = log(16), iters = 400L) {
  n <- 6L
  set.seed(7)
  Lc <- matrix(c(1, 2, 3, 4, 5, 6), iters, n, byrow = TRUE)   # distinct per column, so a wrong row shows
  Lc[, 2] <- log(4); Lc[, 5] <- count_latent_row5             # the two missing count cells
  Lx <- matrix(stats::rnorm(iters * n, 0, 0.5), iters, n)      # z-scale latents, sd 0.5
  Lx[, 3] <- Lx[, 3] + 1; Lx[, 6] <- Lx[, 6] + 3               # the two missing continuous cells
  model <- function(L) structure(list(Liab = L), class = "MCMCglmm")
  d <- data.frame(Species = paste0("s", 1:n), cnt = c(1L, 2L, 3L, 2L, 4L, 5L), x = c(9, 11, 12, 12, 10, 14))
  fits <- list(cnt = model(Lc), x = model(Lx))
  structure(list(all_models = list(fits, fits), all_datasets = list(d, d),
                 scaling = list(list(x = list(x = list(mean = 10, sd = 2))),
                                list(x = list(x = list(mean = 10, sd = 2)))),
                 types = list(cnt = "poisson", x = "gaussian"),
                 miss_dat = data.frame(row = c(2L, 5L, 3L, 6L), colname = c("cnt", "cnt", "x", "x")),
                 phylo_ran = list(cluster = "Species"), n_final = 2L),
            class = "bace_final")
}

test_that("predictive: exact back-transform, Poisson spread, right rows, whole-number bounds", {
  fin <- make_synthetic()
  set.seed(3)
  pre <- get_imputation_intervals(fin, method = "predictive")
  pre <- pre[order(pre$variable, pre$row), ]
  cnt <- pre[pre$variable == "cnt", ]; x <- pre[pre$variable == "x", ]

  # continuous: z-scale latent N(1, 0.5) and N(3, 0.5), scaling mean 10 sd 2 -> data scale
  # means 12 and 16, sd 1, 95% bounds about 10.04-13.96 and 14.04-17.96 (800 draws per cell)
  expect_equal(x$row, c(3L, 6L))
  expect_equal(x$mean, c(12, 16), tolerance = 0.03)
  expect_equal(x$sd,   c(1, 1),   tolerance = 0.1)
  expect_equal(x$lower, c(10.04, 14.04), tolerance = 0.03)
  expect_equal(x$upper, c(13.96, 17.96), tolerance = 0.03)

  # counts: rates 4 and 16, so Poisson means 4 and 16 with sd 2 and 4, bounds near qpois
  expect_equal(cnt$row, c(2L, 5L))
  expect_equal(cnt$mean, c(4, 16), tolerance = 0.05)
  expect_equal(cnt$sd,   c(2, 4),  tolerance = 0.1)
  expect_true(all(cnt$lower == round(cnt$lower)) && all(cnt$upper == round(cnt$upper)))
  expect_true(all(abs(cnt$lower - stats::qpois(0.025, c(4, 16))) <= 1))
  expect_true(all(abs(cnt$upper - stats::qpois(0.975, c(4, 16))) <= 1))
  expect_identical(unique(pre$method), "predictive")
})

test_that("predictive: a count rate above the cap is capped as the imputations were, and says so", {
  old <- getOption("BACE.poisson_rate_cap"); on.exit(options(BACE.poisson_rate_cap = old), add = TRUE)
  # observed counts are 1,3,2,5 (rows never missing) -> cap = 3 * 5 + 5 = 20; row 5's rate is 40
  fin <- make_synthetic(count_latent_row5 = log(40))
  bace_options(poisson_rate_cap = 3)
  set.seed(4)
  expect_warning(pre <- get_imputation_intervals(fin, method = "predictive", variable = "cnt"),
                 "cap = 20 = 3 \\* max")
  expect_equal(pre$mean[pre$row == 5], 20, tolerance = 0.05)   # held at the cap, like the imputations
  expect_equal(pre$mean[pre$row == 2], 4,  tolerance = 0.05)   # the uncapped cell is untouched
  # cap off: the rate comes through at 40, and nothing warns
  bace_options(poisson_rate_cap = Inf)
  set.seed(4)
  expect_silent(off <- get_imputation_intervals(fin, method = "predictive", variable = "cnt"))
  expect_equal(off$mean[off$row == 5], 40, tolerance = 0.05)
})

test_that("an empty result still has every documented column", {
  fin <- make_synthetic()
  fin$miss_dat <- fin$miss_dat[0, ]
  expect_message(e <- get_imputation_intervals(fin), "nothing to interval")
  expect_true("method" %in% names(e)); expect_equal(nrow(e), 0L)
})
