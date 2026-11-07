# bace_validate() and get_imputation_intervals().
#
# The scoring arithmetic is tested on HAND-MADE inputs with answers you can
# check on paper (the scorer is its own function precisely so this is
# possible); the masking rules are tested on plain data frames; one real
# tiny-MCMC test checks the whole pipeline end to end.

# ---------------------------------------------------------------------------
# Masking rules (CRAN-safe, no MCMC)
# ---------------------------------------------------------------------------

testthat::test_that("masking hides only observed cells and never the cluster", {
  set.seed(1)
  d <- data.frame(Species = paste0("sp", 1:40),
                  x = rnorm(40), f = factor(rep(c("a", "b", "c", "d"), 10)))
  d$x[1:5] <- NA   # already missing

  m <- .mask_observed_cells(d, vars = c("x", "f"), prop = 0.2, mask_seed = 7)

  # every hidden cell was observed in the original data
  for (i in seq_len(nrow(m$cells))) {
    testthat::expect_false(is.na(d[m$cells$row[i], m$cells$variable[i]]))
    testthat::expect_true(is.na(m$data[m$cells$row[i], m$cells$variable[i]]))
  }
  # the truth was recorded correctly (stored as character, so compare with
  # tolerance: printing a double to text keeps 15 significant digits)
  xr <- m$cells[m$cells$variable == "x", ]
  testthat::expect_equal(as.numeric(xr$truth), d$x[xr$row])
  # cells that were ALREADY missing are not in the record
  testthat::expect_false(any(m$cells$row[m$cells$variable == "x"] %in% 1:5))
  # the cluster column was never touched
  testthat::expect_identical(m$data$Species, d$Species)
  # roughly the requested share (x: 35 observed, 20% ~ 7)
  testthat::expect_equal(sum(m$cells$variable == "x"), round(0.2 * 35))
})

testthat::test_that("masking never empties a factor class and skips tiny variables", {
  d <- data.frame(Species = paste0("sp", 1:20),
                  f = factor(c(rep("a", 18), "b", "b")),   # class b: 2 cells
                  tiny = c(rnorm(8), rep(NA, 12)))          # 8 observed < 10

  testthat::expect_warning(
    m <- .mask_observed_cells(d, vars = c("f", "tiny"), prop = 0.5,
                              mask_seed = 1),
    "tiny")
  # 'tiny' skipped entirely
  testthat::expect_false("tiny" %in% m$cells$variable)
  # class b still observed at least once after masking
  testthat::expect_true("b" %in% m$data$f[!is.na(m$data$f)])

  # same seed -> same mask; different seed -> different mask
  m2 <- suppressWarnings(.mask_observed_cells(d, c("f", "tiny"), 0.5, 1))
  m3 <- suppressWarnings(.mask_observed_cells(d, c("f", "tiny"), 0.5, 2))
  testthat::expect_identical(m$cells, m2$cells)
  testthat::expect_false(identical(m$cells, m3$cells))
})

testthat::test_that("masking restores the caller's RNG state", {
  d <- data.frame(Species = paste0("sp", 1:20), x = rnorm(20))
  set.seed(99); before <- .Random.seed   # capture AFTER all other RNG use
  invisible(.mask_observed_cells(d, "x", 0.2, mask_seed = 5))
  testthat::expect_identical(.Random.seed, before)
})

# ---------------------------------------------------------------------------
# Scoring arithmetic on hand-made inputs (CRAN-safe)
# ---------------------------------------------------------------------------

# A small world we can score by hand: 6 hidden cells of a gaussian 'x',
# 20 imputations.
make_numeric_fixture <- function(imputed_value_fun, n_cells = 6L, M = 20L) {
  n <- 30L
  truth_data <- data.frame(Species = paste0("sp", 1:n), x = seq(1, 30))
  rows_v <- 1:n_cells
  masked_data <- truth_data
  masked_data$x[rows_v] <- NA
  datasets <- lapply(seq_len(M), function(m) {
    d <- masked_data
    d$x[rows_v] <- imputed_value_fun(truth_data$x[rows_v], m)
    d
  })
  list(truth_data = truth_data, masked_data = masked_data,
       datasets = datasets, rows_v = rows_v, M = M,
       cells = data.frame(row = rows_v, variable = "x",
                          truth = as.character(truth_data$x[rows_v])))
}

testthat::test_that("zRMSE is 0 for a perfect imputer and ~1 for mean-filling", {
  # perfect: every imputation equals the truth
  fx <- make_numeric_fixture(function(truth, m) truth)
  tab <- .score_validation(fx$truth_data, fx$cells, fx$masked_data,
                           fx$datasets, prob_preds = list(),
                           types = list(x = "gaussian"))
  testthat::expect_equal(tab$zrmse, 0)
  testthat::expect_equal(tab$coverage, 1)           # zero-width but exact
  testthat::expect_equal(tab$n_zero_width, 6L)

  # mean-filler: every imputation is the mean of the remaining observed cells
  fx2 <- make_numeric_fixture(function(truth, m)
    rep(mean(7:30), length(truth)))
  tab2 <- .score_validation(fx2$truth_data, fx2$cells, fx2$masked_data,
                            fx2$datasets, prob_preds = list(),
                            types = list(x = "gaussian"))
  # by construction the mean-filler's zRMSE EQUALS the printed baseline
  testthat::expect_equal(tab2$zrmse, tab2$zrmse_baseline)
  # and the baseline is where it should be by hand:
  # sqrt(mean((mean(7:30) - 1:6)^2)) / sd(7:30)
  testthat::expect_equal(tab2$zrmse,
                         sqrt(mean((mean(7:30) - 1:6)^2)) / stats::sd(7:30))
})

testthat::test_that("the t-interval coverage is ~95% for a calibrated imputer", {
  testthat::skip_on_cran()   # a few seconds of simulation
  set.seed(42)
  n_cells <- 400L; M <- 50L
  mu <- rnorm(n_cells, 0, 2); s <- 1
  truth_data <- data.frame(Species = paste0("sp", 1:n_cells),
                           x = rnorm(n_cells, mu, s))
  masked_data <- truth_data; masked_data$x[] <- NA
  # keep 100 extra observed rows so s_obs exists
  truth_data <- rbind(truth_data,
                      data.frame(Species = paste0("obs", 1:100),
                                 x = rnorm(100, 0, 2)))
  masked_data <- rbind(masked_data,
                       data.frame(Species = paste0("obs", 1:100),
                                  x = truth_data$x[(n_cells + 1):(n_cells + 100)]))
  datasets <- lapply(seq_len(M), function(m) {
    d <- masked_data
    d$x[1:n_cells] <- rnorm(n_cells, mu, s)   # same distribution as truth
    d
  })
  cells <- data.frame(row = 1:n_cells, variable = "x",
                      truth = as.character(truth_data$x[1:n_cells]))
  tab <- .score_validation(truth_data, cells, masked_data, datasets,
                           prob_preds = list(), types = list(x = "gaussian"))
  testthat::expect_gt(tab$coverage, 0.92)
  testthat::expect_lt(tab$coverage, 0.98)
})

testthat::test_that("count reference lines match their formulas", {
  fx <- make_numeric_fixture(function(truth, m) truth + (m %% 3))  # some spread
  tab <- .score_validation(fx$truth_data, fx$cells, fx$masked_data,
                           fx$datasets, prob_preds = list(),
                           types = list(x = "poisson"))
  imp <- sapply(fx$datasets, function(d) d$x[fx$rows_v])
  ybar <- rowMeans(imp); s2 <- apply(imp, 1, stats::var)
  s_obs <- stats::sd(7:30)
  testthat::expect_equal(tab$floor_counting, sqrt(mean(ybar)) / s_obs)
  testthat::expect_equal(tab$floor_calibrated,
                         sqrt((1 + 1 / fx$M) * mean(s2)) / s_obs)
  # the coverage expectation is always SIMULATED for counts (a fixed mean-rate
  # cutoff misreports rate mixtures) under a fixed internal seed, so it is
  # near-nominal here and identical across calls
  testthat::expect_gt(tab$coverage_expected, 0.90)
  testthat::expect_lt(tab$coverage_expected, 0.99)
  testthat::expect_identical(.expected_coverage_counts(ybar, fx$M),
                             .expected_coverage_counts(ybar, fx$M))
})

testthat::test_that("logit-scale scoring matches a hand computation", {
  n <- 20L; M <- 5L
  truth_data <- data.frame(Species = paste0("sp", 1:n),
                           p = seq(0.05, 0.95, length.out = n))
  masked_data <- truth_data; masked_data$p[1:4] <- NA
  datasets <- lapply(seq_len(M), function(m) {
    d <- masked_data; d$p[1:4] <- truth_data$p[1:4] * 0.9; d
  })
  cells <- data.frame(row = 1:4, variable = "p",
                      truth = as.character(truth_data$p[1:4]))
  tab <- .score_validation(truth_data, cells, masked_data, datasets,
                           prob_preds = list(), types = list(p = "gaussian"),
                           scale = c(p = "logit"))
  sq <- function(x) stats::qlogis(pmin(pmax(x, 1e-4), 1 - 1e-4))
  zr_hand <- sqrt(mean((sq(truth_data$p[1:4] * 0.9) - sq(truth_data$p[1:4]))^2)) /
    stats::sd(sq(truth_data$p[5:20]))
  testthat::expect_equal(tab$zrmse, zr_hand)
  testthat::expect_identical(tab$scored_on, "logit")
})

testthat::test_that("class scoring: accuracy, baseline, balanced accuracy, Brier by hand", {
  # 4 hidden cells of a 3-class factor; the data's other rows are mostly "a"
  truth_data <- data.frame(
    Species = paste0("sp", 1:14),
    f = factor(c("a", "a", "b", "c", rep("a", 8), "b", "c")))
  masked_data <- truth_data
  masked_data$f[1:4] <- NA
  M <- 4L
  # imputations: cells 1-2 always right ("a"); cell 3 ("b") votes "a";
  # cell 4 ("c") votes "c" (3 of 4 runs)
  imp_per_run <- list(c("a", "a", "a", "c"), c("a", "a", "a", "c"),
                      c("a", "a", "b", "c"), c("a", "a", "a", "b"))
  datasets <- lapply(imp_per_run, function(vals) {
    d <- masked_data
    d$f[1:4] <- factor(vals, levels = levels(truth_data$f))
    d
  })
  # class probabilities per run: rows = the 4 missing cells, cols = classes.
  # Use the same matrix in every run so the average is easy to check.
  pmat <- matrix(c(0.8, 0.1, 0.1,
                   0.7, 0.2, 0.1,
                   0.5, 0.4, 0.1,
                   0.1, 0.2, 0.7), nrow = 4, byrow = TRUE,
                 dimnames = list(paste0("sp", 1:4), c("a", "b", "c")))
  prob_preds <- lapply(1:M, function(m) list(f = pmat))

  cells <- data.frame(row = 1:4, variable = "f",
                      truth = as.character(truth_data$f[1:4]))
  tab <- .score_validation(truth_data, cells, masked_data, datasets,
                           prob_preds = prob_preds,
                           types = list(f = "categorical"))
  # votes: a, a, a, c; truth: a, a, b, c -> accuracy 3/4
  testthat::expect_equal(tab$accuracy, 0.75)
  # per-class: a -> 1, b -> 0, c -> 1; balanced = 2/3
  testthat::expect_equal(tab$balanced_accuracy, 2 / 3)
  # baseline: commonest observed class is "a"; truth has two "a" -> 0.5
  testthat::expect_equal(tab$accuracy_baseline, 0.5)
  # Brier by hand from pmat rows vs truth indicators
  ind <- rbind(c(1, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  testthat::expect_equal(tab$brier, mean(rowSums((pmat - ind)^2)))
})

testthat::test_that("Brier matching survives extra original NAs and replicated species", {
  # variable f is missing at rows 2 and 5 ALREADY; we hide rows 3 and 7.
  # The probability matrix covers all four missing rows in row order; the
  # scorer must pick out rows 3 and 7 by position (2nd and 4th).
  truth_data <- data.frame(
    Species = c("s1", "s1", "s2", "s3", "s4", "s5", "s2", rep("z", 7)),
    f = factor(c("a", "b", "a", "b", "a", "b", "b", rep("a", 7))))
  masked_data <- truth_data
  masked_data$f[c(2, 5)] <- NA           # originally missing
  masked_data$f[c(3, 7)] <- NA           # hidden by the masker
  M <- 2L
  datasets <- lapply(1:M, function(m) {
    d <- masked_data
    d$f[c(2, 3, 5, 7)] <- factor(c("a", "a", "a", "b"),
                                 levels = levels(truth_data$f))
    d
  })
  pmat <- matrix(c(0.9, 0.1,    # row 2 (original NA)
                   0.6, 0.4,    # row 3 (hidden; truth a)
                   0.8, 0.2,    # row 5 (original NA)
                   0.3, 0.7),   # row 7 (hidden; truth b)
                 nrow = 4, byrow = TRUE,
                 dimnames = list(make.unique(c("s1", "s2", "s4", "s2")),
                                 c("a", "b")))
  prob_preds <- lapply(1:M, function(m) list(f = pmat))
  cells <- data.frame(row = c(3, 7), variable = "f", truth = c("a", "b"))
  tab <- .score_validation(truth_data, cells, masked_data, datasets,
                           prob_preds = prob_preds,
                           types = list(f = "threshold"))
  brier_hand <- mean(c(sum((c(0.6, 0.4) - c(1, 0))^2),
                       sum((c(0.3, 0.7) - c(0, 1))^2)))
  testthat::expect_equal(tab$brier, brier_hand)
})

testthat::test_that("the repetition summary averages scores and reports a standard error", {
  t1 <- data.frame(variable = "x", type = "gaussian", scored_on = "raw",
                   n_masked = 6L, zrmse = 0.4, zrmse_baseline = 1.0, r = 0.8,
                   coverage = 0.95, coverage_expected = 0.95,
                   n_zero_width = 0L, floor_counting = NA_real_,
                   floor_calibrated = NA_real_, accuracy = NA_real_,
                   accuracy_baseline = NA_real_, balanced_accuracy = NA_real_,
                   brier = NA_real_, rep = 1L)
  t2 <- t1; t2$zrmse <- 0.6; t2$rep <- 2L
  s <- .summarise_validation(rbind(t1, t2))
  testthat::expect_equal(s$zrmse, 0.5)
  testthat::expect_equal(s$zrmse_se, stats::sd(c(0.4, 0.6)) / sqrt(2))
  testthat::expect_equal(s$reps, 2L)
})

# ---------------------------------------------------------------------------
# get_imputation_intervals (CRAN-safe, fake object)
# ---------------------------------------------------------------------------

testthat::test_that("get_imputation_intervals matches the hand formula and skips factors", {
  M <- 20L
  base <- data.frame(Species = paste0("sp", 1:10),
                     x = c(NA, NA, 3:10),
                     f = factor(c(NA, rep(c("a", "b"), length.out = 9))))
  datasets <- lapply(seq_len(M), function(m) {
    d <- base
    d$x[1:2] <- c(1 + m / 10, 2 - m / 10)
    d$f[1]   <- factor("a", levels = levels(base$f))
    d
  })
  miss <- data.frame(row = c(1L, 2L, 1L), col = c(2L, 2L, 3L),
                     colname = c("x", "x", "f"))
  fake <- structure(list(all_datasets = datasets, miss_dat = miss,
                         phylo_ran = list(cluster = "Species")),
                    class = "bace_final")
  iv <- get_imputation_intervals(fake)
  testthat::expect_identical(unique(iv$variable), "x")   # factor skipped
  vals <- sapply(datasets, function(d) d$x[1])
  hw <- stats::qt(0.975, M - 1) * sqrt(1 + 1 / M) * stats::sd(vals)
  testthat::expect_equal(iv$mean[1], mean(vals))
  testthat::expect_equal(iv$lower[1], mean(vals) - hw)
  testthat::expect_equal(iv$upper[1], mean(vals) + hw)
  testthat::expect_error(get_imputation_intervals(fake, variable = "f"),
                         "numeric")
})

# ---------------------------------------------------------------------------
# End to end with real (tiny) MCMC
# ---------------------------------------------------------------------------

testthat::test_that("bace_validate runs end to end and reproduces its masks", {
  testthat::skip_on_cran()
  set.seed(3)
  n <- 40
  tr <- ape::rcoal(n)
  S <- ape::vcv(tr, corr = TRUE); L <- chol(S)
  d <- data.frame(Species = tr$tip.label)
  d$x1 <- as.numeric(t(L) %*% rnorm(n)) + rnorm(n, sd = .3)
  d$x2 <- 0.6 * d$x1 + rnorm(n, sd = .4)
  d$x1[sample(n, 6)] <- NA     # some genuinely missing cells too

  v <- suppressWarnings(bace_validate(
    "x1 ~ x2", ran_phylo_form = "~ 1 | Species",
    phylo = tr, data = d, prop = 0.25,
    runs = 2, nitt = 800, burnin = 200, thin = 2, n_final = 10,
    skip_conv = TRUE, verbose = FALSE, seed = 11))

  testthat::expect_s3_class(v, "bace_validation")
  testthat::expect_true(all(c("x1", "x2") %in% v$table$variable))
  testthat::expect_true(all(is.finite(v$table$zrmse)))
  testthat::expect_true(all(v$table$coverage >= 0 & v$table$coverage <= 1))
  # masked cells were observed in the input
  for (i in seq_len(nrow(v$masked_cells))) {
    testthat::expect_false(is.na(d[v$masked_cells$row[i],
                                   v$masked_cells$variable[i]]))
  }
  # the print method runs
  testthat::expect_output(print(v), "zRMSE")

  # same seed -> identical masks. (The scores agree only approximately: the
  # convergence phase is deliberately unseeded -- see ?bace -- so the model
  # runs are not bitwise repeatable even at a fixed seed.)
  v2 <- suppressWarnings(bace_validate(
    "x1 ~ x2", ran_phylo_form = "~ 1 | Species",
    phylo = tr, data = d, prop = 0.25,
    runs = 2, nitt = 800, burnin = 200, thin = 2, n_final = 10,
    skip_conv = TRUE, verbose = FALSE, seed = 11))
  testthat::expect_identical(v$masked_cells, v2$masked_cells)
})

# ---------------------------------------------------------------------------
# The bace_complete method: validate a finished run
# ---------------------------------------------------------------------------

testthat::test_that("bace_validate(fit) reuses the finished run's data, tree and settings", {
  # A fake finished run carrying stored settings; the workhorse is mocked so
  # this tests the HANDOVER, with no MCMC.
  tr <- ape::rcoal(8)
  d0 <- data.frame(x1 = c(NA, rnorm(7)), x2 = rnorm(8),
                   Species = tr$tip.label)
  fake <- structure(list(
    final_results = structure(list(
      formulas = list(as.formula("x1 ~ x2"), as.formula("x2 ~ x1")),
      n_final = 7L, sweeps_requested = 4L,
      trees_used = list(tr),
      settings = list(ran_phylo_form = "~ 1 | Species", nitt = 1234,
                      thin = 3, burnin = 111, species = FALSE,
                      ovr_categorical = TRUE, nitt_cat_mult = 2L,
                      priors = NULL)
    ), class = "bace_final"),
    initial_results = list(data = list(Initial_Data = d0, Iter_1 = d0,
                                       Iter_2 = d0))
  ), class = "bace_complete")

  seen <- NULL
  testthat::local_mocked_bindings(
    .bace_validate_run = function(...) { seen <<- list(...); "ran" },
    .package = "BACE")
  out <- bace_validate(fake, prop = 0.3, seed = 9, verbose = FALSE)
  testthat::expect_identical(out, "ran")
  testthat::expect_identical(seen$nitt, 1234)
  testthat::expect_identical(seen$burnin, 111)
  testthat::expect_identical(seen$ovr_categorical, TRUE)
  testthat::expect_identical(seen$nitt_cat_mult, 2L)
  testthat::expect_identical(seen$n_final, 7L)
  testthat::expect_identical(seen$sweeps, 4L)
  testthat::expect_identical(seen$runs, 2L)          # Iter_1, Iter_2
  testthat::expect_identical(seen$data, d0)          # raw data, NAs intact
  testthat::expect_identical(seen$phylo, tr)
  testthat::expect_identical(seen$ran_phylo_form, "~ 1 | Species")
  testthat::expect_identical(seen$fixformula[[1]], "x1 ~ x2")
  testthat::expect_identical(seen$prop, 0.3)

  # a result from an older BACE that stored no settings: clear error
  old <- fake
  old$final_results$settings <- NULL
  testthat::expect_error(bace_validate(old), "store its settings")
})

testthat::test_that("bace_validate(fit) runs end to end on a real tiny fit", {
  testthat::skip_on_cran()
  set.seed(21)
  n <- 30
  tr <- ape::rcoal(n)
  S <- ape::vcv(tr, corr = TRUE); L <- chol(S)
  d <- data.frame(Species = tr$tip.label)
  d$x1 <- as.numeric(t(L) %*% rnorm(n)) + rnorm(n, sd = .3)
  d$x2 <- 0.6 * d$x1 + rnorm(n, sd = .4)
  d$x1[sample(n, 4)] <- NA

  fit <- suppressWarnings(bace("x1 ~ x2", "~ 1 | Species", phylo = tr,
                               data = d, runs = 2, nitt = 800, burnin = 200,
                               thin = 2, n_final = 8, skip_conv = TRUE,
                               verbose = FALSE))
  v <- suppressWarnings(bace_validate(fit, prop = 0.25, seed = 5,
                                      verbose = FALSE))
  testthat::expect_s3_class(v, "bace_validation")
  testthat::expect_true(all(c("x1", "x2") %in% v$table$variable))
  # the reused settings produced the same number of imputations
  testthat::expect_equal(v$n_final, 8)
})

# ---------------------------------------------------------------------------
# Review regressions (adversarial review 2026-11-07)
# ---------------------------------------------------------------------------

testthat::test_that("tied votes get fractional credit, matching the manuscript", {
  # one hidden cell whose M = 4 imputations split 2-2 between 'a' and 'b':
  # the manuscript rule scores 0.5 whichever of the two is true
  truth_data <- data.frame(Species = paste0("sp", 1:12),
                           f = factor(c("a", rep(c("a", "b"), 5), "b")))
  masked_data <- truth_data; masked_data$f[1] <- NA
  mk <- function(vals) lapply(vals, function(vv) {
    d <- masked_data; d$f[1] <- factor(vv, levels = levels(truth_data$f)); d
  })
  datasets <- mk(c("a", "b", "a", "b"))
  cells <- data.frame(row = 1L, variable = "f", truth = "a")
  tab_a <- .score_validation(truth_data, cells, masked_data, datasets,
                             prob_preds = list(), types = list(f = "categorical"))
  testthat::expect_equal(tab_a$accuracy, 0.5)
  cells_b <- cells; cells_b$truth <- "b"
  truth_b <- truth_data; truth_b$f[1] <- "b"
  tab_b <- .score_validation(truth_b, cells_b, masked_data, datasets,
                             prob_preds = list(), types = list(f = "categorical"))
  testthat::expect_equal(tab_b$accuracy, 0.5)

  # the modal baseline uses the same rule: remaining observed = 5 a, 6 b
  # (no tie), truth 'a' -> baseline 0; with a 2-way observed tie it halves
  testthat::expect_equal(tab_a$accuracy_baseline, 0)
  testthat::expect_equal(.vote_credit(list(c("a", "b")), "a"), 0.5)
  testthat::expect_equal(.vote_credit(list(c("a", "b", "b")), "a"), 0)
})

testthat::test_that("mask seeds survive date-style seeds and reproduce", {
  s1 <- .derive_mask_seeds(20261107, 3)   # would overflow seed * 1000L
  s2 <- .derive_mask_seeds(20261107, 3)
  s3 <- .derive_mask_seeds(20261108, 3)
  testthat::expect_identical(s1, s2)
  testthat::expect_false(identical(s1, s3))
  testthat::expect_length(unique(s1), 3L)
  # and the caller's RNG state is untouched
  set.seed(5); before <- .Random.seed
  invisible(.derive_mask_seeds(20261107, 2))
  testthat::expect_identical(.Random.seed, before)
})

testthat::test_that("list formulas mask only imputed (response) variables", {
  # x2 is a complete predictor with NO formula of its own: it must not be
  # masked (masking it used to crash the pipeline downstream)
  tr <- ape::rcoal(12)
  d <- data.frame(Species = tr$tip.label, x1 = rnorm(12), x2 = rnorm(12))
  d$x1[1:2] <- NA

  seen_data <- NULL
  fake_fit <- NULL
  testthat::local_mocked_bindings(
    bace = function(fixformula, ran_phylo_form, phylo, data, ...) {
      seen_data <<- data
      M <- 3L
      imps <- lapply(1:M, function(m) { dd <- data; dd$x1[is.na(dd$x1)] <- 0; dd })
      structure(list(imputed_datasets = imps, prob_preds = list(),
                     final_results = structure(list(types = list(x1 = "gaussian")),
                                               class = "bace_final")),
                class = "bace_complete")
    },
    .package = "BACE")
  v <- suppressWarnings(bace_validate(list("x1 ~ x2"), "~ 1 | Species",
                                      phylo = tr, data = d, prop = 0.3,
                                      seed = 4, verbose = FALSE))
  testthat::expect_identical(unique(v$masked_cells$variable), "x1")
  # x2 reached the pipeline untouched (no new NAs)
  testthat::expect_false(anyNA(seen_data$x2))
})

testthat::test_that("character vectors of formulas are refused; formula objects dispatch", {
  tr <- ape::rcoal(8)
  d <- data.frame(Species = tr$tip.label, x1 = rnorm(8), x2 = rnorm(8))
  testthat::expect_error(
    bace_validate(c("x1 ~ x2", "x2 ~ x1"), "~ 1 | Species", phylo = tr, data = d),
    "list")

  seen <- NULL
  testthat::local_mocked_bindings(
    .bace_validate_run = function(...) { seen <<- list(...); "ran" },
    .package = "BACE")
  out <- bace_validate(x1 ~ x2, ran_phylo_form = "~ 1 | Species",
                       phylo = tr, data = d)
  testthat::expect_identical(out, "ran")
  testthat::expect_identical(seen$fixformula, "x1 ~ x2")
})

testthat::test_that("unknown scale names and values warn instead of silently scoring raw", {
  tr <- ape::rcoal(12)
  d <- data.frame(Species = tr$tip.label, x1 = rnorm(12), x2 = rnorm(12))
  d$x1[1] <- NA
  testthat::local_mocked_bindings(
    bace = function(fixformula, ran_phylo_form, phylo, data, ...) {
      imps <- lapply(1:3, function(m) { dd <- data; dd$x1[is.na(dd$x1)] <- 0; dd })
      structure(list(imputed_datasets = imps, prob_preds = list(),
                     final_results = structure(list(types = list(x1 = "gaussian")),
                                               class = "bace_final")),
                class = "bace_complete")
    },
    .package = "BACE")
  testthat::expect_warning(
    suppressMessages(bace_validate(list("x1 ~ x2"), "~ 1 | Species", phylo = tr,
                                   data = d, prop = 0.3, seed = 1,
                                   verbose = FALSE, scale = c(X = "logit"))),
    "scale")
})

testthat::test_that("bace_complete method refuses silently-dropped setting overrides", {
  tr <- ape::rcoal(8)
  d0 <- data.frame(x1 = c(NA, rnorm(7)), Species = tr$tip.label)
  fake <- structure(list(
    final_results = structure(list(
      formulas = list(as.formula("x1 ~ 1")), n_final = 5L,
      sweeps_requested = 3L, trees_used = list(tr),
      types = list(x1 = "gaussian"),
      settings = list(ran_phylo_form = "~ 1 | Species", nitt = 100, thin = 1,
                      burnin = 10, species = FALSE, ovr_categorical = FALSE,
                      nitt_cat_mult = 1L, priors = NULL)
    ), class = "bace_final"),
    initial_results = list(data = list(Initial_Data = d0, Iter_1 = d0))
  ), class = "bace_complete")
  testthat::local_mocked_bindings(
    .bace_validate_run = function(...) "ran", .package = "BACE")
  testthat::expect_warning(bace_validate(fake, nitt = 99, verbose = FALSE),
                           "taken? from that run|model settings")
})

testthat::test_that("masking never removes the last two observed numeric cells", {
  d <- data.frame(Species = paste0("sp", 1:12),
                  x = c(rnorm(10), NA, NA))   # 10 observed
  m <- .mask_observed_cells(d, "x", prop = 0.95, mask_seed = 3)
  testthat::expect_lte(sum(m$cells$variable == "x"), 8L)
  testthat::expect_gte(sum(!is.na(m$data$x)), 2L)
})

testthat::test_that("summaries drop NA repetitions per metric instead of propagating", {
  t1 <- data.frame(variable = "f", type = "categorical", scored_on = "class",
                   n_masked = 4L, zrmse = NA_real_, zrmse_baseline = NA_real_,
                   r = NA_real_, coverage = NA_real_,
                   coverage_expected = NA_real_, n_zero_width = NA_integer_,
                   floor_counting = NA_real_, floor_calibrated = NA_real_,
                   accuracy = 0.6, accuracy_baseline = 0.4,
                   balanced_accuracy = 0.5, brier = 0.2, rep = 1L)
  t2 <- t1; t2$brier <- NA_real_; t2$rep <- 2L
  s <- .summarise_validation(rbind(t1, t2))
  testthat::expect_equal(s$brier, 0.2)       # the NA rep is dropped
  testthat::expect_equal(s$accuracy, 0.6)
})

testthat::test_that("get_imputation_intervals degrades clearly at the edges", {
  base <- data.frame(Species = paste0("sp", 1:6), x = c(NA, 2:6))
  # M = 1: mean returned, bounds NA, a message, and no qt() warning leaks
  fake1 <- structure(list(all_datasets = list({ d <- base; d$x[1] <- 1; d }),
                          miss_dat = data.frame(row = 1L, col = 2L, colname = "x"),
                          phylo_ran = list(cluster = "Species")),
                     class = "bace_final")
  testthat::expect_message(iv <- get_imputation_intervals(fake1), "one imputation")
  testthat::expect_equal(iv$mean, 1)
  testthat::expect_true(is.na(iv$lower) && is.na(iv$upper))

  # only factor variables missing: zero-row frame with a message, not NULL
  fakef <- structure(list(
    all_datasets = list(data.frame(Species = "s1", f = factor("a"))),
    miss_dat = data.frame(row = 1L, col = 2L, colname = "f"),
    phylo_ran = list(cluster = "Species")), class = "bace_final")
  testthat::expect_message(iv2 <- get_imputation_intervals(fakef),
                           "No numeric")
  testthat::expect_s3_class(iv2, "data.frame")
  testthat::expect_identical(nrow(iv2), 0L)
})
