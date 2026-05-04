# Direct unit tests for the categorical/threshold prediction helpers
# in model_functions.R: .pred_cat, .pred_cat_iter, .pred_threshold,
# .pred_threshold_iter. These iter variants are not exercised through
# the bace_final_imp() pipeline (which always passes formula+data_full
# and so routes through the *_forward_iter variants), so direct mock
# tests are the only way to exercise them.

# ---- Mock builders -------------------------------------------------------

mock_cat_model <- function(n_obs = 6, n_traits = 2, n_iter = 50, seed = 2026) {
  set.seed(seed)
  # Liab: n_iter rows x (n_traits * n_obs) cols
  Liab <- matrix(rnorm(n_iter * n_traits * n_obs, 0, 0.5),
                 nrow = n_iter, ncol = n_traits * n_obs)
  Liab <- coda::mcmc(Liab)
  # Sol: first n_traits cols are intercepts named "trait{Y}.{Level}"
  level_names <- paste0("Level", seq_len(n_traits))
  sol_cols <- paste0("traitY.", level_names)
  Sol <- matrix(rnorm(n_iter * n_traits, 0, 0.5),
                nrow = n_iter, ncol = n_traits)
  colnames(Sol) <- sol_cols
  Sol <- coda::mcmc(Sol)
  # Fixed and Residual structure as MCMCglmm needs
  list(
    Sol      = Sol,
    Liab     = Liab,
    Residual = list(nrl = n_obs),
    Fixed    = list(nll = 0, nfl = n_traits)
  )
}

mock_threshold_model <- function(n_obs = 8, n_iter = 50, n_cuts = 0,
                                  seed = 2026) {
  set.seed(seed)
  Liab <- matrix(rnorm(n_iter * n_obs, 0, 0.5), nrow = n_iter, ncol = n_obs)
  Liab <- coda::mcmc(Liab)
  Sol  <- coda::mcmc(matrix(rnorm(n_iter * 1, 0, 0.5), n_iter, 1))
  colnames(Sol) <- "Intercept"

  if (n_cuts > 0) {
    # CP must be increasing per row (cumulative cut-points).
    # Use rowwise sorted positive increments.
    incs <- matrix(abs(rnorm(n_iter * n_cuts, 0.5, 0.1)),
                   nrow = n_iter, ncol = n_cuts)
    CP <- t(apply(incs, 1, cumsum))
    CP <- coda::mcmc(CP)
  } else {
    CP <- NULL
  }

  list(
    Sol      = Sol,
    Liab     = Liab,
    CP       = CP,
    Residual = list(nrl = n_obs),
    Fixed    = list(nll = 0, nfl = 1)
  )
}

# ---- .pred_cat -----------------------------------------------------------

test_that(".pred_cat returns a data frame with K+1 columns", {
  m <- mock_cat_model(n_obs = 6, n_traits = 3, n_iter = 50)
  out <- BACE:::.pred_cat(m, baseline_name = "Baseline")
  expect_s3_class(out, "data.frame")
  expect_equal(ncol(out), 3 + 1L)  # K traits + baseline
  expect_equal(nrow(out), 6L)
})

test_that(".pred_cat probabilities sum to ~1 per row", {
  m <- mock_cat_model(n_obs = 8, n_traits = 2, n_iter = 100)
  out <- BACE:::.pred_cat(m, baseline_name = "Baseline")
  row_sums <- rowSums(out)
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

test_that(".pred_cat all probabilities in [0, 1]", {
  m <- mock_cat_model(n_obs = 6, n_traits = 3, n_iter = 100)
  out <- BACE:::.pred_cat(m)
  expect_true(all(out >= 0 & out <= 1))
})

test_that(".pred_cat columns are alphabetically ordered", {
  m <- mock_cat_model(n_obs = 5, n_traits = 2, n_iter = 50)
  out <- BACE:::.pred_cat(m, baseline_name = "Apple")  # A precedes Level
  cols <- colnames(out)
  expect_equal(cols, sort(cols))
})

# ---- .pred_cat_iter ------------------------------------------------------

test_that(".pred_cat_iter returns same structure as .pred_cat at one iteration", {
  m <- mock_cat_model(n_obs = 5, n_traits = 2, n_iter = 50)
  out <- BACE:::.pred_cat_iter(m, baseline_name = "Baseline", iteration = 1)
  expect_s3_class(out, "data.frame")
  expect_equal(ncol(out), 2 + 1L)
  expect_equal(nrow(out), 5L)
})

test_that(".pred_cat_iter probabilities sum to ~1 per row", {
  m <- mock_cat_model(n_obs = 6, n_traits = 3, n_iter = 50)
  out <- BACE:::.pred_cat_iter(m, baseline_name = "Baseline", iteration = 5)
  row_sums <- rowSums(out)
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

test_that(".pred_cat_iter at different iterations produces different probabilities", {
  m <- mock_cat_model(n_obs = 4, n_traits = 2, n_iter = 50)
  o1 <- BACE:::.pred_cat_iter(m, iteration = 1)
  o2 <- BACE:::.pred_cat_iter(m, iteration = 25)
  # Different random draws -> different probabilities
  expect_false(isTRUE(all.equal(as.matrix(o1), as.matrix(o2))))
})

test_that(".pred_cat_iter sweep across iterations matches mean of full chain", {
  m <- mock_cat_model(n_obs = 4, n_traits = 2, n_iter = 80)
  full_out <- BACE:::.pred_cat(m)
  # Per-iter probabilities averaged should approximate full posterior mean
  per_iter_sum <- matrix(0, nrow = 4, ncol = 3)
  for (i in seq_len(80)) {
    per_iter_sum <- per_iter_sum +
      as.matrix(BACE:::.pred_cat_iter(m, iteration = i))
  }
  per_iter_mean <- per_iter_sum / 80
  # Should match full_out closely (both are means over the same chain)
  expect_equal(as.numeric(per_iter_mean), as.numeric(as.matrix(full_out)),
               tolerance = 1e-8)
})

# ---- .pred_threshold (binary, no CP) --------------------------------------

test_that(".pred_threshold binary (no CP) returns 2-column probability frame", {
  m <- mock_threshold_model(n_obs = 6, n_iter = 50, n_cuts = 0)
  out <- BACE:::.pred_threshold(m)
  expect_s3_class(out, "data.frame")
  expect_equal(ncol(out), 2L)
  expect_equal(nrow(out), 6L)
})

test_that(".pred_threshold binary probabilities sum to 1", {
  m <- mock_threshold_model(n_obs = 6, n_iter = 100, n_cuts = 0)
  out <- BACE:::.pred_threshold(m)
  row_sums <- rowSums(out)
  expect_true(all(abs(row_sums - 1) < 1e-8))
})

test_that(".pred_threshold default level_names are Level_1, Level_2, ...", {
  # n_cuts cut-points -> ncol(cbind(0, CP)) = n_cuts+1 boundaries -> n_cuts+2 levels
  m <- mock_threshold_model(n_obs = 4, n_iter = 50, n_cuts = 2)
  out <- BACE:::.pred_threshold(m, level_names = NULL)
  expect_equal(colnames(out),
               c("Level_1", "Level_2", "Level_3", "Level_4"))
})

test_that(".pred_threshold with custom level_names labels columns", {
  m <- mock_threshold_model(n_obs = 5, n_iter = 50, n_cuts = 2)
  out <- BACE:::.pred_threshold(m,
                                 level_names = c("low", "med", "high", "top"))
  expect_equal(colnames(out), c("low", "med", "high", "top"))
})

# ---- .pred_threshold_iter ------------------------------------------------

test_that(".pred_threshold_iter binary returns 2-column probability frame", {
  m <- mock_threshold_model(n_obs = 5, n_iter = 50, n_cuts = 0)
  out <- BACE:::.pred_threshold_iter(m, iteration = 3)
  expect_s3_class(out, "data.frame")
  expect_equal(ncol(out), 2L)
  expect_equal(nrow(out), 5L)
})

test_that(".pred_threshold_iter ordinal (CP present) returns K+1 columns", {
  # n_cuts cut-points -> n_cuts+2 levels
  m <- mock_threshold_model(n_obs = 5, n_iter = 50, n_cuts = 3)
  out <- BACE:::.pred_threshold_iter(m, iteration = 10)
  expect_equal(ncol(out), 3 + 2L)
})

test_that(".pred_threshold_iter probabilities sum to 1 per row", {
  m <- mock_threshold_model(n_obs = 6, n_iter = 50, n_cuts = 2)
  out <- BACE:::.pred_threshold_iter(m, iteration = 5)
  row_sums <- rowSums(out)
  expect_true(all(abs(row_sums - 1) < 1e-8))
})

test_that(".pred_threshold_iter sweep matches .pred_threshold mean", {
  # n_cuts=2 -> 4 levels
  m <- mock_threshold_model(n_obs = 4, n_iter = 60, n_cuts = 2)
  full_out <- BACE:::.pred_threshold(m)
  per_iter_sum <- matrix(0, nrow = 4, ncol = 4)
  for (i in seq_len(60)) {
    per_iter_sum <- per_iter_sum +
      as.matrix(BACE:::.pred_threshold_iter(m, iteration = i))
  }
  per_iter_mean <- per_iter_sum / 60
  expect_equal(as.numeric(per_iter_mean), as.numeric(as.matrix(full_out)),
               tolerance = 1e-8)
})

# ---- .impute_levels ------------------------------------------------------

test_that(".impute_levels(sample=FALSE) returns argmax labels", {
  prob <- data.frame(A = c(0.7, 0.1, 0.5),
                     B = c(0.2, 0.6, 0.3),
                     C = c(0.1, 0.3, 0.2))
  out <- BACE:::.impute_levels(prob, levels_var = c("A", "B", "C"),
                                sample = FALSE)
  expect_equal(out, c("A", "B", "A"))
})

test_that(".impute_levels(sample=TRUE) returns valid level labels", {
  set.seed(2026)
  prob <- data.frame(A = c(0.7, 0.1, 0.5),
                     B = c(0.2, 0.6, 0.3),
                     C = c(0.1, 0.3, 0.2))
  out <- BACE:::.impute_levels(prob, levels_var = c("A", "B", "C"),
                                sample = TRUE)
  expect_true(all(out %in% c("A", "B", "C")))
  expect_equal(length(out), 3L)
})
