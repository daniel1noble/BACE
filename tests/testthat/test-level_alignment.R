# Tests for the .align_levels_to_probs helper that reconciles a
# declared factor-level vector with the column count of a probability
# matrix. The helper guards against the failure mode where MCMCglmm
# fits a threshold / categorical model on fewer levels than the
# response factor declares (heavy subsampling, presence-only binary
# indicators, etc.) and the iterator returns a probability matrix
# with ncol(p) < length(levels_var).

test_that(".align_levels_to_probs returns levels_declared when n_cols matches", {
  expect_identical(
    BACE:::.align_levels_to_probs(c("a","b","c"), n_cols = 3L),
    c("a","b","c")
  )
})

test_that(".align_levels_to_probs returns observed levels when those match n_cols", {
  # Declared: 5 levels; observed: 3; n_cols: 3 -> use observed.
  out <- BACE:::.align_levels_to_probs(
    levels_declared = c("1","2","3","4","5"),
    n_cols          = 3L,
    levels_observed = c("2","3","5")
  )
  expect_identical(out, c("2","3","5"))
})

test_that(".align_levels_to_probs falls back to first n_cols declared with a warning", {
  # Declared: 5 levels; observed: 5 (all observed) but model still
  # produced only 2 columns (e.g. MCMCglmm couldn't identify enough
  # cutpoints from a sparse high-K ordinal). Falls back with a warning.
  expect_warning(
    out <- BACE:::.align_levels_to_probs(
      levels_declared = c("1","2","3","4","5"),
      n_cols          = 2L,
      levels_observed = c("1","2","3","4","5"),
      response_var    = "diet_breadth"
    ),
    regexp = "diet_breadth"
  )
  expect_identical(out, c("1","2"))
})

test_that(".align_levels_to_probs errors on impossible n_cols > declared", {
  expect_error(
    BACE:::.align_levels_to_probs(c("a","b"), n_cols = 5L,
                                  response_var = "x"),
    regexp = "5 columns but only 2 declared"
  )
})

test_that(".impute_levels handles ncol(pred_prob) < length(levels_var)", {
  # 4-column declared factor, 2-column predicted probabilities (sparse
  # ordinal fit). Without the alignment, .impute_levels would error;
  # with it, we get the first 2 declared levels and a warning.
  set.seed(1)
  pred_prob <- matrix(c(0.7, 0.3, 0.2, 0.8, 0.5, 0.5),
                       nrow = 3, byrow = TRUE,
                       dimnames = list(NULL, c("Level_1","Level_2")))
  expect_warning(
    out <- BACE:::.impute_levels(pred_prob,
                                  levels_var = c("a","b","c","d"),
                                  sample = FALSE),
    regexp = "Falling back"
  )
  # argmax per row -> column 1 (a), column 2 (b), tie broken to col 1 (a).
  expect_identical(out, c("a","b","a"))
})

test_that(".impute_levels uses level-named columns when categorical-style", {
  # When pred_prob's colnames ARE valid level names (the categorical
  # case), .impute_levels prefers name-based alignment over positional.
  pred_prob <- matrix(c(0.1, 0.9, 0.6, 0.4),
                       nrow = 2, byrow = TRUE,
                       dimnames = list(NULL, c("yes","no")))
  out <- BACE:::.impute_levels(pred_prob,
                               levels_var = c("no","maybe","yes"),
                               sample = FALSE)
  # Argmax row 1 col 2 -> "no". Argmax row 2 col 1 -> "yes".
  expect_identical(out, c("no","yes"))
})
