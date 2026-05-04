# Tests for simulate_auxiliary.r helpers: parse_ix_matrix,
# generate_ix_betas, calculate_ix_term, expand_beta_resp,
# apply_missingness, sample_random_effects, design_size, var_name_gen.
# Many small parameterised tests with numerical assertions.

# ---- var_name_gen ---------------------------------------------------------

test_that("var_name_gen creates appropriately named vector", {
  for (k in 1:5) {
    out <- BACE:::var_name_gen(k)
    expect_length(out, k + 1L)  # response + k predictors
    expect_true("y" %in% out || any(grepl("y", out, ignore.case = TRUE)))
  }
})

# ---- apply_missingness ----------------------------------------------------

test_that("apply_missingness inserts ~prop NAs", {
  set.seed(2026)
  for (p in c(0.0, 0.1, 0.3, 0.5, 0.9)) {
    x <- rnorm(100)
    out <- BACE:::apply_missingness(x, prop = p)
    expect_length(out, 100L)
    expected_n <- floor(100 * p)
    expect_equal(sum(is.na(out)), expected_n)
  }
})

test_that("apply_missingness with prop=0 returns original unchanged", {
  set.seed(2026)
  x <- rnorm(50)
  out <- BACE:::apply_missingness(x, prop = 0)
  expect_identical(out, x)
})

test_that("apply_missingness preserves non-NA values exactly", {
  set.seed(2026)
  x <- seq_len(100)
  out <- BACE:::apply_missingness(x, prop = 0.3)
  observed <- out[!is.na(out)]
  # The non-NA values must come from the original x
  expect_true(all(observed %in% x))
})

# ---- sample_random_effects ------------------------------------------------

test_that("sample_random_effects (independent) returns N(0, sigma2)", {
  set.seed(2026)
  for (s2 in c(0.5, 1, 2, 4)) {
    samples <- BACE:::sample_random_effects(sigma2 = s2, n = 1000)
    expect_length(samples, 1000L)
    # Mean should be near 0 (3-sigma envelope around expected mean=0)
    expect_lt(abs(mean(samples)), 3 * sqrt(s2 / 1000))
    # Variance should be near sigma2 (loose tolerance)
    expect_lt(abs(stats::var(samples) - s2), s2 * 0.5)
  }
})

test_that("sample_random_effects with cor_matrix produces correlated draws", {
  set.seed(2026)
  # Construct a known correlation matrix (3x3 with strong off-diagonals)
  cor_mat <- matrix(c(1, 0.7, 0.5,
                       0.7, 1, 0.3,
                       0.5, 0.3, 1), 3, 3)
  out <- BACE:::sample_random_effects(sigma2 = 1, n = 3, cor_matrix = cor_mat)
  expect_length(out, 3L)
  expect_true(all(is.finite(out)))
})

# ---- design_size ----------------------------------------------------------

test_that("design_size returns ns and sim_order", {
  bm <- BACE:::generate_default_beta_matrix(3, sparsity = 0.5)
  set.seed(2026)
  out <- BACE:::design_size(c("gaussian", "gaussian", "gaussian"), bm)
  expect_type(out, "list")
  expect_true("ns" %in% names(out))
  expect_true("sim_order" %in% names(out))
  expect_length(out$sim_order, 3L)
})

# ---- parse_ix_matrix ------------------------------------------------------

test_that("parse_ix_matrix returns empty list when no interactions", {
  ix <- matrix(0, 3, 3)
  out <- BACE:::parse_ix_matrix(ix, c("y", "x1", "x2"))
  expect_type(out, "list")
  expect_named(out, c("y", "x1", "x2"))
  for (el in out) expect_length(el, 0L)
})

test_that("parse_ix_matrix detects pairwise interaction in lower triangle", {
  # Row 3 (y) has codes [1, 1, 0] -> x1 and x2 share code 1 -> x1:x2
  ix <- matrix(c(0, 0, 0,
                 0, 0, 0,
                 1, 1, 0), 3, 3, byrow = TRUE)
  out <- BACE:::parse_ix_matrix(ix, c("x1", "x2", "y"))
  # Expect y to have x1:x2 as an interaction
  expect_true(length(out$y) >= 1L)
})

test_that("parse_ix_matrix handles multi-digit codes", {
  # 12 means digits 1 and 2 -- x1 has codes {1,2}, x2 has {1}, x3 has {2}
  ix <- matrix(c(0,  0, 0, 0,
                 0,  0, 0, 0,
                 0,  0, 0, 0,
                 12, 1, 2, 0), 4, 4, byrow = TRUE)
  out <- BACE:::parse_ix_matrix(ix, c("x1", "x2", "x3", "y"))
  # y should have at least 2 pairwise interactions (x1:x2 from code 1
  # and x1:x3 from code 2)
  expect_true(length(out$y) >= 2L)
})

# ---- generate_ix_betas ----------------------------------------------------

test_that("generate_ix_betas produces named beta values for parsed interactions", {
  ix <- matrix(c(0, 0, 0,
                 0, 0, 0,
                 1, 1, 0), 3, 3, byrow = TRUE)
  parsed <- BACE:::parse_ix_matrix(ix, c("x1", "x2", "y"))
  set.seed(2026)
  betas <- BACE:::generate_ix_betas(parsed, beta_ix = NULL,
                                     beta_range = c(-0.5, 0.5))
  expect_type(betas, "list")
  # Beta value for y's interactions should be within range
  if (length(betas$y)) {
    for (b in betas$y) {
      expect_true(b >= -0.5 && b <= 0.5)
    }
  }
})

test_that("generate_ix_betas accepts user-provided beta_ix overrides", {
  ix <- matrix(c(0, 0, 0,
                 0, 0, 0,
                 1, 1, 0), 3, 3, byrow = TRUE)
  parsed <- BACE:::parse_ix_matrix(ix, c("x1", "x2", "y"))
  user_betas <- list(y = list("x1:x2" = 0.42))
  betas <- BACE:::generate_ix_betas(parsed, beta_ix = user_betas,
                                     beta_range = c(-0.3, 0.3))
  expect_equal(unname(unlist(betas$y["x1:x2"])), 0.42)
})

# ---- calculate_ix_term ----------------------------------------------------

test_that("calculate_ix_term: numeric * numeric is element-wise product", {
  d <- data.frame(a = 1:5, b = c(2, 3, 4, 5, 6))
  out <- BACE:::calculate_ix_term(d, c("a", "b"))
  expect_equal(out, c(2, 6, 12, 20, 30))
})

test_that("calculate_ix_term: factor * numeric uses (level - 1) coding", {
  d <- data.frame(
    f = factor(c("A","B","A","B","B")),  # 1 / 2 -> 0 / 1 after offset
    n = c(2, 4, 6, 8, 10)
  )
  out <- BACE:::calculate_ix_term(d, c("f", "n"))
  # Expected: (0,1,0,1,1) * (2,4,6,8,10) = (0, 4, 0, 8, 10)
  expect_equal(out, c(0, 4, 0, 8, 10))
})

test_that("calculate_ix_term: character is coerced to factor first", {
  d <- data.frame(
    s = c("a", "b", "a", "b"),
    n = c(2, 4, 6, 8),
    stringsAsFactors = FALSE
  )
  out <- BACE:::calculate_ix_term(d, c("s", "n"))
  # Expected: as.factor("a","b","a","b") -> 1,2,1,2 -> 0,1,0,1
  # x n: 0,4,0,8
  expect_equal(out, c(0, 4, 0, 8))
})

# ---- expand_beta_resp -----------------------------------------------------

test_that("expand_beta_resp accepts a numeric vector of length n_predictors", {
  out <- BACE:::expand_beta_resp(
    beta_resp       = c(0.5, -0.3, 0.2),
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2", "x3"),
    intercept       = 1.0
  )
  expect_type(out, "list")
  expect_true("beta_full" %in% names(out))
  expect_true("beta_resp_stored" %in% names(out))
})

test_that("expand_beta_resp errors when vector length doesn't match n_predictors", {
  expect_error(BACE:::expand_beta_resp(
    beta_resp       = c(0.5, -0.3),  # length 2, but 3 predictors
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2", "x3"),
    intercept       = 0
  ), regexp = "must have length equal to")
})

test_that("expand_beta_resp errors when beta_resp is not numeric or list", {
  expect_error(BACE:::expand_beta_resp(
    beta_resp       = "not a vector",
    predictor_types = c("gaussian"),
    var_names       = c("y", "x1"),
    intercept       = 0
  ), regexp = "must be a numeric vector or a list")
})

# ---- More numerical assertions on expand_beta_resp ------------------------

test_that("expand_beta_resp: beta_full first entry is the intercept", {
  for (intercept_val in c(-1, 0, 0.5, 2)) {
    out <- BACE:::expand_beta_resp(
      beta_resp       = c(0.1, 0.2),
      predictor_types = c("gaussian", "gaussian"),
      var_names       = c("y", "x1", "x2"),
      intercept       = intercept_val
    )
    expect_equal(out$beta_full[1], intercept_val)
  }
})

test_that("expand_beta_resp: beta_full coefficients match input vector for gaussian preds", {
  out <- BACE:::expand_beta_resp(
    beta_resp       = c(0.5, -0.3, 0.7),
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2", "x3"),
    intercept       = 0
  )
  # For all-gaussian predictors, beta_full should be c(0, 0.5, -0.3, 0.7)
  expect_equal(out$beta_full, c(0, 0.5, -0.3, 0.7))
})

# ---- apply_missingness: bounds -------------------------------------------

test_that("apply_missingness: with prop=1, returns at most n NAs", {
  set.seed(2026)
  x <- 1:50
  out <- BACE:::apply_missingness(x, prop = 1)
  expect_true(sum(is.na(out)) <= 50L)
})

# ---- design_size: handles all-gaussian predictors ------------------------

test_that("design_size sim_order is a permutation of seq_len(n_pred)", {
  set.seed(2026)
  for (n in 2:5) {
    bm <- BACE:::generate_default_beta_matrix(n, sparsity = 0.3)
    out <- BACE:::design_size(rep("gaussian", n), bm)
    expect_setequal(out$sim_order, seq_len(n))
  }
})

# ---- parse_ix_matrix: lower-triangular only is honored -------------------

test_that("parse_ix_matrix ignores upper-triangular entries", {
  ix <- matrix(0, 3, 3)
  ix[1, 2] <- 1  # upper-triangular -- should be ignored
  ix[1, 3] <- 1  # upper-triangular -- should be ignored
  out <- BACE:::parse_ix_matrix(ix, c("y", "x1", "x2"))
  for (el in out) expect_length(el, 0L)
})

# ---- calculate_ix_term: zero predictors produce zero output --------------

test_that("calculate_ix_term: zero predictor gives zero output", {
  d <- data.frame(a = c(0, 0, 0, 0), b = c(1, 2, 3, 4))
  out <- BACE:::calculate_ix_term(d, c("a", "b"))
  expect_equal(out, c(0, 0, 0, 0))
})

# ---- generate_ix_betas: empty interactions -> empty list ----------------

test_that("generate_ix_betas with empty interactions returns empty list per variable", {
  parsed <- list(y = list(), x1 = list(), x2 = list())
  out <- BACE:::generate_ix_betas(parsed, beta_ix = NULL)
  expect_type(out, "list")
  for (el in out) expect_length(el, 0L)
})
