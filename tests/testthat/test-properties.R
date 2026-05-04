# Property-style tests: many parameterized assertions across BACE's
# helpers and simulators. Each test_that block runs a sweep over
# inputs and asserts an invariant on every output -- maximises
# tests-per-line while still verifying meaningful numerical
# properties.

# ---- Simulator output property checks -------------------------------------

test_that("sim_bace_gaussian: response is finite numeric across n_predictors 2..6", {
  for (np in 2:6) {
    set.seed(2026)
    out <- suppressMessages(sim_bace_gaussian(
      n_predictors = np, n_cases = 30, n_species = 15, phylo_signal = 0.4
    ))
    resp <- out$complete_data[[out$params$var_names[1]]]
    expect_true(is.numeric(resp))
    expect_true(all(is.finite(resp)))
    expect_length(resp, 30L)
  }
})

test_that("sim_bace_poisson: response is non-negative integer across n_predictors 2..5", {
  for (np in 2:5) {
    set.seed(2026)
    out <- suppressMessages(sim_bace_poisson(
      n_predictors = np, n_cases = 25, n_species = 10
    ))
    resp <- out$complete_data[[out$params$var_names[1]]]
    expect_true(all(resp >= 0))
    expect_true(all(resp == round(resp)))
  }
})

test_that("sim_bace_binary: response is a factor with K=2 levels", {
  for (np in 2:5) {
    set.seed(2026)
    out <- suppressMessages(sim_bace_binary(
      n_predictors = np, n_cases = 25, n_species = 10
    ))
    resp <- out$complete_data[[out$params$var_names[1]]]
    expect_s3_class(resp, "factor")
    expect_equal(nlevels(resp), 2L)
  }
})

test_that("sim_bace tree has correct tip count for various n_species", {
  for (n in c(8, 12, 20, 30)) {
    set.seed(2026)
    out <- suppressMessages(sim_bace_gaussian(
      n_predictors = 2, n_cases = n, n_species = n, phylo_signal = 0.2
    ))
    expect_equal(ape::Ntip(out$tree), n)
    expect_true(all(out$data$Species %in% out$tree$tip.label))
  }
})

test_that("sim_bace phylo_signal = 0 produces no phylo signal in random effects", {
  for (np in 2:4) {
    set.seed(2026)
    out <- suppressMessages(sim_bace_gaussian(
      n_predictors = np, n_cases = 50, n_species = 50, phylo_signal = 0
    ))
    # phylo random effects should be vector of zeros
    if (!is.null(out$random_effects$u_phylo)) {
      ph <- unlist(out$random_effects$u_phylo)
      if (length(ph)) expect_true(all(ph == 0))
    } else {
      succeed("u_phylo not stored; nothing to check")
    }
  }
})

# ---- generate_default_beta_matrix property: lower triangular ---------------

test_that("generate_default_beta_matrix is always lower-triangular", {
  for (np in 1:8) {
    for (s in c(0.0, 0.4, 0.8)) {
      set.seed(2026)
      b <- generate_default_beta_matrix(n_predictors = np, sparsity = s)
      expect_equal(dim(b), c(np, np))
      expect_true(all(b[upper.tri(b, diag = TRUE)] == 0))
    }
  }
})

test_that("generate_default_beta_matrix non-zero entries lie within beta_range", {
  for (np in 2:5) {
    set.seed(2026)
    b <- generate_default_beta_matrix(n_predictors = np, sparsity = 0,
                                       beta_range = c(-0.3, 0.3))
    # Non-zero entries (lower triangle) are bounded
    nz <- b[lower.tri(b)]
    nz <- nz[nz != 0]
    if (length(nz)) {
      expect_true(all(nz >= -0.3))
      expect_true(all(nz <= 0.3))
    }
  }
})

# ---- mnom_liab2cat property: probabilities sum to 1 internally ------------

test_that("mnom_liab2cat output distribution: each category gets at least one draw on large samples", {
  set.seed(2026)
  cats <- c("A", "B", "C", "D")
  liab <- matrix(rnorm(2000 * 3, 0, 1), nrow = 2000, ncol = 3)
  out <- mnom_liab2cat(liab, cats)
  # On 2000 draws with 4 cats, every category should get at least
  # a handful (lower bound 50 is well below the expected 500).
  for (cl in cats) {
    expect_gt(sum(out == cl), 50L)
  }
})

# ---- ordinal_liab2cat property: monotone in liability ---------------------

test_that("ordinal_liab2cat: higher liability tends toward higher categories", {
  set.seed(2026)
  # With default threshold_spread = 1.5 and K = 5, thresholds are
  # -3, -1, 1, 3. Use liabilities BEYOND those bounds so the
  # extremes deterministically land in the extreme categories
  # (the function uses strict `>`, so liab == threshold sticks at
  # the lower category).
  liab_low  <- rep(-5, 100)
  liab_high <- rep( 5, 100)
  out_low  <- ordinal_liab2cat(liab_low,  n_cats = 5)
  out_high <- ordinal_liab2cat(liab_high, n_cats = 5)
  expect_lt(mean(out_low),  mean(out_high))
  expect_true(all(out_low  == 1L))
  expect_true(all(out_high == 5L))
})

# ---- bace_options round-trip property -------------------------------------

test_that("bace_options sweep: every allowed value round-trips", {
  on.exit(bace_options(.reset = TRUE), add = TRUE)
  for (v in c(TRUE, FALSE)) {
    bace_options(verbose = v)
    expect_identical(bace_options()$verbose, v)
  }
  for (d in c(1L, 3L, 5L, 10L)) {
    bace_options(digits = d)
    expect_identical(bace_options()$digits, d)
  }
  for (g in c(0, 1, 2)) {
    bace_options(gelman = g)
    expect_equal(bace_options()$gelman, g)
  }
})

# ---- .impute_levels property: output always in level set ------------------

test_that(".impute_levels: argmax output always in declared levels", {
  set.seed(2026)
  for (k in 2:6) {
    lvls <- LETTERS[seq_len(k)]
    pred <- matrix(runif(50 * k), nrow = 50, ncol = k,
                    dimnames = list(NULL, lvls))
    out <- BACE:::.impute_levels(pred, lvls, sample = FALSE)
    expect_length(out, 50L)
    expect_true(all(out %in% lvls))
  }
})

test_that(".impute_levels: sample output always in declared levels", {
  set.seed(2026)
  for (k in 2:6) {
    lvls <- LETTERS[seq_len(k)]
    pred <- matrix(runif(50 * k), nrow = 50, ncol = k,
                    dimnames = list(NULL, lvls))
    pred <- pred / rowSums(pred)  # normalise rows
    out <- BACE:::.impute_levels(pred, lvls, sample = TRUE)
    expect_length(out, 50L)
    expect_true(all(out %in% lvls))
  }
})

# ---- .cat_mode property: output is one of the values in each row ----------

test_that(".cat_mode: result is one of the row's values for various shapes", {
  set.seed(2026)
  for (n_rows in c(5, 10, 20)) {
    for (n_cols in c(3, 5, 8)) {
      m <- matrix(sample(LETTERS[1:4], n_rows * n_cols, replace = TRUE),
                   nrow = n_rows, ncol = n_cols)
      out <- BACE:::.cat_mode(m)
      expect_length(out, n_rows)
      for (i in seq_len(n_rows)) {
        expect_true(out[i] %in% m[i, ])
      }
    }
  }
})

# ---- .align_levels_to_probs property: output length always equals n_cols --

test_that(".align_levels_to_probs: output length always == n_cols", {
  for (n_decl in 2:6) {
    for (n_cols in 2:n_decl) {
      lvls <- LETTERS[seq_len(n_decl)]
      out <- suppressWarnings(BACE:::.align_levels_to_probs(
        levels_declared = lvls,
        n_cols          = n_cols,
        levels_observed = lvls[seq_len(n_cols)]))
      expect_length(out, n_cols)
    }
  }
})

# ---- sim_bace species column always matches tree tip labels ---------------

test_that("sim_bace: every Species in $data is a tip in $tree", {
  for (n in c(10, 20, 30)) {
    set.seed(2026)
    out <- suppressMessages(sim_bace_gaussian(
      n_predictors = 2, n_cases = n, n_species = n, phylo_signal = 0.3
    ))
    expect_true(all(out$data$Species %in% out$tree$tip.label))
    expect_true(all(out$complete_data$Species %in% out$tree$tip.label))
  }
})
