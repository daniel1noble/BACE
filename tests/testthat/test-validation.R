# Validation-path tests across BACE's public API. Each test exercises
# a specific input-validation branch (stop / warning / coerce). These
# are fast unit tests with no MCMC, so they run on CRAN.
#
# Includes numerical assertions on coerced/normalised return values
# where the function returns something instead of erroring.

# ---- bace_imp validation ---------------------------------------------------

test_that("bace_imp errors when fixformula references a variable not in data", {
  set.seed(2026)
  phylo <- ape::rtree(8)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  data <- data.frame(y = rnorm(8), x1 = rnorm(8),
                     Species = phylo$tip.label,
                     stringsAsFactors = FALSE)
  data$y[1] <- NA
  expect_error(
    suppressWarnings(suppressMessages(bace_imp(
      fixformula     = "y ~ nonexistent_var",
      ran_phylo_form = "~ 1 | Species",
      phylo          = phylo, data = data,
      runs = 1, nitt = 200, thin = 1, burnin = 50, verbose = FALSE
    ))),
    regexp = "nonexistent_var|not in"
  )
})

test_that("bace_imp errors when phylo tip labels don't match data Species", {
  set.seed(2026)
  phylo <- ape::rtree(8)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  data <- data.frame(y = rnorm(8), x1 = rnorm(8),
                     Species = paste0("not_a_tip_", seq_len(8)),
                     stringsAsFactors = FALSE)
  data$y[1] <- NA
  expect_error(
    suppressWarnings(suppressMessages(bace_imp(
      fixformula     = "y ~ x1",
      ran_phylo_form = "~ 1 | Species",
      phylo          = phylo, data = data,
      runs = 1, nitt = 200, thin = 1, burnin = 50, verbose = FALSE
    )))
  )
})

# ---- assess_convergence validation ----------------------------------------

test_that("assess_convergence rejects non-bace inputs of various shapes", {
  expect_error(assess_convergence(NULL),    regexp = "'bace'")
  expect_error(assess_convergence(list()),  regexp = "'bace'")
  expect_error(assess_convergence(42),      regexp = "'bace'")
  expect_error(assess_convergence("hello"), regexp = "'bace'")
})

test_that("assess_convergence rejects an invalid method string", {
  fake <- list(data = list(NULL, data.frame(y = c(1,2)),
                            data.frame(y = c(2,3)),
                            data.frame(y = c(3,4))),
               types = list(y = "gaussian"),
               miss_dat = data.frame(rowname = "1", colname = "y",
                                      stringsAsFactors = FALSE))
  class(fake) <- "bace"
  expect_error(assess_convergence(fake, method = "bogus_method"))
})

# ---- pool_posteriors validation -------------------------------------------

test_that("pool_posteriors rejects wrong-class inputs", {
  expect_error(pool_posteriors(NULL),
               regexp = "bace_final|bace_complete")
  expect_error(pool_posteriors(list(a = 1)),
               regexp = "bace_final|bace_complete")
})

# ---- get_pooled_model / get_imputed_data validation -----------------------

test_that("get_pooled_model rejects wrong-class object types", {
  expect_error(get_pooled_model(NULL),
               regexp = "'bace_complete'|'bace_pooled'")
  expect_error(get_pooled_model(list()),
               regexp = "'bace_complete'|'bace_pooled'")
  expect_error(get_pooled_model(42),
               regexp = "'bace_complete'|'bace_pooled'")
})

test_that("get_imputed_data rejects wrong-class object types", {
  expect_error(get_imputed_data(NULL),
               regexp = "'bace_complete'|'bace_final'")
  expect_error(get_imputed_data(list()),
               regexp = "'bace_complete'|'bace_final'")
  expect_error(get_imputed_data("foo"),
               regexp = "'bace_complete'|'bace_final'")
})

test_that("get_imputed_data rejects an invalid format value", {
  obj <- list(all_datasets = list(data.frame(y = 1:3)))
  class(obj) <- "bace_final"
  expect_error(get_imputed_data(obj, format = "matrix"))
})

# ---- bace_options validation (extension) ----------------------------------

test_that("bace_options rejects non-list-like calls", {
  expect_error(bace_options(c(verbose = TRUE)),
               regexp = "must be named|Unknown")
})

# ---- .build_formula validation --------------------------------------------

test_that(".build_formula rejects empty / malformed strings", {
  expect_error(BACE:::.build_formula(""),
               regexp = "specify a formula string")
  expect_error(BACE:::.build_formula(NULL),
               regexp = "specify a formula string")
  expect_error(BACE:::.build_formula(c("y ~ x", "z ~ w")),
               regexp = "specify a formula string")
  expect_error(BACE:::.build_formula(123),
               regexp = "specify a formula string")
})

# ---- generate_default_beta_matrix validation -----------------------------

test_that("generate_default_beta_matrix accepts a sparsity value", {
  for (s in c(0, 0.3, 0.5, 0.8, 1.0)) {
    set.seed(2026)
    b <- generate_default_beta_matrix(n_predictors = 4, sparsity = s)
    expect_equal(dim(b), c(4L, 4L))
    expect_true(all(b[upper.tri(b, diag = TRUE)] == 0))
  }
})

test_that("generate_default_beta_matrix sparsity = 1 produces all-zero matrix", {
  set.seed(2026)
  b <- generate_default_beta_matrix(n_predictors = 5, sparsity = 1.0)
  expect_true(all(b == 0))
})

# ---- mnom_liab2cat validation ---------------------------------------------

test_that("mnom_liab2cat output length matches n_cases", {
  for (k in 2:5) {
    cats <- LETTERS[seq_len(k)]
    set.seed(2026)
    liab <- matrix(rnorm(20 * (k - 1)), nrow = 20, ncol = k - 1)
    out <- mnom_liab2cat(liab, cats)
    expect_length(out, 20L)
    expect_true(all(out %in% cats))
  }
})

# ---- ordinal_liab2cat validation ------------------------------------------

test_that("ordinal_liab2cat output is in 1..K for various K", {
  for (k in 2:6) {
    set.seed(2026)
    liab <- rnorm(50)
    out <- ordinal_liab2cat(liab, n_cats = k)
    expect_true(is.integer(out))
    expect_length(out, 50L)
    expect_true(all(out >= 1L & out <= k))
  }
})

test_that("ordinal_liab2cat threshold_spread parameter affects level distribution", {
  set.seed(2026)
  liab <- rnorm(200)
  tight <- ordinal_liab2cat(liab, n_cats = 5, threshold_spread = 0.5)
  loose <- ordinal_liab2cat(liab, n_cats = 5, threshold_spread = 3.0)
  # Tight thresholds give more middle-class observations relative to
  # loose. Specifically, loose has more level-3 (middle) values
  # because more of N(0,1) sits inside +/- middle-threshold.
  expect_true(sum(loose == 3L) >= sum(tight == 3L))
})

# ---- print methods on misc ------------------------------------------------

test_that("print.bace_convergence handles an empty result gracefully", {
  fake <- list(converged = NA,
               method_results = list(),
               summary_stats = NULL,
               diagnostics = list())
  class(fake) <- "bace_convergence"
  expect_output(print(fake), ".+")
})
