# Track A hardening tests (2026-07).
#
# A3 — the two-phase invariant at the prediction level:
#   .predict_bace(sample = FALSE) is a deterministic point estimate (consumes
#   no RNG), so the chained-equations convergence diagnostic can read the
#   imputed values settling across runs. .predict_bace(sample = TRUE) draws
#   from the posterior predictive and so varies between calls. If the
#   convergence-chain call in bace_imp() were ever flipped to sample = TRUE,
#   the imputations would stop being deterministic — these tests pin the
#   contract that makes that regression visible.
#
# A5 — pool_posteriors() is a stacking (mixture) combiner, so a pooled
#   coefficient's posterior variance equals the mean within-imputation
#   variance PLUS the between-imputation variance (law of total variance).
#   This is what makes it propagate imputation uncertainty rather than
#   understate it. The test checks the identity numerically.

# ---- Shared fixture: one fitted gaussian phylo model + its dat_prep ---------

make_gaussian_fit <- function(seed = 123, nitt = 1500, burnin = 300, thin = 2) {
  set.seed(seed)
  tr <- ape::rtree(30)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))

  C <- ape::vcv(tr, corr = TRUE)
  x <- as.numeric(MASS::mvrnorm(1, rep(0, 30), C))
  u <- as.numeric(MASS::mvrnorm(1, rep(0, 30), C))     # phylo random effect
  y <- 1 + 2 * x + u + rnorm(30, 0, 0.5)
  d <- data.frame(y = y, x = x, Species = tr$tip.label,
                  stringsAsFactors = FALSE)

  types <- list(y = "gaussian", x = "gaussian")
  dp    <- BACE:::.data_prep(formula = y ~ x, data = d, types = types,
                             ran_cluster = "Species")
  data_i <- dp[[1]]                                    # scaled, as bace_imp uses
  prior  <- BACE:::.make_prior(n_rand = 1, type = "gaussian")
  model  <- BACE:::.model_fit(
    data = data_i, tree = tr, fixformula = y ~ x,
    randformula = ~ Species, type = "gaussian",   # phylo grouping, as bace_imp passes
    prior = prior, nitt = nitt, thin = thin, burnin = burnin
  )
  list(model = model, dat_prep = dp)
}

test_that("A3: .predict_bace(sample = FALSE) is deterministic across calls", {
  testthat::skip_on_cran()
  f <- suppressWarnings(suppressMessages(make_gaussian_fit()))

  # No set.seed between the two calls: sample = FALSE must consume no RNG.
  p1 <- BACE:::.predict_bace(f$model, f$dat_prep, sample = FALSE,
                             response_var = "y", type = "gaussian")
  p2 <- BACE:::.predict_bace(f$model, f$dat_prep, sample = FALSE,
                             response_var = "y", type = "gaussian")
  expect_identical(p1$pred_values, p2$pred_values)
  expect_true(all(is.finite(p1$pred_values)))
})

test_that("A3: .predict_bace(sample = TRUE) draws vary between calls", {
  testthat::skip_on_cran()
  f <- suppressWarnings(suppressMessages(make_gaussian_fit()))

  set.seed(1)
  q1 <- BACE:::.predict_bace(f$model, f$dat_prep, sample = TRUE,
                             response_var = "y", type = "gaussian")
  q2 <- BACE:::.predict_bace(f$model, f$dat_prep, sample = TRUE,
                             response_var = "y", type = "gaussian")
  # Posterior-predictive draws of a continuous trait differ with prob ~1;
  # this proves the sample flag actually changes behaviour (so sample = FALSE
  # in the convergence chain is a meaningful, load-bearing choice).
  expect_false(isTRUE(all.equal(q1$pred_values, q2$pred_values)))
})

test_that("A5: pool_posteriors variance = within + between (stacking identity)", {
  testthat::skip_on_cran()
  set.seed(2027)
  tr <- ape::rtree(35)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  C <- ape::vcv(tr, corr = TRUE)
  x <- as.numeric(MASS::mvrnorm(1, rep(0, 35), C))
  u <- as.numeric(MASS::mvrnorm(1, rep(0, 35), C))
  y <- 1 + 1.5 * x + u + rnorm(35, 0, 0.6)
  d <- data.frame(y = y, x = x, Species = tr$tip.label,
                  stringsAsFactors = FALSE)
  d$y[sample(35, 12)] <- NA                            # ~34% missing on y

  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
    phylo = tr, data = d, runs = 3,
    nitt = 1500, thin = 2, burnin = 300, verbose = FALSE)))
  fin <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = tr, n_final = 5,
    nitt = 1500, thin = 2, burnin = 300, verbose = FALSE)))
  pooled <- suppressWarnings(suppressMessages(pool_posteriors(fin)))

  slope_of <- function(sol) {
    sol <- as.matrix(sol)
    sol[, which(colnames(sol) == "x")]
  }

  var_pooled  <- stats::var(slope_of(pooled$models$y$Sol))
  within_vars <- vapply(fin$all_models,
                        function(m) stats::var(slope_of(m$y$Sol)), numeric(1))
  mean_within <- mean(within_vars)
  between     <- stats::var(vapply(fin$all_models,
                        function(m) mean(slope_of(m$y$Sol)), numeric(1)))

  # Pooling must not lose variance, and must actually carry between-imputation
  # variation, and the mixture identity should hold to within MC tolerance.
  expect_gt(var_pooled, mean_within * 0.98)
  expect_gt(between, 0)
  expect_lt(abs(var_pooled - (mean_within + between)) / (mean_within + between),
            0.4)
})
