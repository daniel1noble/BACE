# Direct unit tests for the prediction helpers in model_functions.R
# (.pred_cont, .pred_count). These are exposed (.-prefixed but
# @export-ed) so we can call them with hand-built mock MCMCglmm
# objects and assert numerical output without paying for a full
# MCMC fit.

# ---- Helper: mock MCMCglmm-like objects -----------------------------------

mock_cont_model <- function(n_obs = 10, n_params = 3, n_iter = 50,
                             seed = 2026) {
  set.seed(seed)
  X <- matrix(rnorm(n_obs * n_params), n_obs, n_params)
  colnames(X) <- paste0("b", seq_len(n_params))
  Sol <- matrix(rnorm(n_iter * n_params), n_iter, n_params)
  colnames(Sol) <- paste0("b", seq_len(n_params))
  Sol <- coda::mcmc(Sol)
  list(X = X, Sol = Sol, Z = NULL)
}

mock_count_model <- function(n_obs = 10, n_iter = 50, seed = 2026) {
  set.seed(seed)
  Liab <- matrix(rnorm(n_iter * n_obs, 0, 0.5), n_iter, n_obs)
  Liab <- coda::mcmc(Liab)
  list(Liab = Liab)
}

# ---- .pred_cont -----------------------------------------------------------

test_that(".pred_cont returns posterior summary with expected columns", {
  m <- mock_cont_model(n_obs = 8, n_params = 2, n_iter = 100)
  out <- BACE:::.pred_cont(m)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("post_mean", "post_sd", "ci_lower", "ci_upper"))
  expect_equal(nrow(out), 8L)
})

test_that(".pred_cont posterior columns are finite", {
  m <- mock_cont_model(n_obs = 12, n_params = 3, n_iter = 200)
  out <- BACE:::.pred_cont(m)
  expect_true(all(is.finite(out$post_mean)))
  expect_true(all(is.finite(out$post_sd)))
  expect_true(all(out$post_sd >= 0))
  # ci_lower <= post_mean <= ci_upper
  expect_true(all(out$ci_lower <= out$post_mean + 1e-10))
  expect_true(all(out$post_mean <= out$ci_upper + 1e-10))
})

test_that(".pred_cont: post_mean equals colMeans(Sol %*% t(X)) when no Z", {
  m <- mock_cont_model(n_obs = 5, n_params = 2, n_iter = 50)
  expected_eta <- as.matrix(m$Sol) %*% t(m$X)
  expected_mean <- colMeans(expected_eta)
  out <- BACE:::.pred_cont(m)
  expect_equal(out$post_mean, as.numeric(expected_mean), tolerance = 1e-8)
})

test_that(".pred_cont errors when model$X is missing", {
  m <- list(Sol = coda::mcmc(matrix(0, 10, 2)))
  expect_error(BACE:::.pred_cont(m),
               regexp = "model\\$X is missing")
})

test_that(".pred_cont errors when model$Sol is missing", {
  m <- list(X = matrix(0, 5, 2))
  expect_error(BACE:::.pred_cont(m),
               regexp = "model\\$Sol is missing")
})

test_that(".pred_cont with Z appends random-effect columns", {
  set.seed(2026)
  X <- matrix(rnorm(8 * 2), 8, 2)
  colnames(X) <- c("b1", "b2")
  Z <- matrix(rnorm(8 * 3), 8, 3)
  colnames(Z) <- c("re1", "re2", "re3")
  Sol <- matrix(rnorm(50 * 5), 50, 5)
  colnames(Sol) <- c("b1", "b2", "re1", "re2", "re3")
  m <- list(X = X, Z = Z, Sol = coda::mcmc(Sol))
  out <- BACE:::.pred_cont(m)
  expect_equal(nrow(out), 8L)
  expect_true(all(is.finite(out$post_mean)))
})

# ---- .pred_count ---------------------------------------------------------

test_that(".pred_count returns posterior summary with expected columns", {
  m <- mock_count_model(n_obs = 10, n_iter = 100)
  out <- BACE:::.pred_count(m)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("post_mean", "post_sd", "ci_lower", "ci_upper"))
  expect_equal(nrow(out), 10L)
})

test_that(".pred_count post_mean is positive (exp() of liability)", {
  m <- mock_count_model(n_obs = 8, n_iter = 100)
  out <- BACE:::.pred_count(m)
  # Posterior mean of mu = E[exp(liab)] is always positive
  expect_true(all(out$post_mean > 0))
  expect_true(all(out$post_sd >= 0))
  expect_true(all(out$ci_lower > 0))
  expect_true(all(out$ci_upper > out$ci_lower))
})

test_that(".pred_count errors when model$Liab is missing", {
  m <- list(X = matrix(0, 10, 2), Sol = coda::mcmc(matrix(0, 10, 2)))
  expect_error(BACE:::.pred_count(m),
               regexp = "Liab is missing")
})

test_that(".pred_count: post_mean approximately equals colMeans(exp(Liab))", {
  m <- mock_count_model(n_obs = 6, n_iter = 200)
  liab <- as.matrix(m$Liab)
  expected_mean <- colMeans(exp(liab))
  out <- BACE:::.pred_count(m)
  expect_equal(out$post_mean, as.numeric(expected_mean), tolerance = 1e-8)
})

# ---- Property tests for prediction helpers -------------------------------

test_that(".pred_cont sweep over n_obs and n_iter", {
  for (n_obs in c(5, 10, 20)) {
    for (n_iter in c(20, 50, 100)) {
      m <- mock_cont_model(n_obs = n_obs, n_params = 2, n_iter = n_iter)
      out <- BACE:::.pred_cont(m)
      expect_equal(nrow(out), n_obs)
      expect_true(all(is.finite(out$post_mean)))
    }
  }
})

test_that(".pred_count sweep over n_obs", {
  for (n_obs in c(3, 5, 8, 15)) {
    m <- mock_count_model(n_obs = n_obs, n_iter = 100)
    out <- BACE:::.pred_count(m)
    expect_equal(nrow(out), n_obs)
    expect_true(all(out$post_mean > 0))
  }
})
