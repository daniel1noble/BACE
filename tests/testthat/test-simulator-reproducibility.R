# Tests for the machine-independence of the trait simulator.
#
# THE BUG THESE EXIST FOR (2026-11-13). sample_random_effects() drew the
# phylogenetic effects with MASS::mvrnorm(), which builds its matrix square root
# from eigen(Sigma, symmetric = TRUE). Eigenvectors are NOT unique: LAPACK fixes
# no sign convention, so any implementation may return a column negated, and for
# repeated eigenvalues the whole subspace is arbitrary. Every such choice is a
# valid square root, so every machine drew a DIFFERENT -- though correctly
# distributed -- set of traits from an IDENTICAL random stream.
#
# Measured across three machines, same seed, same tree, same RNG settings:
#
#   Accelerate (macOS)   x1[1:3] = 7.169  -1.485  0.273   lambda_x1 0.827
#   OpenBLAS   (Linux)   x1[1:3] = 2.047  -1.537  2.636   lambda_x1 0.875
#   AOCL       (Linux)   x1[1:3] = 4.074  -1.417  4.599   lambda_x1 0.850
#
# The tree, tip labels, runif(), sample() and the masked cells were IDENTICAL on
# all three -- only the traits differed, which is what isolates the fault to the
# matrix square root rather than to the RNG.
#
# The fix: draw through the CHOLESKY factor, which is unique for a
# positive-definite matrix, so implementations agree to rounding error instead
# of choosing between valid alternatives.

suppressPackageStartupMessages(library(ape))

make_C <- function(n = 12, seed = 11) {
  set.seed(seed)
  tr <- ape::rcoal(n)
  tr$tip.label <- sprintf("t%02d", seq_len(n))
  list(tree = tr, C = ape::vcv(tr, corr = TRUE))
}

# ---- 1. the golden test: ANY machine running this suite must agree ----------

test_that("the phylogenetic draw reproduces exactly, on any machine", {
  # This is the cross-machine test. The values below were produced on
  # aarch64-apple-darwin with reference BLAS. A Cholesky factor is unique, so a
  # different LAPACK must agree to rounding error -- the tolerance is 1e-8,
  # which passes rounding (~1e-15) and fails a sign flip (order 1).
  #
  # IF THIS FAILS ON YOUR MACHINE, the simulator is not machine-independent and
  # the campaign cannot be pooled across platforms. Do not "fix" it by widening
  # the tolerance or regenerating the numbers.
  z <- make_C()
  expect_equal(sum(z$tree$edge.length), 6.328153991169, tolerance = 1e-10)

  set.seed(2024)
  d <- sample_random_effects(2.5, 12, z$C)
  expect_equal(head(d, 6),
               c(1.552629966241, 0.741103550046, 0.428162137271,
                 0.135288189880, 1.589755566908, 1.831682336927),
               tolerance = 1e-8)
  expect_equal(sum(d), 5.626721268865, tolerance = 1e-8)
})

# ---- 2. it must be the Cholesky route, not something else ------------------

test_that("the draw is exactly t(chol(Sigma)) %*% z", {
  # Pins the method, not just the answer. Reverting to mvrnorm(), or to any
  # eigen- or svd-based square root, fails here even if the distribution is
  # still right.
  z <- make_C()
  sigma2 <- 3.1
  set.seed(7); got <- sample_random_effects(sigma2, 12, z$C)
  set.seed(7); want <- as.vector(t(chol(sigma2 * z$C)) %*% rnorm(12))
  expect_equal(got, want, tolerance = 1e-12)
})

test_that("an eigen square root really is ambiguous (why the above matters)", {
  # Demonstrates the hazard, so this file documents its own reason to exist.
  # Negating eigenvector columns gives an equally valid square root and a
  # completely different draw from the same stream.
  z <- make_C(n = 60, seed = 3)
  eS <- eigen(z$C, symmetric = TRUE)
  set.seed(99); x <- rnorm(60)
  root <- function(V) as.vector(V %*% diag(sqrt(pmax(eS$values, 0))) %*% x)

  a <- root(eS$vectors)
  V <- eS$vectors; flip <- seq(1, 60, 2); V[, flip] <- -V[, flip]
  b <- root(V)

  # both are valid draws from the same distribution ...
  expect_equal(stats::var(a), stats::var(b), tolerance = 0.35)
  # ... but they are not the same numbers, by a wide margin
  expect_gt(max(abs(a - b)), 0.5)

  # the Cholesky factor admits no such choice
  L1 <- t(chol(z$C))
  L2 <- t(chol(z$C))
  expect_identical(L1, L2)
  expect_equal(L1 %*% x, L2 %*% x)
})

# ---- 3. the random stream must stay aligned -------------------------------

test_that("the draw consumes exactly n normals, so nothing downstream shifts", {
  # The masked cells and the tree were identical across machines precisely
  # because RNG consumption matched. If a future implementation drew a
  # different number of values, every downstream draw would move and the
  # campaign's missingness patterns would change silently.
  z <- make_C()
  set.seed(5); invisible(sample_random_effects(2, 12, z$C)); after_fn <- runif(1)
  set.seed(5); invisible(rnorm(12));                         after_ref <- runif(1)
  expect_equal(after_fn, after_ref, tolerance = 1e-12)

  # and the independent branch is plain rnorm, unchanged
  set.seed(8); a <- sample_random_effects(4, 25, NULL)
  set.seed(8); b <- rnorm(25, mean = 0, sd = 2)
  expect_equal(a, b, tolerance = 1e-12)
})

# ---- 4. it still draws the right distribution ------------------------------

test_that("draws have the requested variance and the tree's correlation", {
  skip_on_cran()
  z <- make_C(n = 20, seed = 4)
  sigma2 <- 2.25
  set.seed(1)
  X <- t(replicate(4000, sample_random_effects(sigma2, 20, z$C)))

  # marginal variance of each tip is sigma2 (the diagonal of a corr matrix is 1)
  expect_equal(mean(apply(X, 2, stats::var)), sigma2, tolerance = 0.12)
  # the empirical covariance recovers sigma2 * C
  emp <- stats::cov(X)
  expect_lt(max(abs(emp - sigma2 * z$C)), 0.35)
  # and the mean is zero
  expect_lt(max(abs(colMeans(X))), 0.15)
})

test_that("a near-degenerate matrix gets a deterministic ridge, not an eigen fallback", {
  # The fallback must stay deterministic. Falling back to an eigen square root
  # would quietly reintroduce the very machine-dependence this file guards.
  n <- 6
  S <- matrix(1, n, n)            # rank 1: singular, chol() must fail on it
  expect_error(chol(S), regexp = ".*")

  expect_warning(L <- BACE:::.chol_factor_lower(S), regexp = "positive")
  expect_true(all(dim(L) == c(n, n)))
  expect_true(all(is.finite(L)))
  expect_true(all(L[upper.tri(L)] == 0))        # lower triangular

  # deterministic: the same input gives the identical factor every time
  expect_identical(L, suppressWarnings(BACE:::.chol_factor_lower(S)))

  # a well-behaved matrix takes the plain path and needs no warning
  z <- make_C()
  expect_silent(L2 <- BACE:::.chol_factor_lower(z$C))
  expect_equal(unname(L2 %*% t(L2)), unname(as.matrix(z$C)), tolerance = 1e-10)
})

# ---- 5. end to end through sim_bace ---------------------------------------

test_that("sim_bace is deterministic and carries the requested phylogenetic signal", {
  skip_on_cran()
  skip_if_not_installed("phytools")

  one <- function() {
    set.seed(321)
    suppressWarnings(suppressMessages(
      sim_bace(response_type = "gaussian", predictor_types = "gaussian",
               phylo_signal = c(0.9, 0.9), missingness = c(0, 0),
               n_cases = 80, n_species = 80)))
  }
  a <- one(); b <- one()
  expect_equal(a$complete_data$y, b$complete_data$y, tolerance = 1e-12)
  expect_equal(a$complete_data$x1, b$complete_data$x1, tolerance = 1e-12)

  # high signal must actually come through on a continuous trait, where the
  # latent and observed scales are the same thing
  v <- a$complete_data$x1
  names(v) <- a$complete_data$Species
  lam <- suppressWarnings(
    phytools::phylosig(a$tree, v[a$tree$tip.label], method = "lambda")$lambda)
  expect_gt(lam, 0.5)
})
