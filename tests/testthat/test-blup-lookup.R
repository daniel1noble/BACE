# Regression tests for the 2026-10 BLUP-lookup fix.
#
# Bug: the categorical forward helpers built species-effect column names as
# "trait.level.<tip>" but MCMCglmm names them "trait.level.<REterm>.<tip>",
# so the lookup never matched and every species effect was silently 0 - the
# default categorical path did non-phylogenetic imputation. The threshold
# helpers matched the phylo effect but never added the second species effect
# ("<REterm>2", species = TRUE). All four forward helpers now gather BLUPs
# via .gather_blups (exact prefix + suffix match) with a zero-match tripwire.
#
# These tests pin: (a) hand-computed equality INCLUDING BLUPs, (b) that
# species effects actually change the output (the exact regression), (c) node
# -effect and dotted-tip handling, (d) the tripwire warning, and (e) a real
# MCMCglmm fit, because only a real fit can catch drift between MCMCglmm's
# actual naming and what the helpers assume (mocks could never catch it).

# ---- Mock builders --------------------------------------------------------

# Categorical: response y, baseline "A", non-baseline levels B, C; predictor x;
# 3 species + 1 node decoy + 1 dotted tip; optional Species2 columns.
mock_cat_forward_model <- function(n_iter = 40, species2 = FALSE, seed = 1) {
  set.seed(seed)
  tips <- c("sp1", "sp2", "sp.4")   # includes a dotted tip label
  cols <- c("traity.B", "traity.C",
            "traity.B:x", "traity.C:x",
            paste0("traity.B.Species.", c(tips, "Node1")),
            paste0("traity.C.Species.", c(tips, "Node1")))
  if (species2) {
    cols <- c(cols,
              paste0("traity.B.Species2.", tips),
              paste0("traity.C.Species2.", tips))
  }
  Sol <- matrix(rnorm(n_iter * length(cols), 0, 0.8), n_iter, length(cols))
  colnames(Sol) <- cols
  list(Sol = coda::mcmc(Sol))
}

# Threshold (3 ordered levels -> 1 free cutpoint): intercept + x; phylo and
# (optionally) second species effects; node decoy.
mock_thr_forward_model <- function(n_iter = 40, species2 = FALSE, seed = 2) {
  set.seed(seed)
  tips <- c("sp1", "sp2", "sp.4")
  cols <- c("(Intercept)", "x", paste0("Species.", c(tips, "Node1")))
  if (species2) cols <- c(cols, paste0("Species2.", tips))
  Sol <- matrix(rnorm(n_iter * length(cols), 0, 0.8), n_iter, length(cols))
  colnames(Sol) <- cols
  CP <- coda::mcmc(matrix(abs(rnorm(n_iter, 1, 0.1)), n_iter, 1))
  list(Sol = coda::mcmc(Sol), CP = CP)
}

mock_data <- function() {
  data.frame(y = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
             x = c(-0.5, 0.2, 1.1),
             Species = c("sp1", "sp2", "sp.4"),
             stringsAsFactors = FALSE)
}

# Independent re-implementation of the categorical softmax, taking the BLUP
# values directly (no name lookup), for hand-computed expectations.
cat_softmax_by_hand <- function(beta, data_i, blups_B, blups_C) {
  J  <- 2
  c2 <- (16 * sqrt(3) / (15 * pi))^2
  IJ <- 1 / (J + 1) * (diag(J) + matrix(1, J, J))
  scl <- sqrt(1 + c2 * diag(IJ))
  eta_B <- beta[["traity.B"]] + beta[["traity.B:x"]] * data_i$x + blups_B
  eta_C <- beta[["traity.C"]] + beta[["traity.C:x"]] * data_i$x + blups_C
  eB <- exp(eta_B / scl[1]); eC <- exp(eta_C / scl[2])
  denom <- 1 + eB + eC
  cbind(1 / denom, eB / denom, eC / denom)
}

# ---- Categorical forward: BLUPs included, hand-computed ------------------

test_that("categorical forward predictions include phylo BLUPs (posterior mean)", {
  m <- mock_cat_forward_model()
  d <- mock_data()
  beta <- colMeans(as.matrix(m$Sol))
  out <- BACE:::.pred_cat_forward(m, formula = y ~ x, data_i = d,
                                  cluster_col = "Species")
  expected <- cat_softmax_by_hand(
    beta, d,
    blups_B = beta[paste0("traity.B.Species.", d$Species)],
    blups_C = beta[paste0("traity.C.Species.", d$Species)])
  expect_equal(unname(as.matrix(out)), unname(expected), tolerance = 1e-12)

  # The regression: zeroing the species effects MUST change the output.
  m0 <- m
  m0$Sol[, grep("\\.Species\\.", colnames(m$Sol))] <- 0
  out0 <- BACE:::.pred_cat_forward(m0, formula = y ~ x, data_i = d,
                                   cluster_col = "Species")
  expect_false(isTRUE(all.equal(as.matrix(out), as.matrix(out0))))
})

test_that("categorical forward_iter uses the requested iteration's BLUPs", {
  m <- mock_cat_forward_model()
  d <- mock_data()
  it <- 7L
  beta <- as.matrix(m$Sol)[it, ]
  out <- BACE:::.pred_cat_forward_iter(m, formula = y ~ x, data_i = d,
                                       cluster_col = "Species", iteration = it)
  expected <- cat_softmax_by_hand(
    beta, d,
    blups_B = beta[paste0("traity.B.Species.", d$Species)],
    blups_C = beta[paste0("traity.C.Species.", d$Species)])
  expect_equal(unname(as.matrix(out)), unname(expected), tolerance = 1e-12)
})

test_that("categorical forward adds the second species effect when present", {
  m <- mock_cat_forward_model(species2 = TRUE)
  d <- mock_data()
  beta <- colMeans(as.matrix(m$Sol))
  out <- BACE:::.pred_cat_forward(m, formula = y ~ x, data_i = d,
                                  cluster_col = "Species")
  expected <- cat_softmax_by_hand(
    beta, d,
    blups_B = beta[paste0("traity.B.Species.", d$Species)] +
              beta[paste0("traity.B.Species2.", d$Species)],
    blups_C = beta[paste0("traity.C.Species.", d$Species)] +
              beta[paste0("traity.C.Species2.", d$Species)])
  expect_equal(unname(as.matrix(out)), unname(expected), tolerance = 1e-12)
})

# ---- Threshold forward: both species effects, hand-computed ---------------

test_that("threshold forward includes phylo + second species effects", {
  m <- mock_thr_forward_model(species2 = TRUE)
  d <- mock_data()
  beta <- colMeans(as.matrix(m$Sol))
  cp <- c(0, mean(as.matrix(m$CP)))
  eta <- beta[["(Intercept)"]] + beta[["x"]] * d$x +
    beta[paste0("Species.", d$Species)] +
    beta[paste0("Species2.", d$Species)]
  expected <- cbind(pnorm(cp[1] - eta),
                    pnorm(cp[2] - eta) - pnorm(cp[1] - eta),
                    1 - pnorm(cp[2] - eta))
  out <- BACE:::.pred_threshold_forward(m, formula = y ~ x, data_i = d,
                                        cluster_col = "Species",
                                        level_names = c("A", "B", "C"))
  expect_equal(unname(as.matrix(out)), unname(expected), tolerance = 1e-12)

  it <- 11L
  beta_i <- as.matrix(m$Sol)[it, ]
  cp_i <- c(0, as.matrix(m$CP)[it, ])
  eta_i <- beta_i[["(Intercept)"]] + beta_i[["x"]] * d$x +
    beta_i[paste0("Species.", d$Species)] +
    beta_i[paste0("Species2.", d$Species)]
  expected_i <- cbind(pnorm(cp_i[1] - eta_i),
                      pnorm(cp_i[2] - eta_i) - pnorm(cp_i[1] - eta_i),
                      1 - pnorm(cp_i[2] - eta_i))
  out_i <- BACE:::.pred_threshold_forward_iter(m, formula = y ~ x, data_i = d,
                                               cluster_col = "Species",
                                               level_names = c("A", "B", "C"),
                                               iteration = it)
  expect_equal(unname(as.matrix(out_i)), unname(expected_i), tolerance = 1e-12)
})

# ---- Decoys and tripwire --------------------------------------------------

test_that("node-effect columns are ignored and dotted tip labels resolve", {
  m <- mock_cat_forward_model()
  d <- mock_data()
  beta <- colMeans(as.matrix(m$Sol))
  out <- BACE:::.pred_cat_forward(m, formula = y ~ x, data_i = d,
                                  cluster_col = "Species")
  # Row 3 is the dotted tip "sp.4": its BLUP must be the sp.4 column, not a
  # regex-mangled lookup. Verify via the hand computation restricted to row 3.
  expected3 <- cat_softmax_by_hand(
    beta, d[3, , drop = FALSE],
    blups_B = beta[["traity.B.Species.sp.4"]],
    blups_C = beta[["traity.C.Species.sp.4"]])
  expect_equal(unname(as.matrix(out[3, , drop = FALSE])), unname(expected3),
               tolerance = 1e-12)
})

test_that("zero-match tripwire warns instead of failing silently", {
  m <- mock_cat_forward_model()
  d <- mock_data()
  d$Species <- c("zz1", "zz2", "zz3")   # matches no BLUP column
  expect_warning(
    BACE:::.pred_cat_forward(m, formula = y ~ x, data_i = d,
                             cluster_col = "Species"),
    "naming mismatch")
  m2 <- mock_thr_forward_model()
  expect_warning(
    BACE:::.pred_threshold_forward(m2, formula = y ~ x, data_i = d,
                                   cluster_col = "Species",
                                   level_names = c("A", "B", "C")),
    "naming mismatch")
})

test_that(".gather_blups returns zeros without warning when no RE columns exist", {
  beta <- c(a = 1, b = 2)
  expect_silent(
    out <- BACE:::.gather_blups(beta, names(beta), "Species.", c("sp1", "sp2")))
  expect_identical(out, c(0, 0))
})

# ---- Real fit: the only test class that can catch naming drift ------------

test_that("real categorical fit: forward predictions respond to species effects", {
  testthat::skip_on_cran()
  set.seed(99)
  n <- 25
  tr <- ape::rtree(n)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  C <- ape::vcv(tr, corr = TRUE)
  x <- as.numeric(MASS::mvrnorm(1, rep(0, n), C))
  u <- as.numeric(MASS::mvrnorm(1, rep(0, n), 2 * C))
  lab <- cut(0.5 * x + u + rnorm(n, 0, 0.5),
             breaks = c(-Inf, -0.5, 0.5, Inf), labels = c("A", "B", "C"))
  d <- data.frame(y = factor(lab, levels = c("A", "B", "C")), x = x,
                  Species = tr$tip.label, stringsAsFactors = FALSE)
  types <- list(y = "categorical", x = "gaussian")
  dp <- BACE:::.data_prep(formula = y ~ x, data = d, types = types,
                          ran_cluster = "Species")
  data_i <- dp[[1]]
  prior <- suppressWarnings(
    BACE:::.make_prior(n_rand = 1, n_levels = 3, type = "categorical",
                       fixform = y ~ x, data = data_i, gelman = 2))
  model <- suppressWarnings(BACE:::.model_fit(
    data = data_i, tree = tr, fixformula = y ~ x, randformula = ~ Species,
    type = "categorical", prior = prior, nitt = 800, thin = 2, burnin = 200))

  # MCMCglmm's actual naming must be matched: BLUP columns exist...
  expect_true(any(grepl("\\.Species\\.", colnames(model$Sol))))
  # ...and zeroing them must change the forward predictions (pre-fix: no change).
  p1 <- BACE:::.pred_cat_forward(model, formula = y ~ x, data_i = data_i,
                                 cluster_col = "Species")
  m0 <- model
  m0$Sol[, grep("\\.Species\\.", colnames(model$Sol))] <- 0
  p0 <- BACE:::.pred_cat_forward(m0, formula = y ~ x, data_i = data_i,
                                 cluster_col = "Species")
  expect_gt(max(abs(as.matrix(p1) - as.matrix(p0))), 1e-6)
})
