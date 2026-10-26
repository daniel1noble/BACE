# Regression tests for the 2026-10 BLUP-lookup fix and the posterior-averaged
# forward probabilities.
#
# Bug 1 (lookup): the categorical forward helpers built species-effect column
# names as "trait.level.<tip>" but MCMCglmm names them
# "trait.level.<REterm>.<tip>", so the lookup never matched and every species
# effect was silently 0 - the default categorical path did non-phylogenetic
# imputation. The threshold helpers matched the phylo effect but never added
# the second species effect ("<REterm>2", species = TRUE). All forward
# helpers now gather BLUPs via .gather_blups(_mat) (exact prefix + suffix
# match) with a zero-match tripwire.
#
# Bug 2 (plug-in): the forward probability paths evaluated Phi / softmax at
# posterior-MEAN coefficients, which is overconfident (Jensen). They now
# average the per-draw probabilities over all stored draws; the single-draw
# (_iter) path used for imputation draws is the same core at one iteration.
#
# These tests pin: (a) hand-computed equality INCLUDING BLUPs and averaging,
# (b) that species effects actually change the output, (c) node-effect and
# dotted-tip handling, (d) the tripwire warning, (e) averaged == mean over
# per-iteration sweep, == plug-in exactly when the posterior is degenerate,
# and (f) a real MCMCglmm fit (only a real fit can catch drift between
# MCMCglmm's actual naming and what the helpers assume).

# ---- Mock builders --------------------------------------------------------

# Categorical: response y, baseline "A", non-baseline levels B, C; predictor x;
# 3 species + 1 node decoy + 1 dotted tip; optional Species2 columns.
mock_cat_forward_model <- function(n_iter = 40, species2 = FALSE, seed = 1,
                                   zero_var = FALSE) {
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
  if (zero_var) Sol <- matrix(rep(Sol[1, ], each = n_iter), n_iter)
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

# Independent re-implementations averaging the per-draw probabilities,
# with BLUPs looked up directly by their full column names (no prefix logic).
cat_probs_by_hand <- function(m, d, use_species2 = FALSE, iters = NULL) {
  Sol <- as.matrix(m$Sol)
  if (is.null(iters)) iters <- seq_len(nrow(Sol))
  J <- 2
  c2 <- (16 * sqrt(3) / (15 * pi))^2
  IJ <- 1 / (J + 1) * (diag(J) + matrix(1, J, J))
  scl <- sqrt(1 + c2 * diag(IJ))
  acc <- matrix(0, nrow(d), J + 1)
  for (k in iters) {
    b <- Sol[k, ]
    bl_B <- b[paste0("traity.B.Species.", d$Species)]
    bl_C <- b[paste0("traity.C.Species.", d$Species)]
    if (use_species2) {
      bl_B <- bl_B + b[paste0("traity.B.Species2.", d$Species)]
      bl_C <- bl_C + b[paste0("traity.C.Species2.", d$Species)]
    }
    eB <- exp((b[["traity.B"]] + b[["traity.B:x"]] * d$x + bl_B) / scl[1])
    eC <- exp((b[["traity.C"]] + b[["traity.C:x"]] * d$x + bl_C) / scl[2])
    den <- 1 + eB + eC
    acc <- acc + cbind(1 / den, eB / den, eC / den)
  }
  acc / length(iters)
}

thr_probs_by_hand <- function(m, d, use_species2 = FALSE, iters = NULL) {
  Sol <- as.matrix(m$Sol)
  CP  <- as.matrix(m$CP)
  if (is.null(iters)) iters <- seq_len(nrow(Sol))
  acc <- matrix(0, nrow(d), 3)
  for (k in iters) {
    b <- Sol[k, ]
    cp <- c(0, CP[k, ])
    eta <- b[["(Intercept)"]] + b[["x"]] * d$x +
      b[paste0("Species.", d$Species)]
    if (use_species2) eta <- eta + b[paste0("Species2.", d$Species)]
    acc <- acc + cbind(pnorm(cp[1] - eta),
                       pnorm(cp[2] - eta) - pnorm(cp[1] - eta),
                       1 - pnorm(cp[2] - eta))
  }
  acc / length(iters)
}

# ---- Categorical forward: BLUPs + posterior averaging, hand-computed -------

test_that("categorical forward probabilities average per-draw softmax incl. BLUPs", {
  m <- mock_cat_forward_model()
  d <- mock_data()
  out <- BACE:::.pred_cat_forward(m, formula = y ~ x, data_i = d,
                                  cluster_col = "Species")
  expect_equal(unname(as.matrix(out)), unname(cat_probs_by_hand(m, d)),
               tolerance = 1e-12)

  # The lookup regression: zeroing the species effects MUST change the output.
  m0 <- m
  m0$Sol[, grep("\\.Species\\.", colnames(m$Sol))] <- 0
  out0 <- BACE:::.pred_cat_forward(m0, formula = y ~ x, data_i = d,
                                   cluster_col = "Species")
  expect_false(isTRUE(all.equal(as.matrix(out), as.matrix(out0))))

  # The plug-in regression: averaged probabilities differ from the softmax at
  # posterior-mean coefficients (Jensen), except in degenerate posteriors.
  beta <- colMeans(as.matrix(m$Sol))
  m_pl <- m
  m_pl$Sol <- coda::mcmc(matrix(rep(beta, each = 2), 2,
                                dimnames = list(NULL, names(beta))))
  plugin <- BACE:::.pred_cat_forward(m_pl, formula = y ~ x, data_i = d,
                                     cluster_col = "Species")
  expect_false(isTRUE(all.equal(as.matrix(out), as.matrix(plugin))))
})

test_that("categorical forward_iter equals the hand computation at one draw", {
  m <- mock_cat_forward_model()
  d <- mock_data()
  it <- 7L
  out <- BACE:::.pred_cat_forward_iter(m, formula = y ~ x, data_i = d,
                                       cluster_col = "Species", iteration = it)
  expect_equal(unname(as.matrix(out)),
               unname(cat_probs_by_hand(m, d, iters = it)), tolerance = 1e-12)
})

test_that("averaged forward equals the mean over the per-iteration sweep", {
  m <- mock_cat_forward_model(n_iter = 15)
  d <- mock_data()
  out <- BACE:::.pred_cat_forward(m, formula = y ~ x, data_i = d,
                                  cluster_col = "Species")
  sweep_mean <- Reduce(`+`, lapply(seq_len(15), function(it) {
    as.matrix(BACE:::.pred_cat_forward_iter(m, formula = y ~ x, data_i = d,
                                            cluster_col = "Species",
                                            iteration = it))
  })) / 15
  expect_equal(as.matrix(out), sweep_mean, tolerance = 1e-8)

  # Degenerate posterior: averaging == plug-in exactly.
  mz <- mock_cat_forward_model(zero_var = TRUE)
  outz <- BACE:::.pred_cat_forward(mz, formula = y ~ x, data_i = d,
                                   cluster_col = "Species")
  outz1 <- BACE:::.pred_cat_forward_iter(mz, formula = y ~ x, data_i = d,
                                         cluster_col = "Species", iteration = 1L)
  expect_equal(as.matrix(outz), as.matrix(outz1), tolerance = 1e-12)

  # Rows are probability distributions.
  expect_equal(unname(rowSums(as.matrix(out))), rep(1, 3), tolerance = 1e-10)
})

test_that("categorical forward adds the second species effect when present", {
  m <- mock_cat_forward_model(species2 = TRUE)
  d <- mock_data()
  out <- BACE:::.pred_cat_forward(m, formula = y ~ x, data_i = d,
                                  cluster_col = "Species")
  expect_equal(unname(as.matrix(out)),
               unname(cat_probs_by_hand(m, d, use_species2 = TRUE)),
               tolerance = 1e-12)
})

# ---- Threshold forward: both species effects, averaging, hand-computed -----

test_that("threshold forward includes phylo + second species effects (averaged)", {
  m <- mock_thr_forward_model(species2 = TRUE)
  d <- mock_data()
  out <- BACE:::.pred_threshold_forward(m, formula = y ~ x, data_i = d,
                                        cluster_col = "Species",
                                        level_names = c("A", "B", "C"))
  expect_equal(unname(as.matrix(out)),
               unname(thr_probs_by_hand(m, d, use_species2 = TRUE)),
               tolerance = 1e-12)

  it <- 11L
  out_i <- BACE:::.pred_threshold_forward_iter(m, formula = y ~ x, data_i = d,
                                               cluster_col = "Species",
                                               level_names = c("A", "B", "C"),
                                               iteration = it)
  expect_equal(unname(as.matrix(out_i)),
               unname(thr_probs_by_hand(m, d, use_species2 = TRUE, iters = it)),
               tolerance = 1e-12)
})

# ---- Decoys and tripwire --------------------------------------------------

test_that("node-effect columns are ignored and dotted tip labels resolve", {
  m <- mock_cat_forward_model()
  d <- mock_data()
  out <- BACE:::.pred_cat_forward(m, formula = y ~ x, data_i = d,
                                  cluster_col = "Species")
  # Row 3 is the dotted tip "sp.4": its BLUP must be the sp.4 column, not a
  # regex-mangled lookup. Hand computation restricted to row 3.
  expect_equal(unname(as.matrix(out[3, , drop = FALSE])),
               unname(cat_probs_by_hand(m, d[3, , drop = FALSE])),
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
