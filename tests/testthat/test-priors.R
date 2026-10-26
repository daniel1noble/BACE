# Tests for the user-facing `priors` argument (2026-10), which exposes
# per-variable MCMCglmm prior overrides on bace() / bace_imp() /
# bace_final_imp(). Blocks (G/R/B) replace the default block wholesale;
# R overrides are refused for latent-scale (threshold/ordinal/categorical)
# variables; the final phase reads the same BACE.gelman option as the
# chained phase (it used to hard-code gelman = 2).

types2 <- list(y = "gaussian", x = "threshold")

# ---- .validate_priors ------------------------------------------------------

test_that(".validate_priors accepts NULL and well-formed lists silently", {
  expect_silent(BACE:::.validate_priors(NULL, types2))
  expect_silent(BACE:::.validate_priors(
    list(y = list(R = list(V = 1, nu = 1))), types2, miss_vars = c("y", "x")))
})

test_that(".validate_priors rejects malformed input", {
  expect_error(BACE:::.validate_priors(list(list(G = 1)), types2),
               "named list")
  expect_error(BACE:::.validate_priors(list(zz = list(G = 1)), types2),
               "not in the model formulas")
  expect_error(BACE:::.validate_priors(list(y = list()), types2),
               "non-empty list")
  expect_error(BACE:::.validate_priors(list(y = list(Q = 1)), types2),
               "'G', 'R', 'B'")
})

test_that(".validate_priors refuses R overrides on latent-scale types", {
  expect_error(BACE:::.validate_priors(
    list(x = list(R = list(V = 2))), types2),
    "fixed at 1 on[\\s\\n]*the latent probit scale", perl = TRUE)
  # ...but allows R for gaussian
  expect_silent(BACE:::.validate_priors(
    list(y = list(R = list(V = 2, nu = 1))), types2))
})

test_that(".validate_priors warns when a prior targets a complete variable", {
  expect_warning(BACE:::.validate_priors(
    list(y = list(G = list(list(V = 1, nu = 1)))), types2, miss_vars = "x"),
    "no missing data")
})

# ---- .resolve_prior --------------------------------------------------------

default_g <- BACE:::.make_prior(n_rand = 1, type = "gaussian")

test_that(".resolve_prior returns the default untouched when override is NULL", {
  expect_identical(BACE:::.resolve_prior(default_g, NULL, "gaussian", "y"),
                   default_g)
})

test_that(".resolve_prior replaces blocks wholesale", {
  out <- BACE:::.resolve_prior(default_g,
                               list(R = list(V = 0.5, nu = 3)),
                               "gaussian", "y")
  expect_identical(out$R, list(V = 0.5, nu = 3))
  expect_identical(out$G, default_g$G)   # untouched

  # Half-Cauchy-style G, supplied bare (single random effect): auto-wrapped.
  hc <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25)
  out2 <- BACE:::.resolve_prior(default_g, list(G = hc), "gaussian", "y")
  expect_identical(out2$G[[1]], hc)
  expect_identical(out2$R, default_g$R)
})

test_that(".resolve_prior enforces G count and dimensions", {
  default_2re <- BACE:::.make_prior(n_rand = 2, type = "gaussian")
  expect_error(BACE:::.resolve_prior(
    default_2re, list(G = list(list(V = 1, nu = 1))), "gaussian", "y"),
    "2 element")

  default_cat <- suppressWarnings(BACE:::.make_prior(
    n_rand = 1, n_levels = 4, type = "categorical",
    fixform = y ~ x, data = data.frame(y = factor(c("a", "b", "c", "d")),
                                       x = 1:4), gelman = 2))
  expect_error(BACE:::.resolve_prior(
    default_cat, list(G = list(list(V = 1, nu = 1))), "categorical", "y"),
    "dimension")
  ok <- BACE:::.resolve_prior(
    default_cat, list(G = list(list(V = diag(3), nu = 1))), "categorical", "y")
  expect_identical(ok$G[[1]]$V, diag(3))
})

test_that(".resolve_prior validates B internally and against the default", {
  expect_error(BACE:::.resolve_prior(
    default_g, list(B = list(mu = c(0, 0), V = diag(3))), "gaussian", "y"),
    "must agree")

  default_cat <- suppressWarnings(BACE:::.make_prior(
    n_rand = 1, n_levels = 3, type = "categorical",
    fixform = y ~ x, data = data.frame(y = factor(c("a", "b", "c")),
                                       x = 1:3), gelman = 2))
  n_exp <- length(default_cat$B$mu)
  expect_error(BACE:::.resolve_prior(
    default_cat, list(B = list(mu = 0, V = diag(1))), "categorical", "y"),
    "expects")
  ok <- BACE:::.resolve_prior(
    default_cat, list(B = list(mu = rep(0, n_exp), V = diag(n_exp) * 2)),
    "categorical", "y")
  expect_identical(ok$B$V, diag(n_exp) * 2)
})

# ---- Integration: a supplied prior must actually reach the model -----------

test_that("priors reach .model_fit in both phases; typos warn in bace()", {
  testthat::skip_on_cran()
  set.seed(77)
  n <- 25
  tr <- ape::rtree(n)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  C <- ape::vcv(tr, corr = TRUE)
  x <- as.numeric(MASS::mvrnorm(1, rep(0, n), C))
  u <- as.numeric(MASS::mvrnorm(1, rep(0, n), C))
  y <- 1 + 1.5 * x + u + rnorm(n, 0, 0.6)
  d <- data.frame(y = y, x = x, Species = tr$tip.label,
                  stringsAsFactors = FALSE)
  d$y[sample(n, 6)] <- NA

  # Record (response variable, R$nu) for every fit: the override must reach
  # y's model and ONLY y's model (x keeps the default nu = 0.002 - the
  # chained phase also fits models for complete variables).
  seen <- list()
  real_fit <- BACE:::.model_fit
  testthat::local_mocked_bindings(
    .model_fit = function(data, tree, fixformula, randformula, type, prior, ...) {
      seen[[length(seen) + 1]] <<- c(resp = all.vars(fixformula[[2]]),
                                     nu = prior$R$nu)
      real_fit(data = data, tree = tree, fixformula = fixformula,
               randformula = randformula, type = type, prior = prior, ...)
    },
    .package = "BACE"
  )

  my_priors <- list(y = list(R = list(V = 1, nu = 7)))
  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
    phylo = tr, data = d, runs = 1,
    nitt = 800, thin = 2, burnin = 200, verbose = FALSE,
    priors = my_priors)))
  nus <- do.call(rbind, seen)
  expect_true(all(nus[nus[, "resp"] == "y", "nu"] == "7"))
  expect_true(all(nus[nus[, "resp"] == "x", "nu"] == "0.002"))

  seen <- list()
  fin <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = tr, n_final = 1,
    nitt = 800, thin = 2, burnin = 200, verbose = FALSE,
    priors = my_priors)))
  nus <- do.call(rbind, seen)
  expect_true(all(nus[nus[, "resp"] == "y", "nu"] == "7"))

  # Unknown-argument warning in bace(): arguments that match nothing used to
  # vanish silently into `...`. (Note `prior = ` itself now PARTIALLY matches
  # the `priors` formal under R's argument matching, so that old typo is
  # honoured rather than dropped - the warning covers everything else.)
  # The warning fires before any MCMC, so no tree/data needed.
  w <- character(0)
  withCallingHandlers(
    try(suppressMessages(bace(fixformula = "y ~ x", bogus_arg = 1)),
        silent = TRUE),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("does not use", w)))
})

test_that("final phase honours the BACE.gelman option (was hard-coded 2)", {
  testthat::skip_on_cran()
  set.seed(88)
  n <- 30
  tr <- ape::rtree(n)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  C <- ape::vcv(tr, corr = TRUE)
  x <- as.numeric(MASS::mvrnorm(1, rep(0, n), C))
  u <- as.numeric(MASS::mvrnorm(1, rep(0, n), 2 * C))
  lab <- cut(x + u, breaks = c(-Inf, -0.4, 0.4, Inf), labels = c("a", "b", "c"))
  d <- data.frame(y = factor(lab, levels = c("a", "b", "c")), x = x,
                  Species = tr$tip.label, stringsAsFactors = FALSE)
  d$y[sample(n, 6)] <- NA

  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
    phylo = tr, data = d, runs = 1,
    nitt = 800, thin = 2, burnin = 200, verbose = FALSE)))

  seen_gelman <- numeric(0)
  real_mp <- BACE:::.make_prior
  testthat::local_mocked_bindings(
    .make_prior = function(..., gelman = 0) {
      seen_gelman <<- c(seen_gelman, gelman)
      real_mp(..., gelman = gelman)
    },
    .package = "BACE"
  )

  old <- options(BACE.gelman = 0)
  on.exit(options(old), add = TRUE)
  fin <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = tr, n_final = 1,
    nitt = 800, thin = 2, burnin = 200, verbose = FALSE)))
  # Pre-fix this recorded 2 (hard-coded); it must now follow the option.
  expect_true(all(seen_gelman == 0))
})
