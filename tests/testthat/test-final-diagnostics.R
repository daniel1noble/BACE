# Tests for the 2026-10 diagnostics fix (Ayumi's point 3): ordinal CUTPOINTS
# are now part of the MCMC health check, and the FINAL-phase models (the ones
# pooled for inference) are diagnosed at all - previously only the chained
# phase was, and Ayumi's cutpoint ESS of 5-22 went unreported. Aggregation is
# per-run worst-case across runs; deliberately no cross-run R-hat (the runs
# target different conditional posteriors, so cross-run spread is imputation
# uncertainty, not non-convergence).

# ---- Mock builders ---------------------------------------------------------

# A well-mixing mock chain set: iid normals -> high ESS.
mock_fit <- function(n_iter = 400, n_cuts = 1, seed = 1, random_walk = FALSE) {
  set.seed(seed)
  draw <- function(n) if (random_walk) cumsum(rnorm(n, 0, 1)) else rnorm(n)
  Sol <- cbind(draw(n_iter), draw(n_iter))
  colnames(Sol) <- c("(Intercept)", "x")
  VCV <- cbind(abs(draw(n_iter)) + 0.5, rep(1, n_iter))
  colnames(VCV) <- c("Species", "units")
  out <- list(Sol = coda::mcmc(Sol), VCV = coda::mcmc(VCV),
              Fixed = list(nfl = 2L))
  if (n_cuts > 0) {
    CP <- matrix(abs(draw(n_iter)) + 1, n_iter, n_cuts)
    colnames(CP) <- paste0("cutpoint.traity.", seq_len(n_cuts))
    out$CP <- coda::mcmc(CP)
  }
  out
}

test_that(".diagnose_one_model includes cutpoints and labels origins", {
  m <- mock_fit(n_cuts = 1)
  d <- BACE:::.diagnose_one_model(m, "threshold")
  expect_setequal(unique(d$effect_type), c("Fixed", "Variance", "Cutpoint"))
  expect_identical(d$parameter[d$effect_type == "Cutpoint"],
                   "cutpoint.traity.1")
  # The fixed latent residual is labelled Fixed and skipped...
  expect_identical(d$convergence[d$parameter == "units"], "Fixed")
  # ...but the cutpoint IS tested (iid draws -> healthy ESS).
  expect_gt(d$ess[d$effect_type == "Cutpoint"], 200)
  # Gaussian models: units is a real parameter, not fixed.
  mg <- mock_fit(n_cuts = 0)
  dg <- BACE:::.diagnose_one_model(mg, "gaussian")
  expect_false("Cutpoint" %in% dg$effect_type)
  expect_false(dg$convergence[dg$parameter == "units"] == "Fixed")
})

test_that(".check_final_diagnostics aggregates worst case across runs", {
  # Run 1 healthy; run 2's cutpoint is a random walk -> tiny ESS, flagged.
  final_out <- list(
    all_models = list(
      list(y = mock_fit(n_cuts = 1, seed = 1)),
      list(y = mock_fit(n_cuts = 1, seed = 2, random_walk = TRUE))
    ),
    all_ovr_models = list(list(), list()),
    types = list(y = "threshold")
  )
  diag <- BACE:::.check_final_diagnostics(final_out)
  expect_s3_class(diag, "bace_final_diagnostics")
  bp <- diag$by_param
  cp_row <- bp[bp$effect_type == "Cutpoint", ]
  expect_identical(nrow(cp_row), 1L)
  expect_identical(cp_row$n_runs, 2L)
  # The worst (min-ESS) run must be run 2, the random walk.
  expect_identical(cp_row$min_ess_run, 2L)
  expect_lt(cp_row$min_ess, 200)
  # And it must be in the flagged set.
  expect_true("Cutpoint" %in% diag$flagged$effect_type)
  # Printing lists the flagged cutpoint and returns invisibly.
  expect_output(print(diag), "Cutpoint")
})

test_that(".check_final_diagnostics covers OVR per-level sub-models", {
  final_out <- list(
    all_models = list(list()),
    all_ovr_models = list(
      list(y = list(a = mock_fit(n_cuts = 0, seed = 3),
                    b = mock_fit(n_cuts = 0, seed = 4)))
    ),
    types = list(y = "categorical")
  )
  diag <- BACE:::.check_final_diagnostics(final_out)
  expect_setequal(unique(diag$by_param$response_variable),
                  c("y (a)", "y (b)"))
})

test_that("bace_final objects carry diagnostics end to end", {
  testthat::skip_on_cran()
  set.seed(555)
  n <- 20
  tr <- ape::rtree(n)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  x <- rnorm(n)
  y <- 1 + x + rnorm(n)
  d <- data.frame(y = y, x = x, Species = tr$tip.label,
                  stringsAsFactors = FALSE)
  d$y[sample(n, 5)] <- NA
  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
    phylo = tr, data = d, runs = 1,
    nitt = 600, thin = 2, burnin = 150, verbose = FALSE)))
  fin <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = tr, n_final = 2,
    nitt = 600, thin = 2, burnin = 150, verbose = FALSE)))
  expect_s3_class(fin$diagnostics, "bace_final_diagnostics")
  expect_true(all(c("min_ess", "min_ess_run", "pct_runs_flagged") %in%
                    colnames(fin$diagnostics$by_param)))
  expect_identical(max(fin$diagnostics$by_param$n_runs), 2L)
  expect_output(print(fin), "MCMC health")
})
