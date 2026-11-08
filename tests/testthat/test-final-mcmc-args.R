# final_nitt / final_thin / final_burnin: separate MCMC settings for the
# final imputation runs.
#
# The two phases do different jobs (the chained phase only needs point
# imputations to settle; the final runs supply the draws that are pooled),
# so users routinely want longer final chains. Before this, bace() forced
# one setting on both phases and the case studies had to bypass the wrapper.
# Default NULL must reproduce the old behaviour exactly.

# A fake pipeline that records what each phase was given, so the forwarding
# can be checked without fitting anything.
mock_pipeline <- function(seen, converged) {
  d <- data.frame(y = c(NA, rnorm(7)), x1 = rnorm(8),
                  Species = paste0("s", 1:8))
  testthat::local_mocked_bindings(
    bace_imp = function(..., nitt, thin, burnin, runs) {
      seen$chain <<- list(nitt = nitt, thin = thin, burnin = burnin,
                          runs = runs)
      structure(list(data = list(d)), class = "bace")
    },
    assess_convergence = function(...) list(converged = converged),
    bace_final_imp = function(..., nitt, thin, burnin, n_final) {
      seen$final <<- list(nitt = nitt, thin = thin, burnin = burnin,
                          n_final = n_final)
      structure(list(all_datasets = list(d), all_prob_preds = list(),
                     all_models = list()), class = "bace_final")
    },
    pool_posteriors = function(...) structure(list(models = list()),
                                              class = "bace_pooled"),
    .package = "BACE", .env = parent.frame())
  d
}

testthat::test_that("final_* default to the chained settings (previous behaviour)", {
  seen <- new.env()
  d <- mock_pipeline(seen, converged = TRUE)
  suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", phylo = ape::rcoal(8), data = d,
         nitt = 4321, thin = 7, burnin = 321, verbose = FALSE)))
  testthat::expect_identical(seen$final[c("nitt", "thin", "burnin")],
                             seen$chain[c("nitt", "thin", "burnin")])
  testthat::expect_identical(seen$final$nitt, 4321)
  testthat::expect_identical(seen$final$thin, 7)
  testthat::expect_identical(seen$final$burnin, 321)
})

testthat::test_that("final_* reach the final phase only, leaving the chained phase alone", {
  seen <- new.env()
  d <- mock_pipeline(seen, converged = TRUE)
  suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", phylo = ape::rcoal(8), data = d,
         nitt = 5000, thin = 5, burnin = 1000,
         final_nitt = 50000, final_thin = 45, final_burnin = 5000,
         verbose = FALSE)))
  testthat::expect_identical(seen$chain[c("nitt", "thin", "burnin")],
                             list(nitt = 5000, thin = 5, burnin = 1000))
  testthat::expect_identical(seen$final[c("nitt", "thin", "burnin")],
                             list(nitt = 50000, thin = 45, burnin = 5000))
})

testthat::test_that("final_* are forwarded on the skip_conv / non-converged branch too", {
  # bace() calls bace_final_imp at TWO sites; this is the one the examples
  # and both case studies actually use.
  seen <- new.env()
  d <- mock_pipeline(seen, converged = FALSE)
  suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", phylo = ape::rcoal(8), data = d,
         nitt = 5000, thin = 5, burnin = 1000, final_nitt = 50000,
         skip_conv = TRUE, verbose = FALSE)))
  testthat::expect_identical(seen$final$nitt, 50000)
  # unspecified ones still inherit
  testthat::expect_identical(seen$final$thin, 5)
  testthat::expect_identical(seen$final$burnin, 1000)
})

testthat::test_that("final_* accept per-formula lists, like nitt/thin/burnin", {
  seen <- new.env()
  d <- mock_pipeline(seen, converged = TRUE)
  suppressWarnings(suppressMessages(
    bace(list("y ~ x1", "x1 ~ y"), "~ 1 | Species", phylo = ape::rcoal(8),
         data = d, nitt = 5000,
         final_nitt = list(50000, 20000), final_burnin = list(5000, 2000),
         verbose = FALSE)))
  testthat::expect_identical(seen$final$nitt, list(50000, 20000))
  testthat::expect_identical(seen$final$burnin, list(5000, 2000))
})

testthat::test_that("the two phases really do run different chain lengths", {
  testthat::skip_on_cran()
  set.seed(4)
  n <- 20
  tr <- ape::rcoal(n)
  d <- data.frame(Species = tr$tip.label, y = rnorm(n), x1 = rnorm(n))
  d$y[1:3] <- NA

  # all three settings differ between phases, so an implementation that
  # honoured only final_nitt (or swapped thin/burnin) cannot pass
  fit <- suppressWarnings(bace(
    "y ~ x1", "~ 1 | Species", phylo = tr, data = d,
    runs = 2, nitt = 400, burnin = 100, thin = 2,          # chained -> 150
    final_nitt = 900, final_burnin = 300, final_thin = 3,  # final   -> 200
    n_final = 2, skip_conv = TRUE, verbose = FALSE))

  # stored draws per model = (nitt - burnin) / thin, per phase
  testthat::expect_equal(
    nrow(fit$initial_results$models_last_run$y$Sol), (400 - 100) / 2)
  testthat::expect_equal(
    nrow(fit$final_results$all_models[[1]]$y$Sol), (900 - 300) / 3)

  # both phases' settings are recorded, so bace_validate(fit) can re-run the
  # pipeline as it was actually run (final_results$settings holds the FINAL
  # values; bace() keeps the chained ones alongside them)
  testthat::expect_equal(fit$settings$nitt, 400)
  testthat::expect_equal(fit$settings$final_nitt, 900)
  testthat::expect_equal(fit$final_results$settings$nitt, 900)
})

# ---------------------------------------------------------------------------
# Review regressions (adversarial review 2026-11-08)
# ---------------------------------------------------------------------------

testthat::test_that("the new arguments did not move the existing ones", {
  # They must sit at the END of the formals: inserting them mid-signature
  # silently rebinds positional calls written against earlier versions
  # (bace(..., 600, 2, 100, 3, 4) would read runs/n_final as final_nitt/thin).
  testthat::expect_identical(
    names(formals(bace))[1:12],
    c("fixformula", "ran_phylo_form", "phylo", "data", "nitt", "thin",
      "burnin", "runs", "n_final", "species", "verbose", "plot"))
  testthat::expect_true(all(c("final_nitt", "final_thin", "final_burnin") %in%
                              names(formals(bace))))
})

testthat::test_that("a positional call still means what it meant before", {
  seen <- new.env()
  d <- mock_pipeline(seen, converged = TRUE)
  tr <- ape::rcoal(8)
  suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", tr, d, 600, 2, 100, 3, 4,
         verbose = FALSE)))
  testthat::expect_identical(seen$chain$runs, 3)      # not final_nitt
  testthat::expect_identical(seen$final$n_final, 4)   # not final_thin
})

testthat::test_that("bad final_* settings are caught before the chained phase runs", {
  seen <- new.env()
  d <- mock_pipeline(seen, converged = TRUE)
  tr <- ape::rcoal(8)
  # the error must name the argument the user actually wrote...
  testthat::expect_error(
    suppressMessages(bace("y ~ x1", "~ 1 | Species", tr, d,
                          final_nitt = "5e4", verbose = FALSE)),
    "final_nitt")
  # ...and must fire before any fitting happens
  testthat::expect_null(seen$chain)
})

testthat::test_that("a final chain shorter than the convergence chain warns", {
  seen <- new.env()
  d <- mock_pipeline(seen, converged = TRUE)
  tr <- ape::rcoal(8)
  testthat::expect_warning(
    suppressMessages(bace("y ~ x1", "~ 1 | Species", tr, d,
                          nitt = 5000, final_nitt = 1000, verbose = FALSE)),
    "shorter than")
})

testthat::test_that("the retry ladder still re-runs the chained phase at chained settings", {
  # the second bace_imp() call site: only reached when convergence fails
  seen <- new.env(); seen$chain_calls <- list()
  d <- data.frame(y = c(NA, rnorm(7)), x1 = rnorm(8),
                  Species = paste0("s", 1:8))
  testthat::local_mocked_bindings(
    bace_imp = function(..., nitt, thin, burnin, runs) {
      seen$chain_calls[[length(seen$chain_calls) + 1L]] <<-
        list(nitt = nitt, thin = thin, burnin = burnin, runs = runs)
      structure(list(data = list(d)), class = "bace")
    },
    assess_convergence = function(...) list(converged = FALSE),
    bace_final_imp = function(..., nitt) {
      seen$final_nitt <<- nitt
      structure(list(all_datasets = list(d), all_prob_preds = list(),
                     all_models = list()), class = "bace_final")
    },
    pool_posteriors = function(...) structure(list(models = list()),
                                              class = "bace_pooled"),
    .package = "BACE")
  suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", phylo = ape::rcoal(8), data = d,
         nitt = 500, thin = 5, burnin = 100, runs = 4, max_attempts = 2,
         final_nitt = 9000, verbose = FALSE)))
  testthat::expect_gt(length(seen$chain_calls), 1L)   # ladder was entered
  for (cc in seen$chain_calls) testthat::expect_identical(cc$nitt, 500)
  testthat::expect_identical(seen$final_nitt, 9000)
})

testthat::test_that("bace_validate(fit) re-runs each phase at ITS OWN settings", {
  # bace_final_imp records the settings IT was given, which are the final-phase
  # ones once final_* are used. Without bace() storing both phases,
  # bace_validate(fit) would re-run the convergence phase at the long final
  # chain length -- silently, and far more slowly than the run it validates.
  tr <- ape::rcoal(8)
  d0 <- data.frame(y = c(NA, rnorm(7)), x1 = rnorm(8),
                   Species = tr$tip.label)
  fake <- structure(list(
    final_results = structure(list(
      formulas = list(stats::as.formula("y ~ x1")), n_final = 5L,
      sweeps_requested = 3L, trees_used = list(tr),
      types = list(y = "gaussian"),
      settings = list(ran_phylo_form = "~ 1 | Species",
                      nitt = 900, thin = 3, burnin = 300,   # FINAL phase
                      species = FALSE, ovr_categorical = FALSE,
                      nitt_cat_mult = 1L, priors = NULL)
    ), class = "bace_final"),
    initial_results = list(data = list(Initial_Data = d0, Iter_1 = d0)),
    settings = list(nitt = 400, thin = 2, burnin = 100,     # CHAINED phase
                    final_nitt = 900, final_thin = 3, final_burnin = 300)
  ), class = "bace_complete")

  seen <- NULL
  testthat::local_mocked_bindings(
    .bace_validate_run = function(...) { seen <<- list(...); "ran" },
    .package = "BACE")
  invisible(bace_validate(fake, prop = 0.2, seed = 1, verbose = FALSE))

  testthat::expect_identical(seen$nitt, 400)          # chained, not 900
  testthat::expect_identical(seen$thin, 2)
  testthat::expect_identical(seen$burnin, 100)
  testthat::expect_identical(seen$final_nitt, 900)    # final phase preserved
  testthat::expect_identical(seen$final_thin, 3)
  testthat::expect_identical(seen$final_burnin, 300)

  # objects from before bace() stored its own settings still work: the final
  # settings are used for both phases, as they were before this feature
  old <- fake; old$settings <- NULL
  invisible(bace_validate(old, prop = 0.2, seed = 1, verbose = FALSE))
  testthat::expect_identical(seen$nitt, 900)
})
