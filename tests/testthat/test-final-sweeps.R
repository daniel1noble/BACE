# Tests for the 2026-10 final-phase sweeps fix (Ayumi's point 1).
#
# Bug: each final run made ONE pass over the variables from the same
# chained-phase point-imputed start, so the first-visited variable (the
# analysis response y) conditioned on IDENTICAL point fills of the other
# variables in every run - its between-imputation variance collapsed to MCMC
# noise (observed: SD = 0 across 50 runs) and results depended on formula
# order. Fix: `sweeps` stochastic passes per run (default 3; warm-up passes
# use a short post-burn-in tail; only the last pass is kept), auto-collapsed
# to 1 when at most one variable is missing (then one pass is exact). Plus
# per-run L'Ecuyer-CMRG streams under `seed` (serial == parallel), and the
# parallel fallback now redoes only the failed runs.

# Shared fixture: two phylo-correlated gaussians, BOTH ~30% missing, with a
# cached chained-phase object (built once, reused across tests).
make_init2 <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(303)
    n <- 22
    tr <- ape::rtree(n)
    tr <- ape::compute.brlen(tr, method = "Grafen")
    tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
    C <- ape::vcv(tr, corr = TRUE)
    x <- as.numeric(MASS::mvrnorm(1, rep(0, n), C))
    u <- as.numeric(MASS::mvrnorm(1, rep(0, n), C))
    y <- 1 + 1.5 * x + u + rnorm(n, 0, 0.6)
    d <- data.frame(y = y, x = x, Species = tr$tip.label,
                    stringsAsFactors = FALSE)
    x_miss <- sample(n, 7)
    y_miss <- sample(n, 7)
    d$x[x_miss] <- NA
    d$y[y_miss] <- NA
    init <- suppressWarnings(suppressMessages(bace_imp(
      fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
      phylo = tr, data = d, runs = 2,
      nitt = 600, thin = 2, burnin = 150, verbose = FALSE)))
    cache <<- list(init = init, tr = tr, x_miss = sort(x_miss),
                   y_miss = sort(y_miss))
    cache
  }
})

test_that("sweeps: fit counts, warm-tail nitt, and varying predictor values", {
  testthat::skip_on_cran()
  f <- make_init2()

  log_fits <- list()
  real_fit <- BACE:::.model_fit
  testthat::local_mocked_bindings(
    .model_fit = function(data, tree, fixformula, randformula, type, prior,
                          nitt = 6000, thin = 5, burnin = 1000) {
      log_fits[[length(log_fits) + 1]] <<- list(
        resp = all.vars(fixformula[[2]]), nitt = nitt,
        x_at_miss = data$x[f$x_miss])
      real_fit(data = data, tree = tree, fixformula = fixformula,
               randformula = randformula, type = type, prior = prior,
               nitt = nitt, thin = thin, burnin = burnin)
    },
    .package = "BACE"
  )

  fin <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = f$init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = f$tr, n_final = 3,
    nitt = 600, thin = 2, burnin = 150, verbose = FALSE, sweeps = 2)))

  # 3 runs x 2 sweeps x 2 missing variables = 12 fits.
  expect_length(log_fits, 12L)
  expect_identical(fin$sweeps, 2L)
  expect_identical(fin$sweeps_requested, 2L)

  # Warm-up fits use the shortened tail (burnin + 50 * thin); last-sweep fits
  # use the full nitt. Within each run the first 2 fits are sweep 1 (warm).
  nitts <- vapply(log_fits, `[[`, numeric(1), "nitt")
  warm_expected <- 150 + 50 * 2
  per_run <- split(seq_len(12L), rep(1:3, each = 4L))
  for (idx in per_run) {
    expect_identical(unname(nitts[idx]), c(warm_expected, warm_expected, 600, 600))
  }

  # THE regression (Ayumi's SD = 0): the x values seen by each run's KEPT
  # (last-sweep) y-fit must differ across runs. y is visited first, so the
  # kept y-fit is fit 3 within each run's block of 4.
  kept_y_x <- lapply(per_run, function(idx) log_fits[[idx[3]]]$x_at_miss)
  expect_identical(vapply(log_fits[c(3, 7, 11)], `[[`, character(1), "resp"),
                   rep("y", 3))
  sds <- apply(do.call(rbind, kept_y_x), 2, sd)
  expect_gt(max(sds), 0)
  expect_true(mean(sds > 0) == 1)   # every x-missing cell varies across runs

  # And the completed datasets still contain no NAs.
  for (ds in fin$all_datasets) expect_false(anyNA(ds[, c("y", "x")]))
})

test_that("sweeps auto-collapse to 1 when only one variable is missing", {
  testthat::skip_on_cran()
  set.seed(304)
  n <- 20
  tr <- ape::rtree(n)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  x <- rnorm(n)
  y <- 1 + x + rnorm(n)
  d <- data.frame(y = y, x = x, Species = tr$tip.label,
                  stringsAsFactors = FALSE)
  d$y[sample(n, 6)] <- NA
  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
    phylo = tr, data = d, runs = 1,
    nitt = 600, thin = 2, burnin = 150, verbose = FALSE)))

  n_fits <- 0L
  real_fit <- BACE:::.model_fit
  testthat::local_mocked_bindings(
    .model_fit = function(...) { n_fits <<- n_fits + 1L; real_fit(...) },
    .package = "BACE"
  )
  fin <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = tr, n_final = 2,
    nitt = 600, thin = 2, burnin = 150, verbose = FALSE, sweeps = 3)))

  # One pass is exact here: 2 runs x 1 sweep x 1 missing variable = 2 fits.
  expect_identical(n_fits, 2L)
  expect_identical(fin$sweeps, 1L)
  expect_identical(fin$sweeps_requested, 3L)
})

test_that("seed gives reproducible, serial==parallel final imputations", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")   # mclapply forking
  f <- make_init2()

  run_fin <- function(cores) {
    suppressWarnings(suppressMessages(bace_final_imp(
      bace_object = f$init, fixformula = "y ~ x",
      ran_phylo_form = "~ 1 | Species", phylo = f$tr, n_final = 2,
      nitt = 600, thin = 2, burnin = 150, verbose = FALSE,
      sweeps = 2, seed = 42, n_cores = cores)))
  }

  a <- run_fin(1)
  b <- run_fin(1)
  expect_equal(a$all_datasets, b$all_datasets)

  cc <- run_fin(2)
  expect_equal(a$all_datasets, cc$all_datasets)

  # Different seed, different draws.
  d2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = f$init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = f$tr, n_final = 2,
    nitt = 600, thin = 2, burnin = 150, verbose = FALSE,
    sweeps = 2, seed = 43)))
  expect_false(isTRUE(all.equal(a$all_datasets, d2$all_datasets)))
})

test_that("parallel fallback redoes failed runs on their original streams", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  f <- make_init2()

  # Reference: clean serial run with the seed.
  ref <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = f$init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = f$tr, n_final = 2,
    nitt = 600, thin = 2, burnin = 150, verbose = FALSE,
    sweeps = 2, seed = 7, n_cores = 1)))

  # Make every FORKED worker fail (children have a different pid), so the
  # fallback must redo all runs serially - on their original streams.
  main_pid <- Sys.getpid()
  real_fit <- BACE:::.model_fit
  testthat::local_mocked_bindings(
    .model_fit = function(...) {
      if (Sys.getpid() != main_pid) stop("simulated worker crash")
      real_fit(...)
    },
    .package = "BACE"
  )
  fb <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = f$init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = f$tr, n_final = 2,
    nitt = 600, thin = 2, burnin = 150, verbose = FALSE,
    sweeps = 2, seed = 7, n_cores = 2)))

  expect_false(any(vapply(fb$all_datasets, inherits, logical(1), "try-error")))
  expect_equal(ref$all_datasets, fb$all_datasets)
})

test_that("sweeps argument is validated", {
  testthat::skip_on_cran()
  f <- make_init2()
  expect_error(suppressWarnings(bace_final_imp(
    bace_object = f$init, fixformula = "y ~ x",
    ran_phylo_form = "~ 1 | Species", phylo = f$tr, n_final = 1,
    verbose = FALSE, sweeps = 0)), "sweeps")
})

test_that("a warm-up sweep fit failure is tolerated (skipped Gibbs update + warning)", {
  testthat::skip_on_cran()
  f <- make_init2()
  warm_nitt <- 150 + 50 * 2

  n_calls <- 0L
  real_fit <- BACE:::.model_fit
  testthat::local_mocked_bindings(
    .model_fit = function(data, tree, fixformula, randformula, type, prior,
                          nitt = 6000, thin = 5, burnin = 1000) {
      n_calls <<- n_calls + 1L
      if (n_calls == 1L && nitt == warm_nitt) stop("singular mock (warm-up)")
      real_fit(data = data, tree = tree, fixformula = fixformula,
               randformula = randformula, type = type, prior = prior,
               nitt = nitt, thin = thin, burnin = burnin)
    },
    .package = "BACE"
  )

  warns <- character(0)
  fin <- withCallingHandlers(
    suppressMessages(bace_final_imp(
      bace_object = f$init, fixformula = "y ~ x",
      ran_phylo_form = "~ 1 | Species", phylo = f$tr, n_final = 1,
      nitt = 600, thin = 2, burnin = 150, verbose = FALSE, sweeps = 2)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    })

  expect_true(any(grepl("carrying the previous sweep's values", warns)))
  expect_length(fin$all_datasets, 1L)
  expect_false(anyNA(fin$all_datasets[[1]][, c("y", "x")]))
})

test_that("a LAST-sweep fit failure still aborts the run (kept results never degraded)", {
  testthat::skip_on_cran()
  f <- make_init2()

  real_fit <- BACE:::.model_fit
  testthat::local_mocked_bindings(
    .model_fit = function(data, tree, fixformula, randformula, type, prior,
                          nitt = 6000, thin = 5, burnin = 1000) {
      if (nitt == 600) stop("singular mock (last sweep)")
      real_fit(data = data, tree = tree, fixformula = fixformula,
               randformula = randformula, type = type, prior = prior,
               nitt = nitt, thin = thin, burnin = burnin)
    },
    .package = "BACE"
  )

  expect_error(
    suppressWarnings(suppressMessages(bace_final_imp(
      bace_object = f$init, fixformula = "y ~ x",
      ran_phylo_form = "~ 1 | Species", phylo = f$tr, n_final = 1,
      nitt = 600, thin = 2, burnin = 150, verbose = FALSE, sweeps = 2))),
    regexp = "singular mock \\(last sweep\\)")
})
