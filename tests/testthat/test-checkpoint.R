# Checkpoint/resume.
#
# The promise is narrow and must hold exactly: a resumed run gives the same
# answer as one that was never interrupted, and work produced under different
# settings is never silently mixed. Most of this is tested without fitting
# anything (handcrafted manifests and files); one real fit proves end to end
# that completed work is reloaded rather than recomputed.

# --- helpers ---------------------------------------------------------------

tmpdir_fresh <- function() {
  d <- file.path(tempdir(), paste0("bace_ckpt_", as.integer(runif(1, 1, 1e8))))
  unlink(d, recursive = TRUE); d
}

# --- atomic write ----------------------------------------------------------

testthat::test_that("checkpoint files are written atomically and read back", {
  d <- .ckpt_dir_prepare(tmpdir_fresh())
  f <- file.path(d, "thing.rds")
  .ckpt_write(list(a = 1, b = "x"), f, compress = FALSE)
  testthat::expect_true(file.exists(f))
  testthat::expect_identical(readRDS(f), list(a = 1, b = "x"))
  # no temporary files are left behind
  testthat::expect_length(list.files(d, pattern = "[.]tmp$"), 0L)
  # compressed round-trip works too
  f2 <- file.path(d, "thing2.rds")
  .ckpt_write(list(a = 1), f2, compress = TRUE)
  testthat::expect_identical(readRDS(f2), list(a = 1))
})

testthat::test_that("a failed rename errors instead of silently losing the run", {
  d <- .ckpt_dir_prepare(tmpdir_fresh())
  f <- file.path(d, "thing.rds")
  testthat::local_mocked_bindings(file.rename = function(from, to) FALSE,
                                  .package = "base")
  testthat::expect_error(.ckpt_write(list(a = 1), f), "NOT saved")
  testthat::expect_false(file.exists(f))
})

testthat::test_that("leftover temporary files from a crash are swept on resume", {
  d <- .ckpt_dir_prepare(tmpdir_fresh())
  stale <- file.path(d, "bace_ckpt_leftover.tmp")
  saveRDS(1, stale)
  testthat::expect_message(.ckpt_sweep_stale(d, verbose = TRUE), "leftover")
  testthat::expect_false(file.exists(stale))
})

# --- the manifest: what counts as "the same analysis" ----------------------

testthat::test_that("a matching manifest resumes and a differing one is refused by name", {
  d <- .ckpt_dir_prepare(tmpdir_fresh())
  mf <- file.path(d, "m.rds")
  m1 <- list(data = "hash1", nitt = "hashA", seed = 1, bace_version = "0.2.0")

  # first call writes it and reports "not a resume"
  testthat::expect_false(.ckpt_manifest_check(mf, m1, "test", verbose = FALSE))
  # second identical call is a resume
  testthat::expect_true(.ckpt_manifest_check(mf, m1, "test", verbose = FALSE))

  # a changed setting is refused, and the message names it
  m2 <- m1; m2$nitt <- "hashB"
  testthat::expect_error(.ckpt_manifest_check(mf, m2, "test", verbose = FALSE),
                         "nitt")
  # ... and does not name settings that did not change
  testthat::expect_error(.ckpt_manifest_check(mf, m2, "test", verbose = FALSE),
                         "^(?!.*\\bseed\\b).*$", perl = TRUE)

  # a changed package version gets the reinstall hint
  m3 <- m1; m3$bace_version <- "0.3.0"
  testthat::expect_error(.ckpt_manifest_check(mf, m3, "test", verbose = FALSE),
                         "install_github")
})

testthat::test_that("ladder-mutable fields may differ without refusing the resume", {
  d <- .ckpt_dir_prepare(tmpdir_fresh())
  mf <- file.path(d, "m.rds")
  m1 <- list(data = "h", runs = 10, attempt = 1L)
  .ckpt_manifest_check(mf, m1, "test", mutable = c("runs", "attempt"),
                       verbose = FALSE)
  m2 <- list(data = "h", runs = 15, attempt = 2L)
  testthat::expect_true(.ckpt_manifest_check(mf, m2, "test",
                                             mutable = c("runs", "attempt"),
                                             verbose = FALSE))
  testthat::expect_identical(readRDS(mf)$runs, 15)   # refreshed
})

testthat::test_that("formula environments do not affect the fingerprint", {
  # serialize() follows a formula's environment, so two content-identical
  # formulas built in different places would otherwise look different (and
  # could drag unrelated megabytes into the hash)
  make <- function() { big <- rnorm(1e4); stats::as.formula("y ~ x1") }
  testthat::expect_identical(.ckpt_formula_text(make()),
                             .ckpt_formula_text(make()))
  testthat::expect_identical(.ckpt_formula_text(stats::as.formula("y ~ x1")),
                             .ckpt_formula_text(make()))
  testthat::expect_false(identical(.ckpt_formula_text("y ~ x1"),
                                   .ckpt_formula_text("y ~ x2")))
})

testthat::test_that("hashes are content-based and ignore the R version header", {
  testthat::expect_identical(.hash_obj(list(1, "a")), .hash_obj(list(1, "a")))
  testthat::expect_false(identical(.hash_obj(list(1)), .hash_obj(list(2))))
})

testthat::test_that("the chain file name carries the iteration count", {
  f10 <- .ckpt_chain_files("d", 10)
  f15 <- .ckpt_chain_files("d", 15)
  # each retry attempt saves separately instead of colliding
  testthat::expect_false(identical(f10$object, f15$object))
  testthat::expect_match(basename(f10$object), "chain_runs0010")
})

# --- seed requirement ------------------------------------------------------

testthat::test_that("checkpointing requires a seed, except for the chained phase", {
  tr <- ape::rcoal(8)
  d <- data.frame(Species = tr$tip.label, y = c(NA, rnorm(7)), x1 = rnorm(8))
  dir <- tmpdir_fresh()

  # bace(): refused up front, before any fitting
  testthat::expect_error(
    bace("y ~ x1", "~ 1 | Species", phylo = tr, data = d,
         checkpoint_dir = dir, verbose = FALSE),
    "needs a seed")

  # bace_final_imp(): same requirement
  fake <- structure(list(data = list(d), miss_dat = data.frame(
    row = 1L, col = 2L, colname = "y"), types = list(y = "gaussian"),
    phylo_ran = list(cluster = "Species")), class = "bace")
  testthat::expect_error(
    bace_final_imp(fake, "y ~ x1", "~ 1 | Species", phylo = tr,
                   checkpoint_dir = dir, verbose = FALSE),
    "needs a seed")

  # bace_imp(): no seed needed - this phase is saved, never replayed
  testthat::expect_no_error(
    suppressWarnings(bace_imp("y ~ x1", "~ 1 | Species", phylo = tr, data = d,
                              runs = 1, nitt = 300, burnin = 100, thin = 2,
                              checkpoint_dir = dir, verbose = FALSE)))
})

# --- end to end ------------------------------------------------------------

testthat::test_that("a resumed run reloads finished work and matches an uninterrupted one", {
  testthat::skip_on_cran()
  set.seed(11)
  n <- 16; tr <- ape::rcoal(n)
  d <- data.frame(Species = tr$tip.label, y = rnorm(n), x1 = rnorm(n))
  d$y[1:3] <- NA
  args <- list(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
               phylo = tr, data = d, runs = 2, nitt = 400, burnin = 100,
               thin = 2, n_final = 4, skip_conv = TRUE, verbose = FALSE,
               seed = 42)
  dir <- tmpdir_fresh()
  orig_fit <- .model_fit

  run_counted <- function(arglist) {
    k <- 0
    res <- suppressWarnings(testthat::with_mocked_bindings(
      do.call(bace, arglist),
      .model_fit = function(...) { k <<- k + 1; orig_fit(...) },
      .package = "BACE"))
    list(res = res, fits = k)
  }

  A <- run_counted(c(args, list(checkpoint_dir = dir)))
  testthat::expect_true(file.exists(file.path(dir, "final_run0004.rds")))
  testthat::expect_true(file.exists(file.path(dir, "chain_a01_runs0002.rds")))

  # lose the last two final runs, as a crash would
  file.remove(file.path(dir, c("final_run0003.rds", "final_run0004.rds")))
  B <- run_counted(c(args, list(checkpoint_dir = dir)))

  # only the missing runs were refitted; the chain and runs 1-2 came off disk
  testthat::expect_lt(B$fits, A$fits)
  testthat::expect_identical(B$res$imputed_datasets, A$res$imputed_datasets)
  testthat::expect_identical(B$res$pooled_models$models$y$Sol,
                             A$res$pooled_models$models$y$Sol)
})

testthat::test_that("a resume with changed settings is refused, not silently mixed", {
  testthat::skip_on_cran()
  set.seed(12)
  n <- 14; tr <- ape::rcoal(n)
  d <- data.frame(Species = tr$tip.label, y = rnorm(n), x1 = rnorm(n))
  d$y[1:2] <- NA
  dir <- tmpdir_fresh()
  base_args <- list(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
                    phylo = tr, data = d, runs = 1, nitt = 300, burnin = 100,
                    thin = 2, n_final = 2, skip_conv = TRUE, verbose = FALSE,
                    seed = 7, checkpoint_dir = dir)
  suppressWarnings(do.call(bace, base_args))

  # same directory, different data: must refuse rather than mix
  d2 <- d; d2$x1 <- rnorm(n)
  testthat::expect_error(
    suppressWarnings(do.call(bace, utils::modifyList(base_args,
                                                     list(data = d2)))),
    "different settings|Differing")
})

# ---------------------------------------------------------------------------
# Review regressions (adversarial review 2026-11-08)
# ---------------------------------------------------------------------------

testthat::test_that("a saved run is only reused when the manifest proves it is ours", {
  # The bug this pins: the manifest verdict was computed and then ignored, so
  # a missing or unreadable manifest let ANOTHER analysis's imputations load
  # as if they were ours - silently, and with the wrong values at cells that
  # were missing in this analysis but observed in the other.
  testthat::skip_on_cran()
  set.seed(31)
  n <- 14; tr <- ape::rcoal(n)
  mk <- function(seed) { set.seed(seed)
    d <- data.frame(Species = tr$tip.label, y = rnorm(n), x1 = rnorm(n))
    d$y[1:3] <- NA; d }
  d1 <- mk(1); d2 <- mk(2)
  args <- function(d, dir) list(fixformula = "y ~ x1",
    ran_phylo_form = "~ 1 | Species", phylo = tr, data = d, runs = 1,
    nitt = 300, burnin = 100, thin = 2, n_final = 2, skip_conv = TRUE,
    verbose = FALSE, seed = 5, checkpoint_dir = dir)

  dir <- tmpdir_fresh()
  A <- suppressWarnings(do.call(bace, args(d1, dir)))
  testthat::expect_true(file.exists(file.path(dir, "final_run0001.rds")))

  # lose the manifest, keep the saved runs: must NOT hand analysis 1's
  # imputations to analysis 2
  file.remove(file.path(dir, "final_manifest.rds"))
  chain2 <- suppressWarnings(bace_imp("y ~ x1", "~ 1 | Species", phylo = tr,
                                      data = d2, runs = 1, nitt = 300,
                                      burnin = 100, thin = 2, verbose = FALSE))
  testthat::expect_error(
    suppressWarnings(bace_final_imp(chain2, "y ~ x1", "~ 1 | Species",
                                    phylo = tr, n_final = 2, nitt = 300,
                                    burnin = 100, thin = 2, seed = 5,
                                    verbose = FALSE, checkpoint_dir = dir)),
    "cannot be matched")
})

testthat::test_that("a readable file of the wrong shape is not pooled as an imputation", {
  # Previously any non-NULL object readRDS returned was accepted, silently
  # dropping that imputation from the pooling while n_final still said 4.
  testthat::expect_false(.ckpt_is_run(data.frame(a = 1:3)))
  testthat::expect_false(.ckpt_is_run(list(nonsense = 1)))
  testthat::expect_false(.ckpt_is_run(NULL))
  testthat::expect_true(.ckpt_is_run(list(models = list(), dataset = data.frame(),
                                          prob_preds = list(), ovr_models = list())))

  testthat::skip_on_cran()
  set.seed(32)
  n <- 14; tr <- ape::rcoal(n)
  d <- data.frame(Species = tr$tip.label, y = rnorm(n), x1 = rnorm(n))
  d$y[1:3] <- NA
  a <- list(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
            phylo = tr, data = d, runs = 1, nitt = 300, burnin = 100,
            thin = 2, n_final = 3, skip_conv = TRUE, verbose = FALSE, seed = 9)
  dir <- tmpdir_fresh()
  A <- suppressWarnings(do.call(bace, c(a, list(checkpoint_dir = dir))))
  # plant a readable but wrong file where run 2 lives
  saveRDS(data.frame(a = 1:3), file.path(dir, "final_run0002.rds"))
  B <- suppressWarnings(do.call(bace, c(a, list(checkpoint_dir = dir))))
  # run 2 was recomputed, not swallowed: all three imputations are present
  testthat::expect_length(B$final_results$all_models, 3L)
  testthat::expect_false(any(vapply(B$final_results$all_models,
                                    function(m) is.null(m$y), logical(1))))
  testthat::expect_identical(B$imputed_datasets, A$imputed_datasets)
})

testthat::test_that("the converged dataset is part of the final-phase fingerprint", {
  # Every final run starts from it, and the chained phase is unseeded, so two
  # runs of that phase give different starting points under identical
  # arguments. Without this field a resume could mix runs from both.
  tr <- ape::rcoal(8)
  d <- data.frame(Species = tr$tip.label, y = c(NA, rnorm(7)), x1 = rnorm(8))
  d2 <- d; d2$y[4] <- NA
  testthat::expect_false(identical(.hash_obj(d), .hash_obj(d2)))

  dir <- .ckpt_dir_prepare(tmpdir_fresh())
  mf <- .ckpt_manifest_file(dir)
  m1 <- list(converged_data = .hash_obj(d), seed = 1)
  .ckpt_manifest_check(mf, m1, "final-imputation", verbose = FALSE)
  m2 <- list(converged_data = .hash_obj(d2), seed = 1)
  testthat::expect_error(
    .ckpt_manifest_check(mf, m2, "final-imputation", verbose = FALSE),
    "converged_data")
})

testthat::test_that("a truncated checkpoint is detected and recomputed, not trusted", {
  testthat::skip_on_cran()
  set.seed(33)
  n <- 14; tr <- ape::rcoal(n)
  d <- data.frame(Species = tr$tip.label, y = rnorm(n), x1 = rnorm(n))
  d$y[1:3] <- NA
  a <- list(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
            phylo = tr, data = d, runs = 1, nitt = 300, burnin = 100,
            thin = 2, n_final = 2, skip_conv = TRUE, verbose = FALSE, seed = 4)
  dir <- tmpdir_fresh()
  A <- suppressWarnings(do.call(bace, c(a, list(checkpoint_dir = dir))))
  f <- file.path(dir, "final_run0001.rds")
  writeBin(readBin(f, "raw", 40), f)             # crash mid-write
  B <- suppressWarnings(do.call(bace, c(a, list(checkpoint_dir = dir))))
  testthat::expect_identical(B$imputed_datasets, A$imputed_datasets)
})

testthat::test_that("the pooled subsample is reproducible, so a resume matches", {
  # sample_size subsamples the pooled draws. A resumed run skips the chained
  # phase, so the ambient random stream sits elsewhere; seeding the subsample
  # is what keeps the two identical.
  testthat::skip_on_cran()
  set.seed(34)
  n <- 14; tr <- ape::rcoal(n)
  d <- data.frame(Species = tr$tip.label, y = rnorm(n), x1 = rnorm(n))
  d$y[1:3] <- NA
  a <- list(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
            phylo = tr, data = d, runs = 1, nitt = 400, burnin = 100,
            thin = 2, n_final = 2, skip_conv = TRUE, verbose = FALSE,
            seed = 42, sample_size = 30)
  dir <- tmpdir_fresh()
  A <- suppressWarnings(do.call(bace, c(a, list(checkpoint_dir = dir))))
  B <- suppressWarnings(do.call(bace, c(a, list(checkpoint_dir = dir))))  # full resume
  testthat::expect_identical(B$pooled_models$models$y$Sol,
                             A$pooled_models$models$y$Sol)
})

testthat::test_that("each convergence-retry attempt saves to its own file", {
  # round(runs * 1.5^k) can repeat a value (runs = 1 gives 1, 2, 2), which
  # would make attempt 3 silently replay attempt 2's saved chain.
  f2a <- .ckpt_chain_files("d", 2, tag = "a02")
  f2b <- .ckpt_chain_files("d", 2, tag = "a03")
  testthat::expect_false(identical(f2a$object, f2b$object))
  testthat::expect_match(basename(f2a$object), "a02")
})

testthat::test_that(".with_seed is deterministic and restores the caller's stream", {
  set.seed(1); before <- .Random.seed
  x <- .with_seed(99, runif(3))
  testthat::expect_identical(.Random.seed, before)
  testthat::expect_identical(x, .with_seed(99, runif(3)))
  testthat::expect_false(identical(x, .with_seed(100, runif(3))))
  # NULL seed is a pass-through
  testthat::expect_length(.with_seed(NULL, runif(2)), 2L)
})

testthat::test_that("checkpoint_dir rejects malformed paths clearly", {
  testthat::expect_error(.ckpt_dir_prepare(c("a", "b")), "single directory")
  testthat::expect_error(.ckpt_dir_prepare(NA_character_), "single directory")
  testthat::expect_error(.ckpt_dir_prepare(""), "single directory")
  testthat::expect_error(.ckpt_dir_prepare(42), "single directory")
})
