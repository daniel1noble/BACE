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
  testthat::expect_true(file.exists(file.path(dir, "chain_runs0002.rds")))

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
