# script/tests-rubin/test-bace.R
# Gate G-S3 for the BACE arm in script/rubin_bace.R: fit_bace_mi() + mi_bace_shipped() (BACE's
# own n_final datasets, untouched -- "Current BACE" against the installed, corrected package).
# The 2026-09 mi_bace_resid() tests were removed with that arm (git: eaca7e3).
#
# ONE tiny shared BACE fit for all three expectations (n = 40, seed = 1, lambda = 0.7, rho = 0.5,
# nitt = 3000, burnin = 1000, thin = 10, runs = 2, M = 5) -- BACE itself dominates runtime, so
# everything downstream reuses the single fit rather than re-fitting per expectation.

source(normalizePath(file.path(testthat::test_path(), "..", "campaign_gnn_off_lib.R")))
source(normalizePath(file.path(testthat::test_path(), "..", "rubin_bace.R")))

testthat::skip_if_not_installed("BACE")

suppressPackageStartupMessages({ library(ape) })

cell <- make_cell("types_mixed", 40, seed = 1, lambda = 0.7, rho = 0.5)
df_miss <- cell$df_miss

t_fit0 <- Sys.time()
fit <- suppressWarnings(fit_bace_mi(cell, M = 5L, nitt = 3000, burnin = 1000, thin = 10, runs = 2))
cat(sprintf("[test-bace] BACE tiny fit (n=40, nitt=3000, burnin=1000, thin=10, runs=2, M=5): %.1fs "
           , as.numeric(difftime(Sys.time(), t_fit0, units = "secs"))))
cat(sprintf("(fit$wall_s = %.1fs, converged = %s, ess_med = %.0f)\n",
           fit$wall_s, fit$converged, fit$ess_med %||% NA_real_))

shipped <- mi_bace_shipped(fit$outb, df_miss)
M <- length(fit$outb$imputed_datasets)

testthat::test_that("fit_bace_mi returns the bace object, wall time, and BACE's own diagnostics", {
  testthat::expect_s3_class(fit$outb, "bace_complete")
  testthat::expect_length(fit$outb$imputed_datasets, 5L)
  testthat::expect_true(is.numeric(fit$wall_s) && fit$wall_s > 0)
  testthat::expect_identical(fit$converged, isTRUE(fit$outb$converged))
  bd <- bace_diagnostics(fit$outb)
  testthat::expect_identical(fit$diag[names(bd)], bd)   # plus the input-cleaning record
  testthat::expect_true(is.character(fit$diag$input_fix))
})

testthat::test_that("mi_bace_shipped is outb$imputed_datasets, aligned to df_miss row order/columns only", {
  testthat::expect_type(shipped, "list")
  testthat::expect_length(shipped, M)
  for (i in seq_len(M)) {
    raw <- fit$outb$imputed_datasets[[i]]
    raw <- raw[rownames(df_miss), names(df_miss), drop = FALSE]
    s <- shipped[[i]]
    testthat::expect_identical(rownames(s), rownames(df_miss))
    testthat::expect_identical(names(s), names(df_miss))
    for (v in names(df_miss)) {
      # compare on values only -- alignment is allowed to change storage class (e.g. BACE's
      # rounded-double count predictions vs df_miss's integer class), never the value
      testthat::expect_equal(as.character(s[[v]]), as.character(raw[[v]]),
                             info = sprintf("dataset %d, trait %s", i, v))
    }
  }
})

testthat::test_that("Current BACE's imputations differ across datasets at missing predictor cells (the 2026-10 fix)", {
  # The 2026-09 package left every non-response trait's missing cells identical across the M
  # datasets (one shared point-filled start). The corrected package sweeps each run, so for every
  # trait with missing data the M imputed values at a missing cell must VARY across datasets.
  for (v in names(df_miss)) {
    miss_idx <- which(is.na(df_miss[[v]]))
    if (!length(miss_idx)) next
    draws <- sapply(shipped, function(s) as.character(s[[v]][miss_idx[1]]))
    testthat::expect_gt(length(unique(draws)), 1L,
                        label = sprintf("distinct imputed values across %d datasets for trait %s", M, v))
  }
})
