# Gate G-S3b: the chained BACE arm (mi_bace_chain) and Current BACE both vary the conditioning
# values across imputations. 2026-09 version of this gate asserted the OPPOSITE for the shipped
# arm (identical c1 design matrices across runs) -- that was the improper-MI defect's signature,
# fixed in BACE 0.1.0 (2026-10-26): each final run now sweeps, so c1's predictors differ across
# runs in BOTH arms. The old assertion is inverted and kept as a regression test of the fix.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "script", "campaign_gnn_off_lib.R"))
source(file.path(root, "script", "rubin_bace.R"))

cell <- make_cell("types_mixed", n = 40, seed = 1L, lambda = 0.7, rho = 0.5)
M <- 4L
set.seed(11L)
fb <- fit_bace_mi(cell, M = M, nitt = 3000L, burnin = 1000L, thin = 10L, runs = 2L)
ch <- mi_bace_chain(fb, cell$df_miss, M = M)
message(sprintf("[test-bace-chain] shipped fit %.1fs, chain %.1fs (M = %d)", fb$wall_s, ch$wall_s, M))

testthat::test_that("chained arm returns M datasets with observed cells untouched", {
  testthat::expect_length(ch$datasets, M)
  obs <- !is.na(cell$df_miss)
  for (d in ch$datasets) {
    for (v in names(cell$df_miss)) {
      o <- obs[, v]
      testthat::expect_identical(as.character(d[[v]][o]), as.character(cell$df_miss[[v]][o]))
    }
  }
})

testthat::test_that("both Current BACE and chained runs vary c1's design matrix across runs", {
  X_ship <- lapply(fb$outb$final_results$all_models, function(m) as.matrix(m[["c1"]]$X))
  X_chain <- lapply(ch$models, function(m) as.matrix(m[["c1"]]$X))
  same <- function(L) all(vapply(L[-1], function(x) isTRUE(all.equal(x, L[[1]])), logical(1)))
  testthat::expect_false(same(X_ship))     # regression test of the 2026-10 fix: no shared anchor
  testthat::expect_false(same(X_chain))    # the chain varies c1's predictors across runs
})
