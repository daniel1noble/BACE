# code/tests-rubin/test-campaign-grid.R
#
# WHAT THIS GUARDS: that two different scenarios never write the same file.
#
# Why that matters more than it sounds. The campaign is resumable: rubin_cell.R skips any dataset
# whose .rds already exists. That is what lets a three-week run be stopped and restarted. It also
# means that if two scenarios produce the same filename, the second one is skipped as "already
# done" and NEVER RUNS -- no error, no warning, and the results file for it is simply the other
# scenario's, wearing its name.
#
# The 2026-11 design walks straight into this: the "typical" regime (traits correlated 0.5, no
# covariate) and the "good covariate" regime (traits correlated 0.5, plus a covariate) share
# lambda, rho, n, M and seed. Only the covariate differs. So the covariate has to be in the name.

code_dir <- normalizePath(file.path(testthat::test_path(), ".."))
source(file.path(code_dir, "rubin_lib.R"))

# ---- 1. every scenario of the planned design gets its own filename -----------------------------
testthat::test_that("no two scenarios in the 2026-11 design share a filename", {
  # the design, exactly as rubin_campaign_local.sh writes it
  regimes <- list(list(rho = 0,   n_cov = 0L, rho_cov = 0),
                  list(rho = 0.5, n_cov = 0L, rho_cov = 0),
                  list(rho = 0.5, n_cov = 1L, rho_cov = 0.8))

  tags <- character(0)
  for (n in c(100, 500))
    for (lambda in c(0.3, 0.7, 0.9))
      for (rg in regimes)
        for (miss in c("mcar", "mar", "clade"))
          for (seed in 1:3)
            tags <- c(tags, cell_tag(lambda, rg$rho, miss, 0.30, n, 20, seed,
                                     n_cov = rg$n_cov, rho_cov = rg$rho_cov))

  testthat::expect_equal(length(tags), 2 * 3 * 3 * 3 * 3)   # 54 scenarios x 3 seeds
  testthat::expect_equal(length(unique(tags)), length(tags))
})

# ---- 2. the specific collision this was written to prevent -------------------------------------
testthat::test_that("the two rho = 0.5 regimes are distinguishable by filename alone", {
  typical <- cell_tag(0.3, 0.5, "mcar", 0.30, 100, 20, 1, n_cov = 0L)
  good    <- cell_tag(0.3, 0.5, "mcar", 0.30, 100, 20, 1, n_cov = 1L, rho_cov = 0.8)

  testthat::expect_false(identical(typical, good))
  testthat::expect_true(grepl("_cov1-0.8_", good, fixed = TRUE))
  testthat::expect_false(grepl("cov", typical))
})

# ---- 3. names written before 2026-11 are unchanged, so old results still resume -----------------
testthat::test_that("n_cov = 0 reproduces the pre-2026-11 filename exactly", {
  testthat::expect_equal(
    cell_tag(0.7, 0.5, "mcar", 0.30, 100, 20, 1),
    "rubin_types_mixed_BM_l0.7_r0.5_mcar0.3_n100_M20_s1")
  # a real filename from the 2026-10 campaign, reproduced from its settings
  testthat::expect_equal(
    cell_tag(0.3, 0, "mcar", 0.30, 1000, 20, 7),
    "rubin_types_mixed_BM_l0.3_r0_mcar0.3_n1000_M20_s7")
  # "0.30" as a string and 0.3 as a number must give the same name, or resume breaks on a
  # harmless difference in how the launcher spelled the number
  testthat::expect_equal(cell_tag(0.7, 0.5, "mcar", "0.30", "100", "20", "1"),
                         cell_tag(0.7, 0.5, "mcar", 0.30, 100, 20, 1))
})

# ---- 4. every setting that changes the answer is in the name -----------------------------------
testthat::test_that("changing any result-affecting setting changes the filename", {
  base <- cell_tag(0.7, 0.5, "mcar", 0.30, 100, 20, 1, n_cov = 1L, rho_cov = 0.8)
  variants <- c(
    lambda  = cell_tag(0.9, 0.5, "mcar", 0.30, 100, 20, 1, n_cov = 1L, rho_cov = 0.8),
    rho     = cell_tag(0.7, 0.0, "mcar", 0.30, 100, 20, 1, n_cov = 1L, rho_cov = 0.8),
    miss    = cell_tag(0.7, 0.5, "clade", 0.30, 100, 20, 1, n_cov = 1L, rho_cov = 0.8),
    frac    = cell_tag(0.7, 0.5, "mcar", 0.50, 100, 20, 1, n_cov = 1L, rho_cov = 0.8),
    n       = cell_tag(0.7, 0.5, "mcar", 0.30, 500, 20, 1, n_cov = 1L, rho_cov = 0.8),
    M       = cell_tag(0.7, 0.5, "mcar", 0.30, 100, 50, 1, n_cov = 1L, rho_cov = 0.8),
    seed    = cell_tag(0.7, 0.5, "mcar", 0.30, 100, 20, 2, n_cov = 1L, rho_cov = 0.8),
    n_cov   = cell_tag(0.7, 0.5, "mcar", 0.30, 100, 20, 1, n_cov = 0L, rho_cov = 0.8),
    rho_cov = cell_tag(0.7, 0.5, "mcar", 0.30, 100, 20, 1, n_cov = 1L, rho_cov = 0.5),
    smoke   = cell_tag(0.7, 0.5, "mcar", 0.30, 100, 20, 1, n_cov = 1L, rho_cov = 0.8, smoke = TRUE))

  for (nm in names(variants))
    testthat::expect_false(identical(base, unname(variants[[nm]])),
                           info = paste("changing", nm, "did not change the filename"))
})
