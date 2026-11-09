# code/tests-rubin/test-completecase.R
#
# WHAT THIS GUARDS: the complete-case arm -- the competitor the manuscript's first claim is
# measured against. If this arm is wrong, Claim 1 is wrong, in whichever direction the error runs.
#
# The checks are deliberately boring and hand-checkable:
#   1. the survivor rule is exactly "no NA in any trait", including the edge cases
#   2. it refuses, rather than inventing, when too few species survive
#   3. attrition behaves as the arithmetic says it should under each gap mechanism
#   4. the result table has the same shape as every other arm's, so aggregation cannot break

code_dir <- normalizePath(file.path(testthat::test_path(), ".."))
source(file.path(code_dir, "campaign_gnn_off_lib.R"))
source(file.path(code_dir, "rubin_lib.R"))
source(file.path(code_dir, "rubin_completecase.R"))

# ---- 1. the survivor rule, on a dataset small enough to check by eye ---------------------------
testthat::test_that("complete_species keeps exactly the rows with no gaps", {
  df <- data.frame(
    c1 = c(1,  2, NA,  4,  5),
    c2 = c(1, NA,  3,  4,  5),
    d1 = c(9,  9,  9,  9,  9),            # the MAR driver: fully observed, must not count
    row.names = c("sp1", "sp2", "sp3", "sp4", "sp5"))

  # sp1 complete, sp2 missing c2, sp3 missing c1, sp4 and sp5 complete
  testthat::expect_equal(complete_species(df), c("sp1", "sp4", "sp5"))

  # a column that IS a trait does count against a species
  df2 <- df; df2$c3 <- c(NA, 1, 1, 1, 1)
  testthat::expect_equal(complete_species(df2), c("sp4", "sp5"))

  # if d1 were not excluded, a gap in it would wrongly disqualify a species -- check the exclusion
  # is doing real work rather than being decorative
  df3 <- df; df3$d1[1] <- NA
  testthat::expect_equal(complete_species(df3, ignore = "d1"), c("sp1", "sp4", "sp5"))
  testthat::expect_equal(complete_species(df3, ignore = character(0)), c("sp4", "sp5"))
})

# ---- 2. it refuses rather than producing a meaningless number ----------------------------------
testthat::test_that("score_complete_case declines when too few species survive", {
  cl <- make_cell("types_mixed", 60, 1, lambda = 0.7, rho = 0.5, thresholds = "fixed",
                  driver = TRUE, miss = "mcar", miss_frac = 0.30)

  strict <- score_complete_case(cl$df_miss, cl$tree, truth_rho = 0.5, min_species = 10000L)
  testthat::expect_equal(strict$status, "too_few")
  testthat::expect_null(strict$estimands)
  # the survivor count is still reported when it declines -- that number IS the result
  testthat::expect_true(is.finite(strict$n_complete))
  testthat::expect_equal(strict$frac_kept, strict$n_complete / nrow(cl$df_miss))
})

# ---- 3. attrition matches the arithmetic, and differs by mechanism in the expected direction ---
testthat::test_that("survival rates follow the missingness mechanism", {
  rate <- function(miss, seeds = 1:12) mean(vapply(seeds, function(s) {
    cl <- make_cell("types_mixed", 150, s, lambda = 0.7, rho = 0.5, thresholds = "fixed",
                    driver = TRUE, miss = miss, miss_frac = 0.30)
    length(complete_species(cl$df_miss)) / nrow(cl$df_miss)
  }, numeric(1)))

  mcar <- rate("mcar"); mar <- rate("mar"); clade <- rate("clade")

  # MCAR: 7 traits, each independently 30% missing, so a species survives with probability
  # 0.7^7 = 0.082. This is the number that makes complete-case hopeless at small n, and the
  # manuscript quotes it, so it is pinned here.
  testthat::expect_equal(mcar, 0.7^7, tolerance = 0.03)

  # MAR: gaps are driven by one variable, so they cluster within species and MORE survive.
  testthat::expect_gt(mar, mcar)

  # Clade: whole clades lose everything at once, so the species outside those clades are
  # untouched and most survive -- which is why clade is where the BIAS story lives rather than
  # the attrition story.
  testthat::expect_gt(clade, 0.5)
  testthat::expect_gt(clade, mar)
})

# ---- 4. the output has the same shape as the other arms ----------------------------------------
testthat::test_that("the estimand table matches the shape aggregation expects", {
  cl <- make_cell("types_mixed", 150, 2, lambda = 0.7, rho = 0.5, thresholds = "fixed",
                  driver = TRUE, miss = "clade", miss_frac = 0.30)
  res <- score_complete_case(cl$df_miss, cl$tree, truth_rho = 0.5)

  testthat::expect_equal(res$status, "scored")
  testthat::expect_equal(nrow(res$estimands), 2L)
  testthat::expect_equal(res$estimands$estimand, c("slope", "cor"))
  testthat::expect_true(all(res$estimands$arm == "complete_case"))
  testthat::expect_equal(names(res$estimands),
    c("arm", "estimand", "estimate", "se", "lower", "upper", "df", "fmi", "m_ok",
      "truth", "complete_data", "covered"))

  # the interval must contain its own point estimate and be the right way round
  testthat::expect_true(all(res$estimands$lower < res$estimands$estimate))
  testthat::expect_true(all(res$estimands$estimate < res$estimands$upper))

  # degrees of freedom come from the SURVIVORS, not from the full dataset -- crediting the arm
  # with information it never had would make its intervals too narrow and flatter its coverage
  testthat::expect_equal(res$estimands$df[1], res$n_complete - 2)
  testthat::expect_lt(res$n_complete, nrow(cl$df_miss))
})
