# code/tests-rubin/test-discrete.R
#
# Scoring for the discrete traits, checked on inputs whose right answer is known by hand:
# fill_degenerate on a one-class trait, accuracy and the 95% category set on perfect and on tied
# draws, and the c1 ~ bin downstream estimand.
#
# (The castor Mk tests that used to live here moved with the code they test, to
# archive-2026-09/test-discrete-castor.R, when the frequentist arms left the campaign.)
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages({
  source(file.path(root, "script", "campaign_gnn_off_lib.R")); source(file.path(root, "script", "rubin_lib.R"))
  source(file.path(root, "script", "rubin_discrete.R"))
})
old_kind <- RNGkind()[1]; RNGkind("Mersenne-Twister"); on.exit(RNGkind(old_kind), add = TRUE)
set.seed(11)

cell <- make_cell("types_mixed", 60L, 7L, lambda = 0.7, rho = 0.5, thresholds = "fixed", driver = TRUE)
cell$df_miss <- cell$truth; cell$df_miss[cell$mask] <- NA
for (v in names(cell$truth)) { x <- cell$truth[[v]]; x[cell$mask[, v]] <- NA; cell$df_miss[[v]] <- x }


testthat::test_that("fill_degenerate imputes a single observed class and counts the cells", {
  df <- data.frame(b = factor(c("yes", "yes", NA, NA), levels = c("no", "yes")))
  f <- fill_degenerate(list(df, df), df, "b")
  testthat::expect_true(all(f$sets[[2]]$b == "yes")); testthat::expect_equal(f$n_filled[["b"]], 2L)   # cells, not cell-draws
})

testthat::test_that("score_discrete: perfect draws score 1/0/1/1; a 50:50 binary tie counts half", {
  s <- score_discrete("x", rep(list(cell$truth), 20), cell$truth, cell$mask)
  testthat::expect_equal(s$accuracy, rep(1, 3)); testthat::expect_equal(s$brier, rep(0, 3))
  testthat::expect_equal(s$set_coverage, rep(1, 3)); testthat::expect_equal(s$set_size, rep(1, 3))
  flip <- cell$truth; flip$bin <- factor(ifelse(flip$bin == "yes", "no", "yes"), levels = levels(flip$bin))
  s2 <- score_discrete("x", c(rep(list(cell$truth), 10), rep(list(flip), 10)), cell$truth, cell$mask, "bin")
  testthat::expect_equal(s2$accuracy, 0.5); testthat::expect_equal(s2$brier, 0.5); testthat::expect_equal(s2$set_size, 2)
})

testthat::test_that("est_pgls_xy codes a two-level factor as 0/1 and matches est_pgls_slope_fast", {
  d <- data.frame(row.names = rownames(cell$truth), c1 = as.numeric(cell$truth$bin == "yes"), c2 = cell$truth$c1)
  testthat::expect_equal(est_pgls_xy(cell$truth, cell$tree, y = "c1", x = "bin")$estimate,
                         est_pgls_slope_fast(d, cell$tree)$estimate)
})

testthat::test_that("a trait with one class in the complete data is skipped by both scores", {
  tr <- cell$truth; tr$bin <- factor(rep("no", nrow(tr)))
  s <- score_discrete("x", rep(list(tr), 3), tr, cell$mask)
  testthat::expect_false("bin" %in% s$trait)
  testthat::expect_null(score_discrete_estimand("x", rep(list(tr), 3), tr, cell$tree, NULL, nrow(tr)))
})

testthat::test_that("the c1 ~ bin row carries the complete-data interval (t, n - 2 df)", {
  e <- score_discrete_estimand("x", rep(list(cell$truth), 3), cell$truth, cell$tree, NULL, nrow(cell$truth))
  testthat::expect_equal(e$complete_upper - e$complete_data, stats::qt(0.975, nrow(cell$truth) - 2) * e$complete_se)
  testthat::expect_equal(e$estimate, e$complete_data)
})



testthat::test_that("score_discrete: fractional ties at the edge of the 95% set, one-class traits skipped, all-NA cells counted", {
  tr <- data.frame(row.names = paste0("s", 1:3), k = factor(c("A", "B", "A"), levels = c("A", "B", "C")),
                   o = factor(c("L2", "L2", "L2"), levels = c("L1", "L2", "L3"), ordered = TRUE))
  mk <- matrix(c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE), 3, dimnames = list(rownames(tr), c("k", "o")))
  draw <- function(cls2, cls3) { d <- tr; d$k <- factor(c("A", cls2, cls3), levels = levels(tr$k)); d }
  sets <- c(rep(list(draw("A", NA)), 18), list(draw("B", NA)), list(draw("C", NA)))   # cell 2: 18 A, 1 B, 1 C; cell 3: all NA
  s <- score_discrete("x", sets, tr, mk)
  testthat::expect_equal(nrow(s), 1L)                        # o has one realised class in the complete data: skipped
  testthat::expect_equal(s$n_cells, 1L); testthat::expect_equal(s$n_cells_all_na, 1L)
  testthat::expect_equal(s$set_size, 2); testthat::expect_equal(s$set_coverage, 0.5)   # B tied with C for the second place
  testthat::expect_equal(s$accuracy, 0)
})

testthat::test_that("c1 ~ bin: undefined when bin is constant in every imputed dataset; per-dataset target is rho x GLS slope of L5", {
  one <- cell$truth; one$bin[] <- levels(one$bin)[2L]
  testthat::expect_null(score_discrete_estimand("x", rep(list(one), 3), cell$truth, cell$tree, NULL, nrow(one)))
  eig <- pagel_eigen(cell$tree, rownames(cell$truth)); L5 <- cell$L[rownames(cell$truth), 5]
  e <- score_discrete_estimand("x", rep(list(cell$truth), 3), cell$truth, cell$tree, eig, nrow(one), L5 = L5,
                               lambda = 0.7, rho = 0.5)
  x01 <- as.numeric(cell$truth$bin == "yes"); V <- 0.7 * cov2cor(ape::vcv(cell$tree))[rownames(cell$truth), rownames(cell$truth)] + 0.3 * diag(nrow(one))
  X <- cbind(1, x01); b <- solve(t(X) %*% solve(V, X), t(X) %*% solve(V, L5))[2]
  testthat::expect_equal(e$target_cond, 0.5 * b, tolerance = 1e-8)
  e0 <- score_discrete_estimand("x", rep(list(cell$truth), 3), cell$truth, cell$tree, eig, nrow(one), L5 = L5, lambda = 0.7, rho = 0)
  testthat::expect_equal(e0$target_cond, 0)
})

