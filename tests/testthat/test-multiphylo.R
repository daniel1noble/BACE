# Tree-per-imputation (multiPhylo) and the per-tree inverseA cache.
#
# Design contract (see .agents/plans/next-features-2026-11.md, feature 1):
# tree selection is DETERMINISTIC in the run index (no RNG consumed), the
# convergence chain conditions on trees[[1]] only, final run m is fit under
# tree tree_index[m], and pool_posteriors aligns tip BLUPs by name / drops
# node BLUPs when the runs' trees differ. The identity test at the bottom is
# the load-bearing one: it can only pass if the per-run tree actually reaches
# the fits AND no extra RNG is consumed anywhere on the path.

# ---------------------------------------------------------------------------
# CRAN-safe: selection arithmetic
# ---------------------------------------------------------------------------

testthat::test_that(".select_tree_indices: subsequence, recycling, edge cases", {
  # more trees than runs: evenly spaced subsequence, no duplicates,
  # spanning the whole file
  idx <- .select_tree_indices(100L, 10L)
  testthat::expect_length(idx, 10L)
  testthat::expect_identical(anyDuplicated(idx), 0L)
  testthat::expect_identical(idx[1], 1L)
  testthat::expect_identical(idx[10], 100L)

  # fewer trees than runs: balanced recycling
  testthat::expect_identical(.select_tree_indices(2L, 4L), c(1L, 2L, 1L, 2L))
  idx2 <- .select_tree_indices(20L, 50L)
  testthat::expect_true(all(table(idx2) %in% c(2L, 3L)))  # balanced

  # equal: identity
  testthat::expect_identical(.select_tree_indices(5L, 5L), 1:5)
  # single tree: all ones
  testthat::expect_identical(.select_tree_indices(1L, 7L), rep(1L, 7L))
  # no duplicates for a sweep of sizes (subsequence branch)
  for (K in c(3L, 7L, 50L, 123L)) {
    for (M in c(2L, 3L, K)) {
      if (K >= M) testthat::expect_identical(
        anyDuplicated(.select_tree_indices(K, M)), 0L)
    }
  }
})

# ---------------------------------------------------------------------------
# CRAN-safe: validator (runs inverseA on tiny trees -- no MCMC)
# ---------------------------------------------------------------------------

testthat::test_that(".validate_trees: classes, tip coverage, labels, branch lengths", {
  set.seed(1)
  trA <- ape::rcoal(8)
  trB <- ape::rcoal(8); trB$tip.label <- trA$tip.label
  sp <- trA$tip.label

  # single phylo, multiPhylo, and plain list are all accepted
  testthat::expect_length(.validate_trees(trA, sp, verbose = FALSE), 1L)
  testthat::expect_length(
    suppressWarnings(.validate_trees(c(trA, trB), sp, verbose = FALSE)), 2L)
  testthat::expect_length(
    suppressWarnings(.validate_trees(list(trA, trB), sp, verbose = FALSE)), 2L)
  # junk is refused
  testthat::expect_error(.validate_trees(list(trA, 1), sp),
                         "phylo")

  # tip-mismatch errors NAME the offending tree
  trC <- trB; trC$tip.label[2] <- "not_in_data"
  testthat::expect_error(
    suppressWarnings(.validate_trees(c(trA, trC), sp, verbose = FALSE)),
    "tree 2 of 2")

  # node labels (e.g. BEAST support values, possibly duplicated) are stripped
  # with a message; duplicated node labels must NOT reach inverseA
  trD <- trA; trD$node.label <- rep("0.95", trD$Nnode)
  testthat::expect_message(
    out <- .validate_trees(trD, sp, verbose = TRUE), "node labels stripped")
  testthat::expect_null(out[[1]]$node.label)

  # negative branch lengths error; zero-length branches error with the remedy
  # (inverseA rejects them unconditionally, so this breaks nothing that works
  # today -- it just fails at entry, named, instead of deep inside a fit)
  trE <- trA; trE$edge.length[1] <- -0.01
  testthat::expect_error(.validate_trees(trE, sp, verbose = FALSE), "negative")
  trF <- trA
  e_int <- which(trF$edge[, 2] > ape::Ntip(trF))[1]
  child <- trF$edge[e_int, 2]; len <- trF$edge.length[e_int]
  trF$edge.length[e_int] <- 0
  kids <- which(trF$edge[, 1] == child)
  trF$edge.length[kids] <- trF$edge.length[kids] + len  # stays ultrametric
  testthat::expect_true(ape::is.ultrametric(trF))
  testthat::expect_error(.validate_trees(trF, sp, verbose = FALSE), "1e-8")
})

testthat::test_that(".validate_trees: inverseA rejections surface at entry, named", {
  set.seed(2)
  trA <- ape::rcoal(8)
  # a non-ultrametric tree (inverseA's default scale = TRUE hard-errors);
  # rtree() without ultrametricization is the simplest reproduction
  trNU <- ape::rtree(8); trNU$tip.label <- trA$tip.label
  testthat::expect_false(ape::is.ultrametric(trNU))
  testthat::expect_error(
    .validate_trees(c(trA, trNU), trA$tip.label, verbose = FALSE),
    "tree 2 of 2.*inverseA")
})

testthat::test_that(".validate_trees: warns on inconsistent branch-length scales", {
  set.seed(3)
  trA <- ape::rcoal(8)
  trG <- trA; trG$edge.length <- trG$edge.length * 3  # a rescaled copy
  testthat::expect_warning(
    .validate_trees(c(trA, trG), trA$tip.label, verbose = FALSE),
    "root depths")
})

# ---------------------------------------------------------------------------
# CRAN-safe: the Ainv cache
# ---------------------------------------------------------------------------

testthat::test_that(".get_ainv caches per tree; key distinguishes topology and node labels", {
  set.seed(4)
  trA <- ape::rcoal(8)
  trB <- ape::rcoal(8); trB$tip.label <- trA$tip.label

  .clear_ainv_cache()
  n_computes <- 0L
  testthat::local_mocked_bindings(
    .compute_ainv = function(tree) {
      n_computes <<- n_computes + 1L
      MCMCglmm::inverseA(tree, nodes = "ALL")$Ainv
    },
    .package = "BACE"
  )

  a1 <- .get_ainv(trA)
  a2 <- .get_ainv(trA)            # cache hit
  testthat::expect_identical(n_computes, 1L)
  testthat::expect_identical(a1, a2)

  b1 <- .get_ainv(trB)            # different topology, same tips: new entry
  testthat::expect_identical(n_computes, 2L)
  testthat::expect_false(identical(a1, b1))

  # same tree, different node labels: a DIFFERENT cache entry (inverseA names
  # node rows after node labels, so sharing an entry would corrupt Sol names)
  trA_lab <- trA; trA_lab$node.label <- paste0("clade", seq_len(trA$Nnode))
  l1 <- .get_ainv(trA_lab)
  testthat::expect_identical(n_computes, 3L)
  testthat::expect_false(identical(rownames(a1), rownames(l1)))
  # values are the same matrix, only names differ
  testthat::expect_equal(unname(as.matrix(a1)), unname(as.matrix(l1)))

  .clear_ainv_cache()
  .get_ainv(trA)
  testthat::expect_identical(n_computes, 4L)  # cache really was cleared
  .clear_ainv_cache()
})

# ---------------------------------------------------------------------------
# CRAN-safe: hashing is content-based
# ---------------------------------------------------------------------------

testthat::test_that(".tree_cache_key depends on content, not object identity", {
  set.seed(5)
  trA <- ape::rcoal(8)
  trA2 <- trA  # copy
  testthat::expect_identical(.tree_cache_key(trA), .tree_cache_key(trA2))
  trB <- trA; trB$edge.length[1] <- trB$edge.length[1] * 1.0001
  testthat::expect_false(.tree_cache_key(trA) == .tree_cache_key(trB))
})

# ---------------------------------------------------------------------------
# CRAN-safe: with_imputations per-imputation trees (lm fits, no MCMC)
# ---------------------------------------------------------------------------

testthat::test_that("with_imputations assigns tree m to imputation m, recycled", {
  set.seed(6)
  trA <- ape::rcoal(8); trB <- ape::rcoal(8); trB$tip.label <- trA$tip.label
  imps <- lapply(1:5, function(i) {
    d <- data.frame(x1 = rnorm(20)); d$y <- d$x1 + rnorm(20); d
  })
  depth <- function(tr) round(sum(tr$edge.length), 10)

  seen <- new.env(); seen$d <- numeric(0)
  fits <- with_imputations(imps, function(d, tree) {
    seen$d <- c(seen$d, depth(tree)); stats::lm(y ~ x1, data = d)
  }, tree = c(trA, trB), .progress = FALSE)
  # recycling 1,2,1,2,1
  testthat::expect_identical(seen$d,
    c(depth(trA), depth(trB), depth(trA), depth(trB), depth(trA)))

  # single tree: same tree every fit (previous behaviour)
  seen$d <- numeric(0)
  fits <- with_imputations(imps, function(d, tree) {
    seen$d <- c(seen$d, depth(tree)); stats::lm(y ~ x1, data = d)
  }, tree = trA, .progress = FALSE)
  testthat::expect_identical(unique(seen$d), depth(trA))

  # .f without a tree formal: supplying a tree it cannot receive is an ERROR.
  # Until 1.0.1 this silently dropped the tree, and this test asserted that
  # silence ("no injection, no error") -- it was pinning the defect in GitHub
  # issue #18, where a misnamed argument made every fit use the wrong
  # phylogeny with nothing to warn the user.
  testthat::expect_error(
    with_imputations(imps, function(d) stats::lm(y ~ x1, data = d),
                     tree = c(trA, trB), .progress = FALSE),
    regexp = "cannot receive it")

  # ...and a function that genuinely does not want a tree still works when no
  # tree is supplied.
  fits <- with_imputations(imps, function(d) stats::lm(y ~ x1, data = d),
                           .progress = FALSE)
  testthat::expect_s3_class(fits, "bace_mi_fits")
})

# ---------------------------------------------------------------------------
# CRAN-safe: pool_posteriors under differing trees (handcrafted fake models)
# ---------------------------------------------------------------------------

testthat::test_that("pool_posteriors aligns tip BLUPs by name and drops node BLUPs when trees vary", {
  # Two fake runs whose Sol column NAME SETS are identical but ORDERED
  # differently (the reproduced real-world failure mode), with one dotted tip
  # label. Values are chosen so misalignment is detectable.
  tips <- c("t1", "t2", "G.fortis")
  fixed <- c("(Intercept)", "x1")
  tipc  <- paste0("animal.", tips)
  nodec <- c("animal.Node2", "animal.Node3")
  cols1 <- c(fixed, tipc, nodec)
  cols2 <- c(fixed, tipc[c(3, 1, 2)], nodec[c(2, 1)])  # permuted random cols

  mk_model <- function(cols, base) {
    n <- 10L
    Sol <- matrix(rep(base + seq_along(cols), each = n), nrow = n,
                  dimnames = list(NULL, cols))
    Z <- matrix(0, nrow = 3, ncol = length(cols) - 2L,
                dimnames = list(NULL, setdiff(cols, fixed)))
    VCV <- matrix(1, nrow = n, ncol = 2,
                  dimnames = list(NULL, c("animal", "units")))
    list(Sol = coda::as.mcmc(Sol), VCV = coda::as.mcmc(VCV), Z = Z,
         CP = NULL, Liab = NULL, Deviance = NULL, DIC = NULL)
  }
  run1 <- list(y = mk_model(cols1, base = 0))
  run2 <- list(y = mk_model(cols2, base = 100))

  fake <- structure(list(
    all_models = list(run1, run2),
    tree_index = c(1L, 2L),
    tip_labels = tips
  ), class = "bace_final")

  testthat::expect_message(pl <- pool_posteriors(fake), "different trees")
  m <- pl$models$y
  testthat::expect_true(m$BACE_pooling$trees_vary)
  testthat::expect_identical(m$BACE_pooling$n_node_blup_cols_dropped, 2L)
  testthat::expect_false(any(grepl("\\.Node[0-9]+$", colnames(m$Sol))))
  # by-name alignment: run 2's value for animal.t1 must land in the
  # animal.t1 column (run 2 column value = 100 + position in cols2)
  v_run2_t1 <- 100 + which(cols2 == "animal.t1")
  testthat::expect_identical(unique(m$Sol[11:20, "animal.t1"]), v_run2_t1)
  v_run2_gf <- 100 + which(cols2 == "animal.G.fortis")
  testthat::expect_identical(unique(m$Sol[11:20, "animal.G.fortis"]), v_run2_gf)

  # identical trees (tree_index all 1): previous behaviour, node columns kept
  fake1 <- fake; fake1$tree_index <- c(1L, 1L)
  fake1$all_models <- list(run1, run1)
  pl1 <- pool_posteriors(fake1)
  testthat::expect_false(pl1$models$y$BACE_pooling$trees_vary)
  testthat::expect_true(all(nodec %in% colnames(pl1$models$y$Sol)))

  # legacy objects without tree_index: previous behaviour
  fake0 <- fake; fake0$tree_index <- NULL
  fake0$all_models <- list(run1, run1)
  pl0 <- pool_posteriors(fake0)
  testthat::expect_false(pl0$models$y$BACE_pooling$trees_vary)
})

# ---------------------------------------------------------------------------
# CRAN-safe: bace() forwards the FULL tree set under skip_conv = TRUE
# ---------------------------------------------------------------------------

testthat::test_that("bace(skip_conv = TRUE) forwards the full multiPhylo to bace_final_imp", {
  set.seed(7)
  trA <- ape::rcoal(8); trB <- ape::rcoal(8); trB$tip.label <- trA$tip.label
  trees <- c(trA, trB)
  d <- data.frame(Species = trA$tip.label, y = rnorm(8), x1 = rnorm(8))
  d$y[1:2] <- NA

  seen_phylo <- NULL
  fake_final <- structure(list(all_datasets = list(d), all_prob_preds = list(),
                               all_models = list()), class = "bace_final")
  testthat::local_mocked_bindings(
    bace_imp = function(...) structure(list(data = list(d)), class = "bace"),
    assess_convergence = function(...) list(converged = FALSE),
    bace_final_imp = function(bace_object, ..., phylo) {
      seen_phylo <<- phylo
      fake_final
    },
    pool_posteriors = function(...) structure(list(models = list()),
                                              class = "bace_pooled"),
    .package = "BACE"
  )
  res <- suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", phylo = trees, data = d,
         skip_conv = TRUE, verbose = FALSE)))
  testthat::expect_s3_class(res, "bace_complete")
  # the SECOND call site (skip_conv branch) received the full tree set,
  # not a pre-selected single tree
  testthat::expect_true(inherits(seen_phylo, "multiPhylo") ||
                          (is.list(seen_phylo) && length(seen_phylo) == 2L))
})

# ---------------------------------------------------------------------------
# Real tiny-MCMC: the identity test (the load-bearing one) + single-tree
# equivalence + tree metadata
# ---------------------------------------------------------------------------

testthat::test_that("tree m reaches run m exactly: multiPhylo runs match single-tree runs stream-for-stream", {
  testthat::skip_on_cran()
  set.seed(11)
  n <- 25
  trA <- ape::rcoal(n); trB <- ape::rcoal(n); trB$tip.label <- trA$tip.label
  S <- ape::vcv(trA, corr = TRUE); L <- chol(S)
  d <- data.frame(Species = trA$tip.label)
  d$x1 <- as.numeric(t(L) %*% rnorm(n)) + rnorm(n, sd = .3)
  d$x2 <- 0.6 * d$x1 + as.numeric(t(L) %*% rnorm(n)) * .5 + rnorm(n, sd = .3)
  d$x1[sample(n, 5)] <- NA; d$x2[sample(n, 5)] <- NA

  # ONE shared converged object (the chained phase is unseeded by design, so
  # separate invocations would start the final phases from different states
  # and the equalities below could never hold)
  fit <- suppressWarnings(bace_imp("x1 ~ x2", "~ 1 | Species", phylo = trA,
                                   data = d, runs = 2, nitt = 600,
                                   burnin = 100, thin = 2, verbose = FALSE))

  args <- list(fixformula = "x1 ~ x2", ran_phylo_form = "~ 1 | Species",
               n_final = 4, nitt = 600, burnin = 100, thin = 2,
               verbose = FALSE, seed = 7)
  fA  <- suppressWarnings(do.call(bace_final_imp,
           c(list(bace_object = fit, phylo = trA), args)))
  fB  <- suppressWarnings(do.call(bace_final_imp,
           c(list(bace_object = fit, phylo = trB), args)))
  fAB <- suppressWarnings(do.call(bace_final_imp,
           c(list(bace_object = fit, phylo = c(trA, trB)), args)))

  # metadata
  testthat::expect_identical(fAB$tree_index, c(1L, 2L, 1L, 2L))
  testthat::expect_length(fAB$tree_hashes, 2L)
  testthat::expect_identical(fA$tree_index, rep(1L, 4L))

  # the identity property: run m under the multiPhylo is bitwise the run m
  # under that run's tree alone (same seed). Catches both "tree silently
  # ignored" and any stray RNG consumption on the multiPhylo path.
  testthat::expect_identical(fAB$all_datasets[[1]], fA$all_datasets[[1]])
  testthat::expect_identical(fAB$all_datasets[[3]], fA$all_datasets[[3]])
  testthat::expect_identical(fAB$all_datasets[[2]], fB$all_datasets[[2]])
  testthat::expect_identical(fAB$all_datasets[[4]], fB$all_datasets[[4]])
  testthat::expect_false(identical(fA$all_datasets[[1]], fB$all_datasets[[1]]))

  # end-to-end pooling works under differing trees and reports the drop
  pl <- suppressMessages(pool_posteriors(fAB))
  testthat::expect_true(pl$models$x1$BACE_pooling$trees_vary)
  testthat::expect_false(any(grepl("\\.Node[0-9]+$",
                                   colnames(pl$models$x1$Sol))))
})

testthat::test_that("multiPhylo: serial == parallel under a seed", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  set.seed(12)
  n <- 20
  trA <- ape::rcoal(n); trB <- ape::rcoal(n); trB$tip.label <- trA$tip.label
  d <- data.frame(Species = trA$tip.label, x1 = rnorm(n), x2 = rnorm(n))
  d$x1[sample(n, 4)] <- NA; d$x2[sample(n, 4)] <- NA
  fit <- suppressWarnings(bace_imp("x1 ~ x2", "~ 1 | Species", phylo = trA,
                                   data = d, runs = 2, nitt = 600,
                                   burnin = 100, thin = 2, verbose = FALSE))
  args <- list(bace_object = fit, fixformula = "x1 ~ x2",
               ran_phylo_form = "~ 1 | Species", phylo = c(trA, trB),
               n_final = 4, nitt = 600, burnin = 100, thin = 2,
               verbose = FALSE, seed = 9)
  f_ser <- suppressWarnings(do.call(bace_final_imp, c(args, list(n_cores = 1))))
  f_par <- suppressWarnings(do.call(bace_final_imp, c(args, list(n_cores = 2))))
  testthat::expect_identical(f_ser$all_datasets, f_par$all_datasets)
})
