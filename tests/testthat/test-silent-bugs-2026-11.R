# Regressions for the three SILENT correctness bugs fixed in 1.0.1
# (GitHub issues #16, #17, #18). All three produced ordinary-looking numbers
# that were wrong, with no error and no warning -- the class of defect that
# survives review -- so each test asserts the NUMBER, not just that code runs.

# ---- #18: with_imputations() must never drop the tree in silence ----------

make_imputed_list <- function(M = 4, n = 30, seed = 1) {
  set.seed(seed)
  tr <- ape::rcoal(n)
  tr$tip.label <- paste0("t", seq_len(n))
  base <- data.frame(Species = tr$tip.label, x = rnorm(n))
  base$y <- 0.5 * base$x + rnorm(n)
  list(tree = tr,
       data = lapply(seq_len(M), function(m) {
         d <- base; d$y <- d$y + rnorm(n, 0, 0.3); d
       }))
}

test_that("#18 a tree that cannot be delivered is an error, not a silent drop", {
  s <- make_imputed_list()
  other_tree <- ape::rcoal(30)
  other_tree$tip.label <- paste0("t", 1:30)

  # The silent case: the argument has a DEFAULT, so R never complains and
  # every fit quietly used the wrong phylogeny.
  expect_error(
    with_imputations(s$data,
                     function(d, phy = other_tree) stats::lm(y ~ x, data = d),
                     tree = s$tree, .progress = FALSE),
    regexp = "cannot receive it"
  )
  # the message must name the function's real arguments and the fix
  err <- tryCatch(
    with_imputations(s$data,
                     function(d, phy = other_tree) stats::lm(y ~ x, data = d),
                     tree = s$tree, .progress = FALSE),
    error = function(e) conditionMessage(e))
  expect_match(err, "phy")
  expect_match(err, "wrong phylogeny")
})

test_that("#18 a correctly named `tree` argument still receives the tree", {
  s <- make_imputed_list()
  seen <- list()
  fits <- with_imputations(s$data, function(d, tree) {
    seen[[length(seen) + 1L]] <<- tree
    stats::lm(y ~ x, data = d)
  }, tree = s$tree, .progress = FALSE)
  expect_equal(attr(fits, "n_failed"), 0L)
  expect_length(seen, length(s$data))
  expect_true(all(vapply(seen, function(t) identical(t, s$tree), logical(1))))
})

test_that("#18 a function taking ... receives the tree too", {
  s <- make_imputed_list()
  got <- NULL
  with_imputations(s$data, function(d, ...) {
    got <<- list(...)$tree
    stats::lm(y ~ x, data = d)
  }, tree = s$tree, .progress = FALSE)
  expect_identical(got, s$tree)
})

test_that("#18 no tree supplied is unaffected", {
  s <- make_imputed_list()
  expect_silent(
    fits <- with_imputations(s$data, function(d) stats::lm(y ~ x, data = d),
                             .progress = FALSE))
  expect_equal(attr(fits, "n_failed"), 0L)
})

# ---- #17: pool_mi() must not silently skip the small-sample correction ----

test_that("#17 df_fun returning NULL still yields the Barnard-Rubin correction", {
  skip_if_not_installed("nlme")
  s <- make_imputed_list(M = 5)
  fits <- with_imputations(s$data, function(d, tree)
    nlme::gls(y ~ x, data = d,
              correlation = ape::corBrownian(phy = tree, form = ~Species)),
    tree = s$tree, .progress = FALSE)

  # stats::df.residual() returns NULL for gls -- the natural call that did
  # nothing at all before 1.0.1.
  expect_null(stats::df.residual(fits[[1]]))

  uncorrected <- pool_mi(fits)
  expect_message(corrected <- pool_mi(fits, df_fun = stats::df.residual),
                 regexp = "nobs - n_coef")

  # the whole point: the two must now DIFFER
  expect_false(isTRUE(all.equal(uncorrected$df, corrected$df)))

  # and the correction must match an explicitly supplied correct df_fun
  explicit <- pool_mi(fits, df_fun = function(f)
    stats::nobs(f) - length(stats::coef(f)))
  expect_equal(corrected$df, explicit$df)

  # Barnard-Rubin caps df at the complete-data df (n - p = 30 - 2 = 28)
  expect_true(all(corrected$df < 28))
  # the uncorrected df could exceed it wildly -- that was the defect
  expect_true(any(uncorrected$df > 28))
})

test_that("#17 an underivable complete-data df warns instead of passing silently", {
  fake <- structure(list(), class = "bace_fakefit_test")
  registerS3method("coef", "bace_fakefit_test",
                   function(object, ...) c(a = 1, b = 2))
  registerS3method("vcov", "bace_fakefit_test",
                   function(object, ...) diag(2))
  ff <- structure(list(fake, fake, fake),
                  class = c("bace_mi_fits", "list"),
                  n_fits = 3L, n_failed = 0L, failed = integer(0),
                  is_mcmcglmm = FALSE)
  expect_warning(pool_mi(ff, df_fun = stats::df.residual),
                 regexp = "NOT applied")
})

test_that("#17 .complete_data_df matches df.residual where both are defined", {
  set.seed(2)
  d <- data.frame(x = rnorm(30)); d$y <- d$x + rnorm(30)
  m <- stats::lm(y ~ x, data = d)
  expect_equal(BACE:::.complete_data_df(m, 2L), stats::df.residual(m))
})

test_that("#17 .complete_data_df never partial-matches its way to a wrong n", {
  # `$` partial-matches on lists: list(nitt = 13000)$n is 13000. If the
  # fallback used `fit$n` it would read an MCMC chain length as the sample
  # size and report df = 12998 for a 20-observation fit, which makes the
  # Barnard-Rubin correction a no-op -- the exact failure this fallback
  # exists to prevent, reinstated silently.
  expect_equal(list(nitt = 13000)$n, 13000)          # the hazard is real
  expect_null(list(nitt = 13000)[["n", exact = TRUE]])

  chainy <- structure(list(nitt = 13000L), class = "bace_chainy_test")
  expect_true(is.na(BACE:::.complete_data_df(chainy, 2L)))

  grouped <- structure(list(ngroups = 4L), class = "bace_chainy_test")
  expect_true(is.na(BACE:::.complete_data_df(grouped, 2L)))

  # an exactly named component is still honoured
  exact_n <- structure(list(n = 20L), class = "bace_chainy_test")
  expect_equal(BACE:::.complete_data_df(exact_n, 2L), 18L)

  # and `residuals` is matched exactly too, not via `residualVar`
  rv <- structure(list(residualVar = rep(0, 500)), class = "bace_chainy_test")
  expect_true(is.na(BACE:::.complete_data_df(rv, 2L)))
})

test_that("#17 gls -- the headline case -- gets a correct complete-data df", {
  skip_if_not_installed("nlme")
  set.seed(3)
  d <- data.frame(x = rnorm(20)); d$y <- d$x + rnorm(20)
  g <- nlme::gls(y ~ x, data = d)
  # df.residual() is NULL for gls: that NULL is what used to be swallowed.
  expect_null(stats::df.residual(g))
  expect_equal(BACE:::.complete_data_df(g, 2L), 18)
  expect_equal(BACE:::.complete_data_df(g, 2L),
               stats::df.residual(stats::lm(y ~ x, data = d)))
})

# ---- #16: a multinomial predictor must keep its random effects ------------

test_that("#16 a no-dependency multinomial predictor carries the random effects", {
  # The defect was arithmetic and total: with no dependencies every species
  # got an IDENTICAL liability, so the simulated trait could not carry any
  # phylogenetic, species or residual variation at all. This reproduces the
  # branch's arithmetic directly, so it is fast and exact.
  intercept <- 0.3
  set.seed(1)
  u <- as.numeric(scale(rnorm(150)))         # the species-level random effects
  linear_pred <- intercept + u
  k <- 1; n_liab <- 2
  liab_intercept <- intercept + (k - n_liab / 2) * 0.5

  fixed <- liab_intercept + linear_pred - intercept
  expect_gt(stats::sd(fixed), 0.5)           # pre-fix this was exactly 0
  expect_equal(fixed - mean(fixed), u - mean(u))   # it is the random effects
})

test_that("#16 simulated multinomial predictors respond to phylo_signal", {
  testthat::skip_on_cran()
  skip_if_not_installed("phytools")
  sig <- function(ps, seed) {
    set.seed(seed)
    s <- suppressWarnings(sim_bace(
      response_type = "gaussian", predictor_types = "multinomial3",
      phylo_signal = c(ps, ps), missingness = c(0, 0),
      n_cases = 200, n_species = 200))
    v <- s$data[[setdiff(names(s$data), c("y", "Species"))[1]]]
    x <- as.numeric(as.character(v) == sort(unique(as.character(v)))[1])
    names(x) <- s$data$Species
    suppressWarnings(
      phytools::phylosig(s$tree, x[s$tree$tip.label], method = "lambda")$lambda)
  }
  lo <- vapply(1:4, function(s) sig(0.10, s), numeric(1))
  hi <- vapply(1:4, function(s) sig(0.95, s), numeric(1))
  # Observed-scale lambda for a 3-category trait is far below the latent value
  # phylo_signal sets; what the fix restores is that it responds AT ALL.
  expect_gt(mean(hi), mean(lo) + 0.05)
})

# ---- #18b: the tree-per-imputation PAIRING must match what generated it ----

test_that("#18b with more trees than imputations the two rules used to disagree", {
  # bace_final_imp spreads across the supplied trees; with_imputations used to
  # recycle modularly. They agree only when n_trees <= M, so the common case
  # "100 posterior trees, 50 imputations" mispaired 49 of 50 fits.
  expect_identical(as.integer(BACE:::.select_tree_indices(6, 3)), c(1L, 4L, 6L))
  expect_false(identical(as.integer(BACE:::.select_tree_indices(6, 3)),
                         ((seq_len(3) - 1L) %% 6L) + 1L))
  # ...and they agree where they always did, so old behaviour is preserved
  expect_identical(as.integer(BACE:::.select_tree_indices(2, 5)),
                   ((seq_len(5) - 1L) %% 2L) + 1L)
})

test_that("#18b the recorded $tree_index is what the fits actually receive", {
  set.seed(11)
  n <- 8
  trees <- lapply(1:6, function(i) {
    tr <- ape::rcoal(n); tr$tip.label <- paste0("t", seq_len(n))
    tr$edge.length <- tr$edge.length * i   # makes each tree's depth unique
    tr
  })
  depth <- function(tr) round(sum(tr$edge.length), 8)
  depths <- vapply(trees, depth, numeric(1))
  expect_equal(length(unique(depths)), 6L)   # the probe is well-posed

  M <- 3L
  imps <- lapply(seq_len(M), function(i) {
    d <- data.frame(x1 = rnorm(n)); d$y <- d$x1 + rnorm(n); d
  })
  idx <- as.integer(BACE:::.select_tree_indices(6, M))   # c(1, 4, 6)

  # a bace_final-shaped object carrying the pairing BACE actually used
  obj <- structure(
    list(all_datasets = imps,
         tree_index   = idx,
         tree_hashes  = vapply(trees, BACE:::.tree_cache_key, character(1))),
    class = c("bace_final", "list"))

  seen <- new.env(); seen$d <- numeric(0)
  invisible(with_imputations(obj, function(d, tree) {
    seen$d <- c(seen$d, depth(tree)); stats::lm(y ~ x1, data = d)
  }, tree = structure(trees, class = "multiPhylo"), .progress = FALSE))

  # each fit gets the tree that generated ITS imputation
  expect_equal(seen$d, depths[idx])
  # and specifically NOT the old modular answer
  expect_false(isTRUE(all.equal(seen$d, depths[((seq_len(M) - 1L) %% 6L) + 1L])))
})

test_that("#18b a different tree set is caught by content, not by count", {
  set.seed(12)
  n <- 8
  mk <- function(mult) {
    tr <- ape::rcoal(n); tr$tip.label <- paste0("t", seq_len(n))
    tr$edge.length <- tr$edge.length * mult; tr
  }
  gen_trees  <- lapply(1:4, mk)
  other_trees <- lapply(5:8, mk)          # SAME COUNT, different trees
  imps <- lapply(1:3, function(i) {
    d <- data.frame(x1 = rnorm(n)); d$y <- d$x1 + rnorm(n); d
  })
  obj <- structure(
    list(all_datasets = imps,
         tree_index   = as.integer(BACE:::.select_tree_indices(4, 3)),
         tree_hashes  = vapply(gen_trees, BACE:::.tree_cache_key, character(1))),
    class = c("bace_final", "list"))

  # counts match (4 vs 4), so the pre-1.0.1 guard stayed silent here
  expect_message(
    with_imputations(obj, function(d, tree) stats::lm(y ~ x1, data = d),
                     tree = structure(other_trees, class = "multiPhylo"),
                     .progress = FALSE),
    regexp = "not the 4 tree")

  # the matching set is silent
  expect_silent(
    with_imputations(obj, function(d, tree) stats::lm(y ~ x1, data = d),
                     tree = structure(gen_trees, class = "multiPhylo"),
                     .progress = FALSE))
})

# ---- gaussian back-transform: post_sd must be a STANDARD DEVIATION ---------

test_that("post_sd is scaled, not shifted, when back-transforming to data scale", {
  testthat::skip_on_cran()
  # sd(a + b*X) = |b| * sd(X): the location shift must NOT reach post_sd.
  # Transforming the whole .pred_cont() frame at once added mean(y) to it, so
  # on a variable with mean 1000 and sd 0.01 the reported "SD" was ~1000.003
  # instead of ~0.003, and it went NEGATIVE whenever mean(y) < 0.
  set.seed(101); n <- 25
  tr <- ape::rcoal(n); tr$tip.label <- paste0("sp", seq_len(n))

  run_one <- function(mean_val, sd_val) {
    d <- data.frame(Species = tr$tip.label, x1 = rnorm(n))
    d$y <- mean_val + sd_val * (d$x1 + rnorm(n))
    d$y[1:5] <- NA
    fit <- suppressWarnings(suppressMessages(
      bace_imp("y ~ x1", "~ 1 | Species", phylo = tr, data = d,
               runs = 1, nitt = 1200, burnin = 300, thin = 5, verbose = FALSE)))
    list(fit = fit, d = d, sd_val = stats::sd(d$y, na.rm = TRUE))
  }

  # large positive mean, tiny spread: the bug is enormous here
  r <- run_one(1000, 0.01)
  ps <- r$fit$pred_list_last_run$y$post_sd
  expect_true(all(is.finite(ps)))
  expect_true(all(ps > 0))                      # an SD is non-negative
  expect_lt(max(ps), 10 * r$sd_val)             # ~0.003, NOT ~1000
  expect_lt(max(ps), 1)                         # decisive: bug gives ~1000.003

  # exact invariant: post_sd == z-scale post_sd * sd(y)
  z_sd <- .pred_cont(r$fit$models_last_run$y)$post_sd
  expect_equal(ps, z_sd * r$sd_val)

  # negative mean: the bug produced negative standard deviations
  r2 <- run_one(-5, 0.01)
  ps2 <- r2$fit$pred_list_last_run$y$post_sd
  expect_true(all(ps2 > 0))
  expect_equal(ps2, .pred_cont(r2$fit$models_last_run$y)$post_sd * r2$sd_val)

  # the location columns ARE shifted (we did not over-correct)
  pm <- r$fit$pred_list_last_run$y$post_mean
  expect_gt(min(pm), 900)
  expect_lt(max(pm), 1100)
})

# ---- multinomial H2: the residual trace must be tr(R), not column 1 -------

test_that("multinomial diagonal columns are matched on REAL MCMCglmm names", {
  # MCMCglmm labels traits by factor LEVEL ("traity.B"), never by number.
  # The old regex demanded a digit index, so it never matched a real fit and
  # a silent fallback returned column 1 of R -- one diagonal plus J-1
  # OFF-diagonals -- making tr(R) exactly 1 for every K.
  real3 <- c("traity.B:traity.B.units", "traity.C:traity.B.units",
             "traity.B:traity.C.units", "traity.C:traity.C.units")
  expect_equal(BACE:::.get_multinom_diag_cols(real3, "units", 2L),
               c("traity.B:traity.B.units", "traity.C:traity.C.units"))

  # the fabricated digit form the old test used must still work
  fake <- c("traity.1:traity.1.units", "traity.2:traity.1.units",
            "traity.1:traity.2.units", "traity.2:traity.2.units")
  expect_equal(BACE:::.get_multinom_diag_cols(fake, "units", 2L),
               c("traity.1:traity.1.units", "traity.2:traity.2.units"))

  # factor levels containing dots must not break the split
  dotted <- c("traity.forest.edge:traity.forest.edge.units",
              "traity.open:traity.forest.edge.units",
              "traity.forest.edge:traity.open.units",
              "traity.open:traity.open.units")
  expect_equal(BACE:::.get_multinom_diag_cols(dotted, "units", 2L),
               c("traity.forest.edge:traity.forest.edge.units",
                 "traity.open:traity.open.units"))

  # idh(trait):Species has no ":" and is all diagonal
  idh <- c("traity.B.Species", "traity.C.Species")
  expect_equal(BACE:::.get_multinom_diag_cols(idh, "Species", 2L), idh)

  # ".Species" must not swallow ".Species2"
  mixed <- c("traity.B.Species", "traity.C.Species",
             "traity.B.Species2", "traity.C.Species2")
  expect_equal(BACE:::.get_multinom_diag_cols(mixed, "Species", 2L),
               c("traity.B.Species", "traity.C.Species"))
  expect_equal(BACE:::.get_multinom_diag_cols(mixed, "Species2", 2L),
               c("traity.B.Species2", "traity.C.Species2"))
})

test_that("the picked residual columns give tr(R) = 2(K-1)/K, not 1", {
  # R is FIXED at (I_J + 1 1') / K for the multinomial, so this is exact.
  for (K in 3:6) {
    J  <- K - 1L
    R  <- (diag(J) + matrix(1, J, J)) / K
    nm <- paste0("traity.", LETTERS[seq_len(J)])
    cn <- as.vector(outer(nm, nm, function(a, b) paste0(a, ":", b, ".units")))
    V  <- matrix(as.vector(R), nrow = 1, dimnames = list(NULL, cn))

    diag_cols <- BACE:::.get_multinom_diag_cols(cn, "units", J)
    expect_equal(length(diag_cols), J)
    expect_equal(unname(rowSums(V[, diag_cols, drop = FALSE])),
                 2 * J / K)                       # the truth
    expect_equal(unname(rowSums(V[, diag_cols, drop = FALSE])), sum(diag(R)))

    # what the broken fallback produced: column 1 of R, always exactly 1
    first_J <- cn[seq_len(J)]
    expect_equal(unname(rowSums(V[, first_J, drop = FALSE])), 1)
    if (K > 3L) expect_false(isTRUE(all.equal(2 * J / K, 1)))

    # and the no-columns fallback must be 2J/K, not J^2/K
    expect_false(isTRUE(all.equal((J * J) / K, (2 * J) / K)) && K != 3L)
  }
})

# ---- pool_mi: B = 0 and W = 0 are different limits ------------------------

test_that("pool_mi gives df = Inf and fmi = 0 when between-imputation var is 0", {
  # r -> 0 sends Rubin's df to INFINITY, not to M-1. Sharing the M-1 fallback
  # with the r -> Inf case made BACE report 2/(M+2) of the information missing
  # when in truth none was. mice reaches the same limit (df ~ 9e8, fmi ~ 2e-9)
  # by clamping lambda at 1e-4.
  set.seed(5); n <- 40
  d <- data.frame(x = rnorm(n)); d$y <- 1 + 0.5 * d$x + rnorm(n)
  M <- 10L
  fits <- structure(lapply(seq_len(M), function(i) stats::lm(y ~ x, data = d)),
                    class = c("bace_mi_fits", "list"),
                    n_fits = M, n_failed = 0L, failed = integer(0),
                    is_mcmcglmm = FALSE)
  p <- pool_mi(fits)
  expect_true(all(is.infinite(p$df)))
  expect_equal(p$fmi, rep(0, nrow(p)))
  expect_equal(p$riv, rep(0, nrow(p)))
  expect_false(isTRUE(all.equal(p$df[1], M - 1)))   # the old wrong answer

  # pooled estimate and SE must equal the single fit exactly
  one <- stats::lm(y ~ x, data = d)
  expect_equal(p$estimate, unname(stats::coef(one)))
  expect_equal(p$std.error, unname(sqrt(diag(stats::vcov(one)))))
})

test_that("pool_mi still gives df = M-1 when WITHIN-imputation variance is 0", {
  # The other degenerate case: W = 0 means r = Inf, and (1 + 1/r)^2 -> 1,
  # so df -> M-1. This must NOT have changed.
  M <- 6L
  mk <- function(i) {
    f <- structure(list(i = i), class = "bace_zerovar_test")
    f
  }
  registerS3method("coef", "bace_zerovar_test",
                   function(object, ...) c(a = object$i * 1.0))
  registerS3method("vcov", "bace_zerovar_test",
                   function(object, ...) matrix(0, 1, 1))
  fits <- structure(lapply(seq_len(M), mk),
                    class = c("bace_mi_fits", "list"),
                    n_fits = M, n_failed = 0L, failed = integer(0),
                    is_mcmcglmm = FALSE)
  p <- suppressWarnings(pool_mi(fits))
  expect_equal(p$df[1], M - 1)
})

# ---- a character VECTOR of formulas is the list the user meant -----------

test_that("a character vector of formulas is not silently reduced to its first", {
  expect_equal(BACE:::.normalise_fixformula(c("y ~ x1", "z ~ x1")),
               list("y ~ x1", "z ~ x1"))
  # single string, formula object and list must all pass through untouched
  expect_identical(BACE:::.normalise_fixformula("y ~ x1"), "y ~ x1")
  expect_true(inherits(BACE:::.normalise_fixformula(y ~ x1), "formula"))
  expect_equal(BACE:::.normalise_fixformula(list(y ~ x1)), list(y ~ x1))

  # stats::as.formula() on a length-2 character vector keeps only the FIRST:
  # that is the silent reduction this guards against.
  expect_warning(f <- stats::as.formula(c("y ~ x1", "z ~ x1")))
  expect_equal(all.vars(f[[2]]), "y")       # z is gone
})

test_that("every response in a character vector of formulas gets a model", {
  testthat::skip_on_cran()
  set.seed(3); n <- 20
  tr <- ape::rcoal(n); tr$tip.label <- paste0("s", seq_len(n))
  d <- data.frame(Species = tr$tip.label, x1 = rnorm(n))
  d$y <- d$x1 + rnorm(n); d$z <- d$x1 + rnorm(n)
  d$y[1:3] <- NA; d$z[4:6] <- NA

  fit <- suppressWarnings(suppressMessages(
    bace_imp(c("y ~ x1", "z ~ x1"), "~ 1 | Species", phylo = tr, data = d,
             runs = 1, nitt = 400, burnin = 100, thin = 2, verbose = FALSE)))
  # before the fix only "y" was modelled; "z" vanished with no error
  expect_true(all(c("y", "z") %in% names(fit$models_last_run)))
  expect_false(anyNA(fit$imputed_data$z))
  expect_false(anyNA(fit$imputed_data$y))
})

# ---- an unobserved declared level must not shift the category labels ------

test_that("imputed categories follow the OBSERVED levels, not the first n declared", {
  # When a declared level is never observed in training, MCMCglmm re-indexes
  # the response to 1:K_eff in declaration order, so the prediction has fewer
  # columns than the factor has levels. The right vocabulary is the OBSERVED
  # levels in declaration order. Falling back to "the first n_cols DECLARED"
  # is correct only when the unobserved levels sit at the END -- otherwise
  # every imputation is shifted onto a neighbouring category.
  decl <- c("2", "3", "4")
  pp   <- matrix(c(0.2, 0.8,
                   0.7, 0.3), nrow = 2, byrow = TRUE)

  for (obs in list(c("2", "3"), c("3", "4"), c("2", "4"))) {
    got  <- suppressWarnings(
      BACE:::.impute_levels(pp, decl, sample = FALSE, levels_observed = obs,
                            response_var = "ord"))
    want <- obs[apply(pp, 1, which.max)]
    expect_equal(got, want)
  }

  # the specific regression: observed 3|4 used to come back as 3|2
  shifted <- suppressWarnings(BACE:::.impute_levels(pp, decl, sample = FALSE))
  expect_equal(shifted, c("3", "2"))                       # the old answer
  correct <- suppressWarnings(
    BACE:::.impute_levels(pp, decl, sample = FALSE,
                          levels_observed = c("3", "4")))
  expect_equal(correct, c("4", "3"))                       # the right one
  expect_false(identical(shifted, correct))
})

test_that(".align_levels_to_probs honours the observed set when given one", {
  decl <- c("2", "3", "4")
  expect_equal(BACE:::.align_levels_to_probs(decl, 2, levels_observed = c("3","4")),
               c("3", "4"))
  expect_equal(BACE:::.align_levels_to_probs(decl, 3), decl)
  # with no observed set it can only guess, and must warn while doing so
  expect_warning(g <- BACE:::.align_levels_to_probs(decl, 2),
                 regexp = "Falling back")
  expect_equal(g, c("2", "3"))
  # the warning names the variable rather than "(unknown)"
  expect_warning(BACE:::.align_levels_to_probs(decl, 2, response_var = "ord"),
                 regexp = "'ord'")
})

# ---- heritability: the residual must be READ, not assumed ------------------

test_that("H2 uses the residual the model actually fixed, not a literal 1", {
  # The probit latent scale has no size of its own, so the residual is FIXED to
  # identify it. Everything on that scale -- V_A included -- scales with the
  # pin. Writing 1 into the denominator by hand left the numerator scaled and
  # the denominator not, so the reported H2 climbed with an arbitrary constant.
  mk <- function(V_A, V_R, n = 400) {
    m <- cbind(Species = rep(V_A, n), units = rep(V_R, n))
    structure(list(VCV = m), class = "bace_fakefit_h2")
  }
  for (cc in c(1, 2, 4)) {
    out <- BACE:::.compute_H2_stats(mk(2 * cc, cc), "threshold",
                                    species_col = "Species", species = FALSE)
    # the ratio is scale free: 2c / (2c + c) = 2/3 whatever c is
    expect_equal(unname(out$mean), 2 / 3, tolerance = 1e-8)
  }

  # the old hard-coded denominator would have given a different answer for
  # every c; confirm those answers really do differ, so this test can fail
  old <- vapply(c(1, 2, 4), function(cc) (2 * cc) / (2 * cc + 1), numeric(1))
  expect_equal(length(unique(round(old, 6))), 3L)
})

test_that("a FREE residual prior is refused for threshold-scale traits", {
  testthat::skip_on_cran()
  skip_if_not_installed("ape")
  set.seed(11); n <- 60
  tr <- ape::rcoal(n); tr$tip.label <- sprintf("s%03d", seq_len(n))
  C <- ape::vcv(tr); C <- C / max(C); L <- t(chol(C))
  lia <- as.numeric(L %*% rnorm(n)) + rnorm(n)
  d <- data.frame(Species = tr$tip.label, y = factor(as.integer(lia > 0)))
  G <- list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000))

  grab <- function(pr) {
    w <- NULL
    o <- withCallingHandlers(
      suppressMessages(phylo_signal_summary(
        d, tr, species_col = "Species", variables = "y", methods = "H2",
        prior = pr, nitt = 3000, burnin = 600, thin = 5)),
      warning = function(z) { w <<- c(w, conditionMessage(z))
                              invokeRestart("muffleWarning") })
    list(h2 = o$table$H2_mean[1], refused = any(grepl("must keep the residual FIXED", w)))
  }

  free <- grab(list(R = list(V = 1, nu = 0.002), G = G))
  expect_true(free$refused)
  expect_true(is.na(free$h2))        # NA, never the old confident 1.0000

  fixed <- grab(list(R = list(V = 1, fix = 1), G = G))
  expect_false(fixed$refused)
  expect_true(is.finite(fixed$h2))
})

test_that("phylo_signal_summary and bace() classify every column the same way", {
  bt <- function(x) { d <- data.frame(v = x)
                      BACE:::.get_type(BACE:::.summarise_var_types(d), "v", d) }
  st <- function(x) { d <- data.frame(v = x)
                      BACE:::.detect_types_for_signal(d, "v")[["v"]] }
  cases <- list(c(1L, 2L, 3L, 7L, 2L),          # integer
                c(1L, 2L, NA, 7L, 2L),          # integer with NA
                c(1, 2, 3, 7, 2),               # DOUBLE whole numbers
                c(1, 2, NA, 7, 2),              # ...with one NA
                c(1.5, 2.25, 3.1, 7.7, 2.2))    # fractional
  for (x in cases) expect_identical(bt(x), st(x))

  # the two that used to diverge, pinned explicitly
  expect_identical(st(c(1, 2, 3, 7, 2)), "gaussian")   # was "poisson"
  # and one NA must not change the family
  expect_identical(st(c(1, 2, 3, 7, 2)), st(c(1, 2, NA, 7, 2)))
  expect_identical(st(c(1L, 2L, 3L)), st(c(1L, NA, 3L)))

  # factors still behave
  expect_identical(bt(factor(c("a","b","a"))), st(factor(c("a","b","a"))))
  expect_identical(bt(factor(c("a","b","c"))), st(factor(c("a","b","c"))))
})

# ---- the z-scoring constants must be kept, per run AND per model ----------

test_that("scaling constants are recorded and make the back-transform exact", {
  testthat::skip_on_cran()
  skip_if_not_installed("ape")
  set.seed(5); n <- 50
  tr <- ape::rcoal(n); tr$tip.label <- paste0("sp", seq_len(n))
  d <- data.frame(Species = tr$tip.label, x1 = rnorm(n, 500, 20))
  d$y <- 100 + 2 * d$x1 + rnorm(n, 0, 30)
  d$y[sample(n, 12)] <- NA; d$x1[sample(n, 10)] <- NA

  fit <- suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", phylo = tr, data = d, runs = 3,
         nitt = 6000, burnin = 1000, thin = 5, n_final = 5,
         verbose = FALSE, skip_conv = TRUE, seed = 5)))

  # they exist, one entry per final run
  expect_false(is.null(fit$scaling))
  expect_equal(length(fit$scaling), 5L)
  expect_false(is.null(fit$final_results$scaling))

  own  <- get_scaling(fit, "x1")              # x1 in its own model
  pred <- get_scaling(fit, "x1", model = "y") # x1 as a predictor in y's model
  expect_equal(nrow(own), 5L)
  expect_equal(nrow(pred), 5L)
  expect_true(all(is.finite(own$sd)) && all(own$sd > 0))
  expect_true(all(is.finite(pred$sd)) && all(pred$sd > 0))

  # THE invariant: a variable's divisor as a PREDICTOR is not its observed-only
  # one, because as a predictor it carries that run's imputed values. Storing
  # only the response's constants (as the first cut of this fix did) would make
  # every cross-variable coefficient unconvertible.
  expect_equal(length(unique(round(own$sd, 10))), 1L)   # response: stable
  expect_gt(length(unique(round(pred$sd, 10))), 1L)     # predictor: moves
  expect_gt(max(abs(pred$sd - own$sd[1])) / own$sd[1], 0.01)

  # and the recorded constants recover the true slope
  sy <- get_scaling(fit, "y", model = "y")
  bz <- mean(as.matrix(get_pooled_model(fit, "y")$Sol)[, "x1"])
  expect_lt(abs(bz * mean(sy$sd) / mean(pred$sd) - 2) / 2, 0.05)

  # non-gaussian / unknown names are refused, not silently NA
  expect_error(get_scaling(fit, "not_a_variable"), regexp = "No model recorded")
  expect_error(get_scaling(fit, "y", model = "nope"), regexp = "No model recorded")
  expect_error(get_scaling(structure(list(), class = "lm")),
               regexp = "must be a 'bace_complete'")
})

test_that("a checkpoint without scaling constants is not accepted", {
  # Accepting a pre-1.0.2 checkpoint would leave that run's constants NULL and
  # its coefficients unconvertible -- silently. It must be recomputed instead.
  old <- list(models = list(), dataset = data.frame(), prob_preds = list(),
              ovr_models = list())
  new <- c(old, list(scaling = list()))
  expect_false(BACE:::.ckpt_is_run(old))
  expect_true(BACE:::.ckpt_is_run(new))
})

# ---- pool_posteriors: every run must be in the SAME units before stacking --

test_that("pool_posteriors rescales runs to a common divisor before stacking", {
  testthat::skip_on_cran()
  skip_if_not_installed("ape")
  set.seed(5); n <- 60
  tr <- ape::rcoal(n); tr$tip.label <- paste0("sp", seq_len(n))
  d <- data.frame(Species = tr$tip.label, x1 = rnorm(n, 500, 20))
  d$y <- 100 + 2 * d$x1 + rnorm(n, 0, 40)
  d$y[sample(n, 30)] <- NA; d$x1[sample(n, 30)] <- NA   # heavy, so divisors drift

  fit <- suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", phylo = tr, data = d, runs = 3,
         nitt = 8000, burnin = 1500, thin = 5, n_final = 15,
         verbose = FALSE, skip_conv = TRUE, seed = 5)))

  sy   <- get_scaling(fit, "y",  model = "y")$sd
  sx   <- get_scaling(fit, "x1", model = "y")$sd
  mods <- fit$final_results$all_models

  # the probe only means something if the divisors actually moved
  expect_gt((max(sx) - min(sx)) / mean(sx), 0.10)

  # reference: convert each run with ITS OWN constants, then stack
  ref <- unlist(lapply(seq_along(mods), function(m)
    as.matrix(mods[[m]]$y$Sol)[, "x1"] * sy[m] / sx[m]))

  # what pool_posteriors gives, put on the data scale with one constant
  pp  <- suppressMessages(pool_posteriors(fit$final_results, variable = "y"))
  got <- as.matrix(pp$models$y$Sol)[, "x1"] * mean(sy) / mean(sx)

  w <- function(x) diff(unname(stats::quantile(x, c(0.025, 0.975))))
  expect_equal(mean(got), mean(ref), tolerance = 1e-6)
  expect_equal(w(got),    w(ref),    tolerance = 1e-6)

  # and without the correction it would NOT match -- so this test can fail
  raw <- unlist(lapply(mods, function(m) as.matrix(m$y$Sol)[, "x1"])) *
           mean(sy) / mean(sx)
  expect_gt(abs(w(raw) - w(ref)) / w(ref), 0.01)
})

test_that("pool_posteriors still works when no scaling was recorded", {
  # A bace_final from an older version has no $scaling. Pooling must still
  # run -- uncorrected, as it always did -- rather than erroring.
  testthat::skip_on_cran()
  skip_if_not_installed("ape")
  set.seed(9); n <- 30
  tr <- ape::rcoal(n); tr$tip.label <- paste0("sp", seq_len(n))
  d <- data.frame(Species = tr$tip.label, x1 = rnorm(n))
  d$y <- 0.5 * d$x1 + rnorm(n); d$y[1:8] <- NA
  fit <- suppressWarnings(suppressMessages(
    bace("y ~ x1", "~ 1 | Species", phylo = tr, data = d, runs = 2,
         nitt = 2500, burnin = 500, thin = 5, n_final = 4,
         verbose = FALSE, skip_conv = TRUE, seed = 9)))

  legacy <- fit$final_results
  legacy$scaling <- NULL
  expect_s3_class(suppressMessages(pool_posteriors(legacy, variable = "y")),
                  "bace_pooled")
})

# ---- the step plot must USE the check's yardstick, not invent one ---------

test_that("the check records the yardstick it reached its verdict with", {
  a <- BACE:::.conv_assess_series(rnorm(10, 5, 0.02), rep(1, 10),
                                  n_cells = 20L, kind = "numeric")
  expect_equal(unname(a$stats["spread"]), 1)        # the imputation spread

  b <- BACE:::.conv_assess_series(rep(c(0.50, 0.525), 5), rep(0, 10),
                                  n_cells = 40L, kind = "proportion")
  # a proportion is measured in PROPORTION units; sqrt(p(1-p)) is a function
  # of the mean and must never be used as a spread
  expect_equal(unname(b$stats["spread"]), 1)
  expect_false(isTRUE(all.equal(unname(b$stats["spread"]),
                                mean(sqrt(0.5 * 0.5)))))

  # an unassessable series still carries the field, as NA
  u <- BACE:::.conv_assess_series(c(1, 2), c(1, 1), n_cells = 5L)
  expect_true("spread" %in% names(u$stats))
  expect_true(is.na(unname(u$stats["spread"])))
})

test_that("the plotted step equals the step the check reports", {
  # The panel header prints the check's mean step. The bars used to be drawn
  # on a different yardstick entirely -- sqrt(p(1-p)) -- so the picture and
  # its own caption disagreed, by 2x on an even binary and 18x on a
  # near-certain one, enough to put a settled variable above the alarm line.
  probe <- function(m, n_cells) {
    a <- BACE:::.conv_assess_series(m, sqrt(m * (1 - m)), n_cells = n_cells,
                                    kind = "proportion")
    list(check = unname(a$stats["step"]),
         new   = mean(abs(diff(m))) / unname(a$stats["spread"]),
         old   = mean(abs(diff(m))) / mean(sqrt(m * (1 - m))))
  }

  even <- probe(rep(c(0.50, 0.525), 5), 40L)
  expect_equal(even$new, even$check, tolerance = 1e-10)
  expect_gt(even$old / even$check, 1.5)             # the old panel overstated

  p <- rep(0.999, 10); p[seq(2, 10, 2)] <- 0.994
  near <- probe(p, 40L)
  expect_equal(near$new, near$check, tolerance = 1e-10)
  expect_gt(near$old / near$check, 10)              # ~18x on a near-certain binary
})
