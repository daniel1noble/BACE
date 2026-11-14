# Tests for the ordinal trait generator, and a golden test for every trait type.
#
# THE BUG THESE EXIST FOR (2026-11-14). ordinal_liab2cat() cut the latent
# liability at FIXED, evenly spaced constants (-2.25 / 0 / +2.25 for four
# categories) whatever the liability's scale. In the 2026-11 campaign that
# liability is 0.5 * x1 plus the ordinal's own species, phylogenetic and
# residual effects, with population SD 1.00 / 1.53 / 2.65 at signal 0.3 / 0.7 /
# 0.9 (0.53 / 0.82 / 1.41 for an ordinal with no parent), so at signal 0.3 about
# two species in a hundred reached the outer cutpoints, and with no parent none.
#
# Measured on the 1.0.3 generator (n = 100, 40 seeds per setting, traits related
# at slope 0.5, whole clades hidden):
#
#   signal 0.3   fewer than 4 levels in 27/40; a level unobserved after masking in 31/40
#   signal 0.7                        7/40                                       9/40
#   signal 0.9                        2/40                                       4/40
#
# Seed 1 at signal 0.3 came out 0/47/52/1 (observed 0/30/37/1). One dataset at
# 0.7 (truth 1/43/55/1, observed 0/29/39/0) ran for an hour where its neighbours
# took nineteen minutes: MCMCglmm's cutpoint sampler has nothing to anchor a
# cutpoint between two empty categories.
#
# THE FIX, AND WHY IT IS NOT SAMPLE QUANTILES. sim_bace() now cuts at POPULATION
# thresholds: by default the k/K quantiles of the variable's own latent
# distribution, or the cutpoints given in `ordinal_thresholds` (the campaign
# passes the exact quartiles of its full latent scale). Cutting at the sample
# quantiles would populate every level by construction, but it makes the hidden
# cells' level counts the exact complement of the observed cells' (correlation
# -1.00 across seeds in every setting tried), so an imputer that estimates its
# cutpoints from observed cells -- BACE's threshold model does -- is pushed the
# wrong way, and the observed-majority floor scores 0.11-0.21 against a chance
# level of 0.25. Under population thresholds the correlation is +0.2 to +0.9 and
# the floor 0.27-0.43. The price: at signal 0.9 and n = 100 about 1 dataset in
# 20 has an empty level in truth (a whole tree shifting together); the campaign
# records level counts per dataset. No rule consumes random numbers, so trees,
# masks and every other trait are unchanged by position (verified on 200
# datasets); the quantile rule remains available for direct calls.

# ---- 1. every category is populated, whatever the latent scale or mean -------

test_that("quantile cutpoints populate every category at any latent scale", {
  set.seed(14)
  for (sd in c(0.3, 1.2, 3.2, 10)) {
    for (mu in c(0, 3, -5)) {
      liab <- rnorm(100, mu, sd)
      out  <- ordinal_liab2cat(liab, n_cats = 4)
      expect_identical(as.vector(table(factor(out, levels = 1:4))), c(25L, 25L, 25L, 25L),
                       info = sprintf("sd %.1f, mean %.0f", sd, mu))
    }
  }
  # n not divisible by K: shares differ by at most one case
  liab <- rnorm(100, 0, 1.2)
  for (k in c(3, 5, 6, 7)) {
    tab <- table(factor(ordinal_liab2cat(liab, n_cats = k), levels = seq_len(k)))
    expect_true(all(tab >= floor(100 / k)) && all(tab <= ceiling(100 / k)),
                info = paste("K =", k, ":", paste(tab, collapse = "/")))
  }
  # and at the campaign's larger tree size
  liab <- rnorm(500, 0, 3.2)
  expect_identical(as.vector(table(ordinal_liab2cat(liab, n_cats = 4))), rep(125L, 4))
})

test_that("the fixed cutpoints reproduce the collapse (why the default changed)", {
  # The campaign's ordinal at signal 0.3 has a latent SD of 1.0 (typical
  # regime). With cutpoints at +/-2.25 the two extreme categories share the
  # ~2.4 per cent of species beyond them.
  set.seed(3)
  liab <- rnorm(100, 0, 1.0)
  old  <- suppressWarnings(ordinal_liab2cat(liab, n_cats = 4, method = "fixed"))
  expect_lt(sum(old == 1L) + sum(old == 4L), 15L)
  new  <- ordinal_liab2cat(liab, n_cats = 4)
  expect_identical(sum(new == 1L) + sum(new == 4L), 50L)

  # Tighter still (SD 0.5) and the extremes are EMPTY: the tripwire must say so.
  liab <- rnorm(100, 0, 0.5)
  expect_warning(ordinal_liab2cat(liab, n_cats = 4, method = "fixed"),
                 "2 of the 4 requested categories")
})

# ---- 2. the contract of the new cutter ---------------------------------------

test_that("probs sets unequal shares exactly", {
  set.seed(21)
  liab <- rnorm(100)
  out  <- ordinal_liab2cat(liab, n_cats = 4, probs = c(0.1, 0.2, 0.3, 0.4))
  expect_identical(as.vector(table(factor(out, levels = 1:4))), c(10L, 20L, 30L, 40L))
})

test_that("method = 'thresholds' cuts exactly where it is told", {
  liab <- c(-3, -1.1, -0.9, 0.1, 0.9, 1.1, 3)
  expect_identical(ordinal_liab2cat(liab, n_cats = 4, method = "thresholds", thresholds = c(-1, 0, 1)),
                   c(1L, 1L, 2L, 3L, 3L, 4L, 4L))
  # a value exactly on a threshold stays in the lower category (strict >)
  expect_identical(ordinal_liab2cat(c(-1, 0, 1, 2), n_cats = 4, method = "thresholds", thresholds = c(-1, 0, 1)),
                   c(1L, 2L, 3L, 4L))
  # validation
  expect_error(ordinal_liab2cat(liab, n_cats = 4, method = "thresholds"), "non-decreasing")
  expect_error(ordinal_liab2cat(liab, n_cats = 4, method = "thresholds", thresholds = c(-1, 0)), "3 finite")
  expect_error(ordinal_liab2cat(liab, n_cats = 4, method = "thresholds", thresholds = c(1, 0, -1)), "non-decreasing")
  expect_error(ordinal_liab2cat(liab, n_cats = 4, method = "thresholds", thresholds = c(-1, NA, 1)), "finite")
  # a threshold nobody reaches leaves a category empty: the tripwire names the cause
  expect_warning(ordinal_liab2cat(liab, n_cats = 4, method = "thresholds", thresholds = c(-1, 0, 10)),
                 "whole tree")
})

test_that("the cut is monotone in the liability and returns integer codes", {
  set.seed(22)
  liab <- rnorm(250, 2, 4)
  out  <- ordinal_liab2cat(liab, n_cats = 5)
  expect_true(is.integer(out))
  expect_true(all(diff(out[order(liab)]) >= 0L))
  # the highest liabilities are in the top category, the lowest in the bottom
  expect_identical(out[which.max(liab)], 5L)
  expect_identical(out[which.min(liab)], 1L)
})

test_that("cutting consumes no random numbers (so nothing downstream shifts)", {
  # The campaign's masks, trees and other traits were unchanged by the fix
  # precisely because the cutter never touches the RNG. Keep it that way.
  liab <- seq(-3, 3, length.out = 60)
  set.seed(5); invisible(ordinal_liab2cat(liab, n_cats = 4));                   a <- runif(1)
  set.seed(5); invisible(ordinal_liab2cat(liab, n_cats = 4, method = "fixed")); b <- runif(1)
  set.seed(5); invisible(ordinal_liab2cat(liab, n_cats = 4, method = "thresholds", thresholds = c(-1, 0, 1))); d <- runif(1)
  set.seed(5);                                                                  r <- runif(1)
  expect_identical(a, r)
  expect_identical(b, r)
  expect_identical(d, r)
})

test_that("rounding-level differences in the liability cannot move a code", {
  # Machines agree on the liability to ~1e-15 (Cholesky draw, 1.0.3). The
  # cutpoints sit between order statistics, far from any value, so noise many
  # orders larger than rounding leaves every code where it was.
  set.seed(23)
  liab <- rnorm(100, 0, 1.2)
  base <- ordinal_liab2cat(liab, n_cats = 4)
  for (eps in c(1e-13, 1e-11, 1e-9)) {
    jit <- liab + eps * rnorm(100)
    expect_identical(ordinal_liab2cat(jit, n_cats = 4), base, info = paste("eps", eps))
  }
})

test_that("method = 'fixed' still gives the old evenly spaced constants", {
  # -2.25 / 0 / +2.25 for four categories at the default spread of 1.5
  liab <- c(-3, -2.26, -2.24, -0.1, 0.1, 2.24, 2.26, 3)
  expect_identical(ordinal_liab2cat(liab, n_cats = 4, method = "fixed"),
                   c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L))
  # and the default is NOT the fixed rule
  expect_false(identical(ordinal_liab2cat(c(liab, rnorm(40, 0, 0.1)), n_cats = 4),
                         suppressWarnings(ordinal_liab2cat(c(liab, rnorm(40, 0, 0.1)),
                                                           n_cats = 4, method = "fixed"))))
})

test_that("the tripwire and the argument checks fire", {
  # ties at a cutpoint are the only way the quantile rule can leave a category empty
  expect_warning(ordinal_liab2cat(rep(1, 10), n_cats = 3), "ties")
  expect_error(ordinal_liab2cat(rnorm(10), n_cats = 4, probs = c(0.5, 0.5)), "probs")
  expect_error(ordinal_liab2cat(rnorm(10), n_cats = 4, probs = c(0.1, 0.2, 0.3, 0.3)), "summing to 1")
  expect_error(ordinal_liab2cat(rnorm(3),  n_cats = 4), "at least n_cats")
  expect_error(ordinal_liab2cat(rnorm(10), n_cats = 1), "at least 2")
})

# ---- 3. end to end through sim_bace: the campaign's symptom cannot recur -----

test_that("sim_bace's default thresholds are the variable's own population quantiles", {
  # intercept + qnorm(k/K) * sqrt(sigma2_species + sigma2_phylo + sigma2_residual),
  # with the non-gaussian defaults sigma2 = 0.1 each and
  # sigma2_phylo = signal * 0.2 / (1 - signal). Checked on a response and a predictor.
  set.seed(3)
  s <- suppressWarnings(suppressMessages(
    sim_bace(response_type = "threshold3", predictor_types = c("gaussian", "threshold4"),
             var_names = c("y", "x1", "ord"), beta_matrix = matrix(0, 2, 2),
             beta_resp = list(x1 = 0, ord = 0),
             intercepts = list(predictors = c(0, 1.5), response = -0.5),
             phylo_signal = c(0.5, 0.2, 0.8), n_cases = 300, n_species = 100,
             missingness = NULL)))
  own_sd <- function(ps) sqrt(0.1 + 0.1 + 0.2 * ps / (1 - ps))
  expect_equal(s$params$ordinal_thresholds$y,   -0.5 + qnorm(c(1/3, 2/3)) * own_sd(0.5), tolerance = 1e-12)
  expect_equal(s$params$ordinal_thresholds$ord,  1.5 + qnorm(c(.25, .5, .75)) * own_sd(0.8), tolerance = 1e-12)
  # and with no parent the shares are roughly equal (300 cases on 100 species)
  tab <- table(factor(s$complete_data$ord, levels = 1:4))
  expect_true(all(tab >= 30))
})

test_that("explicit ordinal_thresholds reach the variable, and are checked", {
  # The ord codes must equal the cut of the liability rebuilt from the pieces
  # sim_bace returns: 0.5 * x1 + species + phylogenetic + residual effects.
  thr <- c(-1, 0.2, 1.1)
  set.seed(4)
  s <- suppressWarnings(suppressMessages(
    sim_bace(response_type = "gaussian", predictor_types = c("gaussian", "threshold4"),
             var_names = c("y", "x1", "ord"),
             beta_matrix = matrix(c(0, 0, 0.5, 0), 2, 2, byrow = TRUE),
             beta_resp = list(x1 = 0.3, ord = 0.2), phylo_signal = c(0.5, 0.5, 0.5),
             n_cases = 120, n_species = 120, missingness = NULL,
             ordinal_thresholds = list(ord = thr))))
  idx  <- match(s$complete_data$Species, s$tree$tip.label)
  liab <- 0.5 * s$complete_data$x1 + s$random_effects$u_species$ord[idx] +
    s$random_effects$u_phylo$ord[idx] + s$random_effects$residuals$ord
  cut_at <- function(l, t) { r <- rep(1L, length(l)); for (k in seq_along(t)) r[l > t[k]] <- k + 1L; r }
  expect_identical(as.integer(as.character(s$complete_data$ord)), cut_at(liab, thr))
  expect_identical(s$params$ordinal_thresholds$ord, thr)

  bad <- function(...) suppressWarnings(suppressMessages(
    sim_bace(response_type = "gaussian", predictor_types = c("gaussian", "threshold4"),
             var_names = c("y", "x1", "ord"), beta_matrix = matrix(0, 2, 2),
             beta_resp = list(x1 = 0, ord = 0), n_cases = 30, n_species = 30,
             missingness = NULL, ...)))
  expect_error(bad(ordinal_thresholds = list(x1 = c(-1, 0, 1))), "not thresholdK")
  expect_error(bad(ordinal_thresholds = list(ord = c(-1, 0))), "3 finite")
  expect_error(bad(ordinal_thresholds = c(-1, 0, 1)), "named list")
})

test_that("sim_bace's ordinal predictor has every level at low signal", {
  skip_on_cran()
  # Own-scale population thresholds: every level present, none thin, over five seeds.
  for (seed in 1:5) {
    set.seed(seed)
    sim <- suppressWarnings(suppressMessages(
      sim_bace(response_type = "gaussian",
               predictor_types = c("gaussian", "threshold4"),
               var_names = c("y", "x1", "ord"),
               phylo_signal = c(0.3, 0.3, 0.3), missingness = NULL,
               n_cases = 100, n_species = 100)))
    tab <- table(factor(sim$complete_data$ord, levels = 1:4))
    expect_true(all(tab >= 5L), info = paste("seed", seed, ":", paste(tab, collapse = "/")))
  }
})

# ---- 4. THE GOLDEN TEST FOR EVERY TRAIT TYPE --------------------------------
#
# The campaign's exact recipe (dev/simulation_2026-11/R/02_simulate.R) for the
# scenario n100_ps0.7_typical, seed 1: a 3-level categorical response and six
# predictors -- continuous, continuous, count, binary, 4-level ordinal,
# continuous -- all with phylogenetic signal 0.7 on a 100-species tree, the
# ordinal cut at the quartiles of its population latent scale,
# sqrt((2 * 0.5^2 + 0.2) / (1 - 0.7)). The values below were produced on
# aarch64-apple-darwin (reference BLAS) with BACE 1.0.4 and must be reproduced
# on ANY machine: the continuous values to 1e-8, the integer and factor codes
# exactly.
#
# IF THIS FAILS ON YOUR MACHINE the same seed is not giving the same dataset,
# and simulation results cannot be pooled across machines. Do not "fix" it by
# widening the tolerance or regenerating the numbers; find the step that
# diverged (the tree, a continuous draw, a count, a code) -- that is the point
# of pinning each one separately.

test_that("the same seed gives the identical dataset for every trait type", {
  old <- RNGkind(); on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")

  recipe <- function() {
    predictors <- c("x1", "x2", "cnt", "bin", "ord", "driver")
    types      <- c("gaussian", "gaussian", "poisson", "binary", "threshold4", "gaussian")
    B <- matrix(0, 6, 6, dimnames = list(predictors, predictors))
    for (v in c("x2", "cnt", "bin", "ord")) B[v, "x1"] <- 0.5
    B["driver", "x1"] <- 0.35
    set.seed(1)
    suppressWarnings(suppressMessages(sim_bace(
      response_type = "multinomial3", predictor_types = types,
      var_names = c("cat3", predictors), beta_matrix = B,
      beta_resp = setNames(as.list(rep(0.4, 6)), predictors),
      phylo_signal = rep(0.7, 7), n_cases = 100, n_species = 100,
      missingness = NULL,
      ordinal_thresholds = list(ord = qnorm(c(0.25, 0.5, 0.75)) * sqrt(0.7 / 0.3)))))
  }
  sim <- recipe()
  d   <- sim$complete_data

  # the tree
  expect_equal(sum(sim$tree$edge.length), 230.381009662824, tolerance = 1e-10)
  # labels are sorted in C order on every machine (1.0.4), so these are the
  # same in en_US, C and any other locale
  expect_identical(head(sim$tree$tip.label, 3), c("0i6Rd", "1WsN2", "1wIAO"))
  expect_identical(head(d$Species, 3), c("yjMQ9", "SxiO7", "nSPob"))

  # continuous (Cholesky draw: unique, so machines agree to rounding)
  expect_equal(d$x1[1:3], c(2.952970359363, -3.089835456770, 2.359827961157), tolerance = 1e-8)
  expect_equal(sum(d$x1), -1.547452678902, tolerance = 1e-8)
  expect_equal(d$x2[1:3], c(1.561641072645, -3.704230300603, 4.544254216906), tolerance = 1e-8)
  expect_equal(sum(d$x2), 9.479624531244, tolerance = 1e-8)
  expect_equal(d$driver[1:3], c(0.893563614932, 0.597644010307, -0.974400388109), tolerance = 1e-8)
  expect_equal(sum(d$driver), 29.539815577529, tolerance = 1e-8)

  # count
  expect_identical(d$cnt[1:10], c(1L, 0L, 3L, 3L, 1L, 8L, 1L, 0L, 0L, 18L))
  expect_identical(sum(d$cnt), 398L)

  # binary
  expect_identical(paste(as.character(d$bin[1:20]), collapse = ""), "10100100010010101101")
  expect_identical(as.vector(table(d$bin)), c(52L, 48L))

  # ordinal: population thresholds, so the shares vary -- these exact codes
  expect_equal(sim$params$ordinal_thresholds$ord, qnorm(c(0.25, 0.5, 0.75)) * sqrt(0.7 / 0.3), tolerance = 1e-12)
  expect_identical(paste(as.character(d$ord[1:20]), collapse = ""), "41441422242243333233")
  expect_identical(as.vector(table(d$ord)), c(11L, 35L, 29L, 25L))

  # categorical response
  expect_identical(paste(as.character(d$cat3[1:20]), collapse = ""), "CABCCCCBBCCCCBCCCACB")
  expect_identical(as.vector(table(d$cat3)), c(19L, 19L, 62L))

  # and the whole thing is deterministic in-process
  expect_identical(recipe()$complete_data, d)
})

test_that("the dataset does not depend on the machine's locale", {
  # R CMD check runs tests in the C locale; a laptop usually runs en_US.UTF-8,
  # and the two collate mixed-case strings differently ("1wIAO" < "1WsN2" in
  # en_US, the reverse in C). Found 2026-11-14: the golden test above passed on
  # every value under check but failed on the species names, because the
  # generator sorted names with the session's collation. Worse, a first fix
  # that sorted only sim_tree()'s names in C order silently CHANGED THE VALUES
  # in en_US -- a second, locale-collated ordering downstream had to agree with
  # the first for each case to receive its own species' effect. Names are now
  # sorted in C order everywhere and cases are matched to tips by name, so
  # neither the labels nor the values may move between locales.
  with_collate <- function(loc, code) {
    old <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
    if (!nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", loc))))
      skip(paste("locale", loc, "not available"))
    force(code)
  }
  # A second locale that happens to collate like C would make this test pass
  # on broken code. Find one that really differs on the strings at issue.
  differs <- function(loc) isTRUE(tryCatch(
    with_collate(loc, !identical(sort(c("1wIAO", "1WsN2")), sort(c("1wIAO", "1WsN2"), method = "radix"))),
    skip = function(e) FALSE))
  cands <- c("en_US.UTF-8", "en_GB.UTF-8", "en_AU.UTF-8", "C.UTF-8")
  other <- cands[vapply(cands, differs, logical(1))][1]
  skip_if(is.na(other), "no locale available that collates differently from C")
  one <- function() {
    set.seed(77)
    suppressWarnings(suppressMessages(
      sim_bace(response_type = "gaussian",
               predictor_types = c("gaussian", "threshold4", "binary", "poisson"),
               var_names = c("y", "x1", "ord", "bin", "cnt"),
               phylo_signal = rep(0.7, 5), missingness = NULL,
               n_cases = 60, n_species = 60)))
  }
  a <- with_collate("C", one())
  b <- with_collate(other, one())
  expect_identical(a$tree$tip.label, b$tree$tip.label)
  expect_identical(a$complete_data, b$complete_data)
  # the labels really are in C order, whatever the session's locale
  expect_identical(b$tree$tip.label, sort(b$tree$tip.label, method = "radix"))
})

# ---- 5. the other pairings that used to depend on a sort ----------------

test_that("interaction term names and coefficients do not depend on the locale", {
  # parse_ix_matrix() named "Body:age" in C and "age:Body" in en_US, so a
  # user's beta_ix matched in one locale and was silently replaced by a random
  # draw in the other -- and that draw shifted the whole random stream.
  ix <- matrix(0, 3, 3); ix[3, 1:2] <- 1
  with_collate <- function(loc, code) {
    old <- Sys.getlocale("LC_COLLATE"); on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
    if (!nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", loc)))) skip(paste("locale", loc, "not available"))
    force(code)
  }
  nm <- function() names(parse_ix_matrix(ix, c("Body", "age", "y"))$y)
  expect_identical(with_collate("C", nm()), with_collate("en_US.UTF-8", nm()))
  expect_identical(with_collate("C", nm()), "Body:age")

  # either order of a supplied name reaches the term; a name that matches
  # nothing is an error, never a silent random coefficient
  one <- function(bix) {
    set.seed(5)
    suppressWarnings(suppressMessages(sim_bace(
      response_type = "gaussian", predictor_types = c("gaussian", "gaussian"),
      var_names = c("y", "Body", "age"), ix_matrix = ix, beta_ix = bix,
      n_cases = 50, n_species = 25, missingness = NULL)))
  }
  a <- one(list(y = list("Body:age" = 0.3)))
  b <- one(list(y = list("age:Body" = 0.3)))
  expect_identical(unname(a$params$ix_betas$y), 0.3)
  expect_identical(a$complete_data, b$complete_data)
  expect_identical(with_collate("C", one(list(y = list("age:Body" = 0.3))))$complete_data,
                   with_collate("en_US.UTF-8", one(list(y = list("age:Body" = 0.3))))$complete_data)
  expect_error(one(list(y = list("Body:agee" = 0.3))), "does not exist")
  expect_error(one(list(Body = list("Body:age" = 0.3))), "no interactions")
})

test_that("a one-species draw is an error, not a hang", {
  # ape::vcv() on a one-tip tree never returns. With 2 species and 3 cases drawn
  # with replacement, seed 1 puts every case on one species.
  set.seed(1)
  expect_error(sim_tree(2, 0.8, 0.4, 3), "at least 2")
})
