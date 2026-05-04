# Targeted tests for R/model_functions.R internals.
# Existing tests exercise the threshold + binary + ordinal paths and
# the OVR categorical fallback (ovr_categorical = TRUE). The
# multinomial-probit path through .pred_cat / .pred_cat_iter /
# .pred_cat_forward / .pred_cat_forward_iter is only reachable with
# ovr_categorical = FALSE on a K >= 3 categorical trait, and is the
# largest uncovered chunk of model_functions.R. These tests close
# that gap, plus add direct unit coverage of the small helpers
# .cat_mode and .impute_levels (sample = TRUE path).

# ---- Shared mixed-type fixture ---------------------------------------------

build_mixed_type_fixture <- function(seed = 2026, n = 30) {
  set.seed(seed)
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y_cont = rnorm(n, 0, 1),
    y_cat  = factor(sample(c("A", "B", "C"), n, replace = TRUE),
                    levels = c("A", "B", "C")),
    x1     = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y_cont[c(2, 5, 9)]  <- NA
  data$y_cat[c(3, 7, 12)]  <- NA
  list(phylo = phylo, data = data)
}

# ---- bace_imp with multinomial-probit categorical (.pred_cat path) ---------

test_that("bace_imp(ovr_categorical=FALSE) imputes K=3 categorical end-to-end", {
  testthat::skip_on_cran()
  fix <- build_mixed_type_fixture(seed = 2026, n = 30)
  res <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = list("y_cont ~ y_cat + x1",
                          "y_cat ~ y_cont + x1"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = 2,
    nitt = 600, thin = 5, burnin = 100,
    ovr_categorical = FALSE,
    verbose = FALSE
  )))
  expect_s3_class(res, "bace")
  # Imputed cells should be valid factor levels of y_cat
  imp <- res$data[[length(res$data)]]
  expect_true(all(!is.na(imp$y_cat)))
  expect_true(all(as.character(imp$y_cat) %in% c("A", "B", "C")))
})

# ---- bace_imp with ordinal trait (.pred_threshold_forward_iter etc.) -------

test_that("bace_imp imputes K=4 ordinal trait end-to-end", {
  testthat::skip_on_cran()
  set.seed(2026)
  n <- 30
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y_ord = factor(sample(c("1","2","3","4"), n, replace = TRUE),
                   levels = c("1","2","3","4"), ordered = TRUE),
    x1    = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y_ord[c(2, 5, 9)] <- NA

  res <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_ord ~ x1",
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo,
    data           = data,
    runs           = 2,
    nitt = 600, thin = 5, burnin = 100,
    verbose = FALSE
  )))
  expect_s3_class(res, "bace")
  imp <- res$data[[length(res$data)]]
  expect_true(all(!is.na(imp$y_ord)))
  expect_true(all(as.character(imp$y_ord) %in% c("1","2","3","4")))
})

# ---- .cat_mode: edge cases -------------------------------------------------

test_that(".cat_mode returns most-frequent value per row", {
  m <- matrix(c("A","A","B",
                "B","B","B",
                "A","B","C"), nrow = 3, byrow = TRUE)
  out <- BACE:::.cat_mode(m)
  expect_equal(out[1], "A")  # 2 As, 1 B
  expect_equal(out[2], "B")  # 3 Bs
  # Row 3 is a 3-way tie; must be one of A/B/C
  expect_true(out[3] %in% c("A","B","C"))
})

test_that(".cat_mode handles NA-only rows", {
  m <- matrix(c(NA_character_, NA_character_,
                "A", NA_character_),
              nrow = 2, byrow = TRUE)
  out <- BACE:::.cat_mode(m)
  expect_true(is.na(out[1]))
  expect_equal(out[2], "A")
})

# ---- .impute_levels: sample = TRUE path ------------------------------------

test_that(".impute_levels samples from posterior probability matrix", {
  set.seed(2026)
  pred_prob <- matrix(c(0.1, 0.9,
                         0.5, 0.5,
                         0.9, 0.1),
                       nrow = 3, byrow = TRUE,
                       dimnames = list(NULL, c("yes", "no")))
  out <- BACE:::.impute_levels(pred_prob,
                               levels_var = c("yes", "no"),
                               sample = TRUE)
  expect_length(out, 3L)
  expect_true(all(out %in% c("yes", "no")))
})

# ---- .make_prior: gelman options + each type ------------------------------

test_that(".make_prior produces a valid prior list per type (gelman = 0)", {
  for (type in c("gaussian", "poisson", "categorical", "threshold")) {
    p <- BACE:::.make_prior(n_rand = 1, n_levels = 3, type = type,
                            gelman = 0)
    expect_type(p, "list")
    expect_true("R" %in% names(p))
    expect_true("G" %in% names(p))
  }
})

test_that(".make_prior gelman = 2 for categorical needs fixform + data", {
  # gelman=2 (pseudo-Gelman) for categorical requires the formula and
  # data so it can count the J*nfixef expansion. Pass a minimal valid
  # combo and verify it returns without error.
  d <- data.frame(y = factor(c("A","B","C","A","B"),
                              levels = c("A","B","C")),
                  x = rnorm(5))
  p <- suppressWarnings(BACE:::.make_prior(
    n_rand    = 1,
    n_levels  = 3,
    type      = "categorical",
    gelman    = 2,
    fixform   = y ~ x,
    data      = d))
  expect_type(p, "list")
  expect_true("B" %in% names(p) || "R" %in% names(p))
})

test_that(".make_prior n_rand validation rejects > 2", {
  expect_error(BACE:::.make_prior(n_rand = 3, type = "gaussian",
                                  gelman = 0),
               regexp = "n_rand must be 1.*or 2")
})

# ---- .align_levels_to_probs: extra branches not covered in test-level_alignment

test_that(".align_levels_to_probs: identity when n_cols == declared", {
  out <- BACE:::.align_levels_to_probs(c("a","b","c"), n_cols = 3L)
  expect_identical(out, c("a","b","c"))
})

test_that(".align_levels_to_probs prefers observed when length matches", {
  out <- BACE:::.align_levels_to_probs(
    levels_declared = c("a","b","c","d","e"),
    n_cols = 3L,
    levels_observed = c("b","c","d"))
  expect_identical(out, c("b","c","d"))
})
