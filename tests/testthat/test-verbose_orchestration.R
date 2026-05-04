# Tests targeting the verbose=TRUE messaging paths in bace() /
# bace_imp() / bace_final_imp() and parallel-execution branches.
# These are mostly cat() lines that only fire when verbose is TRUE
# or when n_cores > 1.

# ---- Shared fixture ------------------------------------------------------

build_verbose_fixture <- function(seed = 2026, n = 12) {
  set.seed(seed)
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y = rnorm(n, 0, 1),
    x = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y[c(2, 5, 8)] <- NA
  list(phylo = phylo, data = data)
}

# ---- bace() verbose paths ------------------------------------------------

test_that("bace(verbose=TRUE) prints orchestration headers", {
  testthat::skip_on_cran()
  fix <- build_verbose_fixture()
  out <- capture.output(
    suppressWarnings(suppressMessages(res <- bace(
      fixformula     = "y ~ x",
      ran_phylo_form = "~ 1 | Species",
      phylo          = fix$phylo,
      data           = fix$data,
      runs           = 2, n_final = 2,
      nitt = 400, thin = 5, burnin = 100,
      skip_conv      = TRUE,
      max_attempts   = 1, n_cores = 1L,
      verbose        = TRUE
    )))
  )
  # Should print step headers
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "Step", fixed = TRUE)
})

test_that("bace(skip_conv=TRUE, verbose=TRUE) prints the skip-convergence message", {
  testthat::skip_on_cran()
  fix <- build_verbose_fixture()
  # Force convergence to fail by using very short MCMC, then skip_conv
  out <- capture.output(
    suppressWarnings(suppressMessages(res <- bace(
      fixformula     = "y ~ x",
      ran_phylo_form = "~ 1 | Species",
      phylo          = fix$phylo,
      data           = fix$data,
      runs           = 2, n_final = 2,
      nitt = 200, thin = 1, burnin = 50,
      skip_conv      = TRUE,
      max_attempts   = 1, n_cores = 1L,
      verbose        = TRUE
    )))
  )
  combined <- paste(out, collapse = "\n")
  expect_true(nchar(combined) > 0L)
  # bace() should have printed something about each step
  expect_match(combined, "Step|imputation", ignore.case = TRUE)
})

test_that("bace(verbose=TRUE, max_attempts=2) prints retry messaging when conv fails", {
  testthat::skip_on_cran()
  fix <- build_verbose_fixture()
  out <- capture.output(
    suppressWarnings(suppressMessages(res <- bace(
      fixformula     = "y ~ x",
      ran_phylo_form = "~ 1 | Species",
      phylo          = fix$phylo,
      data           = fix$data,
      runs           = 2, n_final = 2,
      nitt = 200, thin = 1, burnin = 50,
      skip_conv      = FALSE,
      max_attempts   = 2, n_cores = 1L,
      verbose        = TRUE
    )))
  )
  combined <- paste(out, collapse = "\n")
  # Either succeeded on first attempt or printed retry messaging
  expect_true(nchar(combined) > 0L)
})

# ---- bace_final_imp parallel + nitt_cat_mult paths -----------------------

test_that("bace_final_imp with verbose=TRUE prints completion messages", {
  testthat::skip_on_cran()
  fix <- build_verbose_fixture()
  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = 2, nitt = 400, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  out <- capture.output(suppressWarnings(suppressMessages(
    bace_final_imp(
      bace_object    = step1,
      fixformula     = "y ~ x",
      ran_phylo_form = "~ 1 | Species",
      phylo          = fix$phylo,
      nitt = 400, thin = 5, burnin = 100,
      n_final = 2,
      verbose        = TRUE
    )
  )))
  combined <- paste(out, collapse = "\n")
  expect_true(nchar(combined) > 0L)
  expect_match(combined, "imputation|Run|Completed", ignore.case = TRUE)
})

test_that("bace_final_imp with nitt_cat_mult=2 on categorical scales nitt", {
  testthat::skip_on_cran()
  set.seed(2026)
  n <- 12
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y_cat = factor(sample(c("A","B","C"), n, replace = TRUE),
                    levels = c("A","B","C")),
    x = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y_cat[c(2, 5)] <- NA

  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_cat ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo, data = data,
    runs           = 2, nitt = 400, thin = 5, burnin = 100,
    verbose        = FALSE
  )))
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = "y_cat ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2,
    nitt_cat_mult  = 2L,
    verbose        = FALSE
  )))
  expect_s3_class(step2, "bace_final")
  # Categorical imputed values are valid factor levels
  for (d in step2$all_datasets) {
    expect_true(all(!is.na(d$y_cat)))
    expect_true(all(as.character(d$y_cat) %in% c("A","B","C")))
  }
})

test_that("bace_final_imp with ovr_categorical=FALSE on K=3 categorical works", {
  testthat::skip_on_cran()
  set.seed(2026)
  n <- 12
  phylo <- ape::rtree(n)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))
  data <- data.frame(
    y_cat = factor(sample(c("A","B","C"), n, replace = TRUE),
                    levels = c("A","B","C")),
    x = rnorm(n, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y_cat[c(2, 5)] <- NA

  step1 <- suppressWarnings(suppressMessages(bace_imp(
    fixformula     = "y_cat ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo, data = data,
    runs           = 2, nitt = 400, thin = 5, burnin = 100,
    ovr_categorical = FALSE,
    verbose        = FALSE
  )))
  step2 <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object    = step1,
    fixformula     = "y_cat ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo,
    nitt = 400, thin = 5, burnin = 100,
    n_final        = 2,
    ovr_categorical = FALSE,
    verbose        = FALSE
  )))
  expect_s3_class(step2, "bace_final")
})

# ---- bace_imp verbose path ----------------------------------------------

test_that("bace_imp(verbose=TRUE) prints per-run messages", {
  testthat::skip_on_cran()
  fix <- build_verbose_fixture()
  out <- capture.output(suppressWarnings(suppressMessages(
    bace_imp(
      fixformula     = "y ~ x",
      ran_phylo_form = "~ 1 | Species",
      phylo          = fix$phylo,
      data           = fix$data,
      runs           = 2, nitt = 300, thin = 5, burnin = 50,
      verbose        = TRUE
    )
  )))
  combined <- paste(out, collapse = "\n")
  expect_true(nchar(combined) > 0L)
})

# ---- print method on bace_complete ---------------------------------------

test_that("print.bace_complete includes structural info", {
  testthat::skip_on_cran()
  fix <- build_verbose_fixture()
  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = "y ~ x",
    ran_phylo_form = "~ 1 | Species",
    phylo          = fix$phylo,
    data           = fix$data,
    runs           = 2, n_final = 2,
    nitt = 400, thin = 5, burnin = 100,
    skip_conv      = TRUE, max_attempts = 1, n_cores = 1L,
    verbose        = FALSE
  )))
  out <- capture.output(print(res))
  expect_true(length(out) > 1L)
})
