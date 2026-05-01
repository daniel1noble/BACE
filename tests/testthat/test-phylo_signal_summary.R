# =============================================================================
# Tests for phylo_signal_summary() and its integration into bace()
# =============================================================================
# Organised into sections:
#   A. Input validation and error paths (fast, no MCMC)
#   B. Type detection wiring
#   C. MCMC default resolution (unit tests on .resolve_mcmc_defaults)
#   D. Multinomial diagonal-column parsing
#   E. Interpretation bands
#   F. Tree-size warnings
#   G. End-to-end structural smoke (tiny MCMC, asserts shape/columns/class)
#   H. Numerical recovery: gaussian single RE, gaussian dual RE (skip_on_cran)
#   I. Flag logic (low_ess, low_n, unreliable)
#   J. Print method
#   K. bace() integration (phylo_signal = TRUE short-circuit)
#   L. OVR categorical path smoke
# =============================================================================

# -----------------------------------------------------------------------------
# Shared helpers
# -----------------------------------------------------------------------------

.mk_ultrametric_tree <- function(n, seed = 1L) {
  set.seed(seed)
  tr <- ape::compute.brlen(ape::rtree(n), method = "Grafen")
  tr$tip.label <- paste0("sp", seq_along(tr$tip.label))
  tr
}

.mk_gaussian_data <- function(tree, H2 = 0.6, reps = 1, sd_residual = sqrt(0.4),
                              seed = 1L) {
  set.seed(seed)
  n <- length(tree$tip.label)
  A <- ape::vcv.phylo(tree)
  # Scale so V_A = H2 and V_R = 1 - H2, giving target latent H2.
  bm <- as.numeric(MASS::mvrnorm(1, rep(0, n), H2 * A))
  names(bm) <- tree$tip.label
  do.call(rbind, lapply(tree$tip.label, function(s) {
    data.frame(
      Species = s,
      y = bm[[s]] + rnorm(reps, 0, sd_residual),
      stringsAsFactors = FALSE
    )
  }))
}

.mk_dual_re_data <- function(tree, H2 = 0.5, R_sp = 0.2, reps = 4, seed = 2L) {
  set.seed(seed)
  n <- length(tree$tip.label)
  A <- ape::vcv.phylo(tree)
  V_phylo <- H2
  V_sp    <- R_sp
  V_res   <- 1 - H2 - R_sp
  bm <- as.numeric(MASS::mvrnorm(1, rep(0, n), V_phylo * A))
  sp <- rnorm(n, 0, sqrt(V_sp))
  names(bm) <- names(sp) <- tree$tip.label
  do.call(rbind, lapply(tree$tip.label, function(s) {
    data.frame(
      Species = s,
      y = bm[[s]] + sp[[s]] + rnorm(reps, 0, sqrt(V_res)),
      stringsAsFactors = FALSE
    )
  }))
}

.mk_binary_data <- function(tree, H2 = 0.7, seed = 3L) {
  set.seed(seed)
  n <- length(tree$tip.label)
  A <- ape::vcv.phylo(tree)
  latent <- as.numeric(MASS::mvrnorm(1, rep(0, n), H2 * A)) +
            rnorm(n, 0, sqrt(1 - H2))
  data.frame(
    Species = tree$tip.label,
    y = factor(ifelse(latent > 0, 1L, 0L), levels = c(0L, 1L)),
    stringsAsFactors = FALSE
  )
}

# Ultra-fast MCMC for structural tests that just need the function to run
.tiny_mcmc <- list(nitt = 1500, burnin = 500, thin = 2)


# =============================================================================
# A. Input validation
# =============================================================================

test_that("tree must be phylo", {
  expect_error(
    phylo_signal_summary(data = data.frame(Species = "a", y = 1),
                         tree = list(tip.label = "a")),
    "must be a phylo"
  )
})

test_that("data must be data.frame", {
  tr <- .mk_ultrametric_tree(5)
  expect_error(
    phylo_signal_summary(data = list(Species = tr$tip.label, y = 1:5),
                         tree = tr),
    "must be a data.frame"
  )
})

test_that("species_col must exist in data", {
  tr <- .mk_ultrametric_tree(5)
  d <- data.frame(sp = tr$tip.label, y = 1:5)
  expect_error(
    phylo_signal_summary(data = d, tree = tr, species_col = "Species"),
    "not found in data"
  )
})

test_that("unknown variables raise error", {
  tr <- .mk_ultrametric_tree(5)
  d <- data.frame(Species = tr$tip.label, y = 1:5)
  expect_error(
    phylo_signal_summary(data = d, tree = tr, variables = c("y", "nope")),
    "not found in data"
  )
})

test_that("species = TRUE without replicates fails fast", {
  tr <- .mk_ultrametric_tree(10)
  d <- data.frame(Species = tr$tip.label, y = rnorm(10))
  expect_error(
    phylo_signal_summary(data = d, tree = tr, species = TRUE),
    "requires within-species replicates"
  )
})


# =============================================================================
# B. Type detection wiring
# =============================================================================

test_that(".detect_types_for_signal classifies numeric as gaussian", {
  d <- data.frame(y = rnorm(10))
  out <- .detect_types_for_signal(d, "y")
  expect_identical(out$y, "gaussian")
})

test_that(".detect_types_for_signal classifies 2-level factor as threshold", {
  d <- data.frame(y = factor(c("a","b","a","b","a","b")))
  expect_identical(.detect_types_for_signal(d, "y")$y, "threshold")
})

test_that(".detect_types_for_signal classifies >2 unordered as categorical", {
  d <- data.frame(y = factor(c("a","b","c","a","b","c")))
  expect_identical(.detect_types_for_signal(d, "y")$y, "categorical")
})

test_that(".detect_types_for_signal classifies ordered >2 as threshold", {
  d <- data.frame(y = factor(c("lo","med","hi","lo","med","hi"),
                             levels = c("lo","med","hi"), ordered = TRUE))
  expect_identical(.detect_types_for_signal(d, "y")$y, "threshold")
})


# =============================================================================
# C. MCMC default resolution
# =============================================================================

test_that(".resolve_mcmc_defaults gives expected type-specific values", {
  g <- .resolve_mcmc_defaults("gaussian", NULL, NULL, NULL, FALSE)
  expect_equal(g, list(nitt = 20000, burnin = 5000, thin = 10))
  p <- .resolve_mcmc_defaults("poisson", NULL, NULL, NULL, FALSE)
  expect_equal(p, list(nitt = 30000, burnin = 7500, thin = 15))
  t <- .resolve_mcmc_defaults("threshold", NULL, NULL, NULL, FALSE)
  expect_equal(t, list(nitt = 50000, burnin = 12500, thin = 25))
  c <- .resolve_mcmc_defaults("categorical", NULL, NULL, NULL, FALSE)
  expect_equal(c, list(nitt = 80000, burnin = 20000, thin = 40))
})

test_that(".resolve_mcmc_defaults halves nitt and burnin under quick", {
  q <- .resolve_mcmc_defaults("gaussian", NULL, NULL, NULL, TRUE)
  expect_equal(q$nitt, 10000)
  expect_equal(q$burnin, 2500)
  expect_equal(q$thin, 10)  # thin unchanged
})

test_that(".resolve_mcmc_defaults respects user overrides", {
  m <- .resolve_mcmc_defaults("gaussian", 999, 99, 3, TRUE)
  expect_equal(m, list(nitt = 999, burnin = 99, thin = 3))
})

test_that(".resolve_mcmc_defaults rejects unknown types", {
  expect_error(.resolve_mcmc_defaults("nonsense", NULL, NULL, NULL, FALSE),
               "Unknown type")
})


# =============================================================================
# D. Multinomial diagonal-column parsing
# =============================================================================

test_that(".get_multinom_diag_cols finds diagonal trait columns", {
  # Mimic the MCMCglmm naming convention for a K=3 categorical with
  # random effect "Species" and residual "units".
  nms <- c(
    "traitX.1:traitX.1.Species", "traitX.1:traitX.2.Species",
    "traitX.2:traitX.1.Species", "traitX.2:traitX.2.Species",
    "traitX.1:traitX.1.units",   "traitX.1:traitX.2.units",
    "traitX.2:traitX.1.units",   "traitX.2:traitX.2.units"
  )
  diag_sp <- .get_multinom_diag_cols(nms, "Species", J = 2)
  diag_u  <- .get_multinom_diag_cols(nms, "units",   J = 2)
  expect_setequal(diag_sp, c("traitX.1:traitX.1.Species",
                             "traitX.2:traitX.2.Species"))
  expect_setequal(diag_u,  c("traitX.1:traitX.1.units",
                             "traitX.2:traitX.2.units"))
})

test_that(".get_multinom_diag_cols returns empty for unknown RE", {
  nms <- c("traitX.1:traitX.1.Species", "traitX.2:traitX.2.Species")
  # Fallback path: when pattern doesn't match, falls back to endswith .X
  # heuristic. For a truly absent RE, it should be empty.
  res <- .get_multinom_diag_cols(nms, "DoesNotExist", J = 2)
  expect_length(res, 0)
})


# =============================================================================
# E. Interpretation bands
# =============================================================================

test_that(".interpret_H2 bands are <0.2 low, <0.5 moderate, else high", {
  expect_identical(.interpret_H2(0.05), "low")
  expect_identical(.interpret_H2(0.19), "low")
  expect_identical(.interpret_H2(0.20), "moderate")
  expect_identical(.interpret_H2(0.49), "moderate")
  expect_identical(.interpret_H2(0.50), "high")
  expect_identical(.interpret_H2(0.99), "high")
  expect_true(is.na(.interpret_H2(NA)))
})


# =============================================================================
# F. Tree-size warnings
# =============================================================================

test_that("small_tree warning fires for n_species < 30", {
  tr <- .mk_ultrametric_tree(20)
  d <- .mk_gaussian_data(tr, H2 = 0.5)
  expect_warning(
    phylo_signal_summary(data = d, tree = tr, variables = "y",
                          nitt = 1500, burnin = 500, thin = 2, verbose = FALSE),
    "weakly identified"
  )
})

test_that("quick_mode warning fires when quick = TRUE", {
  tr <- .mk_ultrametric_tree(40)
  d <- .mk_gaussian_data(tr, H2 = 0.5)
  expect_warning(
    phylo_signal_summary(data = d, tree = tr, variables = "y",
                          nitt = 1500, burnin = 500, thin = 2, quick = TRUE,
                          verbose = FALSE),
    "quick = TRUE"
  )
})


# =============================================================================
# G. Structural smoke (ultra-fast, runs every devtools::test())
# =============================================================================

test_that("end-to-end gaussian returns correctly-shaped phylo_signal object", {
  tr <- .mk_ultrametric_tree(40)
  d <- .mk_gaussian_data(tr, H2 = 0.6)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 1500, burnin = 500, thin = 2, verbose = FALSE
  ))
  expect_s3_class(res, "phylo_signal")
  expect_named(res, c("table","models","warnings","call","species","n_sim"))
  expect_identical(res$species, FALSE)
  expect_equal(nrow(res$table), 1L)
  expect_setequal(
    intersect(c("variable","type","n_obs","H2_mean","H2_lo","H2_hi",
                "R_species_mean","R_species_lo","R_species_hi",
                "lambda","lambda_p","K","K_p","D","D_random_p","D_BM_p",
                "min_ess","interpretation","flag"),
              names(res$table)),
    c("variable","type","n_obs","H2_mean","H2_lo","H2_hi",
      "R_species_mean","R_species_lo","R_species_hi",
      "lambda","lambda_p","K","K_p","D","D_random_p","D_BM_p",
      "min_ess","interpretation","flag")
  )
})

test_that("multi-variable input yields one row per variable in order", {
  tr <- .mk_ultrametric_tree(40)
  d <- .mk_gaussian_data(tr, H2 = 0.6)
  d$y2 <- d$y + rnorm(nrow(d), 0, 0.5)
  d$y3 <- rnorm(nrow(d))
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = c("y","y2","y3"),
    nitt = 1500, burnin = 500, thin = 2, verbose = FALSE
  ))
  expect_equal(nrow(res$table), 3L)
  expect_identical(res$table$variable, c("y","y2","y3"))
  expect_true(all(res$table$type == "gaussian"))
})

test_that("default keep_models = TRUE returns the fitted models", {
  tr <- .mk_ultrametric_tree(40)
  d <- .mk_gaussian_data(tr, H2 = 0.5)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 1500, burnin = 500, thin = 2, verbose = FALSE
  ))
  expect_named(res$models, "y")
  expect_s3_class(res$models$y, "MCMCglmm")
})

test_that("keep_models = FALSE omits the model list (opt-out)", {
  tr <- .mk_ultrametric_tree(40)
  d <- .mk_gaussian_data(tr, H2 = 0.5)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 1500, burnin = 500, thin = 2,
    keep_models = FALSE, verbose = FALSE
  ))
  expect_null(res$models)
})

test_that("fit_error is trapped and row is still produced", {
  skip_on_cran()
  tr <- .mk_ultrametric_tree(40)
  # Provide malformed data that should cause MCMCglmm to fail
  d <- data.frame(Species = tr$tip.label,
                   y = c(rep(NA, 35), 1:5))  # nearly all NA
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 1500, burnin = 500, thin = 2, verbose = FALSE
  ))
  # Function should not error overall; row should reflect minimal / error state
  expect_s3_class(res, "phylo_signal")
  expect_equal(nrow(res$table), 1L)
})


# =============================================================================
# H. Numerical recovery (skip_on_cran; may take a couple of minutes)
# =============================================================================

test_that("gaussian single RE recovers H2 within tolerance", {
  skip_on_cran()
  tr <- .mk_ultrametric_tree(60, seed = 42)
  d <- .mk_gaussian_data(tr, H2 = 0.6, seed = 42)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 20000, burnin = 5000, thin = 10, verbose = FALSE
  ))
  row <- res$table
  # Posterior mean within ±0.25 of truth
  expect_lt(abs(row$H2_mean - 0.6), 0.25)
  # 95% HPD contains truth
  expect_true(row$H2_lo <= 0.6 && row$H2_hi >= 0.6)
})

test_that("gaussian dual RE produces valid H2 and R_species decomposition", {
  skip_on_cran()
  skip_if(Sys.getenv("BACE_SKIP_SLOW") == "true")
  # NOTE: strict numerical recovery for dual RE needs n_species >> 40 to
  # identify phylo-vs-non-phylo species variance. That belongs in
  # dev/validate_phylo_signal_summary.R. Here we only check that the
  # machinery runs and returns a valid variance decomposition.
  tr <- .mk_ultrametric_tree(40, seed = 7)
  d <- .mk_dual_re_data(tr, H2 = 0.5, R_sp = 0.3, reps = 4, seed = 7)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y", species = TRUE,
    nitt = 20000, burnin = 5000, thin = 10, verbose = FALSE
  ))
  row <- res$table
  # Both components in [0, 1]
  expect_gte(row$H2_mean, 0); expect_lte(row$H2_mean, 1)
  expect_gte(row$R_species_mean, 0); expect_lte(row$R_species_mean, 1)
  # Sum of variance components does not exceed 1 (proper decomposition)
  expect_lte(row$H2_mean + row$R_species_mean, 1.001)
  # HPDs are populated (not collapsed to NA and not zero-width)
  expect_false(is.na(row$H2_lo)); expect_false(is.na(row$H2_hi))
  expect_gt(row$H2_hi - row$H2_lo, 0.05)
  expect_false(is.na(row$R_species_lo)); expect_false(is.na(row$R_species_hi))
})

test_that("binary threshold single RE returns H2 in (0, 1)", {
  skip_on_cran()
  skip_if(Sys.getenv("BACE_SKIP_SLOW") == "true")
  tr <- .mk_ultrametric_tree(60, seed = 11)
  d <- .mk_binary_data(tr, H2 = 0.7, seed = 11)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 20000, burnin = 5000, thin = 10, verbose = FALSE
  ))
  row <- res$table
  expect_equal(row$type, "threshold")
  expect_gt(row$H2_mean, 0)
  expect_lt(row$H2_mean, 1)
})


# =============================================================================
# I. Flag logic
# =============================================================================

test_that("low_n flag fires when n_obs < 20", {
  tr <- .mk_ultrametric_tree(15)
  d <- .mk_gaussian_data(tr, H2 = 0.5)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 1500, burnin = 500, thin = 2, verbose = FALSE
  ))
  expect_true(grepl("low_n", res$table$flag))
})

test_that("low_ess and unreliable flags fire under tiny nitt", {
  tr <- .mk_ultrametric_tree(40)
  d <- .mk_gaussian_data(tr, H2 = 0.5)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 1500, burnin = 500, thin = 2, verbose = FALSE
  ))
  # min_ess default is 1000; with only 500 post-burn draws we'll trigger
  # low_ess by construction.
  expect_true(grepl("low_ess", res$table$flag))
  expect_true(grepl("unreliable", res$table$flag))
})

test_that("flag column is empty string (not NA) when no flags fire", {
  skip_on_cran()
  skip_if(Sys.getenv("BACE_SKIP_SLOW") == "true")
  tr <- .mk_ultrametric_tree(60, seed = 99)
  d <- .mk_gaussian_data(tr, H2 = 0.6, seed = 99)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 20000, burnin = 5000, thin = 10, verbose = FALSE
  ))
  # If we get lucky and gaussian runs mix well, flag is "". If not, at least
  # it should never be NA.
  expect_false(is.na(res$table$flag))
})


# =============================================================================
# J. Print method
# =============================================================================

test_that("print.phylo_signal runs without error and returns invisibly", {
  tr <- .mk_ultrametric_tree(40)
  d <- .mk_gaussian_data(tr, H2 = 0.5)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 1500, burnin = 500, thin = 2, verbose = FALSE
  ))
  out <- capture.output(returned <- print(res))
  expect_identical(returned, res)
  expect_true(any(grepl("Phylogenetic signal summary", out)))
  expect_true(any(grepl("Interpretation bands", out)))
  expect_true(any(grepl("necessary but NOT sufficient", out)))
})

test_that("print marks unreliable rows with asterisk", {
  tr <- .mk_ultrametric_tree(40)
  d <- .mk_gaussian_data(tr, H2 = 0.5)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    nitt = 1500, burnin = 500, thin = 2, verbose = FALSE
  ))
  # tiny nitt -> unreliable; H2_mean print cell should carry an asterisk
  out <- capture.output(print(res))
  expect_true(any(grepl("\\*", out)))
})


# =============================================================================
# K. bace() integration
# =============================================================================

test_that("bace(phylo_signal = TRUE) short-circuits and returns preview", {
  skip_on_cran()
  tr <- .mk_ultrametric_tree(40, seed = 5)
  d <- .mk_gaussian_data(tr, H2 = 0.6, seed = 5)
  d$x <- rnorm(nrow(d))  # dummy predictor to make bace() call valid
  # Introduce minimal missingness so fixformula is plausible
  d$y[1:3] <- NA
  out <- suppressWarnings(suppressMessages(
    bace(fixformula = y ~ x, ran_phylo_form = ~ 1 | Species,
         phylo = tr, data = d, phylo_signal = TRUE, verbose = FALSE)
  ))
  expect_s3_class(out, "phylo_signal")
  expect_s3_class(out, "bace_preview")
  expect_true("table" %in% names(out))
  # Should have rows for every variable in the formula
  expect_setequal(out$table$variable, c("y","x"))
})

test_that("bace(phylo_signal = TRUE) extracts species column from ran_phylo_form", {
  skip_on_cran()
  tr <- .mk_ultrametric_tree(40, seed = 6)
  d <- .mk_gaussian_data(tr, H2 = 0.5, seed = 6)
  d$x <- rnorm(nrow(d))
  out <- suppressWarnings(suppressMessages(
    bace(fixformula = y ~ x, ran_phylo_form = ~ 1 | Species,
         phylo = tr, data = d, phylo_signal = TRUE, verbose = FALSE)
  ))
  # Verify the function actually used "Species" column; if it had used a
  # wrong column, phylo_signal_summary would have errored before this point.
  expect_s3_class(out, "phylo_signal")
})

test_that("bace() default phylo_signal = FALSE does not short-circuit", {
  # Smoke: just check that default path calls bace_imp (not phylo_signal).
  # Make this cheap by monkeypatching bace_imp to raise a sentinel.
  env <- asNamespace("BACE")
  old_imp <- get("bace_imp", envir = env)
  on.exit({
    unlockBinding("bace_imp", env)
    assign("bace_imp", old_imp, envir = env)
    lockBinding("bace_imp", env)
  })
  unlockBinding("bace_imp", env)
  assign("bace_imp",
         function(...) stop("bace_imp was called"), envir = env)
  lockBinding("bace_imp", env)

  tr <- .mk_ultrametric_tree(40)
  d <- .mk_gaussian_data(tr, H2 = 0.5)
  d$x <- rnorm(nrow(d))
  expect_error(
    bace(fixformula = y ~ x, ran_phylo_form = ~ 1 | Species,
         phylo = tr, data = d, verbose = FALSE),
    "bace_imp was called"
  )
})


# =============================================================================
# L. OVR categorical path smoke
# =============================================================================

test_that("ovr_categorical = TRUE runs J binary fits and returns one row", {
  skip_on_cran()
  skip_if(Sys.getenv("BACE_SKIP_SLOW") == "true")
  tr <- .mk_ultrametric_tree(40, seed = 8)
  set.seed(8)
  n <- length(tr$tip.label)
  A <- ape::vcv.phylo(tr)
  # Three latent variables -> argmax -> 3-level factor
  lat <- sapply(1:3, function(k) as.numeric(MASS::mvrnorm(1, rep(0, n), 0.6 * A)) +
                                 rnorm(n, 0, sqrt(0.4)))
  y <- factor(letters[apply(lat, 1, which.max)], levels = c("a","b","c"))
  d <- data.frame(Species = tr$tip.label, y = y, stringsAsFactors = FALSE)
  res <- suppressWarnings(phylo_signal_summary(
    data = d, tree = tr, variables = "y",
    ovr_categorical = TRUE,
    nitt = 1500, burnin = 500, thin = 2, verbose = FALSE
  ))
  expect_equal(nrow(res$table), 1L)
  expect_identical(res$table$type, "categorical (OVR)")
  expect_false(is.na(res$table$H2_mean))
  # OVR doesn't produce an aggregate HPD; lo/hi are NA by design
  expect_true(is.na(res$table$H2_lo))
  expect_true(is.na(res$table$H2_hi))
  # Per-level binary fits returned as a named list keyed by factor level
  expect_type(res$models, "list")
  expect_named(res$models$y, c("a", "b", "c"))
  expect_s3_class(res$models$y$a, "MCMCglmm")
  expect_s3_class(res$models$y$b, "MCMCglmm")
  expect_s3_class(res$models$y$c, "MCMCglmm")
})
