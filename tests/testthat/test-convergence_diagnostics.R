# Tests for R/convergence_diagnostics.R — the chained-equations stability check.
#
# The pre-1.1 suite only smoke-tested this file (expect_s3_class and friends),
# which is how a check that failed settled chains ~97% of the time shipped.
# These tests assert VERDICTS on series whose true status is known by
# construction: settled, drifting, oscillating, collapsing, constant — at a
# trait mean of 0 (where the old %-change criterion exploded) and far from 0.

# ---- the internal rule on constructed series ------------------------------

# A settled chain: iid noise around a constant, no drift, no trend. The noise
# is the Monte-Carlo floor of a mean over n_miss imputed cells.
settled_series <- function(runs, level = 0, n_miss = 30, sd_within = 1, seed = 1) {
  set.seed(seed)
  list(m = level + rnorm(runs, 0, sd_within * sqrt(2 / n_miss)),
       s = rep(sd_within, runs) * exp(rnorm(runs, 0, 0.05)))
}

test_that("a settled chain passes, wherever its mean sits (the old bug)", {
  for (level in c(0, 20)) {
    ch <- settled_series(10, level = level)
    out <- BACE:::.conv_assess_series(ch$m, ch$s, n_cells = 30)
    expect_identical(out$verdict, "ok")
  }
})

test_that("the verdict is location- and scale-invariant", {
  ch <- settled_series(10)
  base <- BACE:::.conv_assess_series(ch$m, ch$s, n_cells = 30)$verdict
  shifted <- BACE:::.conv_assess_series(ch$m * 1000 + 500, ch$s * 1000,
                                        n_cells = 30)$verdict
  expect_identical(base, shifted)
})

test_that("a steady drift is caught at 10 runs, at any level", {
  for (level in c(0, 20)) {
    set.seed(42)
    m <- level + 0.5 * seq_len(10) + rnorm(10, 0, sqrt(2 / 30))
    s <- rep(1, 10)
    out <- BACE:::.conv_assess_series(m, s, n_cells = 30)
    expect_identical(out$verdict, "drift")
  }
})

test_that("a loud flip-flop is caught at 10 runs", {
  set.seed(7)
  m <- 0.75 * (-1)^seq_len(10) + rnorm(10, 0, sqrt(2 / 30))
  s <- rep(1, 10)
  out <- BACE:::.conv_assess_series(m, s, n_cells = 30)
  expect_identical(out$verdict, "oscillation")
})

test_that("a collapsing imputation spread is flagged as stuck, not converged", {
  # The pre-1.1 code special-cased zero variance as 'converged' — exactly
  # backwards: this is van Buuren's locked-up sampler.
  set.seed(8)
  s <- 0.5^(seq_len(10) - 1)
  m <- rnorm(10, 0, sqrt(2 / 30)) * s
  out <- BACE:::.conv_assess_series(m, s, n_cells = 30)
  expect_identical(out$verdict, "stuck")
})

test_that("a genuinely constant chain passes, with a note", {
  out <- BACE:::.conv_assess_series(rep(2.5, 10), rep(0, 10), n_cells = 30)
  expect_identical(out$verdict, "ok")
  expect_match(out$note, "identical imputations")
})

test_that("a geometric approach that settles by mid-chain passes at 10 runs", {
  set.seed(9)
  m <- 2 * 0.5^seq_len(10) + rnorm(10, 0, sqrt(2 / 30))
  s <- rep(1, 10)
  out <- BACE:::.conv_assess_series(m, s, n_cells = 30)
  expect_identical(out$verdict, "ok")
})

test_that("quantised alternation at the noise floor gets a note, not a failure", {
  # A proportion over 30 cells moves in 1/30 jumps and can alternate between
  # two values (lag-1 acf = -1) while moving only at the noise floor. Real
  # chains with exactly this structure produced nominal coverage.
  m <- rep(c(0.033, 0.100), 5)
  s <- sqrt(m * (1 - m))
  out <- BACE:::.conv_assess_series(m, s, n_cells = 30)
  expect_identical(out$verdict, "ok")
  expect_match(out$note, "alternating")
})

test_that("short chains are honest: 5 runs still pass settled chains", {
  ch <- settled_series(5)
  out <- BACE:::.conv_assess_series(ch$m, ch$s, n_cells = 30)
  expect_identical(out$verdict, "ok")
})

# ---- assess_convergence end-to-end on a fabricated bace object ------------
# No MCMC needed: we build the per-run imputed datasets directly, so the test
# pins the full pipeline (series extraction -> rule -> verdicts -> flag).

fake_bace <- function(y_gen, runs = 10, n = 40, miss_rows = 1:10, seed = 3) {
  set.seed(seed)
  base <- data.frame(y = rnorm(n), g = rnorm(n))
  dats <- lapply(seq_len(runs), function(r) {
    d <- base
    d$y[miss_rows] <- y_gen(r, length(miss_rows))
    d
  })
  obj <- list(
    data = c(list(base), dats),   # element 1 is the un-imputed frame
    types = list(y = "gaussian"),
    miss_dat = data.frame(row = miss_rows, colname = "y",
                          stringsAsFactors = FALSE)
  )
  class(obj) <- "bace"
  obj
}

test_that("assess_convergence passes a settled fit and reports structure", {
  bo <- fake_bace(function(r, k) rnorm(k, 5, 1))
  out <- assess_convergence(bo)
  expect_s3_class(out, "bace_convergence")
  expect_true(out$converged)
  expect_identical(unname(out$verdicts["y"]), "ok")
  expect_true(all(c("iteration", "y") %in% names(out$summary_stats)))
  expect_true(all(c("iteration", "y") %in% names(out$sd_stats)))
})

test_that("assess_convergence fails a drifting fit and names the drift", {
  bo <- fake_bace(function(r, k) rnorm(k, r * 1.5, 1))
  out <- assess_convergence(bo)
  expect_false(out$converged)
  expect_identical(unname(out$verdicts["y"]), "drift")
})

test_that("assess_convergence flags a collapsing fit as stuck", {
  bo <- fake_bace(function(r, k) rnorm(k, 0, 0.5^r))
  out <- assess_convergence(bo)
  expect_false(out$converged)
  expect_identical(unname(out$verdicts["y"]), "stuck")
})

test_that("an unrecognised variable type is loud, not silently converged", {
  # The pre-1.1 code fell through silently, leaving an all-zero series that
  # auto-passed every criterion.
  bo <- fake_bace(function(r, k) rnorm(k, 5, 1))
  bo$types$y <- "ordinal"
  expect_error(assess_convergence(bo), regexp = "unrecognised type")
})

# ---- API: defunct and deprecated methods ----------------------------------

test_that("the defective pre-1.1 methods are defunct with an explanation", {
  bo <- fake_bace(function(r, k) rnorm(k, 5, 1))
  for (m in c("summary.acf", "summary.percentage", "summary.trend",
              "summary.geweke", "energy", "wasserstein", "all")) {
    expect_error(assess_convergence(bo, method = m), regexp = "defunct")
  }
})

test_that("method = 'summary' is a deprecated alias for 'stability'", {
  bo <- fake_bace(function(r, k) rnorm(k, 5, 1))
  expect_warning(out <- assess_convergence(bo, method = "summary"),
                 regexp = "deprecated")
  expect_true(out$converged)
})

test_that("the removed tuning arguments fail loudly, not silently", {
  bo <- fake_bace(function(r, k) rnorm(k, 5, 1))
  expect_error(assess_convergence(bo, pct_change_threshold = 0.1),
               regexp = "unused argument")
  expect_error(assess_convergence(bo, lag = 2), regexp = "unused argument")
})

# ---- input validation (kept from the old suite) ---------------------------

test_that("assess_convergence rejects non-bace input", {
  expect_error(assess_convergence(list()),
               regexp = "must be an object of class 'bace'")
})

test_that("assess_convergence warns when too few iterations", {
  fake <- list(
    data     = list(NULL, data.frame(y = c(1, 2)), data.frame(y = c(1, 3))),
    types    = list(y = "gaussian"),
    miss_dat = data.frame(row = 1, colname = "y", stringsAsFactors = FALSE)
  )
  class(fake) <- "bace"
  expect_warning(
    res <- assess_convergence(fake, min_iterations = 5),
    regexp = "Insufficient iterations|Only"
  )
  expect_true(is.na(res$converged))
})

test_that("assess_convergence filters to specified variables and stops when none match", {
  bo <- fake_bace(function(r, k) rnorm(k, 5, 1))
  out <- suppressWarnings(assess_convergence(bo, variables = "y"))
  expect_s3_class(out, "bace_convergence")
  expect_error(
    suppressWarnings(assess_convergence(bo, variables = "nonexistent_var")),
    regexp = "no missing data|None of"
  )
})

# ---- print method ----------------------------------------------------------

test_that("print.bace_convergence prints verdicts and the one-way disclaimer", {
  bo <- fake_bace(function(r, k) rnorm(k, 5, 1))
  out <- assess_convergence(bo)
  expect_output(print(out), "Stability Check")
  expect_output(print(out), "no clear-cut test")
})

test_that("print names the remedy when a pathology is found", {
  bo <- fake_bace(function(r, k) rnorm(k, r * 1.5, 1))
  out <- assess_convergence(bo)
  expect_output(print(out), "more chained runs")
})

# ---- a real (tiny) bace_imp fit still flows through ------------------------

test_that("assess_convergence runs end-to-end on a real tiny fit", {
  testthat::skip_on_cran()
  set.seed(2026)
  phylo <- ape::rtree(12)
  phylo <- ape::compute.brlen(phylo, method = "Grafen")
  phylo$edge.length <- phylo$edge.length /
    max(ape::node.depth.edgelength(phylo))

  data <- data.frame(
    y = rnorm(12, 0, 1),
    x1 = rnorm(12, 0, 1),
    x2 = rnorm(12, 0, 1),
    Species = phylo$tip.label,
    stringsAsFactors = FALSE
  )
  data$y[c(1, 4, 7)] <- NA
  data$x1[c(2, 5)]   <- NA

  bo <- suppressMessages(suppressWarnings(bace_imp(
    fixformula     = list("y ~ x1 + x2", "x1 ~ y + x2"),
    ran_phylo_form = "~ 1 | Species",
    phylo          = phylo,
    data           = data,
    runs           = 5,
    nitt           = 600, thin = 5, burnin = 200,
    verbose        = FALSE
  )))
  out <- suppressWarnings(assess_convergence(bo))
  expect_s3_class(out, "bace_convergence")
  expect_true(all(out$verdicts %in% c("ok", "drift", "oscillation", "stuck")))
  expect_identical(sort(names(out$verdicts)), c("x1", "y"))
  expect_output(print(out), "Stability Check")
})

# ============================================================================
# ACCURACY: measured operating characteristics, pinned as rates.
# 200 constructed chains per truth state, fixed seeds, so the result is
# deterministic: if these rates ever drop, the check got worse and the test
# fails. Calibration provenance: .agents/investigation-convergence-2026-11.md.
# ============================================================================

conv_verdicts <- function(gen, B = 200) {
  vapply(seq_len(B), function(b) {
    ch <- gen(b)
    BACE:::.conv_assess_series(ch$m, ch$s, n_cells = 30)$verdict
  }, character(1))
}
mc_noise <- sqrt(2 / 30)

test_that("accuracy: >= 97% of settled 10-run chains pass", {
  v <- conv_verdicts(function(b) {
    set.seed(b)
    list(m = rnorm(10, 0, mc_noise), s = exp(rnorm(10, 0, 0.05)))
  })
  expect_gte(mean(v == "ok"), 0.97)
})

test_that("accuracy: the verdict distribution is identical when every chain is shifted by +20", {
  # location invariance as a measured property, not a single example;
  # the old %-change criterion failed THIS: its pass rate moved from ~0%
  # at mean 0 to ~100% at mean 20 on identical chains.
  v0 <- conv_verdicts(function(b) {
    set.seed(b)
    list(m = rnorm(10, 0, mc_noise), s = exp(rnorm(10, 0, 0.05)))
  })
  v20 <- conv_verdicts(function(b) {
    set.seed(b)
    list(m = 20 + rnorm(10, 0, mc_noise), s = exp(rnorm(10, 0, 0.05)))
  })
  expect_identical(v0, v20)
})

test_that("accuracy: >= 95% of steady drifts (0.5 spread/run) are caught at 10 runs", {
  v <- conv_verdicts(function(b) {
    set.seed(b)
    list(m = 0.5 * seq_len(10) + rnorm(10, 0, mc_noise),
         s = exp(rnorm(10, 0, 0.05)))
  })
  expect_gte(mean(v == "drift"), 0.95)
})

test_that("accuracy: >= 95% of loud flip-flops are caught at 10 runs", {
  v <- conv_verdicts(function(b) {
    set.seed(b)
    list(m = 0.75 * (-1)^seq_len(10) + rnorm(10, 0, mc_noise),
         s = exp(rnorm(10, 0, 0.05)))
  })
  expect_gte(mean(v == "oscillation"), 0.95)
})

test_that("accuracy: every collapsing chain is caught", {
  v <- conv_verdicts(function(b) {
    set.seed(b)
    s <- 0.5^(seq_len(10) - 1)
    list(m = rnorm(10, 0, mc_noise) * s, s = s)
  })
  expect_identical(unique(v), "stuck")
})

test_that("honesty at 5 runs: settled chains still pass; drift power is limited and known", {
  v_ok <- conv_verdicts(function(b) {
    set.seed(b)
    list(m = rnorm(5, 0, mc_noise), s = exp(rnorm(5, 0, 0.05)))
  })
  expect_gte(mean(v_ok == "ok"), 0.97)
  v_dr <- conv_verdicts(function(b) {
    set.seed(b)
    list(m = 0.5 * seq_len(5) + rnorm(5, 0, mc_noise),
         s = exp(rnorm(5, 0, 0.05)))
  })
  # five points cannot support the full battery; pin a floor AND a ceiling so
  # a future "improvement" that claims near-perfect power at 5 runs fails here.
  expect_gte(mean(v_dr == "drift"), 0.60)
  expect_lte(mean(v_dr == "drift"), 0.999)
})

# ============================================================================
# PRECISION: the reported statistics equal their definitions, hand-computed.
# ============================================================================

test_that("precision: stats match hand computation on a deterministic series", {
  # runs = 6 -> warm-up = first 2, kept = runs 3:6 = 30, 40, 50, 60
  # spread = mean(sd over kept) = 2
  # steps = |10|/2 = 5, 5, 5          -> step = 5
  # net: k = floor(4/3) = 1 -> |60 - 30| / 2 = 15
  # Kendall trend on 4 strictly increasing points: p = 2/4! = 1/12 < 0.10
  # -> verdict drift, and ref = sqrt(2/8) = 0.5
  out <- BACE:::.conv_assess_series(c(10, 20, 30, 40, 50, 60), rep(2, 6),
                                    n_cells = 8)
  expect_identical(out$verdict, "drift")
  expect_equal(unname(out$stats["net"]), 15)
  expect_equal(unname(out$stats["step"]), 5)
  expect_equal(unname(out$stats["ref"]), 0.5)
  expect_equal(unname(out$stats["n_kept"]), 4)
})

test_that("precision: the proportion series and its binomial spread are exact", {
  # 10 imputed cells; run r imputes r cells as "A" and the rest "B"
  # -> proportion of "A" at run r is exactly r/10, sd exactly sqrt(p(1-p))
  runs <- 5
  base <- data.frame(y = factor(rep("B", 20), levels = c("A", "B")))
  base$y[1:10] <- NA
  dats <- lapply(seq_len(runs), function(r) {
    d <- base
    lv <- c(rep("A", r), rep("B", 10 - r))
    d$y[1:10] <- factor(lv, levels = c("A", "B"))
    d
  })
  obj <- list(data = c(list(base), dats),
              types = list(y = "categorical"),
              miss_dat = data.frame(row = 1:10, colname = "y",
                                    stringsAsFactors = FALSE))
  class(obj) <- "bace"
  out <- assess_convergence(obj)
  p_expected <- (1:5) / 10
  expect_equal(out$summary_stats$y, p_expected)
  expect_equal(out$sd_stats$y, sqrt(p_expected * (1 - p_expected)))
})

test_that("precision: sd_ratio for a collapsing chain equals its definition", {
  s <- 0.5^(0:9)                       # first-half median = median(s[1:5])
  m <- rep(0, 10)
  out <- BACE:::.conv_assess_series(m, s, n_cells = 30)
  expect_identical(out$verdict, "stuck")
  expect_equal(unname(out$stats["sd_ratio"]),
               mean(tail(s[max(1, floor(10 / 2)):10], 2)) / median(s[1:5]))
})

# ============================================================================
# BREAKS: inputs that used to crash, silently pass, or silently lie.
# ============================================================================

test_that("break: NA runs do not crash the rule, and too few finite points is loud", {
  # all-NA series: the pre-fix code hit `if (NA)` and crashed
  out <- BACE:::.conv_assess_series(rep(NA_real_, 10), rep(NA_real_, 10),
                                    n_cells = 30)
  expect_identical(out$verdict, "unassessable")
  # a few NA runs: dropped pairwise, the rest assessed
  m <- c(NA, rnorm(9, 5, mc_noise))
  out2 <- BACE:::.conv_assess_series(m, c(NA, rep(1, 9)), n_cells = 30)
  expect_identical(out2$verdict, "ok")
})

test_that("break: Inf is treated as non-finite, not as a huge step", {
  m <- c(Inf, rnorm(9, 0, mc_noise))
  out <- BACE:::.conv_assess_series(m, rep(1, 10), n_cells = 30)
  expect_identical(out$verdict, "ok")
})

test_that("break: a single imputed cell (n_miss = 1) flows through without error", {
  # sd of one value is NA -> the builder records 0 -> the rule must cope
  bo <- local({
    base <- data.frame(y = rnorm(20), g = rnorm(20))
    base$y[1] <- NA
    dats <- lapply(1:10, function(r) { d <- base; d$y[1] <- 5 + 0.01 * r; d })
    obj <- list(data = c(list(base), dats),
                types = list(y = "gaussian"),
                miss_dat = data.frame(row = 1, colname = "y",
                                      stringsAsFactors = FALSE))
    class(obj) <- "bace"
    obj
  })
  out <- assess_convergence(bo)
  expect_s3_class(out, "bace_convergence")
  expect_true(out$verdicts[["y"]] %in%
                c("ok", "drift", "oscillation", "stuck", "unassessable"))
})

test_that("break: one bad variable fails the dataset and is named; the good one stays ok", {
  base <- data.frame(y = rnorm(40), g = rnorm(40))
  base$y[1:10] <- NA
  base$g[11:20] <- NA
  set.seed(4)
  dats <- lapply(1:10, function(r) {
    d <- base
    d$y[1:10]  <- rnorm(10, r * 1.5, 1)   # drifting
    d$g[11:20] <- rnorm(10, 3, 1)         # settled
    d
  })
  obj <- list(data = c(list(base), dats),
              types = list(y = "gaussian", g = "gaussian"),
              miss_dat = data.frame(row = c(1:10, 11:20),
                                    colname = rep(c("y", "g"), each = 10),
                                    stringsAsFactors = FALSE))
  class(obj) <- "bace"
  out <- assess_convergence(obj)
  expect_false(out$converged)
  expect_identical(unname(out$verdicts["y"]), "drift")
  expect_identical(unname(out$verdicts["g"]), "ok")
  # and the old non-monotone 70% rule would have PASSED this dataset at
  # 6 variables; the repaired rule must fail it at any count: verify the
  # aggregation is all(), not a proportion
  expect_identical(out$converged, all(out$verdicts == "ok"))
})

test_that("break: binary trait where every run imputes all-one-level, flipping levels", {
  # Within-run spread is 0 at every run while the mean flips 0 <-> 1.
  # A binary trait is monitored as a PROPORTION by .assess_summary_convergence,
  # which carries its own 0-1 scale, so the flip is caught there.
  m <- rep(c(0, 1), 5)
  s <- rep(0, 10)
  expect_identical(
    BACE:::.conv_assess_series(m, s, n_cells = 10, kind = "proportion")$verdict,
    "oscillation")

  # The numeric path catches it too, as long as it is given a scale that does
  # not come from the series itself (here the SD of a 0/1 variable, ~0.5).
  expect_identical(
    BACE:::.conv_assess_series(m, s, n_cells = 10, kind = "numeric",
                               data_scale = 0.5)$verdict,
    "oscillation")

  # With NO independent scale the honest answer is "unassessable". This test
  # used to assert the old behaviour, where `spread` fell back to sd(m) -- the
  # between-run SD of the series under test. That made the statistic
  # scale-free (mean|diff|/sd -> 2/sqrt(pi) = 1.128 for ANY iid wobble), so a
  # perfectly settled variable was called "drift" 81-89% of the time.
  expect_identical(
    BACE:::.conv_assess_series(m, s, n_cells = 10, kind = "numeric")$verdict,
    "unassessable")
})

test_that("a settled numeric series with no within-run spread is not 'drift'", {
  # The regression: one missing cell (sd of length 1 is NA -> recorded as 0),
  # or rounded counts all landing on the same integer, gives s == 0 at every
  # run. Judged against the variable's OWN observed spread, a chain that moves
  # 0.02 on a variable whose data SD is 1.6 is plainly settled.
  set.seed(1)
  verdicts <- vapply(1:400, function(i)
    BACE:::.conv_assess_series(rnorm(10, 5, 0.02), rep(0, 10), n_cells = 1L,
                               kind = "numeric", data_scale = 1.6)$verdict,
    character(1))
  expect_true(mean(verdicts == "ok") > 0.98)

  # more runs must not make a settled chain look worse (it used to: 0.80 at
  # runs = 10, 0.88 at runs = 50)
  v50 <- vapply(1:200, function(i)
    BACE:::.conv_assess_series(rnorm(50, 5, 0.02), rep(0, 50), n_cells = 1L,
                               kind = "numeric", data_scale = 1.6)$verdict,
    character(1))
  expect_true(mean(v50 == "ok") > 0.98)

  # and a REAL drift on the same scale is still caught
  ramp <- vapply(1:200, function(i)
    BACE:::.conv_assess_series(seq(0, 3 * 1.6, length.out = 10) + rnorm(10, 0, 0.02),
                               rep(0, 10), n_cells = 1L, kind = "numeric",
                               data_scale = 1.6)$verdict, character(1))
  expect_true(all(ramp == "drift"))
})

test_that("break: near-constant series just above the old 1e-10 variance cutoff passes", {
  # the pre-1.1 ACF criterion rejected sd ~ 1e-4 'constant' chains because
  # only var < 1e-10 took the shortcut; the repaired rule must pass them
  set.seed(11)
  m <- 2.5 + rnorm(10, 0, 1e-4)
  s <- rep(1, 10)   # imputation spread is healthy; the MEAN is just stable
  out <- BACE:::.conv_assess_series(m, s, n_cells = 30)
  expect_identical(out$verdict, "ok")
})

# ============================================================================
# REGRESSIONS from the 2026-11-09 adversarial review. Every one of these was a
# demonstrated defect in the first version of this check; each test fails if
# the defect returns. Grouped by what the attacker broke.
# ============================================================================

test_that("regression: one NA imputed cell must not read as a collapsed spread", {
  # sd() without na.rm returned NA -> coerced to 0 -> false 'stuck', which
  # bace() then treats as uncurable and refuses to retry.
  set.seed(5)
  out <- BACE:::.conv_assess_series(c(rnorm(9, 5, 0.2), NA),
                                    c(rep(1, 9), NA), n_cells = 30)
  expect_identical(out$verdict, "ok")
})

test_that("regression: a categorical flipping between its NON-first levels is caught", {
  # Monitoring only the first level's proportion is blind to this: level 1
  # stays flat at 0 while every imputed cell flips between levels 2 and 3.
  runs <- 10
  base <- data.frame(y = factor(rep("a", 20), levels = c("a", "b", "c")))
  base$y[1:4] <- NA
  dats <- lapply(seq_len(runs), function(r) {
    d <- base
    d$y[1:4] <- factor(rep(if (r %% 2) "b" else "c", 4), levels = c("a", "b", "c"))
    d
  })
  obj <- list(data = c(list(base), dats), types = list(y = "categorical"),
              miss_dat = data.frame(row = 1:4, colname = "y",
                                    stringsAsFactors = FALSE))
  class(obj) <- "bace"
  out <- assess_convergence(obj)
  expect_false(out$converged)
  expect_identical(unname(out$verdicts["y"]), "oscillation")
  expect_match(out$notes[["y"]], "level")
})

test_that("regression: a drifting proportion is drift, never 'stuck'", {
  # sqrt(p(1-p)) shrinks as p -> extreme, so reading it as an imputation
  # spread turned a clean drift into a variance collapse and sent the user
  # the wrong remedy (and suppressed the retry that would have helped).
  out <- BACE:::.conv_assess_series(seq(0.1, 0.9, length.out = 10), NULL,
                                    n_cells = 30, kind = "proportion")
  expect_identical(out$verdict, "drift")
})

test_that("regression: a settled proportion over few cells does not false-fire", {
  # sqrt(p(1-p)) is exactly 0 whenever one run's cells all land in one level,
  # which used to trip the collapse gate ~16% of the time at n_cells = 5.
  v <- vapply(1:400, function(b) {
    set.seed(b)
    BACE:::.conv_assess_series(rbinom(10, 5, 0.1) / 5, NULL,
                               n_cells = 5, kind = "proportion")$verdict
  }, character(1))
  expect_gte(mean(v == "ok"), 0.97)
  expect_identical(sum(v == "stuck"), 0L)
})

test_that("regression: huge movement cannot be vetoed by a tie-ridden trend test", {
  # net = 500,000 spreads passed because Kendall's p was 0.14 under heavy ties.
  out <- BACE:::.conv_assess_series(c(rep(0, 9), 1e6), rep(1, 10), n_cells = 20)
  expect_identical(out$verdict, "drift")
})

test_that("regression: a large wave that returns home is still 'not settled'", {
  out <- BACE:::.conv_assess_series(c(0, 0, 0, 0, 0, 3, 6, 6, 3, 0),
                                    rep(1, 10), n_cells = 20)
  expect_identical(out$verdict, "drift")
})

test_that("regression: an exploding imputation spread is flagged, not passed", {
  set.seed(6)
  out <- BACE:::.conv_assess_series(rnorm(10), 1.5^(0:9), n_cells = 30)
  expect_identical(out$verdict, "diverging")
})

test_that("regression: a drift at a tiny absolute scale is not swallowed as constant", {
  # the constant branch used an ABSOLUTE 1e-12 tolerance, so a clean
  # 0 -> 1e-13 ramp was reported as 'identical imputations in every run'.
  out <- BACE:::.conv_assess_series(seq(0, 1e-13, length.out = 10),
                                    rep(1e-14, 10), n_cells = 20)
  expect_identical(out$verdict, "drift")
})

test_that("regression: false oscillation stays rare when many cells are imputed", {
  # ref = sqrt(2/n_cells) shrinks with n_cells but real refit movement does
  # not, so without an absolute amplitude floor the gate stopped screening
  # and a 5-pair correlation decided alone (~25% false at n_cells = 1000).
  for (nc in c(150, 500, 1000)) {
    v <- vapply(1:500, function(b) {
      set.seed(b * 7 + nc)
      e <- rnorm(11, 0, 0.20); m <- e[-1] - 0.3 * e[-11]   # settled, mild alternation
      BACE:::.conv_assess_series(m, rep(1, 10), n_cells = nc)$verdict
    }, character(1))
    expect_lte(mean(v == "oscillation"), 0.05)
    expect_gte(mean(v == "ok"), 0.93)
  }
})

test_that("regression: true oscillation is still caught at large n_cells", {
  for (nc in c(150, 500)) {
    v <- vapply(1:300, function(b) {
      set.seed(b * 11 + nc)
      BACE:::.conv_assess_series(0.75 * (-1)^(1:10) + rnorm(10, 0, sqrt(2/nc)),
                                 rep(1, 10), n_cells = nc)$verdict
    }, character(1))
    expect_gte(mean(v == "oscillation"), 0.95)
  }
})

test_that("regression: 'unassessable' is not a pathology — converged is NA, not FALSE", {
  bo <- local({
    base <- data.frame(y = rnorm(20))
    base$y[1:5] <- NA
    dats <- lapply(1:10, function(r) { d <- base; d$y[1:5] <- NA_real_; d })
    o <- list(data = c(list(base), dats), types = list(y = "gaussian"),
              miss_dat = data.frame(row = 1:5, colname = "y",
                                    stringsAsFactors = FALSE))
    class(o) <- "bace"; o
  })
  out <- assess_convergence(bo)
  expect_identical(unname(out$verdicts["y"]), "unassessable")
  expect_true(is.na(out$converged))
  expect_output(print(out), "Nothing could be assessed")
})

test_that("regression: short chains SAY they could not run the full battery", {
  out <- BACE:::.conv_assess_series(c(1, 1.1, 0.9, 1.05), rep(1, 4), n_cells = 20)
  expect_match(out$note, "too short to test for flip-flopping")
  out3 <- BACE:::.conv_assess_series(c(1, 1.1, 0.9), rep(1, 3), n_cells = 20)
  expect_match(out3$note, "too short to test for drift")
})

test_that("regression: a bace object with nothing missing returns NA, not a crash", {
  obj <- list(data = list(data.frame(y = 1:5), data.frame(y = 1:5),
                          data.frame(y = 1:5), data.frame(y = 1:5)),
              types = list(y = "gaussian"),
              miss_dat = data.frame(row = integer(0), colname = character(0),
                                    stringsAsFactors = FALSE))
  class(obj) <- "bace"
  out <- assess_convergence(obj)
  expect_true(is.na(out$converged))
  expect_output(print(out), "Nothing could be assessed")
})

test_that("regression: min_iterations below 3 is refused, not run", {
  obj <- list(data = list(data.frame(y = 1)), types = list(y = "gaussian"),
              miss_dat = data.frame(row = 1, colname = "y",
                                    stringsAsFactors = FALSE))
  class(obj) <- "bace"
  expect_error(assess_convergence(obj, min_iterations = 0),
               regexp = "min_iterations must be")
})

test_that("regression: a non-finite n_cells does not crash the rule", {
  set.seed(12)
  out <- BACE:::.conv_assess_series(rnorm(10, 0, 0.2), rep(1, 10),
                                    n_cells = NA_real_)
  expect_true(out$verdict %in% c("ok", "drift", "oscillation"))
})
