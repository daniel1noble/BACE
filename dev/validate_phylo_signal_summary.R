# =============================================================================
# Full Monte Carlo validation for phylo_signal_summary()
# =============================================================================
#
# Purpose
# -------
# The fast smoke test in tests/testthat/test-phylo_signal_summary.R runs a
# single replicate per cell at realistic nitt to catch structural / formula
# bugs. This script runs the *full* Monte Carlo grid to assess calibration:
# posterior-mean bias, 95% HPD coverage, and RMSE across replicates.
#
# Design (matches the approved plan)
# ----------------------------------
#   - types:   gaussian, binary, categorical-3, poisson
#   - RE:      single (species=FALSE), dual (species=TRUE, 3 reps/sp)
#   - signals: {0.2, 0.7}
#   - B:       20 Monte Carlo replicates per cell
#   - n_species = 80, n_cases = 240 (3 reps/sp for dual, 1 rep/sp for single)
#   - MCMC: function's type-specific defaults (quick = FALSE)
#
# Pass criteria per cell:
#   - |mean(posterior mean - truth)|   <= 0.05
#   - 95% HPD coverage                  >= 0.90
#   - RMSE across replicates            <= 0.15
#
# Runtime: ~4-6 hours on a modern laptop. Meant for manual invocation:
#   Rscript dev/validate_phylo_signal_summary.R
#
# Output:
#   dev/benchmark_results/phylo_signal_validation/run_<tag>/summary.csv
#   dev/benchmark_results/phylo_signal_validation/run_<tag>/validation_report.md
# =============================================================================

suppressMessages({
  devtools::load_all(quiet = TRUE)
  library(MASS)
  library(ape)
})

# ---- Configuration ----
B           <- 20                 # Monte Carlo replicates per cell
N_SPECIES   <- 80
REPS_DUAL   <- 3
SIGNALS     <- c(0.2, 0.7)
R_SP_FOR_DUAL <- 0.2
TAG         <- format(Sys.time(), "%Y%m%d_%H%M")
OUT_DIR     <- file.path("dev", "benchmark_results",
                          "phylo_signal_validation", paste0("run_", TAG))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOG         <- file(file.path(OUT_DIR, "validation.log"), open = "wt")
on.exit(close(LOG), add = TRUE)

.log <- function(...) {
  msg <- paste0(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
  cat(msg, "\n")
  writeLines(msg, LOG); flush(LOG)
}

.log("Starting phylo_signal_summary validation: tag=", TAG, " B=", B)
.log("n_species=", N_SPECIES, "  signals=", paste(SIGNALS, collapse=","))

# ---- Data generators ----
.mk_tree <- function(n, seed) {
  set.seed(seed)
  tr <- ape::compute.brlen(ape::rtree(n), method = "Grafen")
  tr$tip.label <- paste0("sp", seq_along(tr$tip.label))
  tr
}

.sim_gaussian <- function(tr, H2, dual, R_sp = 0, reps = 1, seed) {
  set.seed(seed)
  n <- length(tr$tip.label); A <- ape::vcv.phylo(tr)
  V_A <- H2; V_S <- if (dual) R_sp else 0
  V_R <- 1 - V_A - V_S
  bm <- as.numeric(MASS::mvrnorm(1, rep(0, n), V_A * A))
  names(bm) <- tr$tip.label
  sp <- if (dual) stats::setNames(rnorm(n, 0, sqrt(V_S)), tr$tip.label) else NULL
  do.call(rbind, lapply(tr$tip.label, function(s) {
    data.frame(Species = s,
               y = bm[[s]] + (if (dual) sp[[s]] else 0) +
                   rnorm(reps, 0, sqrt(V_R)))
  }))
}

.sim_binary <- function(tr, H2, dual, R_sp = 0, reps = 1, seed) {
  set.seed(seed)
  n <- length(tr$tip.label); A <- ape::vcv.phylo(tr)
  V_A <- H2; V_S <- if (dual) R_sp else 0
  V_R <- 1 - V_A - V_S  # residual on latent scale (probit "1")
  bm <- as.numeric(MASS::mvrnorm(1, rep(0, n), V_A * A))
  names(bm) <- tr$tip.label
  sp <- if (dual) stats::setNames(rnorm(n, 0, sqrt(V_S)), tr$tip.label) else NULL
  do.call(rbind, lapply(tr$tip.label, function(s) {
    lat <- bm[[s]] + (if (dual) sp[[s]] else 0) + rnorm(reps, 0, sqrt(V_R))
    data.frame(Species = s,
               y = factor(as.integer(lat > 0), levels = c(0L, 1L)))
  }))
}

.sim_categorical3 <- function(tr, H2, dual, R_sp = 0, reps = 1, seed) {
  set.seed(seed)
  n <- length(tr$tip.label); A <- ape::vcv.phylo(tr)
  V_A <- H2; V_S <- if (dual) R_sp else 0
  V_R <- 1 - V_A - V_S
  # Two latent dims -> argmax over {lat1, lat2, 0} -> 3-level factor
  lat_mat <- sapply(1:2, function(k) {
    bm <- as.numeric(MASS::mvrnorm(1, rep(0, n), V_A * A))
    names(bm) <- tr$tip.label
    sp <- if (dual) stats::setNames(rnorm(n, 0, sqrt(V_S)), tr$tip.label) else NULL
    unlist(lapply(tr$tip.label, function(s) {
      bm[[s]] + (if (dual) sp[[s]] else 0) + rnorm(reps, 0, sqrt(V_R))
    }))
  })
  species_col <- rep(tr$tip.label, each = reps)
  # Append zero column for the reference level
  full <- cbind(lat_mat, 0)
  y <- factor(letters[apply(full, 1, which.max)], levels = c("a","b","c"))
  data.frame(Species = species_col, y = y, stringsAsFactors = FALSE)
}

.sim_poisson <- function(tr, H2, dual, R_sp = 0, reps = 1, seed) {
  set.seed(seed)
  n <- length(tr$tip.label); A <- ape::vcv.phylo(tr)
  V_A <- H2; V_S <- if (dual) R_sp else 0
  V_R <- 1 - V_A - V_S
  bm <- as.numeric(MASS::mvrnorm(1, rep(0, n), V_A * A))
  names(bm) <- tr$tip.label
  sp <- if (dual) stats::setNames(rnorm(n, 0, sqrt(V_S)), tr$tip.label) else NULL
  do.call(rbind, lapply(tr$tip.label, function(s) {
    lin <- 1 + bm[[s]] + (if (dual) sp[[s]] else 0) + rnorm(reps, 0, sqrt(V_R))
    data.frame(Species = s, y = rpois(reps, exp(lin)))
  }))
}

.SIM_FUN <- list(
  gaussian       = .sim_gaussian,
  binary         = .sim_binary,
  "categorical-3"= .sim_categorical3,
  poisson        = .sim_poisson
)

# ---- Monte Carlo loop ----
grid <- expand.grid(
  type    = names(.SIM_FUN),
  re      = c("single", "dual"),
  signal  = SIGNALS,
  stringsAsFactors = FALSE
)

all_rows <- list()
cell_id  <- 0L

for (cell in seq_len(nrow(grid))) {
  type_c  <- grid$type[cell]
  re_c    <- grid$re[cell]
  sig_c   <- grid$signal[cell]
  dual_c  <- re_c == "dual"
  reps_c  <- if (dual_c) REPS_DUAL else 1L

  .log(sprintf("Cell %d/%d: type=%s re=%s signal=%.1f",
               cell, nrow(grid), type_c, re_c, sig_c))

  for (b in seq_len(B)) {
    cell_id <- cell_id + 1L
    seed_b  <- 1000L + cell_id
    tr_b    <- .mk_tree(N_SPECIES, seed = seed_b)
    dat_b   <- .SIM_FUN[[type_c]](tr_b, H2 = sig_c, dual = dual_c,
                                  R_sp = R_SP_FOR_DUAL, reps = reps_c,
                                  seed = seed_b + 1L)

    res <- tryCatch(
      suppressWarnings(phylo_signal_summary(
        data = dat_b, tree = tr_b, species_col = "Species",
        variables = "y", species = dual_c, verbose = FALSE
      )),
      error = function(e) { .log("  err rep=", b, ": ", conditionMessage(e)); NULL }
    )
    if (is.null(res)) next

    row <- res$table
    all_rows[[length(all_rows) + 1L]] <- data.frame(
      cell_id     = cell_id,
      type        = type_c,
      re          = re_c,
      signal      = sig_c,
      rep         = b,
      H2_mean     = row$H2_mean,
      H2_lo       = row$H2_lo,
      H2_hi       = row$H2_hi,
      min_ess     = row$min_ess,
      flag        = row$flag,
      stringsAsFactors = FALSE
    )
  }
}

raw <- do.call(rbind, all_rows)
write.csv(raw, file.path(OUT_DIR, "raw_replicates.csv"), row.names = FALSE)
.log("Wrote raw_replicates.csv (", nrow(raw), " rows)")

# ---- Aggregate per cell ----
agg <- do.call(rbind, lapply(split(raw, list(raw$type, raw$re, raw$signal),
                                    drop = TRUE), function(g) {
  truth <- g$signal[1]
  bias  <- mean(g$H2_mean, na.rm = TRUE) - truth
  rmse  <- sqrt(mean((g$H2_mean - truth)^2, na.rm = TRUE))
  cov95 <- mean(g$H2_lo <= truth & g$H2_hi >= truth, na.rm = TRUE)
  data.frame(
    type   = g$type[1],
    re     = g$re[1],
    signal = truth,
    B_fit  = nrow(g),
    bias   = bias,
    rmse   = rmse,
    cov95  = cov95,
    pass_bias  = abs(bias) <= 0.05,
    pass_rmse  = rmse <= 0.15,
    pass_cov95 = cov95 >= 0.90,
    stringsAsFactors = FALSE
  )
}))
write.csv(agg, file.path(OUT_DIR, "summary.csv"), row.names = FALSE)

# ---- Human-readable report ----
rpt <- file.path(OUT_DIR, "validation_report.md")
con <- file(rpt, open = "wt")
writeLines(c(
  "# phylo_signal_summary() validation report",
  paste0("Tag: ", TAG),
  paste0("B = ", B, "  n_species = ", N_SPECIES),
  "",
  "## Pass criteria",
  "- |bias| <= 0.05, RMSE <= 0.15, 95% HPD coverage >= 0.90",
  "",
  "## Per-cell summary",
  "",
  paste(c(names(agg), ""), collapse = "|"),
  paste(c(rep("---", ncol(agg)), ""), collapse = "|")
), con)
for (r in seq_len(nrow(agg))) {
  writeLines(paste(c(format(unlist(agg[r, ])), ""), collapse = "|"), con)
}
writeLines("", con)
failed <- !(agg$pass_bias & agg$pass_rmse & agg$pass_cov95)
writeLines(sprintf("## Cells failing any criterion: %d / %d",
                   sum(failed), nrow(agg)), con)
close(con)
.log("Wrote ", rpt)
.log("Validation complete. Failures: ", sum(failed), "/", nrow(agg))

# Signal non-zero exit if any cell failed
if (any(failed)) {
  quit(status = 1)
}
