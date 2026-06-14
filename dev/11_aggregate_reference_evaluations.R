# =============================================================================
# 11_aggregate_reference_evaluations.R
#
# Walk dev/simulation_results/evaluation_results/{dataset}/eval_rep_NN.rds,
# pull out per-rep metrics, and emit:
#
#   - per_rep_long.csv      one row per (dataset, rep, variable, metric)
#   - per_rep_long.rds      same in RDS
#   - beta_compare_long.csv one row per (dataset, rep, coefficient) with
#                           BACE vs oracle vs dialled comparators
#   - timings.csv           per-rep bace_runtime / oracle_runtime
#   - summary_report.md     headline numbers aggregated across reps
#
# All paths are relative to the repo root; runnable both locally and from
# the GitHub Actions aggregate job (where artifacts land under
# evaluation_results/eval-*/ before being renamed in-place).
# =============================================================================

suppressPackageStartupMessages({
  library(stats)
  library(utils)
})

EVAL_ROOT <- file.path("dev", "simulation_results", "evaluation_results")
DATASETS  <- c("sim_ideal", "sim_typical", "sim_heterogeneous", "sim_hard")

# -----------------------------------------------------------------------------
# Load every per-rep RDS bundle we can find.
# -----------------------------------------------------------------------------
load_evals <- function() {
  files <- character(0)
  for (ds in DATASETS) {
    d <- file.path(EVAL_ROOT, ds)
    if (!dir.exists(d)) next
    files <- c(files, list.files(d, pattern = "^eval_rep_\\d+\\.rds$",
                                 full.names = TRUE))
  }
  cat("Found", length(files), "eval bundles\n")
  lapply(files, function(f) {
    out <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(out)) {
      cat("  [skip] unreadable:", f, "\n")
      return(NULL)
    }
    out$.file <- f
    out
  }) |> Filter(f = function(x) !is.null(x))
}

evals <- load_evals()
if (length(evals) == 0L) {
  message("No evaluation bundles found under ", EVAL_ROOT)
  quit(status = 0)
}

# -----------------------------------------------------------------------------
# Per-rep long form: cell metrics
# -----------------------------------------------------------------------------
cell_rows <- list()
for (ev in evals) {
  if (!identical(ev$status, "ok")) next
  if (is.null(ev$cell_metrics) || nrow(ev$cell_metrics) == 0L) next
  cm <- ev$cell_metrics
  cm$dataset   <- ev$dataset
  cm$rep_id    <- ev$rep_id
  cm$mechanism <- ev$mechanism
  cm$rate      <- ev$rate
  cell_rows[[length(cell_rows) + 1]] <- cm
}
cell_long <- if (length(cell_rows) > 0L) do.call(rbind, cell_rows) else NULL

# -----------------------------------------------------------------------------
# Per-rep long form: beta comparisons (BACE vs oracle vs dialled)
# -----------------------------------------------------------------------------
beta_rows <- list()
for (ev in evals) {
  if (!identical(ev$status, "ok")) next
  if (is.null(ev$beta_compare) || nrow(ev$beta_compare) == 0L) next
  bc <- ev$beta_compare
  bc$dataset   <- ev$dataset
  bc$rep_id    <- ev$rep_id
  bc$mechanism <- ev$mechanism
  bc$rate      <- ev$rate
  beta_rows[[length(beta_rows) + 1]] <- bc
}
beta_long <- if (length(beta_rows) > 0L) do.call(rbind, beta_rows) else NULL

# -----------------------------------------------------------------------------
# Per-rep long form: timings + status
# -----------------------------------------------------------------------------
timing_rows <- lapply(evals, function(ev) {
  data.frame(
    dataset        = ev$dataset,
    rep_id         = ev$rep_id,
    status         = ev$status,
    bace_converged = if (!is.null(ev$bace_converged)) ev$bace_converged else NA,
    bace_attempts  = if (!is.null(ev$bace_attempts))  ev$bace_attempts  else NA_integer_,
    bace_runtime_s = if (!is.null(ev$bace_runtime))   ev$bace_runtime   else NA_real_,
    oracle_runtime_s = if (!is.null(ev$oracle_runtime)) ev$oracle_runtime else NA_real_,
    error          = if (!is.null(ev$error)) ev$error else NA_character_,
    stringsAsFactors = FALSE
  )
})
timings <- do.call(rbind, timing_rows)

# -----------------------------------------------------------------------------
# Write tidy outputs
# -----------------------------------------------------------------------------
if (!is.null(cell_long))
  utils::write.csv(cell_long, file.path(EVAL_ROOT, "per_rep_long.csv"),
                   row.names = FALSE)
if (!is.null(beta_long))
  utils::write.csv(beta_long, file.path(EVAL_ROOT, "beta_compare_long.csv"),
                   row.names = FALSE)
utils::write.csv(timings, file.path(EVAL_ROOT, "timings.csv"),
                 row.names = FALSE)
saveRDS(list(cell = cell_long, beta = beta_long, timings = timings),
        file.path(EVAL_ROOT, "per_rep_long.rds"))

# -----------------------------------------------------------------------------
# Markdown summary report (headline numbers per dataset, suitable for
# pasting into $GITHUB_STEP_SUMMARY).
# -----------------------------------------------------------------------------
report <- c(
  "# Simulated reference benchmark — summary",
  "",
  sprintf("Generated at %s.", format(Sys.time())),
  sprintf("%d successful reps / %d total bundles loaded.",
          sum(timings$status == "ok"), nrow(timings)),
  ""
)

# Status table
report <- c(report,
  "## Replicate status",
  "",
  "| Dataset | OK | Error | Timeout | Other |",
  "|:--------|---:|------:|--------:|------:|"
)
stat_tab <- table(factor(timings$dataset, levels = DATASETS), timings$status)
col_count <- function(tab, ds, nm)
  if (nm %in% colnames(tab)) tab[ds, nm] else 0L
for (ds in DATASETS) {
  ok    <- col_count(stat_tab, ds, "ok")
  err   <- col_count(stat_tab, ds, "error")
  tmo   <- col_count(stat_tab, ds, "timeout")
  total <- sum(stat_tab[ds, ])
  oth   <- total - ok - err - tmo
  report <- c(report, sprintf("| %s | %d | %d | %d | %d |",
                              ds, ok, err, tmo, oth))
}
report <- c(report, "")

# Convergence
report <- c(report,
  "## BACE convergence rate (among status = ok)",
  "",
  "| Dataset | reps | converged | rate |",
  "|:--------|-----:|----------:|-----:|"
)
ok_t <- subset(timings, status == "ok")
for (ds in DATASETS) {
  sub <- subset(ok_t, dataset == ds)
  if (nrow(sub) == 0L) next
  conv_rate <- mean(sub$bace_converged, na.rm = TRUE)
  report <- c(report, sprintf("| %s | %d | %d | %.2f |",
                              ds, nrow(sub), sum(sub$bace_converged, na.rm = TRUE),
                              conv_rate))
}
report <- c(report, "")

# Timings
report <- c(report,
  "## Per-rep wallclock (median seconds)",
  "",
  "| Dataset | BACE | Oracle |",
  "|:--------|-----:|-------:|"
)
for (ds in DATASETS) {
  sub <- subset(ok_t, dataset == ds)
  if (nrow(sub) == 0L) next
  bm <- median(sub$bace_runtime_s, na.rm = TRUE)
  om <- median(sub$oracle_runtime_s, na.rm = TRUE)
  report <- c(report, sprintf("| %s | %.0f | %.0f |", ds, bm, om))
}
report <- c(report, "")

# Headline cell metrics (mean across reps)
if (!is.null(cell_long)) {
  report <- c(report,
    "## Cell-level imputation accuracy (mean across reps)",
    ""
  )
  for (ds in DATASETS) {
    sub <- subset(cell_long, dataset == ds)
    if (nrow(sub) == 0L) next
    agg <- aggregate(value ~ variable + metric, data = sub, FUN = mean)
    report <- c(report, sprintf("### %s (%d reps)", ds,
                                length(unique(sub$rep_id))), "")
    metrics_of_interest <- c("nrmse", "mae", "accuracy", "ordinal_mae",
                             "coverage95", "brier")
    agg <- subset(agg, metric %in% metrics_of_interest)
    agg <- agg[order(agg$variable, agg$metric), ]
    report <- c(report,
                "| variable | metric | mean |",
                "|:---------|:-------|-----:|")
    for (i in seq_len(nrow(agg))) {
      report <- c(report, sprintf("| %s | %s | %.3f |",
                                   agg$variable[i], agg$metric[i],
                                   agg$value[i]))
    }
    report <- c(report, "")
  }
}

# Headline beta-coverage metrics
if (!is.null(beta_long)) {
  report <- c(report,
    "## Beta coverage on the response equation (y ~ x1 + x2 + x3 + x4)",
    "",
    "Coverage of the **oracle** posterior mean (column `cover_oracle`) is",
    "the standard MI-calibration target; should approach 0.95 with good",
    "calibration. Coverage of the **dialled-in** true coefficient",
    "(column `cover_dialled`) is reported where the coding maps cleanly",
    "(intercept, x1, x2.B, x2.C, x3); x4 is omitted because of the",
    "ordered-factor poly-contrast mismatch.",
    ""
  )
  for (ds in DATASETS) {
    sub <- subset(beta_long, dataset == ds)
    if (nrow(sub) == 0L) next
    agg <- aggregate(cbind(cover_oracle, cover_dialled,
                           bias_vs_oracle, bias_vs_dialled,
                           ci_width_ratio) ~ coef,
                     data = sub, FUN = function(z) mean(z, na.rm = TRUE))
    report <- c(report, sprintf("### %s", ds), "",
                "| coef | cover_oracle | cover_dialled | bias_oracle | bias_dialled | CI_width_ratio |",
                "|:-----|-----:|-----:|-----:|-----:|-----:|")
    for (i in seq_len(nrow(agg))) {
      report <- c(report, sprintf("| %s | %.2f | %s | %+.3f | %s | %.2f |",
        agg$coef[i],
        agg$cover_oracle[i],
        if (is.nan(agg$cover_dialled[i])) "—" else sprintf("%.2f", agg$cover_dialled[i]),
        agg$bias_vs_oracle[i],
        if (is.nan(agg$bias_vs_dialled[i])) "—" else sprintf("%+.3f", agg$bias_vs_dialled[i]),
        agg$ci_width_ratio[i]
      ))
    }
    report <- c(report, "")
  }
}

writeLines(report, file.path(EVAL_ROOT, "summary_report.md"))
cat("\nAggregate done. Wrote per_rep_long.csv, beta_compare_long.csv,",
    "timings.csv, summary_report.md to:\n  ", EVAL_ROOT, "\n")
