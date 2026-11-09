# =============================================================================
# Compare BACE benchmark runs across versions
# =============================================================================
#
# Purpose
# -------
# `dev/avonet_benchmark.R` writes tagged output to
# `dev/benchmark_results/avonet/run_<tag>/` on every run. This script
# scans those tagged runs and produces a side-by-side comparison table
# (per-trait metrics, deltas vs a chosen baseline), so you can see at a
# glance whether a BACE change improved or regressed imputation quality.
#
# Usage
# -----
#   # Compare all AVONET runs (no baseline):
#   Rscript dev/compare_benchmark_runs.R
#
#   # Compare all runs against a specific baseline tag:
#   Rscript dev/compare_benchmark_runs.R --baseline=<run_tag>
#
# Output
# ------
#   dev/benchmark_results/avonet/comparison.csv   wide table of metrics
#   dev/benchmark_results/avonet/comparison.txt   printable summary
# =============================================================================

BENCH_DIR <- file.path("dev", "benchmark_results", "avonet")

# ---- Parse args -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
baseline_tag <- NA_character_
for (a in args) {
  if (startsWith(a, "--baseline=")) {
    baseline_tag <- sub("^--baseline=", "", a)
  }
}

# ---- Discover tagged runs ---------------------------------------------------
run_dirs <- list.files(BENCH_DIR, pattern = "^run_", full.names = TRUE)
run_dirs <- run_dirs[file.info(run_dirs)$isdir]
if (length(run_dirs) == 0) {
  stop("No tagged runs found in ", BENCH_DIR,
       ". Run dev/avonet_benchmark.R first.")
}

# ---- Load summary per run ---------------------------------------------------
loaded <- lapply(run_dirs, function(d) {
  tag <- sub("^run_", "", basename(d))
  summary_path  <- file.path(d, "summary.csv")
  metadata_path <- file.path(d, "metadata.json")
  if (!file.exists(summary_path)) {
    warning("Missing summary.csv in ", d); return(NULL)
  }
  summary_df <- read.csv(summary_path, stringsAsFactors = FALSE)
  summary_df$run_tag <- tag

  meta <- if (file.exists(metadata_path) &&
              requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::fromJSON(metadata_path)
  } else list()
  list(summary = summary_df, metadata = meta)
})
loaded <- loaded[!vapply(loaded, is.null, logical(1))]
if (length(loaded) == 0) stop("No valid run summaries found.")

all_summary <- do.call(rbind, lapply(loaded, `[[`, "summary"))

# ---- Wide comparison: metrics x (trait, run_tag) ----------------------------
metrics <- c("nrmse", "mae_fit", "mae_raw", "correlation", "coverage95",
             "accuracy", "balanced_accuracy", "brier",
             "lambda", "K", "D")

cat("===========================================================\n")
cat("  AVONET benchmark comparison\n")
cat("  Runs found :", length(loaded), "\n")
for (l in loaded) {
  cat("    ", l$summary$run_tag[1],
      "  (converged =",
      if (!is.null(l$metadata$converged)) l$metadata$converged else "?",
      if (!is.null(l$metadata$runtime_min))
        paste0(" runtime = ", l$metadata$runtime_min, "m") else "",
      ")\n")
}
if (!is.na(baseline_tag)) {
  if (!baseline_tag %in% all_summary$run_tag) {
    stop("Baseline tag not found: ", baseline_tag,
         "\nAvailable: ", paste(unique(all_summary$run_tag), collapse = ", "))
  }
  cat("  Baseline   :", baseline_tag, "\n")
}
cat("===========================================================\n\n")

# ---- Build per-metric tables ------------------------------------------------
comparison_tables <- lapply(metrics, function(m) {
  if (!(m %in% names(all_summary))) return(NULL)   # metric not in these runs
  tab <- reshape(
    all_summary[, c("trait", "run_tag", m)],
    idvar     = "trait", timevar = "run_tag", direction = "wide")
  names(tab) <- sub(paste0("^", m, "\\."), "", names(tab))
  tab$metric <- m
  tab[, c("metric", "trait",
          setdiff(names(tab), c("metric", "trait")))]
})
comparison_tables <- comparison_tables[
  !vapply(comparison_tables, is.null, logical(1))]
wide <- do.call(rbind, comparison_tables)

# Deltas vs baseline
if (!is.na(baseline_tag)) {
  other_runs <- setdiff(unique(all_summary$run_tag), baseline_tag)
  for (r in other_runs) {
    wide[[paste0("delta.", r)]] <- wide[[r]] - wide[[baseline_tag]]
  }
}

write.csv(wide, file.path(BENCH_DIR, "comparison.csv"), row.names = FALSE)

# ---- Printable summary ------------------------------------------------------
sink_file <- file.path(BENCH_DIR, "comparison.txt")
sink(sink_file, split = TRUE)
cat("AVONET benchmark comparison\n")
cat("Generated:", format(Sys.time()), "\n\n")
for (m in metrics) {
  sub <- wide[wide$metric == m, ]
  if (nrow(sub) == 0) next
  cat("\n--- ", m,
      if (m %in% c("nrmse", "mae_fit", "mae_raw", "brier")) " (lower is better)" else
      if (m == "coverage95") " (target = 0.95)" else
      if (m %in% c("lambda", "K")) " (phylo signal; 0-1; BM-like ~ 1)" else
      if (m == "D") " (phylo signal; 0 = BM, 1 = random, <0 = conserved)" else
      " (higher is better)", " ---\n", sep = "")
  print(sub[, setdiff(names(sub), "metric")], row.names = FALSE, digits = 3)
}
sink()

cat("\nWrote:\n")
cat("  ", file.path(BENCH_DIR, "comparison.csv"), "\n")
cat("  ", sink_file, "\n")
