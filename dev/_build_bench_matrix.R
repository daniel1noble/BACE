# =============================================================================
# _build_bench_matrix.R
#
# Emit the GitHub Actions matrix for the cross-dataset (real-data) benchmark
# as JSON, sharding each dataset into independent replicate jobs.
#
# Usage (from the workflow `prepare` job):
#   Rscript dev/_build_bench_matrix.R <mode> <n_reps>
#     mode   : "smoke" (1 rep/dataset, tiny budget via BENCH_SMOKE) or
#              "production" (n_reps per dataset at full budget)
#     n_reps : replicate mask draws per dataset in production (default 10)
#
# Writes `matrix=<json>` to $GITHUB_OUTPUT (and stdout), where <json> is
# {"include":[{dataset,script,rep}, ...]} for
# `strategy.matrix: ${{ fromJSON(needs.prepare.outputs.matrix) }}`.
#
# Each job runs its dataset wrapper with BENCH_REP_ID set to `rep`; the
# benchmark engine turns that into an independent seed (subsample + mask).
# =============================================================================

# Dataset -> wrapper script map (single source of truth for the matrix).
# fishbase is intentionally omitted until data-raw/make_fishbase.R exists.
DATASETS <- list(
  list(dataset = "avonet",    script = "dev/00_benchmark_AVONET.R"),
  list(dataset = "pantheria", script = "dev/03_benchmark_pantheria.R"),
  list(dataset = "amphibio",  script = "dev/04_benchmark_amphibio.R"),
  list(dataset = "bien",      script = "dev/05_benchmark_bien.R"),
  list(dataset = "globtherm", script = "dev/06_benchmark_globtherm.R"),
  list(dataset = "leptraits", script = "dev/07_benchmark_leptraits.R")
)

args   <- commandArgs(trailingOnly = TRUE)
mode   <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else "production"
n_reps <- if (length(args) >= 2L && nzchar(args[[2]])) as.integer(args[[2]]) else 10L
if (is.na(n_reps) || n_reps < 1L) n_reps <- 10L
if (identical(mode, "smoke")) n_reps <- 1L

entries <- character(0)
for (d in DATASETS) {
  for (rep in seq_len(n_reps)) {
    entries <- c(entries, sprintf(
      '{"dataset":"%s","script":"%s","rep":%d}', d$dataset, d$script, rep))
  }
}

json <- paste0('{"include":[', paste(entries, collapse = ","), "]}")

gho  <- Sys.getenv("GITHUB_OUTPUT", unset = "")
line <- paste0("matrix=", json)
if (nzchar(gho)) cat(line, "\n", file = gho, sep = "", append = TRUE)
cat("mode=", mode, " n_reps=", n_reps, " jobs=", length(entries), "\n", sep = "")
cat(line, "\n", sep = "")
