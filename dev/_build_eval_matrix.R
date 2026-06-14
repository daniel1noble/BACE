# =============================================================================
# _build_eval_matrix.R
#
# Emit the GitHub Actions `evaluate` matrix as JSON, sharding each reference
# dataset's replicates into small contiguous chunks so every matrix job
# finishes well under the 6h hosted-runner ceiling.
#
# Usage (from the workflow `prepare` job, AFTER datasets are generated):
#   Rscript dev/_build_eval_matrix.R <mode> <chunk_size>
#     mode       : "smoke"/"calibrate" (1 rep/dataset) or "production"
#                  (all reps, chunked)
#     chunk_size : reps per matrix job in production mode (default 4)
#
# Writes a line `matrix=<json>` to $GITHUB_OUTPUT (or stdout when run
# locally), where <json> is {"include":[{dataset,start,end}, ...]} suitable
# for `strategy.matrix: ${{ fromJSON(needs.prepare.outputs.matrix) }}`.
# =============================================================================

REF_ROOT <- file.path("dev", "simulation_results", "reference_datasets")
DATASETS <- c("sim_ideal", "sim_typical", "sim_heterogeneous", "sim_hard")

args       <- commandArgs(trailingOnly = TRUE)
mode       <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else "production"
chunk_size <- if (length(args) >= 2L && nzchar(args[[2]])) as.integer(args[[2]]) else 4L
if (is.na(chunk_size) || chunk_size < 1L) chunk_size <- 4L

avail_reps <- function(ds) {
  d <- file.path(REF_ROOT, ds)
  if (!dir.exists(d)) return(integer(0))
  reps <- list.files(d, pattern = "^rep_\\d+\\.rds$")
  sort(as.integer(sub("^rep_(\\d+)\\.rds$", "\\1", reps)))
}

# Build {dataset,start,end} chunks ------------------------------------------
entries <- list()
add_entry <- function(ds, start, end) {
  entries[[length(entries) + 1L]] <<- sprintf(
    '{"dataset":"%s","start":%d,"end":%d}', ds, start, end)
}

for (ds in DATASETS) {
  reps <- avail_reps(ds)
  if (length(reps) == 0L) next

  if (mode %in% c("smoke", "calibrate")) {
    # One rep per dataset, each its own job. smoke = tiny budget; calibrate =
    # full budget, used to measure real per-rep wall-clock and convergence.
    add_entry(ds, reps[[1]], reps[[1]])
    next
  }

  # Production: contiguous chunks over the actual rep ids present on disk.
  starts <- seq(1L, length(reps), by = chunk_size)
  for (s in starts) {
    e <- min(s + chunk_size - 1L, length(reps))
    add_entry(ds, reps[[s]], reps[[e]])
  }
}

if (length(entries) == 0L)
  stop("No reference-dataset reps found under ", REF_ROOT,
       " — did dev/09_generate_reference_datasets.R run first?")

json <- paste0('{"include":[', paste(unlist(entries), collapse = ","), "]}")

# Emit ----------------------------------------------------------------------
gho <- Sys.getenv("GITHUB_OUTPUT", unset = "")
line <- paste0("matrix=", json)
if (nzchar(gho)) {
  cat(line, "\n", file = gho, sep = "", append = TRUE)
}
cat("mode=", mode, " chunk_size=", chunk_size,
    " jobs=", length(entries), "\n", sep = "")
cat(line, "\n", sep = "")
