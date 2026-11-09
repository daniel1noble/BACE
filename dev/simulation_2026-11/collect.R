# =================================================================================================
# collect.R -- gather the per-dataset results into tables.
#
#   Rscript collect.R results/
#
# Reads every .rds run_one.R wrote and produces five CSVs. The two "per_dataset" files are the raw
# material; the three summaries are what the manuscript reads.
# =================================================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x
a <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(a)) a[1] else "results"
files <- list.files(res_dir, pattern = "[.]rds$", full.names = TRUE)
if (!length(files)) stop("no results in ", res_dir, call. = FALSE)
message("reading ", length(files), " results from ", res_dir)

r <- lapply(files, readRDS)

## STEP 0. Separate the datasets that failed. They are a RESULT -- the failure rate belongs in the
## paper -- so they are counted and reported here rather than quietly skipped.
failed <- Filter(function(x) !is.na(x$failed_at %||% NA), r)
r      <- Filter(function(x)  is.na(x$failed_at %||% NA), r)
if (length(failed)) {
  ft <- table(vapply(failed, function(x) paste0(x$failed_at, ": ", x$error), character(1)))
  message("\n", length(failed), " of ", length(files), " datasets failed (",
          round(100 * length(failed) / length(files), 1), "%):")
  for (i in seq_along(ft)) message("  ", ft[i], "x  ", names(ft)[i])
  utils::write.csv(do.call(rbind, lapply(failed, function(x)
    data.frame(id = x$id, seed = x$seed, failed_at = x$failed_at, error = x$error))),
    file.path(res_dir, "failures.csv"), row.names = FALSE)
}
if (!length(r)) stop("every dataset failed", call. = FALSE)

## STEP 1. The scenario each dataset belongs to, so every table can be grouped by it.
key <- function(x) data.frame(id = x$id, seed = x$seed, n_species = x$scenario$n_species,
                              phylo_signal = x$scenario$phylo_signal, regime = x$scenario$regime,
                              mechanism = x$scenario$mechanism)

## STEP 2. One row per dataset per arm: the slope, and whether its interval caught the truth.
slopes <- do.call(rbind, lapply(r, function(x) if (is.null(x$slopes)) NULL else
                                cbind(key(x), x$slopes)))

## STEP 3. One row per dataset per trait: how well the hidden values came back.
recovery <- do.call(rbind, lapply(r, function(x) if (is.null(x$recovery)) NULL else
                                  cbind(key(x), x$recovery)))

## STEP 4. One row per dataset: everything else worth knowing about how the run went.
runs <- do.call(rbind, lapply(r, function(x) cbind(key(x),
  data.frame(realised_frac  = x$realised_frac,
             lambda_x1      = x$realised$lambda_x1,
             lambda_bin     = x$realised$lambda_bin,
             cor_x1_x2      = x$realised$cor_x1_x2,
             cor_x1_aux     = x$realised$cor_x1_aux,
             cc_n_complete  = x$complete_case$n_complete,
             cc_frac_kept   = x$complete_case$frac_kept,
             cc_status      = x$complete_case$status,
             bace_converged = x$bace_converged,
             wall_s         = x$wall_s,
             bace_version   = x$bace_version))))

## STEP 5. Summaries, averaged over replicates within each scenario. `n` is how many replicates
##         contributed -- always report it, because a scenario where the complete-case arm often
##         declined has fewer, and that is itself a result.
by <- c("n_species", "phylo_signal", "regime", "mechanism")

summary_slopes <- stats::aggregate(
  cbind(bias, covered, width = upper - lower, n = 1) ~ ., data = slopes[, c(by, "arm", "bias", "covered", "lower", "upper")],
  FUN = mean, na.rm = TRUE)

summary_recovery <- stats::aggregate(
  cbind(zrmse, accuracy) ~ ., data = recovery[, c(by, "trait", "kind", "zrmse", "accuracy")],
  FUN = function(z) mean(z, na.rm = TRUE), na.action = stats::na.pass)

summary_runs <- stats::aggregate(
  cbind(bace_converged, cc_frac_kept, lambda_x1, cor_x1_x2, wall_s) ~ ., 
  data = runs[, c(by, "bace_converged", "cc_frac_kept", "lambda_x1", "cor_x1_x2", "wall_s")],
  FUN = mean, na.rm = TRUE)

## STEP 6. Write them.
out <- function(x, name) { f <- file.path(res_dir, name); utils::write.csv(x, f, row.names = FALSE)
                           message("  ", name, "  (", nrow(x), " rows)") }
out(slopes,           "per_dataset_slopes.csv")
out(recovery,         "per_dataset_recovery.csv")
out(summary_slopes,   "summary_slopes.csv")
out(summary_recovery, "summary_recovery.csv")
out(summary_runs,     "summary_runs.csv")

## STEP 7. The one check worth doing every time: the ORACLE arm sees the complete truth, so its
##         coverage must be about 0.95. If it is not, the assumed true slope is wrong and every
##         other arm's coverage is meaningless.
orc <- summary_slopes[summary_slopes$arm == "oracle", ]
message("\noracle coverage (should be ~0.95): ",
        paste(round(range(orc$covered), 3), collapse = " to "))
if (any(orc$covered < 0.90 | orc$covered > 0.99))
  warning("oracle coverage is off 0.95 in some scenarios -- check the assumed true slope before ",
          "reading anything into the other arms.", call. = FALSE)
