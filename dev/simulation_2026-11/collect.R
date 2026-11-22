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

## Files written by an older run_one.R can lack a column that a newer one writes (the scoring
## gained two log-scale columns on 2026-11-16, mid-campaign). Base rbind() refuses to stack data
## frames with different columns, so this fills a missing column with NA before stacking.
stack_rows <- function(pieces) {
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0) return(NULL)
  all_columns <- unique(unlist(lapply(pieces, names)))
  for (i in seq_along(pieces)) {
    for (column in setdiff(all_columns, names(pieces[[i]]))) pieces[[i]][[column]] <- NA
    pieces[[i]] <- pieces[[i]][, all_columns, drop = FALSE]
  }
  do.call(rbind, pieces)
}

## STEP 1. The scenario each dataset belongs to, so every table can be grouped by it.
key <- function(x) data.frame(id = x$id, seed = x$seed, n_species = x$scenario$n_species,
                              phylo_signal = x$scenario$phylo_signal, regime = x$scenario$regime,
                              mechanism = x$scenario$mechanism)

## STEP 2. One row per dataset per arm: the slope, and whether its interval caught the truth.
slopes <- do.call(rbind, lapply(r, function(x) if (is.null(x$slopes)) NULL else
                                cbind(key(x), x$slopes)))

## STEP 3. One row per dataset per trait: how well the hidden values came back.
recovery <- stack_rows(lapply(r, function(x) if (is.null(x$recovery)) NULL else
                                cbind(key(x), x$recovery)))

## STEP 3b. One row per dataset per numeric trait: did the 95% interval for each hidden value
##          contain the truth, and how wide was it (BACE's get_imputation_intervals()). Shinichi
##          read this from the .rds files by hand for the first campaign; now it is a table.
coverage <- stack_rows(lapply(r, function(x) if (is.null(x$cell_coverage)) NULL else
                                cbind(key(x), x$cell_coverage)))
## and the model-based interval (BACE >= 1.0.6), written to its own tables when present
coverage_predictive <- stack_rows(lapply(r, function(x) if (is.null(x$cell_coverage_predictive)) NULL else
                                           cbind(key(x), x$cell_coverage_predictive)))

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
             ## BACE 1.0.1 names what the stability check found, per variable; a dataset
             ## with any pathology records which. "" = every variable ok.
             bace_pathology = paste(unique(setdiff(
               unlist(x$bace_conv_detail$verdicts), "ok")), collapse = "+"),
             wall_s         = x$wall_s,
             ## how often BACE's poisson rate cap truncated an imputed rate (see 03_fit.R):
             ## total rates capped in this dataset, and which variables it hit. 0 and "" = never.
             n_rates_capped = sum(vapply(x$rate_cap %||% list(), function(e) e$n_rates, integer(1))),
             capped_vars    = paste(names(x$rate_cap %||% list()), collapse = "+"),
             bace_version   = x$bace_version))))

## STEP 5. Summaries, averaged over replicates within each scenario. `n` is how many replicates
##         contributed -- always report it, because a scenario where the complete-case arm often
##         declined has fewer, and that is itself a result.
by <- c("n_species", "phylo_signal", "regime", "mechanism")

summary_slopes <- stats::aggregate(
  cbind(bias, covered, width = upper - lower) ~ ., data = slopes[, c(by, "arm", "bias", "covered", "lower", "upper")],
  FUN = mean, na.rm = TRUE)
## `n` = how many replicates contributed to each row. (Until 2026-11-22 this column averaged a
## column of ones, so it read 1 everywhere; Shinichi spotted it.)
replicates <- stats::aggregate(cbind(n = seed) ~ ., data = slopes[, c(by, "arm", "seed")], FUN = length)
summary_slopes <- merge(summary_slopes, replicates, by = c(by, "arm"))

summary_recovery <- stats::aggregate(
  cbind(zrmse, zrmse_log, bias_log, zrmse_logdraws, bias_logdraws, accuracy) ~ .,
  data = recovery[, c(by, "trait", "kind", "zrmse", "zrmse_log", "bias_log",
                      "zrmse_logdraws", "bias_logdraws", "accuracy")],
  FUN = function(z) mean(z, na.rm = TRUE), na.action = stats::na.pass)

summary_coverage <- if (is.null(coverage)) NULL else {
  cov_mean <- stats::aggregate(cbind(coverage, width) ~ ., data = coverage[, c(by, "variable", "coverage", "width")],
                               FUN = mean, na.rm = TRUE)
  cov_n <- stats::aggregate(cbind(n = seed) ~ ., data = coverage[, c(by, "variable", "seed")], FUN = length)
  merge(cov_mean, cov_n, by = c(by, "variable"))
}

summary_runs <- stats::aggregate(
  cbind(bace_converged, drift = bace_pathology == "drift",
        cc_frac_kept, lambda_x1, cor_x1_x2, wall_s, n_rates_capped,
        any_rate_capped = n_rates_capped > 0) ~ .,
  data = cbind(runs[, c(by, "bace_converged", "cc_frac_kept", "lambda_x1",
                        "cor_x1_x2", "wall_s", "n_rates_capped")],
               bace_pathology = runs$bace_pathology),
  FUN = mean, na.rm = TRUE)

## STEP 6. Write them.
out <- function(x, name) { f <- file.path(res_dir, name); utils::write.csv(x, f, row.names = FALSE)
                           message("  ", name, "  (", nrow(x), " rows)") }
out(slopes,           "per_dataset_slopes.csv")
out(recovery,         "per_dataset_recovery.csv")
out(summary_slopes,   "summary_slopes.csv")
out(summary_recovery, "summary_recovery.csv")
if (!is.null(coverage)) { out(coverage, "per_dataset_coverage.csv"); out(summary_coverage, "summary_coverage.csv") }
if (!is.null(coverage_predictive)) {
  cp_mean <- stats::aggregate(cbind(coverage, width) ~ ., data = coverage_predictive[, c(by, "variable", "coverage", "width")],
                              FUN = mean, na.rm = TRUE)
  out(coverage_predictive, "per_dataset_coverage_predictive.csv"); out(cp_mean, "summary_coverage_predictive.csv")
}
out(summary_runs,     "summary_runs.csv")

## STEP 7. The one check worth doing every time: the ORACLE arm sees the complete truth, so its
##         coverage must be about 0.95. If it is not, the assumed true slope is wrong and every
##         other arm's coverage is meaningless.
##         COUNT EACH DATASET ONCE. The truth data do not depend on the missingness mechanism, so
##         the oracle fit for a given seed is the same in the mcar, mar and clade files; and the
##         floor and typical regimes share their random draws (the typical slope is the floor
##         slope plus 0.5, exactly), so those are one fit too. Judged on the files, one unlucky
##         seed is counted up to 18 times. On 2026-11-15 the first 301 files read 0.89 this way
##         while the same seeds, counted once, were ordinary and 150 seeds gave 0.94-0.97.
orc_rows <- slopes[slopes$arm == "oracle", ]
orc_rows$draw <- paste(orc_rows$n_species, orc_rows$phylo_signal,
                       ifelse(orc_rows$regime == "floor", "typical", orc_rows$regime), orc_rows$seed)
orc_rows <- orc_rows[!duplicated(orc_rows$draw), ]
orc  <- stats::aggregate(covered ~ n_species + phylo_signal, orc_rows, mean)
reps <- min(table(paste(orc_rows$n_species, orc_rows$phylo_signal)))
message("\noracle coverage (should be ~0.95), one count per distinct dataset: ",
        paste(round(range(orc$covered), 3), collapse = " to "),
        "   [", nrow(orc_rows), " distinct datasets; ", reps, " per n x signal]")

## Only worth checking once there are enough replicates to measure a proportion. At 1 replicate
## coverage can only be 0 or 1, so the check would cry wolf on every rehearsal and train everyone
## to ignore it. Monte Carlo SE of a 0.95 proportion is sqrt(.95 * .05 / reps), so 20 replicates
## is about where the 0.90-0.99 band starts to mean anything.
if (reps >= 20 && any(orc$covered < 0.90 | orc$covered > 0.99))
  warning("oracle coverage is off 0.95 in some scenarios -- the assumed true slope is suspect, ",
          "and until that is resolved no other arm's coverage means anything.", call. = FALSE)
