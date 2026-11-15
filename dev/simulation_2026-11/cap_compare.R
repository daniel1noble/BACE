# =================================================================================================
# cap_compare.R -- compare the cap-on and cap-off runs written by cap_experiment.sh.
#
#   Rscript cap_compare.R results_cap/
#
# WHAT THIS ANSWERS. BACE caps an imputed count's rate at 3 x the largest count it can see + 5.
# When the missing counts are the large ones, that cap pulls the imputations down. How much?
# Every dataset in the experiment was imputed twice from the same seed -- once with the cap
# (_capon) and once without (_capoff) -- so the only difference between a pair of files is the cap.
#
# HOW TO READ THE OUTPUT. "bias on the log scale" is the average of imputed-minus-true over the
# count trait's hidden cells, on the log(1 + count) scale, in units of the observed spread.
# Negative means the imputed counts run low. If it is clearly negative with the cap on and near
# zero with the cap off, the cap is doing damage; if the two agree, it is not. The "paired
# difference" is cap on minus cap off, averaged over datasets, with its standard error.
# =================================================================================================

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) > 0) args[1] else "results_cap"

## STEP 1. Find the result files. run_one.R named them <scenario>_s<seed>_capon.rds or _capoff.rds.
files <- list.files(results_dir, pattern = "_cap(on|off)[.]rds$", full.names = TRUE)
if (length(files) == 0) stop("no _capon / _capoff result files in ", results_dir, call. = FALSE)

## STEP 2. Pull the numbers we need out of ONE result file.
read_one <- function(file) {
  x   <- readRDS(file)
  arm <- if (grepl("_capon[.]rds$", file)) "cap_on" else "cap_off"

  ## A dataset that failed has no numbers. It still gets a full row -- every column, all NA --
  ## so the rows stack; it is reported and then left out. (Stacking a short row crashed this
  ## script on 2026-11-16; Shinichi found it.)
  if (!is.na(x$failed_at)) {
    return(data.frame(id = x$id, seed = x$seed, arm = arm, failed = TRUE,
                      rates_capped = NA, error_raw = NA, error_log = NA, bias_log = NA,
                      error_logdraws = NA, bias_logdraws = NA, coverage = NA, width = NA,
                      slope_bias = NA, slope_covered = NA, wall_minutes = NA))
  }

  ## The count trait's recovery row (error and bias on the hidden cells).
  count_row <- x$recovery[x$recovery$trait == "cnt", ]

  ## Its 95% interval coverage and width, if the intervals were computed.
  coverage <- NA
  width    <- NA
  if (!is.null(x$cell_coverage)) {
    cover_row <- x$cell_coverage[x$cell_coverage$variable == "cnt", ]
    if (nrow(cover_row) == 1) {
      coverage <- cover_row$coverage
      width    <- cover_row$width
    }
  }

  ## The downstream slope of x2 on x1 from the BACE arm, compared with the truth recorded in
  ## the file (0.5 when the traits are related, 0 in the "floor" regime).
  slope_row <- x$slopes[x$slopes$arm == "bace", ]
  slope_bias    <- if (nrow(slope_row) == 1) slope_row$estimate - slope_row$truth else NA
  slope_covered <- if (nrow(slope_row) == 1) slope_row$covered else NA

  ## The log-scale measures that average the imputations on the log scale first (added to the
  ## scoring on 2026-11-16; NA in files written before that).
  error_logdraws <- if ("zrmse_logdraws" %in% names(count_row)) count_row$zrmse_logdraws else NA
  bias_logdraws  <- if ("bias_logdraws"  %in% names(count_row)) count_row$bias_logdraws  else NA

  ## How many imputed rates the cap truncated in this run (0 when the cap was off).
  rates_capped <- 0
  for (one_variable in x$rate_cap) rates_capped <- rates_capped + one_variable$n_rates

  data.frame(id = x$id, seed = x$seed, arm = arm, failed = FALSE,
             rates_capped   = rates_capped,
             error_raw      = count_row$zrmse,
             error_log      = count_row$zrmse_log,
             bias_log       = count_row$bias_log,
             error_logdraws = error_logdraws,
             bias_logdraws  = bias_logdraws,
             coverage       = coverage,
             width          = width,
             slope_bias     = slope_bias,
             slope_covered  = as.numeric(slope_covered),
             wall_minutes   = x$wall_s / 60)
}

## STEP 3. One row per run.
rows <- list()
for (file in files) rows[[length(rows) + 1]] <- read_one(file)
runs <- do.call(rbind, rows)

if (any(runs$failed)) {
  message(sum(runs$failed), " run(s) failed and are left out: ",
          paste(runs$id[runs$failed], runs$seed[runs$failed], runs$arm[runs$failed], collapse = "; "))
  runs <- runs[!runs$failed, ]
}
utils::write.csv(runs, file.path(results_dir, "cap_per_run.csv"), row.names = FALSE)

## STEP 4. Pair each cap-on run with the cap-off run of the SAME dataset (same scenario and seed).
cap_on  <- runs[runs$arm == "cap_on", ]
cap_off <- runs[runs$arm == "cap_off", ]
pair_key_on  <- paste(cap_on$id,  cap_on$seed)
pair_key_off <- paste(cap_off$id, cap_off$seed)
cap_off <- cap_off[match(pair_key_on, pair_key_off), ]   # reorder cap-off to line up with cap-on
has_pair <- !is.na(cap_off$id)
cap_on  <- cap_on[has_pair, ]
cap_off <- cap_off[has_pair, ]
cat("\n", nrow(cap_on), "datasets with both a cap-on and a cap-off run\n")

## STEP 5. For each measure: the cap-on mean, the cap-off mean, and the paired difference.
report <- function(label, on, off) {
  difference <- on - off
  n  <- sum(!is.na(difference))
  se <- stats::sd(difference, na.rm = TRUE) / sqrt(n)
  cat(sprintf("  %-24s cap on %7.3f   cap off %7.3f   difference %+7.3f  (SE %.3f, n = %d)\n",
              label, mean(on, na.rm = TRUE), mean(off, na.rm = TRUE),
              mean(difference, na.rm = TRUE), se, n))
}
cat(sprintf("  rates capped per dataset with the cap on: median %d, range %d to %d\n",
            as.integer(median(cap_on$rates_capped)), min(cap_on$rates_capped), max(cap_on$rates_capped)))
cat("\n  count trait, hidden cells:\n")
report("error, raw scale",   cap_on$error_raw,  cap_off$error_raw)
report("error, log scale",   cap_on$error_log,  cap_off$error_log)
report("bias, log scale",    cap_on$bias_log,   cap_off$bias_log)
report("error, log of draws", cap_on$error_logdraws, cap_off$error_logdraws)
report("bias, log of draws",  cap_on$bias_logdraws,  cap_off$bias_logdraws)
report("95% interval coverage", cap_on$coverage, cap_off$coverage)
report("interval width",     cap_on$width,      cap_off$width)
cat("\n  downstream slope of x2 on x1, BACE arm (estimate minus the truth recorded for that dataset):\n")
report("slope bias",         cap_on$slope_bias, cap_off$slope_bias)
report("interval covered",   cap_on$slope_covered,  cap_off$slope_covered)
report("wall time, minutes", cap_on$wall_minutes,   cap_off$wall_minutes)

## STEP 6. The same, dataset by dataset, so one odd dataset cannot hide behind the averages.
cat("\n  per dataset:\n")
per_dataset <- data.frame(id = cap_on$id, seed = cap_on$seed,
                          rates_capped  = cap_on$rates_capped,
                          bias_cap_on   = round(cap_on$bias_log, 3),
                          bias_cap_off  = round(cap_off$bias_log, 3),
                          error_cap_on  = round(cap_on$error_log, 3),
                          error_cap_off = round(cap_off$error_log, 3))
print(per_dataset, row.names = FALSE)
