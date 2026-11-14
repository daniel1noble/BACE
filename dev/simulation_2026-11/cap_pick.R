# =================================================================================================
# cap_pick.R -- which campaign datasets did BACE's count cap hit hardest?
#
#   Rscript cap_pick.R results/ 10        # prints "scenario seed", one per line, ten of them
#
# Every campaign result file records how many imputed rates the cap truncated (rate_cap). This
# ranks the finished datasets by that number and prints the top ones, so cap_experiment.sh reruns
# the datasets where the cap did the most, with the cap on and off. Measured on one dataset each
# way: at signal 0.3 the cap fired 39 times, at signal 0.9 not at all -- it bites when the counts
# BACE can see are small, not when the signal is strong -- so guessing the scenario is worse than
# reading the records.
# =================================================================================================

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) > 0) args[1] else "results"
how_many    <- if (length(args) > 1) as.integer(args[2]) else 10L

## STEP 1. Read every finished campaign result and count its cap events.
files <- list.files(results_dir, pattern = "_s[0-9]+[.]rds$", full.names = TRUE)
if (length(files) == 0) stop("no campaign results in ", results_dir, call. = FALSE)
rows <- list()
for (file in files) {
  x <- readRDS(file)
  if (!is.na(x$failed_at)) next                      # a failed dataset has no imputations
  rates_capped <- 0
  for (one_variable in x$rate_cap) rates_capped <- rates_capped + one_variable$n_rates
  rows[[length(rows) + 1]] <- data.frame(id = x$id, seed = x$seed, rates_capped = rates_capped)
}
counts <- do.call(rbind, rows)

## STEP 2. Keep the datasets where the cap fired most, and print them for the launcher.
counts <- counts[order(-counts$rates_capped), ]
top <- head(counts, how_many)
message(nrow(counts), " finished datasets; cap fired in ", sum(counts$rates_capped > 0),
        "; top ", nrow(top), " have ", paste(top$rates_capped, collapse = ", "), " rates capped")
cat(paste(top$id, top$seed), sep = "\n")
