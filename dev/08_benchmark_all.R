# =============================================================================
# 08_benchmark_all.R - Cross-dataset BACE benchmark driver
# =============================================================================
# Runs benchmark_dataset() on every bundled dataset in turn and stacks the
# resulting metric tables into a single tidy frame for cross-dataset
# comparison. Each dataset's individual run also writes its own bundle
# under dev/benchmark_results/<dataset>/.
#
# This script is the answer to "is BACE working everywhere?" — the printed
# table tells you, per dataset and per trait, where BACE is calibrated
# (coverage95 ~ 0.95), where it's accurate (low NRMSE / high accuracy),
# and where it's struggling (low correlation, low balanced_accuracy).
#
# Datasets included:
#   pantheria  : mammals, mixed cont/count/ordinal/binary
#   amphibio   : amphibians, mixed cont/binary/ordinal/categorical
#   bien       : plants, all-continuous
#   globtherm  : thermal limits, all-continuous (raw scale)
#   leptraits  : lepidopterans, mixed cont/binary
#
# AVONET is benchmarked separately by 00_benchmark_AVONET.R; you can
# stack its output here too once that script is migrated to the engine.
# =============================================================================

devtools::load_all(quiet = TRUE)
source("dev/benchmark_engine.R")

# Tiny null-coalescing helper used in the per-dataset config block
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Configuration ---------------------------------------------------------
# Per-dataset config: log_traits + (optionally) MCMC overrides.
# Defaults inherited from benchmark_dataset() unless overridden.
DATASETS <- list(
  pantheria = list(
    log_traits = c("body_mass_g", "head_body_length_mm",
                   "gestation_d", "max_longevity_m")
  ),
  amphibio  = list(
    log_traits = c("body_size_mm", "body_mass_g")
  ),
  bien      = list(
    log_traits = c("height_m", "leaf_area", "sla",
                   "seed_mass", "wood_density")
  ),
  globtherm = list(
    log_traits = character(0)
  ),
  leptraits = list(
    log_traits = c("wingspan_lower", "forewing_length_lower",
                   "n_hostplant_families"),
    pre_hook   = function(traits) {
      # Coerce 0/1 monthly indicators to factor for binary (threshold) fitting.
      months <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")
      for (m in intersect(months, colnames(traits))) {
        traits[[m]] <- factor(traits[[m]], levels = c(0L, 1L))
      }
      traits
    }
  )
)

SUBSET_N     <- 2000L
NITT         <- 20000
THIN         <- 15
BURNIN       <- 4000
RUNS         <- 5
N_FINAL      <- 10
MAX_ATTEMPTS <- 2
N_CORES      <- 4L

# ---- Run all datasets ------------------------------------------------------
results <- list()
for (nm in names(DATASETS)) {
  cfg <- DATASETS[[nm]]

  load(sprintf("dev/testing_data/data/%s_traits.rda", nm))
  load(sprintf("dev/testing_data/data/%s_tree.rda",   nm))
  traits <- get(sprintf("%s_traits", nm))
  tree   <- get(sprintf("%s_tree",   nm))

  if (!is.null(cfg$pre_hook)) traits <- cfg$pre_hook(traits)

  cat(sprintf("\n\n############################################\n"))
  cat(sprintf("###  %s\n", nm))
  cat(sprintf("############################################\n"))

  results[[nm]] <- tryCatch(
    benchmark_dataset(
      traits_df    = traits,
      tree         = tree,
      dataset_name = nm,
      log_traits   = cfg$log_traits %||% character(0),
      subset_n     = SUBSET_N,
      nitt = NITT, thin = THIN, burnin = BURNIN,
      runs = RUNS, n_final = N_FINAL,
      max_attempts = MAX_ATTEMPTS, n_cores = N_CORES,
      verbose = TRUE
    ),
    error = function(e) {
      cat(sprintf("\n[%s] FAILED: %s\n", nm, conditionMessage(e)))
      NULL
    }
  )
}

# ---- Stack metrics + write cross-dataset table -----------------------------
ok <- !vapply(results, is.null, logical(1))
if (!any(ok)) stop("All benchmarks failed.")

all_metrics <- do.call(rbind, lapply(results[ok], `[[`, "metrics"))
all_signal  <- do.call(rbind, lapply(seq_along(results)[ok], function(i) {
  cbind(dataset = names(results)[i], results[[i]]$phylo_signal)
}))
runtimes <- data.frame(
  dataset     = names(results)[ok],
  runtime_min = vapply(results[ok], `[[`, numeric(1), "runtime_min"),
  converged   = vapply(results[ok], function(r) r$res$converged, logical(1)),
  stringsAsFactors = FALSE
)

date_str <- format(Sys.Date(), "%Y%m%d")
out_dir  <- file.path("dev", "benchmark_results", "cross_dataset",
                       paste0("run_", date_str))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

utils::write.csv(all_metrics, file.path(out_dir, "all_metrics.csv"),
                 row.names = FALSE)
utils::write.csv(all_signal,  file.path(out_dir, "all_phylo_signal.csv"),
                 row.names = FALSE)
utils::write.csv(runtimes,    file.path(out_dir, "runtimes.csv"),
                 row.names = FALSE)

cat("\n\n============================================\n")
cat("  Cross-dataset summary\n")
cat("============================================\n\n")
cat("Continuous / count traits (NRMSE / cor / coverage):\n")
cont <- all_metrics[all_metrics$type %in% c("continuous", "count"), ,
                    drop = FALSE]
print(cont[, c("dataset", "trait", "type", "n_hidden",
               "nrmse", "correlation", "coverage95")],
      digits = 3, row.names = FALSE)

cat("\nCategorical / binary / ordinal traits (acc / bal_acc / Brier):\n")
cat_m <- all_metrics[all_metrics$type %in% c("categorical","binary","ordinal"), ,
                     drop = FALSE]
print(cat_m[, c("dataset", "trait", "type", "n_hidden",
                "accuracy", "balanced_accuracy", "brier", "mae_level")],
      digits = 3, row.names = FALSE)

cat("\nRuntimes (min):\n")
print(runtimes, row.names = FALSE)

cat(sprintf("\nFull tables written to: %s\n", out_dir))
