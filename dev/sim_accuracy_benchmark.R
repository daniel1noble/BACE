# =============================================================================
# sim_accuracy_benchmark.R - Per-type imputation accuracy with known truth
# =============================================================================
#
# For each response type {gaussian, poisson, binary, categorical (K=3),
# ordinal (K=4)} this script:
#   1. simulates a complete dataset via sim_bace() with high phylogenetic
#      signal (lambda = 0.9), 60 species, 120 cases;
#   2. masks ~20% of the response (MCAR);
#   3. runs bace() with conservative MCMC settings + skip_conv = TRUE;
#   4. compares the imputed values at the masked cells against truth using
#      the type-aware metrics from dev/benchmark_engine.R;
#   5. compares against a column-mean / modal-class BASELINE — BACE has to
#      beat this to demonstrate the phylogenetic structure is being used.
#
# Output:
#   - dev/sim_accuracy_results/<TIMESTAMP>/summary.csv
#   - dev/sim_accuracy_results/<TIMESTAMP>/per_type/<type>_metrics.csv
#   - console table summarising "BACE > baseline?" per type
#
# References
#   Rubin 1987          - Multiple Imputation for Nonresponse in Surveys.
#   van Buuren 2018 §3  - imputation accuracy diagnostics on simulated data.
#   Gneiting & Raftery 2007, JASA 102:359 - Brier score for class forecasts.
#   Münkemüller et al 2012, MEE 3:743 - reporting phylo signal alongside.
# =============================================================================

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(ape)
})

# Pull metric helpers from the benchmark engine so the simulated and
# empirical benchmarks share one source of truth for metric definitions.
ENGINE <- file.path("dev", "benchmark_engine.R")
if (!file.exists(ENGINE)) {
  stop("benchmark_engine.R not found at ", ENGINE,
       "; run from the package root.")
}
source_metric_helpers <- function(path) {
  src <- readLines(path)
  for (fn in c(".metrics_continuous", ".majority_vote",
               ".metrics_categorical", ".metrics_ordinal")) {
    start_i <- grep(paste0("^", fn, " <- function"), src, fixed = FALSE)[1]
    if (is.na(start_i)) stop("Function not found: ", fn)
    depth <- 0; end_i <- NA
    for (i in start_i:length(src)) {
      depth <- depth + sum(gregexpr("\\{", src[i])[[1]] > 0) -
                       sum(gregexpr("\\}", src[i])[[1]] > 0)
      if (i > start_i && depth == 0) { end_i <- i; break }
    }
    eval(parse(text = src[start_i:end_i]), envir = globalenv())
  }
}
source_metric_helpers(ENGINE)

# ---- Configuration -------------------------------------------------------

CFG <- list(
  n_species   = 60,
  n_cases     = 120,
  missing_pct = 0.20,
  phylo_sig   = 0.90,
  runs        = 4,
  n_final     = 5,
  nitt        = 2500,
  thin        = 10,
  burnin      = 500,
  seed        = 2026
)

# Default output dir is `singlerep_latest` — same name overwrites in
# place each run so we don't pollute the repo with timestamped dirs.
# Override with env var BACE_SIM_OUT to write somewhere else.
OUT_ROOT <- Sys.getenv("BACE_SIM_OUT", unset = "")
if (!nzchar(OUT_ROOT)) {
  OUT_ROOT <- file.path("dev", "sim_accuracy_results", "singlerep_latest")
}
dir.create(file.path(OUT_ROOT, "per_type"), recursive = TRUE,
           showWarnings = FALSE)

# Per-type response specs. Predictors are two complete gaussian columns
# in every case so the only thing varying across the sweep is the
# response type and family.
TYPE_SPECS <- list(
  list(label = "gaussian", response_type = "gaussian",
       metric = "nrmse", lower_is_better = TRUE),
  list(label = "poisson",  response_type = "poisson",
       metric = "nrmse", lower_is_better = TRUE),
  list(label = "binary",   response_type = "binary",
       metric = "balanced_accuracy", lower_is_better = FALSE),
  list(label = "categorical_K3", response_type = "multinomial3",
       metric = "balanced_accuracy", lower_is_better = FALSE),
  list(label = "ordinal_K4",     response_type = "threshold4",
       metric = "mae_level", lower_is_better = TRUE)
)

# ---- Per-type simulation + imputation harness ----------------------------

run_one_type <- function(spec, cfg) {
  t0 <- Sys.time()
  set.seed(cfg$seed)

  # 1. simulate complete data (no NAs) with high phylo signal
  sim <- sim_bace(
    response_type   = spec$response_type,
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    phylo_signal    = rep(cfg$phylo_sig, 3),
    n_cases         = cfg$n_cases,
    n_species       = cfg$n_species,
    missingness     = c(0, 0, 0)
  )
  complete_data <- sim$complete_data
  tree <- sim$tree
  names(complete_data)[names(complete_data) == "species"] <- "Species"

  # 2. apply MCAR mask to response
  set.seed(cfg$seed + 1)
  miss_idx <- sample.int(nrow(complete_data),
                         size = floor(nrow(complete_data) * cfg$missing_pct))
  data_masked <- complete_data
  data_masked$y[miss_idx] <- NA

  # 3. run BACE
  res <- suppressWarnings(suppressMessages(bace(
    fixformula     = "y ~ x1 + x2",
    ran_phylo_form = "~ 1 | Species",
    phylo          = tree,
    data           = data_masked,
    runs           = cfg$runs,
    n_final        = cfg$n_final,
    nitt           = cfg$nitt,
    thin           = cfg$thin,
    burnin         = cfg$burnin,
    skip_conv      = TRUE,
    max_attempts   = 1,
    n_cores        = 1L,
    verbose        = FALSE
  )))

  # 4. extract n_final imputed datasets, build (n_missing x n_final) matrix
  imp_datasets <- res$imputed_datasets
  truth_vals <- complete_data$y[miss_idx]

  # 5. type-aware metrics
  if (spec$response_type %in% c("gaussian", "poisson")) {
    imp_mat <- vapply(imp_datasets,
                       function(d) as.numeric(d$y[miss_idx]),
                       numeric(length(miss_idx)))
    if (is.vector(imp_mat)) imp_mat <- matrix(imp_mat, ncol = 1)
    sd_full <- stats::sd(complete_data$y, na.rm = TRUE)
    bace_metrics <- .metrics_continuous(imp_mat, as.numeric(truth_vals),
                                         sd_full)
    # Baseline: column-mean of OBSERVED values
    base_pred <- mean(data_masked$y, na.rm = TRUE)
    base_imp_mat <- matrix(base_pred, nrow = length(miss_idx),
                            ncol = cfg$n_final)
    base_metrics <- .metrics_continuous(base_imp_mat, as.numeric(truth_vals),
                                         sd_full)
  } else {
    # categorical / binary / ordinal
    truth_chr <- as.character(truth_vals)
    levels_all <- if (is.factor(complete_data$y)) levels(complete_data$y)
                   else sort(unique(as.character(complete_data$y)))
    im <- vapply(imp_datasets,
                  function(d) as.character(d$y[miss_idx]),
                  character(length(miss_idx)))
    if (is.vector(im)) im <- matrix(im, ncol = 1)
    if (spec$response_type %in% c("binary", "multinomial3")) {
      bace_metrics <- .metrics_categorical(im, truth_chr, levels_all)
    } else {
      # threshold/ordinal: levels are ordered by factor declaration
      bace_metrics <- .metrics_ordinal(im, truth_chr, levels_all)
    }
    # Baseline: modal class of OBSERVED response
    obs_y <- as.character(data_masked$y[!is.na(data_masked$y)])
    modal <- names(sort(table(obs_y), decreasing = TRUE))[1]
    base_im <- matrix(modal, nrow = length(miss_idx), ncol = cfg$n_final)
    if (spec$response_type %in% c("binary", "multinomial3")) {
      base_metrics <- .metrics_categorical(base_im, truth_chr, levels_all)
    } else {
      base_metrics <- .metrics_ordinal(base_im, truth_chr, levels_all)
    }
  }

  # 6. compose summary row
  bace_value <- bace_metrics[[spec$metric]]
  base_value <- base_metrics[[spec$metric]]
  if (spec$lower_is_better) {
    bace_better <- bace_value < base_value
  } else {
    bace_better <- bace_value > base_value
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  [%s] %s: BACE=%.3f  baseline=%.3f  better=%s  (%.1fs)\n",
              spec$label, spec$metric, bace_value, base_value,
              bace_better, elapsed))

  # write per-type detail
  detail <- data.frame(
    type            = spec$label,
    response_type   = spec$response_type,
    n_complete      = nrow(complete_data),
    n_missing       = length(miss_idx),
    metric          = spec$metric,
    BACE_value      = bace_value,
    baseline_value  = base_value,
    BACE_better     = bace_better,
    elapsed_seconds = elapsed,
    stringsAsFactors = FALSE
  )
  utils::write.csv(detail,
                   file.path(OUT_ROOT, "per_type",
                             sprintf("%s_metrics.csv", spec$label)),
                   row.names = FALSE)

  detail
}

# ---- Sweep ---------------------------------------------------------------

cat("\n========================================================\n")
cat("BACE per-type imputation accuracy benchmark\n")
cat("--------------------------------------------------------\n")
cat(sprintf("Config: n_species=%d  n_cases=%d  missing=%.0f%%  signal=%.2f\n",
            CFG$n_species, CFG$n_cases, 100 * CFG$missing_pct,
            CFG$phylo_sig))
cat(sprintf("MCMC:   runs=%d  n_final=%d  nitt=%d  thin=%d  burnin=%d\n",
            CFG$runs, CFG$n_final, CFG$nitt, CFG$thin, CFG$burnin))
cat(sprintf("Output: %s\n", OUT_ROOT))
cat("========================================================\n\n")

# Optional: run a single type via env var BACE_SIM_TYPE=<label>; default
# is the full sweep.
target <- Sys.getenv("BACE_SIM_TYPE", unset = "")
specs <- if (nzchar(target)) {
  hit <- TYPE_SPECS[vapply(TYPE_SPECS, function(s) s$label == target,
                           logical(1))]
  if (length(hit) == 0) stop("Unknown BACE_SIM_TYPE: ", target,
                              " — must be one of ",
                              paste(vapply(TYPE_SPECS, `[[`, "", "label"),
                                    collapse = ", "))
  hit
} else {
  TYPE_SPECS
}

results <- lapply(specs, function(s) run_one_type(s, CFG))
summary_df <- do.call(rbind, results)

utils::write.csv(summary_df, file.path(OUT_ROOT, "summary.csv"),
                 row.names = FALSE)

cat("\n========================================================\n")
cat("SUMMARY\n")
cat("--------------------------------------------------------\n")
print(summary_df, row.names = FALSE)
cat("\nResults written to: ", OUT_ROOT, "\n", sep = "")
n_better <- sum(summary_df$BACE_better)
cat(sprintf("BACE beat baseline on %d / %d types.\n",
            n_better, nrow(summary_df)))
