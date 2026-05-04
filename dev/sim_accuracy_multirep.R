# =============================================================================
# sim_accuracy_multirep.R - Multi-replicate per-type accuracy with paired
#                           BACE-vs-baseline test + extra diagnostics for
#                           the poisson NRMSE pathology flagged by the
#                           single-rep run.
# =============================================================================
#
# What this script answers
# ------------------------
#   1. Across N_REPS independent simulations per response type, does BACE
#      reliably beat a column-mean / modal-class baseline? Reports paired
#      Wilcoxon signed-rank p-values so a single bad draw doesn't decide
#      the verdict.
#   2. For continuous types (gaussian, poisson), reports THREE accuracy
#      metrics so we can tell sampling-variance pathologies from real
#      signal-recovery failures:
#         - nrmse_perimp     mean over imputations of NRMSE-per-imputation
#                            (current single-rep headline; penalises BACE's
#                            posterior-predictive sampling variance)
#         - nrmse_pointest   NRMSE of the per-cell posterior MEAN against
#                            truth (integrates out the per-iteration
#                            sampling variance — this is the fair head-to-
#                            head against the column-mean baseline)
#         - spearman         Spearman rank correlation of the per-cell
#                            posterior mean against truth (rank-based, so
#                            scale-free)
#         - coverage95       fraction of truth values inside BACE's per-cell
#                            95% posterior predictive interval (calibration,
#                            not accuracy)
#
# Runtime
# -------
#   ~5 reps x 5 types x ~3-5 sec/run ~= 1.5 min. Set N_REPS via env var
#   BACE_SIM_REPS (default 20). For final manuscript numbers run with
#   N_REPS = 50 or 100.
# =============================================================================

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({ library(ape) })

ENGINE <- file.path("dev", "benchmark_engine.R")
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

N_REPS <- as.integer(Sys.getenv("BACE_SIM_REPS", unset = "20"))

CFG <- list(
  n_species   = 60,
  n_cases     = 120,
  missing_pct = 0.20,
  phylo_sig   = 0.90,
  runs        = 4,
  # n_final = 50 is required for honest coverage95 estimation. With
  # only 5 imputations the empirical [2.5%, 97.5%] quantile is ~ [min,
  # max], which approximates a 67% PI not 95% — coverage looks bad
  # even when the posterior is correctly calibrated. Verified
  # empirically: gaussian coverage ~0.61 at n_final=5 vs ~0.92 at
  # n_final=50. See van Buuren 2018 §2.8 on choosing n_imputations.
  n_final     = 50,
  nitt        = 2500,
  thin        = 10,
  burnin      = 500,
  base_seed   = 2026
)

# Default output dir is `multirep_latest_n<REPS>` — same name overwrites
# in place each run so we don't pollute the repo with timestamped dirs.
# Override with env var BACE_SIM_OUT to write somewhere else.
OUT_ROOT <- Sys.getenv("BACE_SIM_OUT", unset = "")
if (!nzchar(OUT_ROOT)) {
  OUT_ROOT <- file.path("dev", "sim_accuracy_results",
                        sprintf("multirep_latest_n%d", N_REPS))
}
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

TYPE_SPECS <- list(
  list(label = "gaussian", response_type = "gaussian", family = "continuous"),
  list(label = "poisson",  response_type = "poisson",  family = "continuous"),
  list(label = "binary",   response_type = "binary",   family = "categorical"),
  list(label = "categorical_K3", response_type = "multinomial3",
       family = "categorical"),
  list(label = "ordinal_K4", response_type = "threshold4",
       family = "ordinal")
)

# ---- Per-rep harness -----------------------------------------------------
#
# Returns a one-row data.frame with all relevant metrics for one rep
# of one type.

run_one_rep <- function(spec, cfg, rep_id) {
  set.seed(cfg$base_seed + rep_id * 13)

  # 1. simulate complete data
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
  set.seed(cfg$base_seed + rep_id * 13 + 1)
  miss_idx <- sample.int(nrow(complete_data),
                         size = floor(nrow(complete_data) * cfg$missing_pct))
  data_masked <- complete_data
  data_masked$y[miss_idx] <- NA

  # 3. run BACE. n_cores parallelises bace_final_imp (mclapply, so
  # macOS / Linux only — set to 1 on Windows). Runtime savings are
  # large at n_final=50 (~10x speedup with 4 cores).
  n_cores_use <- if (.Platform$OS.type == "windows") 1L
                  else as.integer(Sys.getenv("BACE_SIM_CORES", unset = "4"))
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
    n_cores        = n_cores_use,
    verbose        = FALSE
  )))

  imp_datasets <- res$imputed_datasets
  truth_vals <- complete_data$y[miss_idx]

  # 4. compute metrics — different sets per family
  if (spec$family == "continuous") {
    imp_mat <- vapply(imp_datasets,
                       function(d) as.numeric(d$y[miss_idx]),
                       numeric(length(miss_idx)))
    if (is.vector(imp_mat)) imp_mat <- matrix(imp_mat, ncol = 1)
    sd_full <- stats::sd(complete_data$y, na.rm = TRUE)

    # BACE: per-imputation NRMSE (current headline; penalises sampling)
    bm <- .metrics_continuous(imp_mat, as.numeric(truth_vals), sd_full)

    # BACE: NRMSE on the per-cell posterior mean (integrates out sampling)
    pe_pred <- rowMeans(imp_mat)
    nrmse_pointest <- sqrt(mean((pe_pred - truth_vals)^2)) / sd_full
    spearman <- if (length(unique(truth_vals)) > 1 &&
                    length(unique(pe_pred)) > 1) {
      suppressWarnings(stats::cor(pe_pred, as.numeric(truth_vals),
                                  method = "spearman"))
    } else NA_real_

    # baseline: column-mean of observed
    base_pred <- mean(data_masked$y, na.rm = TRUE)
    base_imp_mat <- matrix(base_pred, nrow = length(miss_idx),
                            ncol = cfg$n_final)
    bb <- .metrics_continuous(base_imp_mat, as.numeric(truth_vals), sd_full)
    base_nrmse_pointest <- sqrt(mean((base_pred - truth_vals)^2)) / sd_full
    base_spearman <- NA_real_  # constant prediction has no rank

    data.frame(
      type             = spec$label,
      family           = spec$family,
      rep              = rep_id,
      n_missing        = length(miss_idx),
      bace_nrmse_perimp     = bm$nrmse,
      bace_nrmse_pointest   = nrmse_pointest,
      bace_spearman         = spearman,
      bace_coverage95       = bm$coverage95,
      base_nrmse_perimp     = bb$nrmse,
      base_nrmse_pointest   = base_nrmse_pointest,
      base_spearman         = base_spearman,
      base_coverage95       = bb$coverage95,
      stringsAsFactors = FALSE
    )
  } else {
    truth_chr <- as.character(truth_vals)
    levels_all <- if (is.factor(complete_data$y)) levels(complete_data$y)
                   else sort(unique(as.character(complete_data$y)))
    im <- vapply(imp_datasets,
                  function(d) as.character(d$y[miss_idx]),
                  character(length(miss_idx)))
    if (is.vector(im)) im <- matrix(im, ncol = 1)

    if (spec$family == "ordinal") {
      bm <- .metrics_ordinal(im, truth_chr, levels_all)
    } else {
      bm <- .metrics_categorical(im, truth_chr, levels_all)
    }

    obs_y <- as.character(data_masked$y[!is.na(data_masked$y)])
    modal <- names(sort(table(obs_y), decreasing = TRUE))[1]
    base_im <- matrix(modal, nrow = length(miss_idx), ncol = cfg$n_final)
    if (spec$family == "ordinal") {
      bb <- .metrics_ordinal(base_im, truth_chr, levels_all)
    } else {
      bb <- .metrics_categorical(base_im, truth_chr, levels_all)
    }

    data.frame(
      type             = spec$label,
      family           = spec$family,
      rep              = rep_id,
      n_missing        = length(miss_idx),
      bace_accuracy          = bm$accuracy,
      bace_balanced_accuracy = bm$balanced_accuracy,
      bace_brier             = bm$brier,
      bace_mae_level         = if (spec$family == "ordinal") bm$mae_level
                               else NA_real_,
      base_accuracy          = bb$accuracy,
      base_balanced_accuracy = bb$balanced_accuracy,
      base_brier             = bb$brier,
      base_mae_level         = if (spec$family == "ordinal") bb$mae_level
                               else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

# ---- Per-type aggregation --------------------------------------------------

aggregate_continuous <- function(df) {
  with(df, {
    # paired Wilcoxon signed-rank: BACE point-estimate vs. baseline
    wt <- tryCatch(
      stats::wilcox.test(bace_nrmse_pointest, base_nrmse_pointest,
                         paired = TRUE, exact = FALSE),
      error = function(e) list(p.value = NA_real_)
    )
    data.frame(
      type = unique(df$type),
      family = unique(df$family),
      n_reps = nrow(df),
      bace_nrmse_perimp_mean   = mean(bace_nrmse_perimp,   na.rm = TRUE),
      bace_nrmse_perimp_sd     = stats::sd(bace_nrmse_perimp,   na.rm = TRUE),
      bace_nrmse_pointest_mean = mean(bace_nrmse_pointest, na.rm = TRUE),
      bace_nrmse_pointest_sd   = stats::sd(bace_nrmse_pointest, na.rm = TRUE),
      bace_spearman_mean       = mean(bace_spearman,       na.rm = TRUE),
      bace_spearman_sd         = stats::sd(bace_spearman,       na.rm = TRUE),
      bace_coverage95_mean     = mean(bace_coverage95,     na.rm = TRUE),
      base_nrmse_perimp_mean   = mean(base_nrmse_perimp,   na.rm = TRUE),
      base_nrmse_pointest_mean = mean(base_nrmse_pointest, na.rm = TRUE),
      paired_wilcox_p_pointest = wt$p.value,
      bace_better_pointest     = mean(bace_nrmse_pointest <
                                       base_nrmse_pointest, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
}

aggregate_categorical <- function(df) {
  is_ord <- unique(df$family) == "ordinal"
  metric <- if (is_ord) "mae_level" else "balanced_accuracy"
  bace_vec <- df[[paste0("bace_", metric)]]
  base_vec <- df[[paste0("base_", metric)]]
  wt <- tryCatch(
    stats::wilcox.test(bace_vec, base_vec, paired = TRUE, exact = FALSE),
    error = function(e) list(p.value = NA_real_)
  )
  bace_better <- if (is_ord) mean(bace_vec < base_vec, na.rm = TRUE)
                 else mean(bace_vec > base_vec, na.rm = TRUE)
  data.frame(
    type = unique(df$type),
    family = unique(df$family),
    n_reps = nrow(df),
    bace_metric = metric,
    bace_metric_mean = mean(bace_vec, na.rm = TRUE),
    bace_metric_sd   = stats::sd(bace_vec, na.rm = TRUE),
    base_metric_mean = mean(base_vec, na.rm = TRUE),
    base_metric_sd   = stats::sd(base_vec, na.rm = TRUE),
    bace_acc_mean    = mean(df$bace_accuracy, na.rm = TRUE),
    base_acc_mean    = mean(df$base_accuracy, na.rm = TRUE),
    bace_brier_mean  = mean(df$bace_brier, na.rm = TRUE),
    base_brier_mean  = mean(df$base_brier, na.rm = TRUE),
    paired_wilcox_p  = wt$p.value,
    bace_better      = bace_better,
    stringsAsFactors = FALSE
  )
}

# ---- Sweep ---------------------------------------------------------------

cat("\n========================================================\n")
cat("BACE multi-replicate per-type accuracy benchmark\n")
cat("--------------------------------------------------------\n")
cat(sprintf("Reps:   %d per type\n", N_REPS))
cat(sprintf("Config: n_species=%d  n_cases=%d  missing=%.0f%%  signal=%.2f\n",
            CFG$n_species, CFG$n_cases, 100 * CFG$missing_pct,
            CFG$phylo_sig))
cat(sprintf("Output: %s\n", OUT_ROOT))
cat("========================================================\n\n")

all_rows <- list()
for (spec in TYPE_SPECS) {
  cat(sprintf("--- %s -----------------------------------\n", spec$label))
  type_rows <- list()
  for (rep_id in seq_len(N_REPS)) {
    t0 <- Sys.time()
    one <- run_one_rep(spec, CFG, rep_id)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    one$elapsed_seconds <- elapsed
    type_rows[[rep_id]] <- one
    cat(sprintf("  rep %2d/%d  (%.1fs)\n", rep_id, N_REPS, elapsed))
  }
  type_df <- do.call(rbind, type_rows)
  utils::write.csv(type_df,
                   file.path(OUT_ROOT, sprintf("%s_reps.csv", spec$label)),
                   row.names = FALSE)
  all_rows[[spec$label]] <- type_df
}

# ---- Aggregate + report ---------------------------------------------------

cont_summary <- do.call(rbind,
                        lapply(all_rows[c("gaussian", "poisson")],
                               aggregate_continuous))
cat_summary  <- do.call(rbind,
                        lapply(all_rows[c("binary", "categorical_K3",
                                          "ordinal_K4")],
                               aggregate_categorical))

utils::write.csv(cont_summary,
                 file.path(OUT_ROOT, "continuous_summary.csv"),
                 row.names = FALSE)
utils::write.csv(cat_summary,
                 file.path(OUT_ROOT, "categorical_summary.csv"),
                 row.names = FALSE)

cat("\n=========================================================\n")
cat("RESULTS — BACE imputed-vs-truth on held-out missing cells\n")
cat("=========================================================\n")
cat("baseline = no-skill imputer (column-mean / modal class).\n")
cat("BACE has to beat this floor to demonstrate phylogenetic value.\n\n")

cat("--- Continuous types (gaussian, poisson) ---\n")
clean_cont <- data.frame(
  type             = cont_summary$type,
  reps             = cont_summary$n_reps,
  bace_NRMSE       = round(cont_summary$bace_nrmse_pointest_mean, 3),
  baseline_NRMSE   = round(cont_summary$base_nrmse_pointest_mean, 3),
  bace_spearman    = round(cont_summary$bace_spearman_mean, 3),
  bace_coverage95  = round(cont_summary$bace_coverage95_mean, 3),
  paired_p         = signif(cont_summary$paired_wilcox_p_pointest, 3),
  bace_better_pct  = round(100 * cont_summary$bace_better_pointest),
  stringsAsFactors = FALSE
)
print(clean_cont, row.names = FALSE)

cat("\n--- Categorical / ordinal types (binary, K=3 cat, K=4 ord) ---\n")
clean_cat <- data.frame(
  type             = cat_summary$type,
  reps             = cat_summary$n_reps,
  metric           = cat_summary$bace_metric,
  bace             = round(cat_summary$bace_metric_mean, 3),
  baseline         = round(cat_summary$base_metric_mean, 3),
  bace_accuracy    = round(cat_summary$bace_acc_mean, 3),
  bace_brier       = round(cat_summary$bace_brier_mean, 3),
  paired_p         = signif(cat_summary$paired_wilcox_p, 3),
  bace_better_pct  = round(100 * cat_summary$bace_better),
  stringsAsFactors = FALSE
)
print(clean_cat, row.names = FALSE)

cat("\nMetric definitions:\n")
cat("  NRMSE          = sqrt(mean((imputed - truth)^2)) / sd(complete y)\n")
cat("                   imputed = per-cell posterior mean across n_final draws\n")
cat("  spearman       = rank correlation of posterior-mean imputation vs truth\n")
cat("  coverage95     = fraction of truth values inside BACE's 95% PI\n")
cat("                   (well-calibrated -> ~0.95)\n")
cat("  bace_accuracy  = fraction of held-out cells where BACE's modal\n")
cat("                   imputed class matches truth\n")
cat("  brier          = mean squared error of class-probability vector vs\n")
cat("                   one-hot truth (lower = better)\n")
cat("  paired_p       = paired Wilcoxon BACE-vs-baseline across reps\n")
cat("  bace_better_pct= fraction of reps where BACE strictly beat baseline\n")

cat("\nResults written to: ", OUT_ROOT, "\n", sep = "")
