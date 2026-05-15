# =============================================================================
# 10_evaluate_reference_datasets.R
#
# Run BACE end-to-end on each replicate of the four reference datasets
# produced by dev/09_generate_reference_datasets.R, and score three benchmark
# axes against the saved truth:
#
#   1. Cell-level imputation accuracy on hidden cells (NRMSE / MAE / cor /
#      accuracy / balanced_accuracy / ordinal_mae / coverage95 / brier).
#      Reuses evaluate_imputation_ensemble() from 02_benchmark_simulated_full.R.
#
#   2. Beta-coverage of the response equation, scored against BOTH:
#        (a) an Oracle MCMCglmm fit on the complete (un-masked) data, with
#            the same formula/family BACE uses, and the same MCMC budget.
#            This is the standard MI-calibration comparator
#            (van Buuren 2018 FIMD chapter 9).
#        (b) the dialled-in true_beta_resp_full from the simulation truth.
#            We skip x4 here because of the contrast-coding mismatch
#            (sim_bace treats x4 ordered as a single integer coef;
#            BACE's model.matrix gives two poly-contrast columns).
#
# Phylogenetic signal recovery (axis 3 in the qmd) is intentionally NOT
# computed here. It is a post-hoc step (separate script) that runs
# phylo_signal_summary() once per (dataset, rep) on the first imputed
# dataset and compares to true_phylo_signal. Bundling it into this loop
# would add ~3 hr of MCMC and is independent of bace() itself.
#
# Output layout
#   dev/simulation_results/evaluation_results/
#     {sim_ideal,sim_typical,sim_heterogeneous,sim_hard}/eval_rep_NN.rds
#     evaluation_summary.csv  (built at the end from per-rep RDS files)
#
# Resumable: an eval_rep_NN.rds whose `status` is "ok" is skipped.
# =============================================================================

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape)
  library(MASS)
  library(parallel)
  library(coda)
  library(MCMCglmm)
})

# -----------------------------------------------------------------------------
# Reuse evaluate_imputation_ensemble() from 02_benchmark_simulated_full.R
# (same helper used in 01 and 02 — single source of truth for cell metrics).
# -----------------------------------------------------------------------------
source_funcs <- function(path, fn_names) {
  src <- readLines(path)
  for (fn in fn_names) {
    start_i <- grep(paste0("^", fn, " <- function"), src)[1]
    if (is.na(start_i)) stop("Function not found in ", path, ": ", fn)
    depth <- 0; end_i <- NA
    for (i in start_i:length(src)) {
      depth <- depth +
        sum(gregexpr("\\{", src[i])[[1]] > 0) -
        sum(gregexpr("\\}", src[i])[[1]] > 0)
      if (i > start_i && depth == 0) { end_i <- i; break }
    }
    eval(parse(text = src[start_i:end_i]), envir = globalenv())
  }
}
source_funcs(
  file.path("dev", "02_benchmark_simulated_full.R"),
  c("evaluate_imputation_ensemble")
)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
REF_ROOT  <- file.path("dev", "simulation_results", "reference_datasets")
EVAL_ROOT <- file.path("dev", "simulation_results", "evaluation_results")
if (!dir.exists(EVAL_ROOT)) dir.create(EVAL_ROOT, recursive = TRUE)

# Production MCMC budget (user choice).
# Individual fields can be overridden via env vars BACE_NITT / BACE_BURNIN
# / BACE_THIN / BACE_RUNS / BACE_NFINAL — useful for fast smoke tests.
.envint <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name, unset = NA_character_)))
  if (is.na(v)) default else v
}
MCMC <- list(
  nitt    = .envint("BACE_NITT",   50000L),
  burnin  = .envint("BACE_BURNIN", 10000L),
  thin    = .envint("BACE_THIN",   25L),
  runs    = .envint("BACE_RUNS",   10L),
  n_final = .envint("BACE_NFINAL", 50L),
  max_attempts = 2L
)

# Outer parallelism: rep-level via mclapply. Each bace() call runs single-
# threaded (n_cores = 1L) to avoid nested parallelism contention.
# Override on the command line / in CI via env BACE_EVAL_CORES.
N_CORES_OUTER <- suppressWarnings(as.integer(
  Sys.getenv("BACE_EVAL_CORES", unset = "6")
))
if (is.na(N_CORES_OUTER) || N_CORES_OUTER < 1L) N_CORES_OUTER <- 6L

# Datasets and reps to evaluate. NULL = all reps in each dataset folder.
DATASETS <- c("sim_ideal", "sim_typical", "sim_heterogeneous", "sim_hard")

# -----------------------------------------------------------------------------
# Per-rep evaluator
# -----------------------------------------------------------------------------
evaluate_one_rep <- function(ds_name, rep_id, mcmc = MCMC) {

  rep_file  <- file.path(REF_ROOT, ds_name, sprintf("rep_%02d.rds", rep_id))
  meta_file <- file.path(REF_ROOT, ds_name, "meta.rds")
  out_dir   <- file.path(EVAL_ROOT, ds_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_file  <- file.path(out_dir, sprintf("eval_rep_%02d.rds", rep_id))

  # Resume: skip if a successful output already exists.
  if (file.exists(out_file)) {
    prev <- tryCatch(readRDS(out_file), error = function(e) NULL)
    if (!is.null(prev) && identical(prev$status, "ok")) return(invisible(prev))
  }

  rb   <- readRDS(rep_file)
  meta <- readRDS(meta_file)

  result <- list(
    dataset   = ds_name,
    rep_id    = rep_id,
    seed      = rb$seed,
    mechanism = rb$mechanism,
    rate      = rb$rate,
    n_cases   = nrow(rb$complete_data),
    settings  = mcmc,
    status    = "starting",
    error     = NA_character_,
    started_at = Sys.time()
  )

  ok <- tryCatch({

    # ------------------------------------------------------------------
    # 1. BACE on the masked data
    # ------------------------------------------------------------------
    t_bace <- Sys.time()
    bace_res <- bace(
      fixformula     = list(
        "y  ~ x1 + x2 + x3 + x4",
        "x1 ~ y  + x2 + x3 + x4",
        "x2 ~ y  + x1 + x3 + x4",
        "x3 ~ y  + x1 + x2 + x4",
        "x4 ~ y  + x1 + x2 + x3"
      ),
      ran_phylo_form = "~1|Species",
      phylo          = rb$tree,
      data           = rb$miss_data,
      nitt           = mcmc$nitt, thin = mcmc$thin, burnin = mcmc$burnin,
      runs           = mcmc$runs, n_final = mcmc$n_final,
      species        = FALSE,
      verbose        = FALSE,
      skip_conv      = FALSE,
      max_attempts   = mcmc$max_attempts,
      n_cores        = 1L
    )
    bace_runtime <- as.numeric(difftime(Sys.time(), t_bace, units = "secs"))

    # ------------------------------------------------------------------
    # 2. Cell-level metrics on hidden cells
    # ------------------------------------------------------------------
    var_types_eval <- c(y = "gaussian", x1 = "categorical", x2 = "categorical",
                        x3 = "count",   x4 = "ordered")
    cell_df <- evaluate_imputation_ensemble(
      true_data = rb$complete_data,
      imp_list  = bace_res$imputed_datasets,
      miss_mask = rb$miss_mask,
      var_types = var_types_eval
    )

    # ------------------------------------------------------------------
    # 3. Oracle MCMCglmm on complete_data (same formula/family/budget)
    # ------------------------------------------------------------------
    t_oracle <- Sys.time()
    cd <- rb$complete_data
    cor_mat <- ape::vcv(rb$tree, corr = TRUE)
    Ainv    <- solve(cor_mat)

    prior_oracle <- list(
      R = list(V = 1, nu = 0.002),
      G = list(G1 = list(V = 1, nu = 0.002))
    )
    oracle_fit <- tryCatch(
      MCMCglmm::MCMCglmm(
        fixed   = y ~ x1 + x2 + x3 + x4,
        random  = ~ Species,
        ginverse = list(Species = as(Ainv, "dgCMatrix")),
        family  = "gaussian",
        data    = cd,
        prior   = prior_oracle,
        nitt    = mcmc$nitt, burnin = mcmc$burnin, thin = mcmc$thin,
        verbose = FALSE
      ),
      error = function(e) {
        message("Oracle fit failed for ", ds_name, " rep ", rep_id, ": ", e$message)
        NULL
      }
    )
    oracle_runtime <- as.numeric(difftime(Sys.time(), t_oracle, units = "secs"))

    # ------------------------------------------------------------------
    # 4. Beta-coverage metrics (y equation only)
    #
    # BACE fits MCMCglmm with pr = TRUE, so Sol's columns are:
    #   [fixed effects ... | random-effect BLUPs (Species.XXXXX)]
    # We only want the fixed effects for beta-coverage. The split point
    # is model$Fixed$nfl when available; otherwise drop "Species." cols.
    # ------------------------------------------------------------------
    fixef_sol <- function(model) {
      sol <- model$Sol
      nfl <- tryCatch(model$Fixed$nfl, error = function(e) NULL)
      if (!is.null(nfl) && is.finite(nfl) && nfl > 0L && nfl <= ncol(sol)) {
        return(sol[, seq_len(nfl), drop = FALSE])
      }
      keep <- !startsWith(colnames(sol), "Species.")
      sol[, keep, drop = FALSE]
    }
    summarise_sol <- function(sol) {
      data.frame(
        coef  = colnames(sol),
        mean  = apply(sol, 2, mean),
        sd    = apply(sol, 2, sd),
        lo025 = apply(sol, 2, quantile, probs = 0.025),
        hi975 = apply(sol, 2, quantile, probs = 0.975),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
    bace_summary <- summarise_sol(fixef_sol(bace_res$pooled_models$models$y))
    oracle_summary <- if (!is.null(oracle_fit)) {
      summarise_sol(fixef_sol(oracle_fit))
    } else {
      data.frame(coef = character(0), mean = numeric(0), sd = numeric(0),
                 lo025 = numeric(0), hi975 = numeric(0),
                 stringsAsFactors = FALSE)
    }

    # Merge: BACE vs oracle (coverage of oracle mean by BACE's 95% CI).
    beta_oracle <- merge(bace_summary, oracle_summary,
                         by = "coef", suffixes = c("_bace", "_oracle"),
                         all = TRUE)
    beta_oracle$bias_vs_oracle <- beta_oracle$mean_bace - beta_oracle$mean_oracle
    beta_oracle$cover_oracle   <- with(beta_oracle,
      mean_oracle >= lo025_bace & mean_oracle <= hi975_bace)
    beta_oracle$ci_width_ratio <- with(beta_oracle,
      (hi975_bace - lo025_bace) / (hi975_oracle - lo025_oracle))

    # Merge: BACE vs dialled-in true coefficient.
    # The dialled vector matches BACE's column names only for the
    # coefficients with unambiguous coding:
    #   (Intercept), x2 dummies (treatment-coded), x3 (numeric).
    # x1 and x4 are mapped through R's contr.poly when they enter as
    # ordered factors (sim_bace returns them ordered, factor() inherits),
    # so the BACE column is x1.L / x4.L / x4.Q. Matching these to a
    # treatment-style dialled value requires a sqrt(K^2-1)/sqrt(K)-style
    # rescaling that we deliberately defer to the oracle comparator
    # rather than hardcode here.
    true_full <- meta$spec$beta_resp
    dial_map <- list(
      "(Intercept)" = meta$true_intercepts$response,
      "x11"         = true_full$x1,       # only if x1 ended up unordered factor
      "x2B"         = true_full$x2[1],
      "x2C"         = true_full$x2[2],
      "x3"          = true_full$x3
    )
    beta_oracle$true_dialled <- vapply(beta_oracle$coef, function(nm) {
      v <- dial_map[[nm]]
      if (is.null(v)) NA_real_ else as.numeric(v)
    }, numeric(1))
    beta_oracle$cover_dialled <- with(beta_oracle,
      !is.na(true_dialled) & true_dialled >= lo025_bace & true_dialled <= hi975_bace)
    beta_oracle$bias_vs_dialled <- with(beta_oracle, mean_bace - true_dialled)

    result$bace_runtime    <- bace_runtime
    result$oracle_runtime  <- oracle_runtime
    result$bace_converged  <- bace_res$converged
    result$bace_attempts   <- bace_res$n_attempts
    result$cell_metrics    <- cell_df
    result$bace_y_summary  <- bace_summary
    result$oracle_y_summary <- oracle_summary
    result$beta_compare    <- beta_oracle
    result$status          <- "ok"
    result$completed_at    <- Sys.time()
    TRUE
  }, error = function(e) {
    result$status <<- "error"
    result$error  <<- conditionMessage(e)
    result$completed_at <<- Sys.time()
    FALSE
  })

  saveRDS(result, out_file)
  invisible(result)
}

# -----------------------------------------------------------------------------
# Main: build the (dataset, rep_id) work list and run in parallel.
# -----------------------------------------------------------------------------
build_worklist <- function() {
  do.call(rbind, lapply(DATASETS, function(ds) {
    reps <- list.files(file.path(REF_ROOT, ds), pattern = "^rep_\\d+\\.rds$")
    if (length(reps) == 0L) return(NULL)
    rep_ids <- as.integer(sub("^rep_(\\d+)\\.rds$", "\\1", reps))
    data.frame(dataset = ds, rep_id = sort(rep_ids),
               stringsAsFactors = FALSE)
  }))
}

# CLI:
#   Rscript dev/10_*.R                    # all reps, all datasets
#   Rscript dev/10_*.R sim_hard           # all reps of sim_hard (CI uses this)
#   Rscript dev/10_*.R sim_hard 1         # one rep, smoke test
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L) {
  worklist <- data.frame(dataset = args[[1]],
                         rep_id  = as.integer(args[[2]]),
                         stringsAsFactors = FALSE)
} else if (length(args) == 1L) {
  reps <- list.files(file.path(REF_ROOT, args[[1]]),
                     pattern = "^rep_\\d+\\.rds$")
  if (length(reps) == 0L)
    stop("No reps found for dataset: ", args[[1]])
  rep_ids <- sort(as.integer(sub("^rep_(\\d+)\\.rds$", "\\1", reps)))
  worklist <- data.frame(dataset = args[[1]], rep_id = rep_ids,
                         stringsAsFactors = FALSE)
} else {
  worklist <- build_worklist()
}

cat("Worklist:", nrow(worklist), "rep(s).\n")
cat("Cores   :", N_CORES_OUTER, "\n")
cat("Budget  : nitt=", MCMC$nitt,
    " burnin=", MCMC$burnin,
    " thin=",   MCMC$thin,
    " runs=",   MCMC$runs,
    " n_final=", MCMC$n_final, "\n", sep = "")
cat("Started :", format(Sys.time()), "\n\n")

t_all <- Sys.time()

if (nrow(worklist) == 1L) {
  # Single rep — run in-process for easier debugging during smoke test.
  ev <- evaluate_one_rep(worklist$dataset[1], worklist$rep_id[1])
  cat("\n=== Smoke-test result ===\n")
  cat("Status:", ev$status, "\n")
  if (ev$status == "ok") {
    cat("BACE runtime  :", round(ev$bace_runtime, 1), "s\n")
    cat("Oracle runtime:", round(ev$oracle_runtime, 1), "s\n")
    cat("BACE converged:", ev$bace_converged,
        "  (attempts =", ev$bace_attempts, ")\n")
    cat("\nCell metrics:\n"); print(ev$cell_metrics, row.names = FALSE)
    cat("\nBeta comparison (y equation):\n"); print(ev$beta_compare, row.names = FALSE)
  } else {
    cat("Error:", ev$error, "\n")
  }
} else {
  res <- mclapply(seq_len(nrow(worklist)), function(i) {
    ds  <- worklist$dataset[i]
    rid <- worklist$rep_id[i]
    t0 <- Sys.time()
    out <- tryCatch(
      evaluate_one_rep(ds, rid),
      error = function(e) list(status = "outer_error", error = conditionMessage(e))
    )
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%s] %s rep %02d -> %s (%.1f min)\n",
                format(Sys.time(), "%H:%M:%S"),
                ds, rid, out$status, elapsed))
    out$status
  }, mc.cores = N_CORES_OUTER)

  cat("\n=== Status counts ===\n")
  print(table(unlist(res)))
}

cat("\nTotal wallclock:",
    round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1),
    "min\n")
