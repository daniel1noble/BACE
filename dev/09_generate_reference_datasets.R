# =============================================================================
# 09_generate_reference_datasets.R
#   Generate the four reference simulated datasets specified in
#   notes/simulated_benchmark_datasets.qmd. Each dataset is a spec under
#   which N_REPS replicates are simulated. The dataset-level "truth"
#   (beta_matrix, beta_resp, intercepts, phylo signal vector, sigmas) is
#   FIXED across replicates so that posterior coverage can be aggregated
#   against a single ground truth. The replicate-level objects (tree,
#   complete_data, miss_mask, random_effects, residuals) vary.
#
#   Output layout
#     dev/simulation_results/reference_datasets/
#       sim_ideal/
#         meta.rds                  # dataset-level truth (shared by all reps)
#         rep_01.rds, rep_02.rds...
#       sim_typical/
#       sim_heterogeneous/
#       sim_hard/
#
#   This script only GENERATES data — it does not run bace(). Evaluation
#   is a separate step.
# =============================================================================

devtools::load_all(quiet = TRUE)
library(ape)
library(MASS)

# -----------------------------------------------------------------------------
# Source inject_missingness() and friends from 02_benchmark_simulated_full.R
# (same pattern used by 01_benchmark_simulated.R). Pulls the MAR/MNAR helpers
# without triggering the main benchmark loop.
# -----------------------------------------------------------------------------
source_funcs <- function(path, fn_names) {
  src <- readLines(path)
  for (fn in fn_names) {
    start_i <- grep(paste0("^", fn, " <- function"), src)[1]
    if (is.na(start_i)) stop("Function not found in ", path, ": ", fn)
    depth <- 0
    end_i <- NA
    for (i in start_i:length(src)) {
      depth <- depth +
        sum(gregexpr("\\{", src[i])[[1]] > 0) -
        sum(gregexpr("\\}", src[i])[[1]] > 0)
      if (i > start_i && depth == 0) {
        end_i <- i
        break
      }
    }
    eval(parse(text = src[start_i:end_i]), envir = globalenv())
  }
}

# Pull the trait-to-numeric coercion and intercept-calibration helpers from
# 02_benchmark_simulated_full.R. We do NOT source inject_missingness from
# there: its species-coverage safeguard un-masks every row when n_cases ==
# n_species (one obs/species), which is the regime all four reference
# datasets use. BACE's chained-equation fitter fills placeholders for NA
# predictors and never drops rows, so the safeguard is unnecessary in this
# pipeline. We re-implement the mechanism logic locally without it.
source_funcs(
  file.path("dev", "02_benchmark_simulated_full.R"),
  c(".trait_to_numeric", ".calibrate_intercept")
)

#' Inject missingness under a specified mechanism, without the
#' species-coverage safeguard. Matches the mechanism semantics of
#' 02_benchmark_simulated_full.R::inject_missingness() otherwise.
#' DEP_STRENGTH is taken as an argument here (not from globalenv) so the
#' caller has explicit control per dataset.
inject_missingness_no_safeguard <- function(complete_data, tree, mechanism,
                                            rate, dep_strength,
                                            vars = c("y","x1","x2","x3","x4")) {
  n <- nrow(complete_data)
  miss_mask  <- as.data.frame(matrix(FALSE, nrow = n, ncol = length(vars),
                                     dimnames = list(NULL, vars)))
  miss_data  <- complete_data
  miss_probs <- as.data.frame(matrix(NA_real_, nrow = n, ncol = length(vars),
                                     dimnames = list(NULL, vars)))

  Sigma_phylo <- if (mechanism == "phylo_MAR") ape::vcv(tree, corr = TRUE) else NULL

  trait_mar_covariate <- c(y = "x3", x1 = "y", x2 = "y",
                           x3 = "x1", x4 = "y")

  # Robust standardiser: zero-variance inputs would otherwise produce
  # all-NA linpreds and crash .calibrate_intercept. Fall back to zero
  # (i.e., MCAR for that variable).
  safe_scale <- function(val) {
    v_sd <- stats::sd(val, na.rm = TRUE)
    if (!is.finite(v_sd) || v_sd == 0) return(rep(0, length(val)))
    as.numeric(scale(val))
  }

  for (v in vars) {
    linpred <- switch(mechanism,
      MCAR = rep(0, n),
      phylo_MAR = {
        z_sp <- MASS::mvrnorm(1, mu = rep(0, nrow(Sigma_phylo)),
                              Sigma = Sigma_phylo)
        names(z_sp) <- rownames(Sigma_phylo)
        -dep_strength * z_sp[as.character(complete_data$Species)]
      },
      trait_MAR = {
        cov_v <- trait_mar_covariate[[v]]
        val   <- .trait_to_numeric(complete_data[[cov_v]])
        -dep_strength * safe_scale(val)
      },
      trait_MNAR = {
        val <- .trait_to_numeric(complete_data[[v]])
        -dep_strength * safe_scale(val)
      },
      stop("Unknown mechanism: ", mechanism)
    )

    c_hat <- .calibrate_intercept(linpred, rate)
    p     <- plogis(c_hat + linpred)
    miss  <- rbinom(n, 1, p) == 1

    miss_mask[[v]]  <- miss
    miss_data[[v]][miss] <- NA
    miss_probs[[v]] <- p
  }

  list(miss_data = miss_data, miss_mask = miss_mask, miss_probs = miss_probs)
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

#' Build a fixed 4x4 beta_matrix for the predictor system. Predictors are
#' x1 (binary), x2 (multinomial3), x3 (poisson), x4 (threshold3).
#' Entry (i, j) = effect of predictor j on predictor i; diagonal = 0.
#' Off-diagonal entries are zeroed with prob = sparsity; non-zero values
#' are drawn uniformly with random sign from the specified magnitude range.
#'
#' Constraint: sim_bace's design_size() requires a DAG over predictors so
#' it can simulate them in topological order. We fix the order
#' x1 -> x2 -> x3 -> x4 and only allow strictly-lower-triangular entries
#' (i.e., predictor i can depend only on x1, ..., x{i-1}). This makes x1
#' the root and gives every dataset the same dependency *direction*;
#' density and magnitudes still vary by spec.
make_beta_matrix <- function(sparsity, mag_range, seed) {
  set.seed(seed)
  p <- 4
  pred_names <- c("x1", "x2", "x3", "x4")
  M <- matrix(0, nrow = p, ncol = p,
              dimnames = list(pred_names, pred_names))
  for (i in seq_len(p)) {
    if (i == 1) next                       # root: no incoming edges
    for (j in seq_len(i - 1)) {             # strict lower triangle only
      if (runif(1) >= sparsity) {
        M[i, j] <- runif(1, mag_range[1], mag_range[2]) * sample(c(-1, 1), 1)
      }
    }
  }
  M
}

#' Convert sim_bace output column types so that .detect_type() in BACE
#' assigns the right family:
#'   x1 (integer 0/1) -> factor with 2 levels -> binary -> threshold (K=2)
#'   x2 (character A/B/C) -> factor with 3 levels -> categorical
#'   x4 (integer 1..3) -> ordered factor 1<2<3 -> threshold (ordinal)
post_process_types <- function(df) {
  df$x1 <- factor(df$x1, levels = c(0, 1))
  df$x2 <- factor(df$x2, levels = c("A", "B", "C"))
  df$x4 <- factor(df$x4, levels = c(1, 2, 3), ordered = TRUE)
  df
}

# -----------------------------------------------------------------------------
# Dataset specifications (mirror notes/simulated_benchmark_datasets.qmd)
# -----------------------------------------------------------------------------

SPECS <- list(
  sim_ideal = list(
    n_cases       = 750,
    n_species     = 750,
    phylo_signal  = c(y = 0.85, x1 = 0.85, x2 = 0.85, x3 = 0.85, x4 = 0.85),
    beta_sparsity = 0.2,
    mag_range     = c(0.3, 0.7),
    beta_resp     = list(x1 = 0.5, x2 = c(0.4, -0.2), x3 = 0.3, x4 = -0.3),
    intercepts    = list(predictors = c(0, 0, 0.5, 0), response = 0),
    mechanism     = "MCAR",
    rate          = 0.15,
    dep_strength  = NA_real_,
    n_reps        = 30,
    base_seed     = 2026
  ),
  sim_typical = list(
    n_cases       = 500,
    n_species     = 500,
    phylo_signal  = c(y = 0.60, x1 = 0.60, x2 = 0.60, x3 = 0.60, x4 = 0.60),
    beta_sparsity = 0.4,
    mag_range     = c(0.2, 0.5),
    beta_resp     = list(x1 = 0.4, x2 = c(0.3, -0.1), x3 = 0.2, x4 = -0.2),
    intercepts    = list(predictors = c(0, 0, 0.5, 0), response = 0),
    mechanism     = "MCAR",
    rate          = 0.30,
    dep_strength  = NA_real_,
    n_reps        = 50,
    base_seed     = 2027
  ),
  sim_heterogeneous = list(
    n_cases       = 600,
    n_species     = 600,
    phylo_signal  = c(y = 0.85, x1 = 0.20, x2 = 0.60, x3 = 0.85, x4 = 0.20),
    beta_sparsity = 0.4,
    mag_range     = c(0.2, 0.5),
    beta_resp     = list(x1 = 0.5, x2 = c(0.3, -0.15), x3 = 0.3, x4 = -0.2),
    intercepts    = list(predictors = c(0, 0, 0.5, 0), response = 0),
    mechanism     = "trait_MAR",
    rate          = 0.35,
    dep_strength  = 1.5,
    n_reps        = 50,
    base_seed     = 2028
  ),
  sim_hard = list(
    n_cases       = 300,
    n_species     = 300,
    phylo_signal  = c(y = 0.20, x1 = 0.20, x2 = 0.20, x3 = 0.20, x4 = 0.20),
    beta_sparsity = 0.7,
    mag_range     = c(0.1, 0.3),
    beta_resp     = list(x1 = 0.3, x2 = c(0.2, -0.1), x3 = 0.1, x4 = -0.1),
    intercepts    = list(predictors = c(0, 0, 0.5, 0), response = 0),
    mechanism     = "trait_MNAR",
    rate          = 0.50,
    dep_strength  = 2.5,
    n_reps        = 50,
    base_seed     = 2029
  )
)

OUT_ROOT <- file.path("dev", "simulation_results", "reference_datasets")
if (!dir.exists(OUT_ROOT)) dir.create(OUT_ROOT, recursive = TRUE)

# -----------------------------------------------------------------------------
# Main generation loop
# -----------------------------------------------------------------------------
overall_t0 <- Sys.time()
summary_rows <- list()

for (ds_name in names(SPECS)) {

  spec <- SPECS[[ds_name]]
  ds_dir <- file.path(OUT_ROOT, ds_name)
  if (!dir.exists(ds_dir)) dir.create(ds_dir, recursive = TRUE)

  cat(sprintf("\n=== %s ===\n", ds_name))
  cat(sprintf("  n_cases=%d  n_species=%d  rate=%.2f  mech=%s  reps=%d\n",
              spec$n_cases, spec$n_species, spec$rate,
              spec$mechanism, spec$n_reps))

  # Resume: skip dataset if all rep files + meta already exist.
  rep_files <- list.files(ds_dir, pattern = "^rep_\\d+\\.rds$",
                          full.names = TRUE)
  if (file.exists(file.path(ds_dir, "meta.rds")) &&
      length(rep_files) >= spec$n_reps) {
    cat(sprintf("  [skip] %d reps already present\n", length(rep_files)))
    # Rebuild summary rows from on-disk files so the final report covers
    # every dataset, not just the freshly generated ones.
    for (rf in rep_files) {
      rb <- readRDS(rf)
      obs_rates <- sapply(c("y", "x1", "x2", "x3", "x4"),
                         function(v) mean(rb$miss_mask[[v]]))
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        dataset     = ds_name,
        rep_id      = rb$rep_id,
        seed        = rb$seed,
        mechanism   = rb$mechanism,
        target_rate = rb$rate,
        y_rate      = obs_rates["y"],
        x1_rate     = obs_rates["x1"],
        x2_rate     = obs_rates["x2"],
        x3_rate     = obs_rates["x3"],
        x4_rate     = obs_rates["x4"]
      )
    }
    next
  }

  # --- dataset-level truth (fixed across replicates) ----------------------
  beta_matrix <- make_beta_matrix(
    sparsity   = spec$beta_sparsity,
    mag_range  = spec$mag_range,
    seed       = spec$base_seed
  )

  meta <- list(
    dataset             = ds_name,
    spec                = spec,
    true_beta_matrix    = beta_matrix,
    true_beta_resp      = spec$beta_resp,
    true_intercepts     = spec$intercepts,
    true_phylo_signal   = spec$phylo_signal,
    skeleton            = list(
      response_type   = "gaussian",
      predictor_types = c("binary", "multinomial3", "poisson", "threshold3"),
      var_names       = c("y", "x1", "x2", "x3", "x4")
    ),
    generated_at        = Sys.time()
  )
  saveRDS(meta, file.path(ds_dir, "meta.rds"))

  dep_strength <- if (is.na(spec$dep_strength)) 1.5 else spec$dep_strength

  ds_t0 <- Sys.time()
  for (rep_id in seq_len(spec$n_reps)) {
    rep_seed <- spec$base_seed + 1000L * rep_id
    set.seed(rep_seed)

    # Always simulate complete data first (no NAs from sim_bace), then
    # apply the chosen mechanism in one shared path.
    sim <- suppressMessages(sim_bace(
      response_type   = "gaussian",
      predictor_types = c("binary", "multinomial3", "poisson", "threshold3"),
      var_names       = c("y", "x1", "x2", "x3", "x4"),
      beta_matrix     = beta_matrix,
      beta_resp       = spec$beta_resp,
      intercepts      = spec$intercepts,
      phylo_signal    = spec$phylo_signal,
      n_cases         = spec$n_cases,
      n_species       = spec$n_species,
      missingness     = rep(0, 5)
    ))

    complete_data <- post_process_types(sim$complete_data)

    mi <- inject_missingness_no_safeguard(
      complete_data = complete_data,
      tree          = sim$tree,
      mechanism     = spec$mechanism,
      rate          = spec$rate,
      dep_strength  = dep_strength
    )
    miss_data  <- mi$miss_data
    miss_mask  <- mi$miss_mask
    miss_probs <- mi$miss_probs

    rep_bundle <- list(
      complete_data        = complete_data,
      miss_data            = miss_data,
      miss_mask            = miss_mask,
      miss_probs           = miss_probs,
      tree                 = sim$tree,
      true_random_effects  = sim$random_effects,
      true_sigmas          = sim$params$sigmas,
      true_beta_resp_full  = sim$params$beta_resp_full,
      seed                 = rep_seed,
      rep_id               = rep_id,
      mechanism            = spec$mechanism,
      rate                 = spec$rate
    )

    saveRDS(rep_bundle,
            file.path(ds_dir, sprintf("rep_%02d.rds", rep_id)))

    # Audit: realised per-variable missingness rates this rep
    obs_rates <- sapply(c("y", "x1", "x2", "x3", "x4"),
                       function(v) mean(miss_mask[[v]]))
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      dataset   = ds_name,
      rep_id    = rep_id,
      seed      = rep_seed,
      mechanism = spec$mechanism,
      target_rate = spec$rate,
      y_rate    = obs_rates["y"],
      x1_rate   = obs_rates["x1"],
      x2_rate   = obs_rates["x2"],
      x3_rate   = obs_rates["x3"],
      x4_rate   = obs_rates["x4"]
    )
  }

  ds_elapsed <- difftime(Sys.time(), ds_t0, units = "secs")
  cat(sprintf("  %d reps generated in %.1fs\n",
              spec$n_reps, as.numeric(ds_elapsed)))
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
summary_df <- do.call(rbind, summary_rows)
saveRDS(summary_df, file.path(OUT_ROOT, "generation_summary.rds"))
utils::write.csv(summary_df,
                 file.path(OUT_ROOT, "generation_summary.csv"),
                 row.names = FALSE)

cat("\n=== Per-dataset realised missingness rates (mean across reps) ===\n")
agg <- aggregate(
  cbind(y_rate, x1_rate, x2_rate, x3_rate, x4_rate) ~ dataset + mechanism + target_rate,
  data = summary_df,
  FUN  = function(z) round(mean(z), 3)
)
print(agg, row.names = FALSE)

cat(sprintf("\nTotal wallclock: %.1fs\n",
            as.numeric(difftime(Sys.time(), overall_t0, units = "secs"))))
cat("Output root:", OUT_ROOT, "\n")
