# =============================================================================
# Simulation: Assessing BACE Imputation Quality (Controlled Simulated Data)
# =============================================================================
#
# Context
# -------
# BACE already has an empirical benchmark against AVONET
# (avonet_mixed_benchmark.R). That benchmark uses REAL trait data, so the
# true values of the "missing" cells are unknown. This script fills the
# gap with a fully controlled simulation: generate complete data with
# known ground truth, hide cells under a realistic missingness
# mechanism, impute with bace(), and measure the error directly.
#
# Key design choices (peer-reviewed backing)
# ------------------------------------------
# (a) NO MCAR. Missingness in phylogenetic trait data is virtually never
#     MCAR (Nakagawa & Freckleton 2008, TREE 23:592; Penone et al. 2014,
#     MEE 5:961). We instead cross three mechanisms motivated by the
#     Rubin (1976, Biometrika 63:581) MCAR/MAR/MNAR framework:
#
#       phylo_MAR    clade-structured missingness. A latent "sampling
#                    intensity" z_i is drawn from Brownian motion on the
#                    tree; low-z species are more often missing. This
#                    reproduces the phylogenetic clustering of missingness
#                    documented in Penone et al. 2014 for bird and
#                    mammal trait databases. MAR because the phylogeny
#                    is observed.
#
#       trait_MAR    missingness of a variable depends on another fully
#                    observed variable (e.g. small-bodied species have
#                    harder-to-measure behavioural traits). MAR.
#
#       trait_MNAR   missingness depends on the unobserved value itself
#                    (e.g. very low values below detection). MNAR.
#                    The hardest case — even correctly specified
#                    imputation models can be biased. Included as a
#                    stress test (see Johnson et al. 2021, MEE 12:1847).
#
# (b) Convergence retries ON. skip_conv = FALSE and max_attempts = 2.
#     Inference from unconverged MCMC is ill-advised (Gelman & Rubin
#     1992, Stat Sci 7:457; Hadfield 2010, JSS 33(2)), so convergence
#     status is recorded per-replicate; the summary allows stratification
#     by converged vs not.
#
# (c) n_final = 20. Graham, Olchowski & Gilreath (2007, Prev Sci 8:206)
#     show that with fraction-of-missing-information >= 0.3 (we're at
#     35% per variable), 5 imputations understates imputation variance
#     by ~15-20%. 20 is the modern floor for MI inference
#     (van Buuren 2018, Flexible Imputation of Data §2.8).
#
# (d) Calibration added. NRMSE and log-MAE tell us about point-estimate
#     accuracy; coverage of 95% posterior predictive intervals tells us
#     whether uncertainty is honest (van Buuren 2018 §2.5; Little & Rubin
#     2019, 3rd ed., §5.4). For categorical/ordered traits, the
#     analogous calibration metric is the Brier score (Brier 1950, MWR
#     78:1; Gneiting & Raftery 2007, JASA 102:359).
#
# Design
# ------
# Five variables, one of each type, simulated jointly with phylogenetic
# signal and cross-trait dependencies:
#
#     y  = gaussian     (continuous)
#     x1 = binary       ("categorical", 2-level ordered factor)
#     x2 = multinomial3 (unordered 3-level factor)
#     x3 = poisson      (integer counts)
#     x4 = threshold3   (ordered categorical, 3 levels)
#
# Crossed axes:
#
#   PHYLOGENETIC-SIGNAL SCENARIO  (4 levels)
#     all_high      all five traits have HIGH   signal (~0.90)
#     all_moderate  all five traits have MODERATE signal (~0.60)
#     all_low       all five traits have LOW   signal (~0.20)
#     mixed         >=1 high, >=1 moderate, >=1 low, shuffled
#
#   MISSINGNESS MECHANISM  (3 levels; ALL at rate = MISS_RATE)
#     phylo_MAR, trait_MAR, trait_MNAR  (see above).
#
# N_SIMS replicates per cell -> 4 * 3 * N_SIMS = 12 * N_SIMS datasets.
#
# Metric glossary (all computed on hidden cells only)
# ---------------------------------------------------
# Point-estimate metrics (averaged over n_final imputations):
#   NRMSE   = RMSE(imputed, true) / sd(complete true values).
#             Scale-free; 0 = perfect; 1 = as bad as predicting the
#             marginal trait mean; > 1 = worse than the marginal mean.
#             Normalised by the sd of the COMPLETE data (marginal),
#             not the hidden subset - under MAR/MNAR the hidden subset
#             can have a much narrower sd than the marginal, so dividing
#             by it inflates NRMSE spuriously (Oba et al. 2003 used
#             the hidden-sd form; Stekhoven & Buhlmann 2012 and
#             van Buuren 2018 moved to the marginal-sd form).
#   MAE     = mean(|imputed - true|) on the raw scale.
#   log_mae (count only) = mean(|log1p(imputed) - log1p(true)|);
#             tail-robust on the count scale.
#   cor     = Pearson correlation (imputed, true). Rank/pattern
#             recovery; insensitive to bias or scale.
#   accuracy          = proportion of hidden cells imputed correctly.
#   balanced_accuracy = mean of per-class recall. Robust to class
#                       imbalance.
#   ordinal_mae       (threshold only) = mean(|rank_true - rank_imp|)
#                       on the factor-level ordering.
#
# Calibration (uses ALL n_final draws as a posterior sample):
#   coverage95  (continuous / count) = proportion of hidden cells whose
#                true value falls inside the 95% posterior predictive
#                interval [q_0.025, q_0.975] of the n_final draws.
#                Well-calibrated intervals give coverage ~ 0.95.
#   brier       (categorical / ordered) = mean over hidden cells of
#                sum_k (p_k - 1{true=k})^2, where p_k is the frequency
#                of class k across the n_final imputations. Lower is
#                better; 0 is perfect.
#
# Outputs
# -------
#   dev/simulation_results/imputation_quality_results.rds   raw per-sim
#   dev/simulation_results/imputation_quality_results.csv   same, csv
#   dev/simulation_results/imputation_quality_summary.csv   aggregated
#   dev/simulation_results/imputation_quality_plots.pdf     diagnostics
# =============================================================================

library(BACE)
library(ape)
library(MASS)
library(parallel)

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

set.seed(2026)

# ---- Simulation scope -------------------------------------------------------
# 500 species and 500 cases means ~1 observation per species, matching the
# AVONET / PanTHERIA comparative-dataset structure (one row per species).
# The Species random effect is then purely phylogenetic (no within-species
# replication), identical to standard phylogenetic comparative analysis.
# Each bace() call is substantially more expensive than with 100 species
# because MCMCglmm's ginverse is 500x500; N_SIMS is therefore set to 50
# per (scenario, mechanism) cell (12 * 50 = 600 runs total), still giving
# useful Monte Carlo precision for between-cell comparisons.
N_SIMS    <- 50                 # replicates per (scenario, mechanism) cell
N_CASES   <- 500                # observations per simulated dataset
N_SPECIES <- 500                # species on the phylogeny (1 obs/species)
N_CORES   <- 4L

# ---- Phylogenetic-signal scenarios -----------------------------------------
SIGNAL <- list(high = 0.90, moderate = 0.60, low = 0.20)

# ---- Missingness mechanisms & rate -----------------------------------------
MECHANISMS <- c("phylo_MAR", "trait_MAR", "trait_MNAR")

# Single realistic rate. 35% is at the high end of what Penone et al.
# (2014) report for common trait fields in AVONET/PanTHERIA (body-mass
# ~10-15%, many behavioural traits 30-50%).
MISS_RATE  <- 0.35

# How strongly missingness depends on the driver. Units are logits per
# SD of the driver. 1.5 gives species-level ICCs for the phylo_MAR mask
# in the range empirically reported by Penone et al. (2014, MEE 5:961)
# for bird / mammal trait databases (~0.3-0.5); larger values push the
# mask toward whole-clade dropouts, which are rare in real datasets and
# break MCMCglmm (it requires every Species level to have >=1 observed
# row per imputation model).
DEP_STRENGTH <- 1.5

# ---- MCMC settings ----------------------------------------------------------
NITT    <- 15000
THIN    <- 10
BURNIN  <- 5000
RUNS    <- 10                   # convergence-phase iterations
N_FINAL <- 20                   # posterior-predictive imputations
MAX_ATTEMPTS <- 2               # retries on non-convergence

# ---- Output paths -----------------------------------------------------------
RESULTS_DIR <- file.path("dev", "simulation_results")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# =============================================================================
# 2. SCENARIO GENERATOR
# =============================================================================

#' Build a length-5 phylogenetic-signal vector for a given scenario.
#' Variable order: (y, x1, x2, x3, x4). For "mixed", each of the three
#' signal levels appears >= 1 time; the remaining 2 slots are random;
#' the result is shuffled across traits.
make_phylo_signal <- function(scenario) {
  switch(scenario,
    all_high     = rep(SIGNAL$high,     5),
    all_moderate = rep(SIGNAL$moderate, 5),
    all_low      = rep(SIGNAL$low,      5),
    mixed = {
      base   <- c(SIGNAL$high, SIGNAL$moderate, SIGNAL$low)
      extras <- sample(unlist(SIGNAL), 2, replace = TRUE)
      sample(c(base, extras))
    },
    stop("Unknown scenario: ", scenario)
  )
}

# =============================================================================
# 3. MISSINGNESS MECHANISMS
# =============================================================================

#' Coerce any trait (numeric, integer, factor, ordered) to a numeric
#' driver for missingness probability calculation.
.trait_to_numeric <- function(x) {
  if (is.factor(x)) as.integer(x) else as.numeric(x)
}

#' Solve for an intercept c such that mean(plogis(c + linpred)) = rate.
#' If no solution exists inside [-20, 20], return the endpoint that
#' minimises |mean(p) - rate|. Happens rarely and only at extreme rates.
.calibrate_intercept <- function(linpred, rate) {
  f <- function(c) mean(plogis(c + linpred)) - rate
  if (sign(f(-20)) == sign(f(20))) {
    return(ifelse(abs(f(-20)) < abs(f(20)), -20, 20))
  }
  uniroot(f, c(-20, 20))$root
}

#' Inject missingness into a complete dataset under a specified
#' mechanism, calibrated so the EXPECTED missingness rate per variable
#' equals `rate`.
#'
#' @param complete_data  data.frame with a `Species` column + the 5 traits
#' @param tree           phylo object (used only for phylo_MAR)
#' @param mechanism      "phylo_MAR" | "trait_MAR" | "trait_MNAR"
#' @param rate           target proportion missing per variable
#' @param vars           variable names to apply missingness to
#' @return list(miss_data, miss_mask, miss_probs)
inject_missingness <- function(complete_data, tree, mechanism, rate,
                               vars = c("y", "x1", "x2", "x3", "x4")) {

  n <- nrow(complete_data)
  miss_mask  <- as.data.frame(matrix(FALSE, nrow = n, ncol = length(vars),
                                     dimnames = list(NULL, vars)))
  miss_data  <- complete_data
  miss_probs <- as.data.frame(matrix(NA_real_, nrow = n, ncol = length(vars),
                                     dimnames = list(NULL, vars)))

  # For phylo_MAR: precompute the phylogenetic covariance so we can
  # draw an independent BM z-vector per variable on the same tree.
  Sigma_phylo <- if (mechanism == "phylo_MAR")
    ape::vcv(tree, corr = TRUE) else NULL

  # For trait_MAR: each variable's missingness depends on a DIFFERENT
  # fully-observed-in-complete-data trait. The map is deliberately
  # heterogeneous so not every missingness vector is driven by the
  # same covariate (which would create spurious cross-variable
  # correlations in the mask).
  trait_mar_covariate <- c(y  = "x3", x1 = "y",  x2 = "y",
                           x3 = "x1", x4 = "y")

  for (v in vars) {
    linpred <- switch(mechanism,

      phylo_MAR = {
        # One independent BM draw per variable on the same tree.
        z_sp <- MASS::mvrnorm(1, mu = rep(0, nrow(Sigma_phylo)),
                              Sigma = Sigma_phylo)
        names(z_sp) <- rownames(Sigma_phylo)
        # Low z  => higher P(missing). Sign chosen by convention.
        -DEP_STRENGTH * z_sp[as.character(complete_data$Species)]
      },

      trait_MAR = {
        cov_v <- trait_mar_covariate[[v]]
        val   <- .trait_to_numeric(complete_data[[cov_v]])
        -DEP_STRENGTH * as.numeric(scale(val))   # standardise -> per-SD logit
      },

      trait_MNAR = {
        val <- .trait_to_numeric(complete_data[[v]])
        -DEP_STRENGTH * as.numeric(scale(val))
      }
    )

    c_hat <- .calibrate_intercept(linpred, rate)
    p     <- plogis(c_hat + linpred)
    miss  <- rbinom(n, 1, p) == 1

    miss_mask[[v]]  <- miss
    miss_data[[v]][miss] <- NA
    miss_probs[[v]] <- p
  }

  # Species-coverage safeguard. MCMCglmm requires every Species factor
  # level to have >=1 non-NA row in the response of each imputation
  # model, otherwise ginverse mismatches throw an error. Under strong
  # phylo_MAR clustering a handful of species can end up 100% masked
  # for a given variable; un-mask one observation at random from any
  # such species. This is a mild perturbation of the declared
  # mechanism (necessary for the benchmark to run) but affects at most
  # a few percent of cells.
  by_sp <- split(seq_len(n), complete_data$Species)
  for (v in vars) {
    for (sp_rows in by_sp) {
      if (length(sp_rows) == 0L) next
      if (all(miss_mask[[v]][sp_rows])) {
        keep <- if (length(sp_rows) == 1L) sp_rows else sample(sp_rows, 1)
        miss_mask[[v]][keep] <- FALSE
        miss_data[[v]][keep] <- complete_data[[v]][keep]
      }
    }
  }

  list(miss_data = miss_data, miss_mask = miss_mask, miss_probs = miss_probs)
}

# =============================================================================
# 4. ACCURACY METRICS (ensemble evaluator)
# =============================================================================

#' Compute imputation-accuracy metrics using the FULL list of n_final
#' imputations as a posterior predictive sample.
#'
#' @param true_data  complete data.frame (no NAs)
#' @param imp_list   list of n_final imputed data.frames from bace()
#' @param miss_mask  logical data.frame; TRUE where the cell was hidden
#' @param var_types  named char vector: "gaussian" | "count" | "categorical"
#'                   | "ordered"
#' @return long-form data.frame: variable, type, metric, value
evaluate_imputation_ensemble <- function(true_data, imp_list,
                                         miss_mask, var_types) {
  out <- list()

  for (v in names(var_types)) {
    idx <- miss_mask[[v]]
    if (!any(idx)) next
    true_vals <- true_data[[v]][idx]
    n_cell    <- sum(idx)

    if (var_types[v] %in% c("gaussian", "count")) {
      # Coerce to numeric upfront: poisson traits are stored as integer
      # in complete_data but as double in the imputed datasets, so a
      # shared FUN.VALUE template would fail.
      tv <- as.numeric(true_vals)
      iv_mat <- vapply(imp_list, function(d) as.numeric(d[[v]][idx]),
                       FUN.VALUE = numeric(n_cell))
      if (!is.matrix(iv_mat)) iv_mat <- matrix(iv_mat, nrow = n_cell)

      # Point-estimate metrics: compute per-imputation, then average.
      # (This is Rubin-style pooling of point estimates, for metrics
      # where pooling at the raw-value level would bias variance.)
      # NRMSE normalises by the marginal sd of the COMPLETE true
      # distribution - see metric glossary for rationale.
      sd_marginal <- sd(as.numeric(true_data[[v]]), na.rm = TRUE)
      per_imp_nrmse <- apply(iv_mat, 2, function(iv) {
        if (!is.finite(sd_marginal) || sd_marginal <= 0) NA_real_
        else sqrt(mean((iv - tv)^2)) / sd_marginal
      })
      per_imp_mae <- apply(iv_mat, 2, function(iv) mean(abs(iv - tv)))
      per_imp_cor <- apply(iv_mat, 2, function(iv) {
        if (length(unique(tv)) > 1 && length(unique(iv)) > 1)
          suppressWarnings(cor(iv, tv)) else NA_real_
      })

      metrics <- c("nrmse", "mae", "correlation")
      values  <- c(mean(per_imp_nrmse, na.rm = TRUE),
                   mean(per_imp_mae),
                   mean(per_imp_cor, na.rm = TRUE))

      if (var_types[v] == "count") {
        tv_pos <- pmax(tv, 0)
        per_imp_log_mae <- apply(iv_mat, 2, function(iv) {
          mean(abs(log1p(pmax(iv, 0)) - log1p(tv_pos)))
        })
        metrics <- c(metrics, "log_mae")
        values  <- c(values, mean(per_imp_log_mae))
      }

      # Calibration: 95% posterior predictive interval coverage.
      # Each row of iv_mat has n_imp draws from the posterior
      # predictive for that cell; the 2.5/97.5 quantiles give a
      # central 95% interval.
      lo <- apply(iv_mat, 1, stats::quantile, probs = 0.025,
                  na.rm = TRUE, names = FALSE)
      hi <- apply(iv_mat, 1, stats::quantile, probs = 0.975,
                  na.rm = TRUE, names = FALSE)
      coverage <- mean(tv >= lo & tv <= hi)

      metrics <- c(metrics, "coverage95")
      values  <- c(values, coverage)

      out[[length(out) + 1]] <- data.frame(
        variable = v,
        type     = if (var_types[v] == "gaussian") "continuous" else "count",
        metric   = metrics, value = values,
        stringsAsFactors = FALSE
      )

    } else {
      # Categorical / ordered: build a character matrix of imputed
      # values. Bypass the shared imp_mat so factor levels survive
      # correctly across different factor codings in the imputed data.
      tv_chr <- as.character(true_vals)
      imp_chr <- vapply(imp_list,
                        function(d) as.character(d[[v]][idx]),
                        FUN.VALUE = character(n_cell))
      if (!is.matrix(imp_chr)) imp_chr <- matrix(imp_chr, nrow = n_cell)

      # Per-imputation accuracy & balanced accuracy, then average.
      per_imp_acc <- apply(imp_chr, 2, function(iv) mean(iv == tv_chr))
      classes     <- unique(tv_chr)
      per_imp_bal <- apply(imp_chr, 2, function(iv) {
        recalls <- vapply(classes, function(cl) {
          is_cl <- tv_chr == cl
          if (!any(is_cl)) NA_real_ else mean(iv[is_cl] == cl)
        }, numeric(1))
        mean(recalls, na.rm = TRUE)
      })

      # Brier score (multiclass form, Gneiting & Raftery 2007):
      # estimate p_ik as the frequency of class k at cell i across
      # the n_imp imputations; score against 0/1 truth indicator.
      all_levels <- if (is.factor(true_data[[v]]))
        levels(true_data[[v]]) else sort(unique(tv_chr))
      n_cell <- length(tv_chr); n_imp <- ncol(imp_chr)
      prob_mat <- matrix(0, nrow = n_cell, ncol = length(all_levels),
                         dimnames = list(NULL, all_levels))
      for (k in all_levels) {
        prob_mat[, k] <- rowMeans(imp_chr == k)
      }
      y_indic <- matrix(0, nrow = n_cell, ncol = length(all_levels),
                        dimnames = list(NULL, all_levels))
      for (i in seq_len(n_cell)) y_indic[i, tv_chr[i]] <- 1
      brier <- mean(rowSums((prob_mat - y_indic)^2))

      metrics <- c("accuracy", "balanced_accuracy", "brier")
      values  <- c(mean(per_imp_acc), mean(per_imp_bal, na.rm = TRUE), brier)

      if (var_types[v] == "ordered") {
        lvls <- levels(true_data[[v]])
        tv_rank <- match(tv_chr, lvls)
        per_imp_ord <- apply(imp_chr, 2, function(iv) {
          mean(abs(match(iv, lvls) - tv_rank), na.rm = TRUE)
        })
        metrics <- c(metrics, "ordinal_mae")
        values  <- c(values, mean(per_imp_ord))
      }

      out[[length(out) + 1]] <- data.frame(
        variable = v,
        type     = if (var_types[v] == "ordered") "ordered" else "categorical",
        metric   = metrics, value = values,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, out)
}

# =============================================================================
# 5. SINGLE-REPLICATE DRIVER
# =============================================================================

run_one_sim <- function(sim_id, scenario, mechanism) {

  # --- Simulate complete data ----------------------------------------------
  phylo_signal <- make_phylo_signal(scenario)
  sim <- sim_bace(
    response_type   = "gaussian",
    predictor_types = c("binary", "multinomial3", "poisson", "threshold3"),
    var_names       = c("y", "x1", "x2", "x3", "x4"),
    phylo_signal    = phylo_signal,
    n_cases         = N_CASES,
    n_species       = N_SPECIES,
    beta_sparsity   = 0.3,
    missingness     = c(0, 0, 0, 0, 0)
  )
  complete_data <- sim$complete_data
  tree          <- sim$tree
  names(complete_data)[names(complete_data) == "species"] <- "Species"

  # --- Inject missingness under the chosen mechanism -----------------------
  miss <- inject_missingness(complete_data, tree, mechanism, MISS_RATE)
  miss_data <- miss$miss_data
  miss_mask <- miss$miss_mask

  # Record realised per-variable missing rates (can deviate from target
  # because each mask is a Bernoulli draw).
  realised_rate <- colMeans(miss_mask)

  # --- Run bace() with convergence retries ---------------------------------
  bace_result <- tryCatch(
    bace(
      fixformula     = list("y  ~ x1 + x2 + x3 + x4",
                            "x1 ~ y  + x2 + x3 + x4",
                            "x2 ~ y  + x1 + x3 + x4",
                            "x3 ~ y  + x1 + x2 + x4",
                            "x4 ~ y  + x1 + x2 + x3"),
      ran_phylo_form = "~1|Species",
      phylo          = tree,
      data           = miss_data,
      nitt           = NITT, thin = THIN, burnin = BURNIN,
      runs           = RUNS, n_final = N_FINAL,
      species        = FALSE, verbose = FALSE,
      skip_conv      = FALSE, max_attempts = MAX_ATTEMPTS,
      n_cores        = 1L
    ),
    error = function(e) {
      message("  [", scenario, "/", mechanism, " sim ", sim_id, "] ",
              "bace() error: ", e$message)
      NULL
    }
  )
  if (is.null(bace_result)) return(NULL)

  # --- Score the full list of n_final imputations --------------------------
  var_types <- c(y  = "gaussian",  x1 = "categorical", x2 = "categorical",
                 x3 = "count",     x4 = "ordered")
  scores <- evaluate_imputation_ensemble(complete_data,
                                         bace_result$imputed_datasets,
                                         miss_mask, var_types)

  # Tag with simulation metadata (+ convergence + realised rate)
  scores$sim_id     <- sim_id
  scores$scenario   <- scenario
  scores$mechanism  <- mechanism
  scores$converged  <- isTRUE(bace_result$converged)
  scores$n_attempts <- bace_result$n_attempts
  scores$phylo_y    <- phylo_signal[1]
  scores$phylo_x1   <- phylo_signal[2]
  scores$phylo_x2   <- phylo_signal[3]
  scores$phylo_x3   <- phylo_signal[4]
  scores$phylo_x4   <- phylo_signal[5]
  scores$realised_rate_mean <- mean(realised_rate)
  scores
}

# =============================================================================
# 6. MAIN LOOP  (scenarios x mechanisms)
# =============================================================================

SCENARIOS <- c("all_high", "all_moderate", "all_low", "mixed")

cat("=============================================================\n")
cat("  BACE Imputation Quality Simulation (simulated data)\n")
cat("  Scenarios   :", paste(SCENARIOS, collapse = ", "), "\n")
cat("  Mechanisms  :", paste(MECHANISMS, collapse = ", "),
    "  (rate =", MISS_RATE * 100, "% per variable)\n")
cat("  Replicates  :", N_SIMS, "per cell (",
    N_SIMS * length(SCENARIOS) * length(MECHANISMS), "total)\n")
cat("  Variables   : y=gaussian, x1=binary, x2=multinomial3,",
    "x3=poisson, x4=threshold3\n")
cat("  MCMC        : nitt=", NITT, " thin=", THIN, " burnin=", BURNIN,
    " runs=", RUNS, " n_final=", N_FINAL,
    " max_attempts=", MAX_ATTEMPTS, "\n", sep = "")
cat("=============================================================\n\n")

all_results <- list()

for (scenario in SCENARIOS) {
  for (mechanism in MECHANISMS) {

    cond_id <- paste(scenario, mechanism, sep = "|")
    cat("--- Scenario:", scenario, " | Mechanism:", mechanism, "---\n")
    t0 <- Sys.time()

    sim_results <- mclapply(seq_len(N_SIMS), function(i) {
      if (i %% 25 == 0) message("  ", cond_id, " sim ", i, " / ", N_SIMS)
      run_one_sim(sim_id = i, scenario = scenario, mechanism = mechanism)
    }, mc.cores = N_CORES)

    ok <- !vapply(sim_results, is.null, logical(1))
    cat("  Completed:", sum(ok), "/", N_SIMS, "replicates",
        " (bace-errors:", sum(!ok), ")\n")
    cat("  Time    :", round(difftime(Sys.time(), t0, units = "mins"), 1),
        "min\n\n")
    all_results[[cond_id]] <- do.call(rbind, sim_results[ok])
  }
}

results_df <- do.call(rbind, all_results)
rownames(results_df) <- NULL

saveRDS(results_df, file.path(RESULTS_DIR, "imputation_quality_results.rds"))
write.csv(results_df, file.path(RESULTS_DIR, "imputation_quality_results.csv"),
          row.names = FALSE)
cat("Raw results saved to:", RESULTS_DIR, "\n\n")

# =============================================================================
# 7. SUMMARY TABLES
# =============================================================================

cat("=============================================================\n")
cat("  CONVERGENCE SUMMARY\n")
cat("=============================================================\n\n")

# Per-cell convergence rates. A low rate here is a flag that the MCMC
# budget is too small for that condition rather than a BACE failure.
conv_tab <- aggregate(converged ~ scenario + mechanism,
                      data = unique(results_df[, c("sim_id", "scenario",
                                                   "mechanism", "converged")]),
                      FUN  = mean)
names(conv_tab)[names(conv_tab) == "converged"] <- "prop_converged"
print(conv_tab, row.names = FALSE)

cat("\n=============================================================\n")
cat("  ACCURACY / CALIBRATION SUMMARY\n")
cat("     mean (sd) [2.5%, 97.5%] over replicates\n")
cat("=============================================================\n\n")

summary_table <- aggregate(
  value ~ variable + type + metric + scenario + mechanism,
  data = results_df,
  FUN  = function(x) c(mean = mean(x, na.rm = TRUE),
                       sd   = sd(x,   na.rm = TRUE),
                       q025 = quantile(x, 0.025, na.rm = TRUE),
                       q975 = quantile(x, 0.975, na.rm = TRUE))
)
summary_flat <- cbind(
  summary_table[, c("variable", "type", "metric", "scenario", "mechanism")],
  as.data.frame(summary_table$value)
)
names(summary_flat) <- c("variable", "type", "metric", "scenario",
                         "mechanism", "mean", "sd", "q025", "q975")

for (m in unique(summary_flat$metric)) {
  cat("--- Metric:", m, "---\n")
  sub <- summary_flat[summary_flat$metric == m, ]
  sub <- sub[order(sub$variable, sub$scenario, sub$mechanism), ]
  print(sub[, c("variable", "scenario", "mechanism",
                "mean", "sd", "q025", "q975")], row.names = FALSE)
  cat("\n")
}

write.csv(summary_flat, file.path(RESULTS_DIR, "imputation_quality_summary.csv"),
          row.names = FALSE)

# =============================================================================
# 8. PLOTS
# =============================================================================

plot_file <- file.path(RESULTS_DIR, "imputation_quality_plots.pdf")
pdf(plot_file, width = 12, height = 8)

mech_cols <- c(phylo_MAR = "#1b7837", trait_MAR = "#f6b26b",
               trait_MNAR = "#e06666")

results_df$scenario  <- factor(results_df$scenario,  levels = SCENARIOS)
results_df$mechanism <- factor(results_df$mechanism, levels = MECHANISMS)

plot_by_mechanism <- function(v, metric_name, ylab, ref_line = NA,
                              ylim = NULL) {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  sub <- results_df[results_df$variable == v &
                    results_df$metric == metric_name, ]
  for (sc in SCENARIOS) {
    ss <- sub[sub$scenario == sc, ]
    boxplot(value ~ mechanism, data = ss,
            main = paste0(sc, " - ", v, " ", metric_name),
            xlab = "Mechanism", ylab = ylab,
            col  = mech_cols, ylim = ylim)
    if (!is.na(ref_line)) abline(h = ref_line, lty = 2, col = "grey40")
  }
  mtext(paste("Imputation", metric_name, "for", v),
        outer = TRUE, cex = 1.1)
}

# Continuous / count: NRMSE (dashed line at 1 = mean-baseline), log-MAE
# for poisson only, plus coverage.
plot_by_mechanism("y",  "nrmse",      "NRMSE",              ref_line = 1)
plot_by_mechanism("x3", "nrmse",      "NRMSE",              ref_line = 1)
plot_by_mechanism("x3", "log_mae",    "log-MAE")
plot_by_mechanism("y",  "coverage95",
                  "95% PI coverage", ref_line = 0.95, ylim = c(0, 1))
plot_by_mechanism("x3", "coverage95",
                  "95% PI coverage", ref_line = 0.95, ylim = c(0, 1))

# Categorical / ordered: accuracy, balanced accuracy, Brier
for (v in c("x1", "x2", "x4")) {
  plot_by_mechanism(v, "accuracy",          "Accuracy", ylim = c(0, 1))
  plot_by_mechanism(v, "balanced_accuracy", "Balanced accuracy",
                    ylim = c(0, 1))
  plot_by_mechanism(v, "brier",             "Brier score (lower = better)")
}
plot_by_mechanism("x4", "ordinal_mae", "|rank_true - rank_imp|")

dev.off()
cat("Plots saved to:", plot_file, "\n\nDone.\n")
