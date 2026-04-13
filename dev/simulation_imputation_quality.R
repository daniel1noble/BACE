# =============================================================================
# Simulation: Assessing BACE Imputation Quality
# =============================================================================
#
# Goal: Evaluate how well bace() recovers true values for missing data across
#       500 simulated datasets, 3 levels of missingness (20%, 40%, 70%), and a
#       mixture of continuous + categorical variables with phylogenetic signal.
#
# Design:
#   - Response:  gaussian (continuous)
#   - Predictors: x1 = gaussian, x2 = binary, x3 = multinomial (3 levels)
#   - Phylogenetic signal in every variable (moderate, lambda ~ 0.4–0.7)
#   - 100 species, 200 observations
#   - 500 simulation replicates per missingness level (1500 total)
#
# Metrics (missing values only):
#   Continuous: RMSE, correlation, coverage of 95% credible intervals
#   Categorical: classification accuracy, Brier score
# =============================================================================

library(BACE)
library(ape)
library(parallel)

# ---- Configuration ----------------------------------------------------------

set.seed(2026)

N_SIMS        <- 1000
MISS_LEVELS   <- c(0.20, 0.40, 0.70)
N_CASES       <- 200
N_SPECIES     <- 100
N_CORES       <- 4L       # for bace() final runs AND outer parallel loop
RESULTS_DIR   <- file.path("dev", "simulation_results")

# MCMC settings (light – increase for production)
NITT    <- 15000
THIN    <- 10
BURNIN  <- 5000
RUNS    <- 10
N_FINAL <- 5

# Create results directory
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ---- Helper: Performance Metrics -------------------------------------------

#' Compute imputation quality for one simulation replicate
#' @param true_data  Complete data.frame (no NAs)
#' @param imp_data   Imputed data.frame (from bace)
#' @param miss_mask  Logical data.frame — TRUE where values were set to NA
#' @param var_types  Named character vector: variable name -> "continuous" or "categorical"
#' @return data.frame with one row per variable, columns: variable, type, metric, value
evaluate_imputation <- function(true_data, imp_data, miss_mask, var_types) {
  results <- list()

  for (v in names(var_types)) {
    idx <- miss_mask[[v]]
    if (sum(idx) == 0) next

    true_vals <- true_data[[v]][idx]
    imp_vals  <- imp_data[[v]][idx]

    if (var_types[v] == "continuous") {
      # RMSE
      rmse <- sqrt(mean((as.numeric(imp_vals) - as.numeric(true_vals))^2))
      # Correlation
      cor_val <- if (length(unique(true_vals)) > 1) {
        cor(as.numeric(imp_vals), as.numeric(true_vals))
      } else {
        NA_real_
      }
      # Mean absolute error
      mae <- mean(abs(as.numeric(imp_vals) - as.numeric(true_vals)))

      results[[length(results) + 1]] <- data.frame(
        variable = v, type = "continuous",
        metric = c("rmse", "correlation", "mae"),
        value  = c(rmse, cor_val, mae),
        stringsAsFactors = FALSE
      )

    } else {
      # Classification accuracy
      accuracy <- mean(as.character(imp_vals) == as.character(true_vals))

      results[[length(results) + 1]] <- data.frame(
        variable = v, type = "categorical",
        metric = "accuracy",
        value  = accuracy,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results)
}

# ---- Helper: Run One Simulation Replicate -----------------------------------

#' Simulate data, apply missingness, run bace, evaluate
#' @param sim_id     Integer simulation ID
#' @param miss_prop  Proportion of missingness (applied to all variables with missing data)
#' @return data.frame of metrics, or NULL on failure
run_one_sim <- function(sim_id, miss_prop) {

  # --- 1. Simulate complete data with phylogenetic signal ---
  sim <- sim_bace(
    response_type  = "gaussian",
    predictor_types = c("gaussian", "binary", "multinomial3"),
    var_names       = c("y", "x1", "x2", "x3"),
    phylo_signal    = c(0.7, 0.6, 0.8, 0.65),  # moderate signal in all variables
    n_cases         = N_CASES,
    n_species       = N_SPECIES,
    beta_sparsity   = 0.3,
    missingness     = c(0, 0, 0, 0)           # no missingness yet
  )

  complete_data <- sim$complete_data
  tree          <- sim$tree

  # --- 2. Apply missingness to x1, x2, x3 (not y) ---
  miss_data <- complete_data
  miss_mask <- data.frame(
    x1 = logical(nrow(complete_data)),
    x2 = logical(nrow(complete_data)),
    x3 = logical(nrow(complete_data))
  )

  for (v in c("x1", "x2", "x3")) {
    n      <- nrow(miss_data)
    n_miss <- floor(n * miss_prop)
    idx    <- sample(seq_len(n), n_miss, replace = FALSE)
    miss_mask[[v]][idx] <- TRUE
    miss_data[[v]][idx] <- NA
  }

  # --- 3. Run bace imputation ---
  bace_result <- tryCatch({
    bace(
      fixformula    = list("y ~ x1 + x2 + x3",
                           "x1 ~ y + x2 + x3",
                           "x2 ~ y + x1 + x3",
                           "x3 ~ y + x1 + x2"),
      ran_phylo_form = "~1|Species",
      phylo          = tree,
      data           = miss_data,
      nitt           = NITT,
      thin           = THIN,
      burnin         = BURNIN,
      runs           = RUNS,
      n_final        = N_FINAL,
      species        = FALSE,
      verbose        = FALSE,
      skip_conv      = TRUE,
      n_cores        = 1L       # inner parallelism off (outer loop is parallel)
    )
  }, error = function(e) {
    message("  [sim ", sim_id, ", miss=", miss_prop, "] bace() error: ", e$message)
    return(NULL)
  })

  if (is.null(bace_result)) return(NULL)

  # --- 4. Evaluate: average across the N_FINAL imputed datasets ---
  var_types <- c(x1 = "continuous", x2 = "categorical", x3 = "categorical")

  all_evals <- lapply(bace_result$imputed_datasets, function(imp_ds) {
    evaluate_imputation(complete_data, imp_ds, miss_mask, var_types)
  })

  # Pool metrics across imputations (mean)
  pooled <- do.call(rbind, all_evals)
  pooled <- aggregate(value ~ variable + type + metric, data = pooled, FUN = mean)

  pooled$sim_id    <- sim_id
  pooled$miss_prop <- miss_prop

  return(pooled)
}

# ---- Main Simulation Loop ---------------------------------------------------

cat("=============================================================\n")
cat("  BACE Imputation Quality Simulation\n")
cat("  Replicates :", N_SIMS, "\n")
cat("  Missingness :", paste0(MISS_LEVELS * 100, "%", collapse = ", "), "\n")
cat("  Variables   : y (gaussian), x1 (gaussian), x2 (binary), x3 (multinomial3)\n")
cat("  Phylo signal: y=0.7, x1=0.6, x2=0.8, x3=0.65\n")
cat("  MCMC        : nitt=", NITT, " thin=", THIN, " burnin=", BURNIN, "\n")
cat("=============================================================\n\n")

all_results <- list()

for (miss_prop in MISS_LEVELS) {

  cat("--- Missingness:", miss_prop * 100, "% ---\n")
  t0 <- Sys.time()

  # Run in parallel across simulation replicates
  sim_results <- mclapply(seq_len(N_SIMS), function(i) {
    if (i %% 50 == 0) message("  sim ", i, " / ", N_SIMS)
    run_one_sim(sim_id = i, miss_prop = miss_prop)
  }, mc.cores = N_CORES)

  # Collect successful results
  sim_results <- sim_results[!sapply(sim_results, is.null)]
  cat("  Completed:", length(sim_results), "/", N_SIMS, "replicates\n")
  cat("  Time:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n\n")

  all_results[[as.character(miss_prop)]] <- do.call(rbind, sim_results)
}

# Combine everything
results_df <- do.call(rbind, all_results)
rownames(results_df) <- NULL

# Save raw results
saveRDS(results_df, file.path(RESULTS_DIR, "imputation_quality_results.rds"))
write.csv(results_df, file.path(RESULTS_DIR, "imputation_quality_results.csv"),
          row.names = FALSE)

cat("Results saved to:", RESULTS_DIR, "\n\n")

# ---- Summarise Results ------------------------------------------------------

cat("=============================================================\n")
cat("  SUMMARY\n")
cat("=============================================================\n\n")

summary_table <- aggregate(
  value ~ variable + type + metric + miss_prop,
  data = results_df,
  FUN  = function(x) {
    c(mean = mean(x, na.rm = TRUE),
      sd   = sd(x, na.rm = TRUE),
      q025 = quantile(x, 0.025, na.rm = TRUE),
      q975 = quantile(x, 0.975, na.rm = TRUE))
  }
)

# Flatten the matrix column
summary_flat <- cbind(
  summary_table[, c("variable", "type", "metric", "miss_prop")],
  as.data.frame(summary_table$value)
)

# Print nicely
for (m in unique(summary_flat$metric)) {
  cat("--- Metric:", m, "---\n")
  sub <- summary_flat[summary_flat$metric == m, ]
  sub$miss_pct <- paste0(sub$miss_prop * 100, "%")
  print(sub[, c("variable", "miss_pct", "mean", "sd", "q025", "q975")], row.names = FALSE)
  cat("\n")
}

# Save summary
write.csv(summary_flat, file.path(RESULTS_DIR, "imputation_quality_summary.csv"),
          row.names = FALSE)

# ---- Plot Results -----------------------------------------------------------

plot_file <- file.path(RESULTS_DIR, "imputation_quality_plots.pdf")
pdf(plot_file, width = 10, height = 8)

# --- Panel 1: Continuous variable (x1) performance ---
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# RMSE by missingness
x1_rmse <- results_df[results_df$variable == "x1" & results_df$metric == "rmse", ]
boxplot(value ~ miss_prop, data = x1_rmse,
        main = "x1 (continuous): RMSE",
        xlab = "Missingness", ylab = "RMSE",
        names = paste0(sort(unique(x1_rmse$miss_prop)) * 100, "%"),
        col = c("#93c47d", "#f6b26b", "#e06666"))

# Correlation by missingness
x1_cor <- results_df[results_df$variable == "x1" & results_df$metric == "correlation", ]
boxplot(value ~ miss_prop, data = x1_cor,
        main = "x1 (continuous): Correlation",
        xlab = "Missingness", ylab = "Correlation (imputed vs true)",
        names = paste0(sort(unique(x1_cor$miss_prop)) * 100, "%"),
        col = c("#93c47d", "#f6b26b", "#e06666"))
abline(h = 1, lty = 2, col = "grey50")

# --- Panel 2: Categorical variable accuracy ---
x2_acc <- results_df[results_df$variable == "x2" & results_df$metric == "accuracy", ]
boxplot(value ~ miss_prop, data = x2_acc,
        main = "x2 (binary): Accuracy",
        xlab = "Missingness", ylab = "Accuracy",
        names = paste0(sort(unique(x2_acc$miss_prop)) * 100, "%"),
        col = c("#93c47d", "#f6b26b", "#e06666"))

x3_acc <- results_df[results_df$variable == "x3" & results_df$metric == "accuracy", ]
boxplot(value ~ miss_prop, data = x3_acc,
        main = "x3 (multinomial3): Accuracy",
        xlab = "Missingness", ylab = "Accuracy",
        names = paste0(sort(unique(x3_acc$miss_prop)) * 100, "%"),
        col = c("#93c47d", "#f6b26b", "#e06666"))

# --- Panel 3: Overlaid density plots ---
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

x1_rmse_split <- split(x1_rmse$value, x1_rmse$miss_prop)
xlim <- range(unlist(x1_rmse_split), na.rm = TRUE)
cols <- c("#93c47d", "#f6b26b", "#e06666")

plot(density(x1_rmse_split[[1]], na.rm = TRUE), col = cols[1], lwd = 2,
     xlim = xlim, main = "x1 RMSE Distribution", xlab = "RMSE")
lines(density(x1_rmse_split[[2]], na.rm = TRUE), col = cols[2], lwd = 2)
lines(density(x1_rmse_split[[3]], na.rm = TRUE), col = cols[3], lwd = 2)
legend("topright", legend = paste0(MISS_LEVELS * 100, "%"),
       col = cols, lwd = 2, title = "Missingness")

# Accuracy distributions for x2
x2_split <- split(x2_acc$value, x2_acc$miss_prop)
xlim2 <- range(unlist(x2_split), na.rm = TRUE)
plot(density(x2_split[[1]], na.rm = TRUE), col = cols[1], lwd = 2,
     xlim = xlim2, main = "x2 (binary) Accuracy", xlab = "Accuracy")
lines(density(x2_split[[2]], na.rm = TRUE), col = cols[2], lwd = 2)
lines(density(x2_split[[3]], na.rm = TRUE), col = cols[3], lwd = 2)
legend("topleft", legend = paste0(MISS_LEVELS * 100, "%"),
       col = cols, lwd = 2, title = "Missingness")

# Accuracy distributions for x3
x3_split <- split(x3_acc$value, x3_acc$miss_prop)
xlim3 <- range(unlist(x3_split), na.rm = TRUE)
plot(density(x3_split[[1]], na.rm = TRUE), col = cols[1], lwd = 2,
     xlim = xlim3, main = "x3 (multinomial3) Accuracy", xlab = "Accuracy")
lines(density(x3_split[[2]], na.rm = TRUE), col = cols[2], lwd = 2)
lines(density(x3_split[[3]], na.rm = TRUE), col = cols[3], lwd = 2)
legend("topleft", legend = paste0(MISS_LEVELS * 100, "%"),
       col = cols, lwd = 2, title = "Missingness")

dev.off()
cat("Plots saved to:", plot_file, "\n")

cat("\nDone.\n")
