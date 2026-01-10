# Example: Using the linear_models component from sim_bace

library(ape)  # For phylogenetic trees
source('R/simulate_auxiliary.r')
source('R/simulate_simBACE.R')

# Simulate data
my_fixeff <- make_fixeff(
  var_names = c('y', 'x1', 'x2'),
  predictor_types = c('binary', 'gaussian'),
  formulas = list(
    y ~ x1 + x2,
    x2 ~ x1
  ),
  betas = list(
    y = c(x1 = 0.6, x2 = 0.4),
    x2 = c(x1 = 0.5)
  )
)

my_raneff <- make_raneff(
  var_names = c('y', 'x1', 'x2'),
  phylo_frac = c(0.4, 0.3, 0.2),
  species_frac = c(0.2, 0.2, 0.3)
)

set.seed(42)
sim <- sim_bace(
  response_type = 'gaussian',
  predictor_types = c('binary', 'gaussian'),
  fixeff = my_fixeff,
  raneff = my_raneff,
  n_cases = 100,
  n_species = 25,
  missingness = c(0.1, 0.15, 0.1)
)

# Access linear model components
cat("\n=== Linear Models Available ===\n")
cat("Variables:", names(sim$linear_models), "\n\n")

# Example 1: Inspect the response variable linear model
cat("=== Response (y) Linear Model ===\n")
cat("Link function:", sim$linear_models$y$link, "\n")
cat("Design matrix X dimensions:", dim(sim$linear_models$y$X), "\n")
cat("  - Rows (observations):", nrow(sim$linear_models$y$X), "\n")
cat("  - Columns (predictors + intercept):", ncol(sim$linear_models$y$X), "\n")
cat("  - Column names:", colnames(sim$linear_models$y$X), "\n")
cat("\nBeta coefficients:", sim$linear_models$y$beta, "\n")
cat("  - Interpretation: intercept=", sim$linear_models$y$beta[1],
    ", x1=", sim$linear_models$y$beta[2],
    ", x2=", sim$linear_models$y$beta[3], "\n", sep="")

cat("\nRandom effects design matrix Z dimensions:", dim(sim$linear_models$y$Z), "\n")
cat("Random effects u (phylo + species) length:", length(sim$linear_models$y$u), "\n")
cat("Residuals e length:", length(sim$linear_models$y$e), "\n")

# Example 2: Inspect a predictor variable
cat("\n=== Predictor (x2) Linear Model ===\n")
cat("Link function:", sim$linear_models$x2$link, "\n")
cat("Beta coefficients:", sim$linear_models$x2$beta, "\n")
cat("  - intercept=", sim$linear_models$x2$beta[1],
    ", x1=", sim$linear_models$x2$beta[2], "\n", sep="")

# Example 3: Categorical variable
cat("\n=== Predictor (x1 - binary) Linear Model ===\n")
cat("Link function:", sim$linear_models$x1$link, "\n")
cat("Design matrix dimensions:", dim(sim$linear_models$x1$X), "\n")
cat("Beta (intercept only):", sim$linear_models$x1$beta, "\n")
cat("Note: This is an independent predictor with no dependencies\n")

# Example 4: Reconstruct linear predictor for response
cat("\n=== Reconstructing Linear Predictor ===\n")
cat("For the response 'y', the true generative model is:\n")
cat("  y = X * beta + Z * u + e\n")
cat("where:\n")
cat("  X = design matrix with intercept, x1, x2\n")
cat("  beta = [", paste(round(sim$linear_models$y$beta, 3), collapse=", "), "]\n", sep="")
cat("  Z = random effects design matrix\n")
cat("  u = combined phylogenetic + species random effects\n")
cat("  e = residual error\n\n")

# Reconstruct (before centering)
linear_pred_raw <- as.numeric(sim$linear_models$y$X %*% sim$linear_models$y$beta) +
                   as.numeric(sim$linear_models$y$Z %*% sim$linear_models$y$u) +
                   sim$linear_models$y$e

cat("Raw linear predictor (before centering):\n")
cat("  Mean:", round(mean(linear_pred_raw), 3), "\n")
cat("  SD:", round(sd(linear_pred_raw), 3), "\n")
cat("  Range: [", round(min(linear_pred_raw), 3), ",", round(max(linear_pred_raw), 3), "]\n")

# The observed data is centered
cat("\nObserved y (after centering):\n")
cat("  Mean:", round(mean(sim$data$y, na.rm=TRUE), 3), "\n")
cat("  SD:", round(sd(sim$data$y, na.rm=TRUE), 3), "\n")
cat("  Range: [", round(min(sim$data$y, na.rm=TRUE), 3), ",",
    round(max(sim$data$y, na.rm=TRUE), 3), "]\n")

cat("\nNote: For gaussian variables, sim_bace centers the linear predictor to mean=0\n")
cat("      This preserves the variance structure for parameter recovery.\n")
cat("      The stored linear_models components show the TRUE generative process.\n")

# Example 5: Using components for custom analysis
cat("\n=== Using Components for Analysis ===\n")
cat("You can use these components to:\n")
cat("  1. Understand the true data generating process\n")
cat("  2. Compute exact predictions with known parameters\n")
cat("  3. Test parameter recovery in your analysis methods\n")
cat("  4. Decompose variance into phylo, species, and residual\n")
cat("  5. Examine the design matrices used for simulation\n")

# Variance decomposition example
phylo_var <- var(as.numeric(sim$linear_models$y$Z %*% sim$random_effects$u_phylo$y))
species_var <- var(as.numeric(sim$linear_models$y$Z %*% sim$random_effects$u_species$y))
residual_var <- var(sim$linear_models$y$e)
total_var <- phylo_var + species_var + residual_var

cat("\nVariance decomposition for y:\n")
cat("  Phylogenetic:", round(phylo_var, 4), "(", round(100*phylo_var/total_var, 1), "%)\n")
cat("  Species:", round(species_var, 4), "(", round(100*species_var/total_var, 1), "%)\n")
cat("  Residual:", round(residual_var, 4), "(", round(100*residual_var/total_var, 1), "%)\n")
cat("  Total:", round(total_var, 4), "\n")

cat("\n=== Example Complete ===\n")
