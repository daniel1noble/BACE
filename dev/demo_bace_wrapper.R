# ========================================================================
# Demo: Using the BACE Wrapper Function with Convergence Checking
# and Posterior Pooling
# ========================================================================
# 
# This script demonstrates the complete BACE workflow:
# 1. Initial imputation with convergence checking
# 2. Final imputation runs (if convergence is achieved)
# 3. Posterior pooling across multiple imputations
# 4. Extracting and interpreting results

library(BACE)
library(ape)
library(phytools)

# ========================================================================
# Step 1: Generate Example Data
# ========================================================================

set.seed(42)

# Create a phylogenetic tree with 30 species
phylo <- rtree(30)
phylo <- compute.brlen(phylo, method = "Grafen")
phylo <- force.ultrametric(phylo)

# Create dataset with multiple variable types
data <- data.frame(
  # Response variable (count data)
  y = rpois(30, lambda = 5),
  
  # Continuous predictors
  x1 = rnorm(30, 10, 2),
  x2 = rnorm(30, 5, 1),
  
  # Categorical predictor
  x3 = factor(rep(c("A", "B", "C"), length.out = 30)),
  
  # Another continuous
  x4 = rnorm(30, 7, 1.5),
  
  # Species identifier
  Species = phylo$tip.label
)

# Introduce missing data (MCAR - missing completely at random)
data$y[sample(1:30, 8)] <- NA
data$x1[sample(1:30, 5)] <- NA
data$x2[sample(1:30, 6)] <- NA
data$x4[sample(1:30, 4)] <- NA

# Check missingness pattern
cat("Missing data summary:\n")
sapply(data, function(x) sum(is.na(x)))

# ========================================================================
# Step 2: Run Complete BACE Analysis
# ========================================================================

cat("\n\n========== Running Complete BACE Analysis ==========\n\n")

# Option 1: Simple analysis with single formula
# (only if one variable has missing data you care about)
if (FALSE) {
  result_simple <- bace(
    fixformula = "y ~ x1 + x2 + x4",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 10,        # Initial runs for convergence checking
    n_final = 10,     # Final runs for posterior pooling
    nitt = 6000,      # MCMC iterations
    thin = 5,         # Thinning
    burnin = 1000,    # Burn-in
    verbose = TRUE,   # Print progress
    plot = TRUE       # Show convergence plots
  )
}


# Option 2: Analysis with multiple formulas (recommended for multiple missing variables)
result <- bace(
  fixformula = list(
    "y ~ x1 + x2 + x4",      # Model for y
    "x1 ~ x2 + x4",          # Model for x1
    "x2 ~ x1 + x4",          # Model for x2
    "x4 ~ x1 + x2"           # Model for x4
  ),
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  runs = 10,           # Initial runs for convergence
  n_final = 10,        # Final imputation runs
  nitt = 6000,         # MCMC iterations
  thin = 5,
  burnin = 1000,
  verbose = TRUE,
  species = FALSE,     # Include species-level random effects
  plot = FALSE,        # Set to TRUE to see convergence plots
  max_attempts = 3     # Maximum convergence attempts
)

# ========================================================================
# Step 3: Examine Results
# ========================================================================

cat("\n\n========== Examining Results ==========\n\n")

# Print overall summary
print(result)

# Check convergence
cat("\n\nConvergence Status:\n")
cat("Converged:", result$converged, "\n")
cat("Attempts needed:", result$n_attempts, "\n")

# View convergence details
print(result$convergence)

# Plot convergence diagnostics (if not plotted earlier)
if (result$converged) {
  plot(result$convergence)
}

# ========================================================================
# Step 4: Access Pooled Posterior Results
# ========================================================================

cat("\n\n========== Pooled Posterior Results ==========\n\n")

# Access pooled models
print(result$pooled_models)

# Get pooled MCMCglmm model for specific variable
cat("\n\nPooled MCMCglmm model for 'y':\n")
pooled_y <- result$pooled_models$models$y

# These models are TRUE MCMCglmm objects, so standard methods work!
print(pooled_y)

# Standard MCMCglmm summary
cat("\n\nSummary for 'y' (using standard MCMCglmm summary method):\n")
summary(pooled_y)

# You can also directly access components like any MCMCglmm object
cat("\n\nFixed Effects (direct access to Sol matrix, first 5 rows):\n")
print(head(pooled_y$Sol, 5))

cat("\n\nVariance Components (direct access to VCV matrix, first 5 rows):\n")
print(head(pooled_y$VCV, 5))

# ========================================================================
# Step 5: Working with Pooled Models (Just like MCMCglmm!)
# ========================================================================

cat("\n\n========== Working with Pooled MCMCglmm Models ==========\n\n")

# The pooled models ARE MCMCglmm objects, so everything works!
y_model <- result$pooled_models$models$y

# Check that it's an MCMCglmm object
cat("Class of pooled model:\n")
print(class(y_model))

# Posterior samples for fixed effects (direct access)
cat("\n\nFixed effect posterior samples (first 5):\n")
print(head(y_model$Sol, 5))

# Calculate custom posterior summaries
# Example: Probability that x1 effect > 0
if ("x1" %in% colnames(y_model$Sol)) {
  prob_x1_positive <- mean(y_model$Sol[, "x1"] > 0)
  cat("\n\nProbability that x1 effect > 0:", round(prob_x1_positive, 3), "\n")
}

# Standard MCMCglmm diagnostics work too
if (require(coda)) {
  cat("\n\nEffective sample size for fixed effects:\n")
  print(coda::effectiveSize(y_model$Sol))
}

# Plot posterior distributions
if (require(ggplot2)) {
  library(ggplot2)
  
  # Convert to data frame for plotting
  sol_df <- as.data.frame(y_model$Sol)
  
  # Plot x1 effect if present
  if ("x1" %in% colnames(y_model$Sol)) {
    p1 <- ggplot(sol_df, aes(x = x1)) +
      geom_density(fill = "steelblue", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "Posterior Distribution of x1 Effect on y",
           subtitle = paste("Based on", nrow(sol_df), "pooled posterior samples"),
           x = "Effect Size", y = "Density") +
      theme_minimal()
    
    print(p1)
  }
}

# Standard MCMCglmm plots also work!
# Uncomment to see standard MCMCglmm diagnostic plots
# plot(y_model)

# ========================================================================
# Step 6: Access Imputed Datasets
# ========================================================================

cat("\n\n========== Accessing Imputed Datasets ==========\n\n")

# All final imputed datasets
final_datasets <- result$final_results$all_datasets

cat("Number of final imputed datasets:", length(final_datasets), "\n")

# Look at first imputed dataset
imputed_data_1 <- final_datasets[[1]]
cat("\nFirst imputed dataset (first 5 rows):\n")
print(head(imputed_data_1, 5))

# Verify no missing data in imputed datasets
cat("\nMissing values in first imputed dataset:\n")
print(sapply(imputed_data_1, function(x) sum(is.na(x))))

# Compare original and imputed values for one observation
cat("\n\nComparison for observation with originally missing y:\n")
orig_missing_idx <- which(is.na(data$y))[1]
cat("Original y:", data$y[orig_missing_idx], "(missing)\n")
cat("Imputed y (dataset 1):", imputed_data_1$y[orig_missing_idx], "\n")
cat("Imputed y (dataset 5):", final_datasets[[5]]$y[orig_missing_idx], "\n")

# ========================================================================
# Step 7: Model Comparison (Optional)
# ========================================================================

cat("\n\n========== Model Information ==========\n\n")

# DIC for pooled model
if (!is.null(result$pooled_models$models$y$DIC)) {
  cat("Mean DIC for 'y' model:", round(result$pooled_models$models$y$DIC, 2), "\n")
}

# Number of posterior samples
cat("Total posterior samples:", result$pooled_models$models$y$total_samples, "\n")
cat("Samples per imputation:", result$pooled_models$models$y$n_samples_per_imputation, "\n")
cat("Number of imputations:", result$pooled_models$models$y$n_imputations, "\n")

# ========================================================================
# Step 8: Advanced Usage - Using Components Separately
# ========================================================================

cat("\n\n========== Advanced: Accessing Individual Components ==========\n\n")

# If you need more control, you can use functions separately:

if (FALSE) {  # Set to TRUE to run
  
  # Step 1: Initial imputation
  initial <- bace_imp(
    fixformula = list("y ~ x1 + x2 + x4", "x1 ~ x2 + x4", 
                     "x2 ~ x1 + x4", "x4 ~ x1 + x2"),
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data,
    runs = 10,
    nitt = 6000,
    thin = 5,
    burnin = 1000,
    verbose = TRUE
  )
  
  # Step 2: Assess convergence
  conv <- assess_convergence(initial, method = "summary")
  print(conv)
  plot(conv)
  
  # Step 3: Final imputation (only if converged)
  if (conv$converged) {
    final <- bace_final_imp(
      bace_object = initial,
      fixformula = list("y ~ x1 + x2 +  using standard MCMCglmm methods
    cat("\n\nUsing standard MCMCglmm summary:\n")
    summary(pooled$models$y)
    
    # The pooled models are true MCMCglmm objects!
    cat("\nClass of pooled model:", class(pooled$models$y), "\n")
  }
}

cat("\n\n========== Demo Complete ==========\n\n")
cat("Key Points:\n")
cat("- Pooled models are TRUE MCMCglmm objects\n")
cat("- Use standard MCMCglmm methods: summary(), print(), plot()\n")
cat("- Results account for both estimation AND imputation uncertainty
      thin = 5,
      burnin = 1000,
      verbose = TRUE
    )
    
    # Step 4: Pool posteriors
    pooled <- pool_posteriors(final)
    print(pooled)
    
    # Step 5: Analyze specific variable
    summary(pooled$models$y)
  }
}

cat("\n\n========== Demo Complete ==========\n\n")
cat("See ?bace for full documentation\n")
cat("See ?pool_posteriors for posterior pooling details\n")
cat("See ?assess_convergence for convergence diagnostic options\n")
