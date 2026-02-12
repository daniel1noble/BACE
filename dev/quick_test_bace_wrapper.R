# ========================================================================
# Quick Test Example for BACE Wrapper Functions
# ========================================================================
# A minimal example to test the new bace wrapper functionality
# Uses small data and few iterations for quick testing

library(BACE)

# Set seed for reproducibility
set.seed(123)

# Create small phylogenetic tree (10 species for speed)
phylo <- ape::rtree(10)
phylo <- ape::compute.brlen(phylo, method = "Grafen")

# Simple dataset
data <- data.frame(
  y = rpois(10, lambda = 5),
  x1 = rnorm(10, 10, 2),
  x2 = rnorm(10, 5, 1),
  Species = phylo$tip.label
)

# Introduce missing data
data$y[c(1, 5)] <- NA
data$x1[c(2, 7)] <- NA

cat("Original data with missing values:\n")
print(data)

# ========================================================================
# Test 1: bace_imp + bace_final_imp + pool_posteriors (manual workflow)
# ========================================================================

cat("\n\n========== Test 1: Manual Workflow ==========\n\n")

# Step 1: Initial imputation
cat("Running initial imputation...\n")
initial <- bace_imp(
  fixformula = list("y ~ x1 + x2", "x1 ~ x2"),
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  runs = 3,
  nitt = 500,
  thin = 2,
  burnin = 100,
  verbose = TRUE
)

# Step 2: Assess convergence
cat("\n\nAssessing convergence...\n")
conv <- assess_convergence(initial, method = "summary")
print(conv)

# Step 3: Final imputation
cat("\n\nRunning final imputations...\n")
final <- bace_final_imp(
  bace_object = initial,
  fixformula = list("y ~ x1 + x2", "x1 ~ x2"),
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  n_final = 5,
  nitt = 500,
  thin = 2,
  burnin = 100,
  verbose = TRUE
)

print(final)

# Step 4: Pool posteriors
cat("\n\nPooling posteriors...\n")
pooled <- pool_posteriors(final)
print(pooled)

# Examine results for y
cat("\n\nSummary for variable 'y':\n")
summary_y <- summary(pooled$models$y)
print(summary_y)

# ========================================================================
# Test 2: Complete bace wrapper (automated workflow)
# ========================================================================

cat("\n\n========== Test 2: Automated Wrapper Workflow ==========\n\n")

result <- bace(
  fixformula = list("y ~ x1 + x2", "x1 ~ x2", "x2 ~ x1"),
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  runs = 3,
  n_final = 5,
  nitt = 500,
  thin = 2,
  burnin = 100,
  verbose = TRUE,
  plot = FALSE
)

summary(result$pooled_models$models$y)
summary(result)
print(result)

# Access pooled results
cat("\n\nPooled model for 'y':\n")
print(result$pooled_models$models$y)

cat("\n\nSummary (using standard MCMCglmm summary method):\n")
summary(result$pooled_models$models$y)

# Verify it's a true MCMCglmm object
cat("\n\nClass check:\n")
cat("Is MCMCglmm object:", inherits(result$pooled_models$models$y, "MCMCglmm"), "\n")
cat("Full class:", paste(class(result$pooled_models$models$y), collapse = ", "), "\n")

# Check convergence
cat("\n\nConvergence status:", result$converged, "\n")

# ========================================================================
# Test 3: Verify imputed data
# ========================================================================

cat("\n\n========== Test 3: Verify Imputed Data ==========\n\n")

# Get one of the final imputed datasets
imputed <- result$final_results$all_datasets[[1]]

cat("Imputed dataset (no missing values):\n")
print(imputed)

cat("\nMissing values check:\n")
print(sapply(imputed, function(x) sum(is.na(x))))

# Compare original vs imputed for missing observations
cat("\n\nOriginal y:", data$y, "\n")
cat("Imputed y:", imputed$y, "\n")

cat("\n\n========== All Tests Complete ==========\n\n")
