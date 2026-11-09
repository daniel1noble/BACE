# Demonstration of MCMC parameter list functionality
# This script shows how to use different MCMC settings for different models

library(BACE)
library(ape)

set.seed(123)

# Create example data
phylo <- rtree(15)
phylo <- compute.brlen(phylo, method = "Grafen")

data <- data.frame(
  y = rnorm(15, 10, 2),
  x1 = rnorm(15, 5, 1),
  x2 = rnorm(15, 3, 0.5),
  x3 = rnorm(15, 7, 1.5),
  Species = phylo$tip.label
)

# Introduce missing data
data$y[c(1, 3, 5)] <- NA
data$x1[c(2, 4, 6)] <- NA
data$x2[c(7, 9, 11)] <- NA

cat("\n=== Example 1: Single MCMC settings for all models (backward compatible) ===\n")
result1 <- bace_imp(
  fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3", "x2 ~ x1 + x3"),
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  runs = 2,
  nitt = 1000,
  thin = 5,
  burnin = 100,
  verbose = TRUE
)

cat("\nModel 1 (y) MCMC samples:", nrow(result1$models_last_run$y$Sol), "\n")
cat("Model 2 (x1) MCMC samples:", nrow(result1$models_last_run$x1$Sol), "\n")
cat("Model 3 (x2) MCMC samples:", nrow(result1$models_last_run$x2$Sol), "\n")
cat("Expected samples: (1000-100)/5 =", (1000-100)/5, "\n")

cat("\n=== Example 2: Model-specific MCMC settings (new functionality) ===\n")
result2 <- bace_imp(
  fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3", "x2 ~ x1 + x3"),
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  runs = 2,
  nitt = list(1000, 1500, 2000),    # Different nitt for each model
  thin = list(5, 10, 5),             # Different thin for each model
  burnin = list(100, 200, 300),      # Different burnin for each model
  verbose = TRUE
)

cat("\nModel 1 (y) MCMC samples:", nrow(result2$models_last_run$y$Sol), "\n")
cat("Expected: (1000-100)/5 =", (1000-100)/5, "\n")

cat("\nModel 2 (x1) MCMC samples:", nrow(result2$models_last_run$x1$Sol), "\n")
cat("Expected: (1500-200)/10 =", (1500-200)/10, "\n")

cat("\nModel 3 (x2) MCMC samples:", nrow(result2$models_last_run$x2$Sol), "\n")
cat("Expected: (2000-300)/5 =", (2000-300)/5, "\n")

cat("\n=== Example 3: Mixed specification (some list, some single) ===\n")
result3 <- bace_imp(
  fixformula = list("y ~ x1 + x2", "x1 ~ x2 + x3"),
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  runs = 2,
  nitt = list(1000, 1500),  # List for nitt
  thin = 5,                  # Single value for thin (applied to both)
  burnin = list(100, 200),   # List for burnin
  verbose = TRUE
)

cat("\nModel 1 (y) MCMC samples:", nrow(result3$models_last_run$y$Sol), "\n")
cat("Expected: (1000-100)/5 =", (1000-100)/5, "\n")

cat("\nModel 2 (x1) MCMC samples:", nrow(result3$models_last_run$x1$Sol), "\n")
cat("Expected: (1500-200)/5 =", (1500-200)/5, "\n")

cat("\n=== Success! All examples demonstrate the new functionality ===\n")
cat("\nKey features:\n")
cat("1. Backward compatible: Single values still work as before\n")
cat("2. Model-specific MCMC: Provide lists to customize each model\n")
cat("3. Mixed specification: Can mix single values and lists\n")
cat("4. Automatic validation: List lengths are checked against formula count\n")
