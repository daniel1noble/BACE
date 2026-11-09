# Example Usage: Species Random Effects Decomposition in BACE
# This script demonstrates the new species random effects feature

library(BACE)
library(ape)

# ==============================================================================
# Example 1: Create Data with Replicated Species
# ==============================================================================

set.seed(123)

# Create a phylogenetic tree with 10 species
phylo <- rtree(10)
phylo <- compute.brlen(phylo, method = "Grafen")

# Create dataset with 3 replicates per species (30 total observations)
n_species <- 10
n_reps <- 3
n_total <- n_species * n_reps

data <- data.frame(
  # Response variable
  body_mass = rnorm(n_total, mean = 100, sd = 20),
  
  # Predictors
  habitat = factor(rep(c("Forest", "Grassland", "Desert"), length.out = n_total)),
  temperature = rnorm(n_total, mean = 25, sd = 5),
  precipitation = rnorm(n_total, mean = 1000, sd = 300),
  
  # Species identifier (replicated)
  Species = rep(phylo$tip.label, each = n_reps)
)

# Introduce missing data
data$body_mass[sample(1:n_total, 5)] <- NA
data$temperature[sample(1:n_total, 3)] <- NA
data$precipitation[sample(1:n_total, 4)] <- NA

# View replication structure
cat("Species replication:\n")
print(table(data$Species))

# ==============================================================================
# Example 2: Default Behavior (species = FALSE)
# ==============================================================================

cat("\n\n=== Example 2: Default Behavior (Phylogenetic effect only) ===\n")

result_default <- bace_imp(
  fixformula = "body_mass ~ habitat + temperature + precipitation",
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  nitt = 30000,
  thin = 5,
  burnin = 	5000,
  runs = 5,
  species = FALSE,  # Default: only phylogenetic random effect
  verbose = TRUE
)

summary(result_default$models_last_run[[1]])
cat("\nModel structure (default):\n")
cat("- Random effects: Phylogenetic effect only\n")
cat("- Number of random effects: 1\n")

# ==============================================================================
# Example 3: Decompose Phylogenetic and Species Effects (species = TRUE)
# ==============================================================================

cat("\n\n=== Example 3: Decompose Phylo + Species Effects ===\n")

result_decomposed <- bace_imp(
  fixformula = "body_mass ~ habitat + temperature + precipitation",
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  nitt = 30000,
  thin = 5,
  burnin = 	5000,
  runs = 5,
  species = TRUE,  # Decompose into phylo + non-phylo species effects
  verbose = TRUE
)
summary(result_decomposed$models_last_run[[1]])
cat("\nModel structure (decomposed):\n")
cat("- Random effects: Phylogenetic + Non-phylogenetic species effects\n")
cat("- Number of random effects: 2\n")
cat("- Phylo effect: Uses phylogenetic covariance matrix\n")
cat("- Species effect: Uses identity matrix (uncorrelated)\n")

# ==============================================================================
# Example 4: Random Slopes on Phylogenetic Effect
# ==============================================================================

cat("\n\n=== Example 4: Random Slopes (Phylo only) ===\n")

# Note: Random slopes only applied to phylogenetic effect
# Non-phylogenetic species effect always uses random intercept
result_int <- bace_imp(
  fixformula = "body_mass ~ habitat + temperature",
  ran_phylo_form = "~ us(1) | Species",  # Random slope for temperature
  phylo = phylo,
  data = data,
  nitt = 1000,
  thin = 5,
  burnin = 200,
  runs = 5,
  species = FALSE,
  verbose = TRUE
)
summary(result_int$models_last_run[[1]])
cat("\nModel structure (random slopes):\n")
cat("- Phylogenetic effect: Random intercept + random slope for temperature\n")
cat("- Species effect: Random intercept only\n")

# ==============================================================================
# Example 5: Multiple Variables with Missing Data
# ==============================================================================

cat("\n\n=== Example 5: Multiple Variables (Formula List) ===\n")

# When multiple variables have missing data, provide a list of formulas
result_multiple <- bace_imp(
  fixformula = list(
    "body_mass ~ habitat + temperature + precipitation",
    "temperature ~ habitat + precipitation",
    "precipitation ~ habitat + temperature"
  ),
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  nitt = 1000,
  thin = 5,
  burnin = 200,
  runs = 5,
  species = TRUE,
  verbose = TRUE
)

cat("\nMultiple variable imputation completed\n")

# ==============================================================================
# Example 6: Insufficient Replication (Will Throw Error)
# ==============================================================================

cat("\n\n=== Example 6: Insufficient Replication (Error Example) ===\n")

# Create data with NO replication (each species appears once)
data_no_reps <- data.frame(
  body_mass = rnorm(10, mean = 100, sd = 20),
  temperature = rnorm(10, mean = 25, sd = 5),
  Species = phylo$tip.label
)
data_no_reps$body_mass[c(1, 3)] <- NA

# This will throw an error
tryCatch({
  result_error <- bace_imp(
    fixformula = "body_mass ~ temperature",
    ran_phylo_form = "~ 1 | Species",
    phylo = phylo,
    data = data_no_reps,
    runs = 1,
    species = TRUE  # ERROR: No replicated species
  )
}, error = function(e) {
  cat("\nExpected Error:\n")
  cat(paste("  ", e$message, "\n"))
})

# ==============================================================================
# Example 7: Low Replication (Will Show Warning)
# ==============================================================================

cat("\n\n=== Example 7: Low Replication (Warning Example) ===\n")

# Create data where only 2 out of 10 species are replicated (20%)
species_vec <- c(
  rep(phylo$tip.label[1:8], each = 1),  # 8 species, 1 obs each
  rep(phylo$tip.label[9:10], each = 2)  # 2 species, 2 obs each
)

data_low_reps <- data.frame(
  body_mass = rnorm(12, mean = 100, sd = 20),
  temperature = rnorm(12, mean = 25, sd = 5),
  Species = species_vec
)
data_low_reps$body_mass[c(1, 3)] <- NA

# This will show a warning
result_warning <- bace_imp(
  fixformula = "body_mass ~ temperature",
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data_low_reps,
  nitt = 500,
  thin = 2,
  burnin = 100,
  runs = 2,
  species = TRUE,  # WARNING: Low replication
  verbose = TRUE
)

cat("\nCompleted with warning about low replication\n")

# ==============================================================================
# Example 8: Accessing Results
# ==============================================================================

cat("\n\n=== Example 8: Accessing Results ===\n")

# Access imputed data from different iterations
cat("\nAvailable iterations:\n")
print(names(result_decomposed$data))

# Final imputed dataset
final_data <- result_decomposed$data[[length(result_decomposed$data)]]
cat("\nFinal imputed data (first 6 rows):\n")
print(head(final_data))

# Check if all missing values are imputed
cat("\nMissing values remaining:\n")
cat("  body_mass:", sum(is.na(final_data$body_mass)), "\n")
cat("  temperature:", sum(is.na(final_data$temperature)), "\n")
cat("  precipitation:", sum(is.na(final_data$precipitation)), "\n")

# Access model diagnostics
if (!is.null(result_decomposed$diagnostics)) {
  cat("\nDiagnostics available:\n")
  print(result_decomposed$diagnostics)
}

# ==============================================================================
# Summary and Recommendations
# ==============================================================================

cat("\n\n=== Summary and Recommendations ===\n\n")

cat("When to use species = TRUE:\n")
cat("1. You have replicated observations per species (ideally >30% of species)\n")
cat("2. You want to partition variance into phylogenetic vs non-phylogenetic components\n")
cat("3. You have sufficient statistical power (more replicated species = better)\n\n")

cat("When to use species = FALSE (default):\n")
cat("1. Each species appears only once in your dataset\n")
cat("2. You only care about phylogenetic signal (not decomposition)\n")
cat("3. You want simpler models with fewer parameters\n\n")

cat("Random slopes:\n")
cat("- Can be specified in ran_phylo_form (e.g., '~ us(1 + x) | Species')\n")
cat("- Only applied to phylogenetic effect when species = TRUE\n")
cat("- Non-phylogenetic species effect always uses random intercept only\n\n")

cat("Script completed successfully!\n")
