# Test script for Examples 1-3 only (no random slopes)
library(BACE)
library(ape)

set.seed(123)
phylo <- rtree(10)
phylo <- compute.brlen(phylo, method = "Grafen")

n_species <- 10
n_reps <- 3
n_total <- n_species * n_reps

data <- data.frame(
  body_mass = rnorm(n_total, mean = 100, sd = 20),
  habitat = factor(rep(c("Forest", "Grassland", "Desert"), length.out = n_total)),
  temperature = rnorm(n_total, mean = 25, sd = 5),
  precipitation = rnorm(n_total, mean = 1000, sd = 300),
  Species = rep(phylo$tip.label, each = n_reps)
)

data$body_mass[sample(1:n_total, 5)] <- NA
data$temperature[sample(1:n_total, 3)] <- NA
data$precipitation[sample(1:n_total, 4)] <- NA

cat("\n=== Example 2: Default Behavior (species = FALSE) ===\n")
result_default <- bace_imp(
  fixformula = "body_mass ~ habitat + temperature + precipitation",
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  nitt = 2000,
  thin = 5,
  burnin = 500,
  runs = 2,
  species = FALSE,
  verbose = TRUE
)
cat("Example 2 completed successfully!\n")

cat("\n=== Example 3: Decompose Phylo + Species Effects (species = TRUE) ===\n")
result_decomposed <- bace_imp(
  fixformula = "body_mass ~ habitat + temperature + precipitation",
  ran_phylo_form = "~ 1 | Species",
  phylo = phylo,
  data = data,
  nitt = 2000,
  thin = 5,
  burnin = 500,
  runs = 2,
  species = TRUE,
  verbose = TRUE
)
cat("Example 3 completed successfully!\n")
cat("\nAll tests passed!\n")
