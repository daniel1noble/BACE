# =============================================================================
# 06_benchmark_globtherm.R - GlobTherm thermal limits
# =============================================================================
# 5 continuous traits: Tmax (degC), Tmin (degC), lat_max, long_max,
# elevation_max. NONE log-transformed: temperatures are linear and can be
# negative, lat/long are angular, elevation can be sea-level / negative.
# Tree is taxonomic (Class/Order/Family/Genus) with Grafen branch lengths.
# Smallest dataset (1969 spp) — full dataset fits inside the 2000 budget.
# Source: Bennett et al. 2018 Sci. Data 5:180022.
# =============================================================================

devtools::load_all(quiet = TRUE)
source("dev/benchmark_engine.R")

load("dev/testing_data/data/globtherm_traits.rda")
load("dev/testing_data/data/globtherm_tree.rda")

result <- benchmark_dataset(
  traits_df    = globtherm_traits,
  tree         = globtherm_tree,
  dataset_name = "globtherm",
  log_traits   = character(0),
  subset_n     = 2000L,
  nitt = 20000, thin = 15, burnin = 4000,
  runs = 5, n_final = 10,
  # n_cores=2 (not 4) to keep MCMCglmm parallel memory pressure under
  # the standard GHA runner's 7GB. globtherm's first cloud run hit
  # SIGTERM ~56min in, mid-final-imputation, after the convergence
  # phase had already passed -- consistent with peak-memory exhaustion
  # from 4 parallel MCMCglmm processes on 2000 spp + phylo random.
  max_attempts = 2, n_cores = 2L,
  verbose = TRUE
)
