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
  # n_cores=1 (sequential) for taxonomy-built-tree datasets. n_cores>=2
  # caused GHA runner shutdowns mid-final-imputation in two consecutive
  # runs (25287097270, 25294935456) at the same Step-3 phase: parallel
  # R workers don't stream log output to the parent process, GHA's
  # runner detects no stdout heartbeat and kills the runner. Sequential
  # avoids this -- every per-iteration message appears in main stdout.
  # Cost is ~2x final-imputation time, still well inside the 5h45m
  # workflow budget for 5 continuous traits.
  max_attempts = 2, n_cores = 1L,
  verbose = TRUE
)
