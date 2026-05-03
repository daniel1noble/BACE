# =============================================================================
# 05_benchmark_bien.R - BIEN plant traits
# =============================================================================
# 5 continuous traits, all log-transformed (height, leaf area, SLA, seed
# mass, wood density — all conventionally log-scaled in plant ecology).
# Largest bundled dataset: 19,109 species; subset to 2000 for the
# routine-run benchmark.
# Source: BIEN (botanicalinformation.org); tree from V.PhyloMaker2
# scenario.3 (Jin & Qian 2022 Ecography 45:e06141).
# =============================================================================

devtools::load_all(quiet = TRUE)
source("dev/benchmark_engine.R")

load("dev/testing_data/data/bien_traits.rda")
load("dev/testing_data/data/bien_tree.rda")

LOG_TRAITS <- c("height_m", "leaf_area", "sla", "seed_mass", "wood_density")

result <- benchmark_dataset(
  traits_df    = bien_traits,
  tree         = bien_tree,
  dataset_name = "bien",
  log_traits   = LOG_TRAITS,
  subset_n     = 2000L,
  nitt = 20000, thin = 15, burnin = 4000,
  runs = 5, n_final = 10,
  max_attempts = 2, n_cores = 4L,
  verbose = TRUE
)
