# =============================================================================
# 04_benchmark_amphibio.R - AmphiBIO amphibians
# =============================================================================
# 2 continuous (log-transformed), 2 binary (diurnal, nocturnal), 1 ordinal
# (diet_breadth K=5), 1 categorical (habitat K=4).
# Source: Oliveira et al. 2017 Sci. Data 4:170123.
#
# Caveats: diurnal / nocturnal are presence-only in AmphiBIO (NA = no
# record, not necessarily absence). body_mass_g has ~91% missingness.
# The benchmark masks an extra 10% of *observed* cells regardless.
# =============================================================================

devtools::load_all(quiet = TRUE)
source("dev/benchmark_engine.R")

load("dev/testing_data/data/amphibio_traits.rda")
load("dev/testing_data/data/amphibio_tree.rda")

LOG_TRAITS <- c("body_size_mm", "body_mass_g")

result <- benchmark_dataset(
  traits_df    = amphibio_traits,
  tree         = amphibio_tree,
  dataset_name = "amphibio",
  log_traits   = LOG_TRAITS,
  subset_n     = 2000L,
  nitt = 20000, thin = 15, burnin = 4000,
  runs = 5, n_final = 10,
  max_attempts = 2, n_cores = 4L,
  verbose = TRUE
)
