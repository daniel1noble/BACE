# =============================================================================
# 04_benchmark_amphibio.R - AmphiBIO amphibians
# =============================================================================
# Trait set aligned with Shinichi's pigauto cross-dataset bench
# (2026-05-04 spec): 2 continuous traits only (Body_size_mm,
# Body_mass_g). Discrete columns (diurnal / nocturnal binary,
# diet_breadth ordinal, habitat categorical) are intentionally
# skipped here -- they were skipped in pigauto's v1 amphibio bench
# because the threshold-joint baseline hits an Rphylopars
# singular-matrix error on the AmphiBIO taxonomic tree, AND
# diurnal / nocturnal are presence-only ("yes" recorded vs NA = no
# record) which makes MCMCglmm threshold fitting on a single
# observed level fail with "Mixed model equations singular".
# Source: Oliveira et al. 2017 Sci. Data 4:170123.
# =============================================================================

devtools::load_all(quiet = TRUE)
source("dev/benchmark_engine.R")

load("dev/testing_data/data/amphibio_traits.rda")
load("dev/testing_data/data/amphibio_tree.rda")

LOG_TRAITS   <- c("body_size_mm", "body_mass_g")
TRAIT_SUBSET <- c("body_size_mm", "body_mass_g")

result <- benchmark_dataset(
  traits_df    = amphibio_traits,
  tree         = amphibio_tree,
  dataset_name = "amphibio",
  log_traits   = LOG_TRAITS,
  trait_subset = TRAIT_SUBSET,
  subset_n     = 2000L,
  nitt = 20000, thin = 15, burnin = 4000,
  runs = 5, n_final = 10,
  max_attempts = 2, n_cores = 4L,
  verbose = TRUE
)
