# =============================================================================
# 00_benchmark_AVONET.R - AVONET birds (Tobias et al. 2022)
# =============================================================================
# 8 continuous (mass / morphometric / range / centroid) + 2 ordinal
# (habitat_density, migration, K=3 each) + 3 categorical (trophic_level
# K=4, primary_lifestyle K=5, habitat K=11) + tax columns (dropped).
# Phylogeny: Hackett et al. 2008 Science 320:1763.
#
# Migrated to the cross-dataset benchmark engine in commit 487c506.
# Previous monolithic implementation is in git history.
#
# Note on trait scope: the .rda built by data-raw/make_avonet.R adds
# habitat_density / migration / habitat to the eight continuous +
# trophic_level / primary_lifestyle that the historical AVONET
# benchmark used. Including them gives a richer trait-type mix
# (continuous + ordinal + categorical) at the cost of more MCMC
# work per call. Set `trait_subset` below to revert to the historical
# 10-trait scope if you want metric-comparable runs.
# =============================================================================

devtools::load_all(quiet = TRUE)
source("dev/benchmark_engine.R")

load("dev/testing_data/data/avonet_traits.rda")
load("dev/testing_data/data/avonet_tree.rda")

# Trait set aligned with Shinichi's pigauto cross-dataset bench
# (2026-05-04 spec): 4 continuous morphometrics + 2 categorical +
# 1 ordinal. 7 traits total, makes BACE's results directly
# comparable to the pigauto AVONET headline.
#
# Caveat on Primary.Lifestyle: Shinichi's pigauto bench reports
# K=7 levels but the avonet_traits.rda built from this repo's
# AVONET.csv has K=5 (Aerial / Aquatic / Generalist / Insessorial /
# Terrestrial). Different AVONET version. BACE's results will not
# be exactly comparable on this trait until both pipelines align
# on the same source CSV.
LOG_TRAITS <- c("mass_g", "wing_length_mm",
                "beak_length_culmen_mm", "tarsus_length_mm")

TRAIT_SUBSET <- c(LOG_TRAITS,
                  "trophic_level",      # categorical K=4
                  "primary_lifestyle",  # categorical K=5 (pigauto K=7)
                  "migration")          # ordinal K=3

result <- benchmark_dataset(
  traits_df    = avonet_traits,
  tree         = avonet_tree,
  dataset_name = "avonet",
  log_traits   = LOG_TRAITS,
  trait_subset = TRAIT_SUBSET,
  subset_n     = 2000L,
  nitt = 20000, thin = 15, burnin = 4000,
  runs = 5,    n_final = 10,
  max_attempts = 2, n_cores = 4L,
  verbose = TRUE
)
