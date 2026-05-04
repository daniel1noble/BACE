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

# Mass / morphometrics / range_size span 5+ orders of magnitude in
# birds (Tobias et al. 2022; Felsenstein 1985 Am Nat 125:1) and are
# conventionally log-transformed before phylogenetic analysis.
# Latitude / longitude are angular -- not log-transformed.
LOG_TRAITS <- c("mass_g", "wing_length_mm", "beak_length_culmen_mm",
                "tarsus_length_mm", "tail_length_mm", "range_size_km2")

# Use the historical 10-trait benchmark scope (8 continuous + 2
# categorical). The full 13-trait set (adding habitat_density /
# migration ordinals + habitat K=11 categorical) ran 5h45m past the
# GHA workflow timeout in run 25287097270 -- chained equations
# scale O(n_traits) per iteration and the K=11 habitat OVR pulls
# 11 binary threshold fits per chain step. The historical 10 traits
# already cover continuous + categorical paths cleanly.
TRAIT_SUBSET <- c(LOG_TRAITS, "centroid_lat", "centroid_lon",
                  "trophic_level", "primary_lifestyle")

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
