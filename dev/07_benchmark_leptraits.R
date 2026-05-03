# =============================================================================
# 07_benchmark_leptraits.R - LepTraits lepidopteran traits
# =============================================================================
# 4 continuous (wingspan, forewing length, flight duration, n hostplant
# families) + 12 monthly flight indicators (Jan..Dec, 0/1).
# The monthly indicators are stored as integer in the .rda but are really
# binary 0/1 — coerce to factor here so the engine treats them as binary
# (threshold model) rather than count (Poisson).
# Tree is taxonomic (Family/Genus/Species) with Grafen branch lengths
# because LepTraits does not ship a phylogeny.
# Source: Shirey et al. 2022 Sci. Data 9:382.
# =============================================================================

devtools::load_all(quiet = TRUE)
source("dev/benchmark_engine.R")

load("dev/testing_data/data/leptraits_traits.rda")
load("dev/testing_data/data/leptraits_tree.rda")

# Coerce monthly indicators from integer 0/1 to factor "0"/"1" so the
# engine routes them through the binary (threshold) path, not Poisson.
month_cols <- c("Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec")
for (m in intersect(month_cols, colnames(leptraits_traits))) {
  leptraits_traits[[m]] <- factor(leptraits_traits[[m]], levels = c(0L, 1L))
}

# Wingspan and forewing span >2 orders of magnitude across butterflies/moths.
# n_hostplant_families is a small count (0-30) — log1p so the count's
# heavy right tail doesn't dominate the gaussian fit.
LOG_TRAITS <- c("wingspan_lower", "forewing_length_lower",
                "n_hostplant_families")

result <- benchmark_dataset(
  traits_df    = leptraits_traits,
  tree         = leptraits_tree,
  dataset_name = "leptraits",
  log_traits   = LOG_TRAITS,
  subset_n     = 2000L,
  nitt = 20000, thin = 15, burnin = 4000,
  runs = 5, n_final = 10,
  max_attempts = 2, n_cores = 4L,
  verbose = TRUE
)
