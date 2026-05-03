# =============================================================================
# 03_benchmark_pantheria.R - PanTHERIA mammals
# =============================================================================
# 4 continuous (log-transformed), 1 count (litter_size, integer -> Poisson),
# 2 ordinal (diet/habitat breadth), 1 binary (terrestriality).
# Source: Jones et al. 2009 Ecology 90:2648.
# =============================================================================

devtools::load_all(quiet = TRUE)
source("dev/benchmark_engine.R")

load("dev/testing_data/data/pantheria_traits.rda")
load("dev/testing_data/data/pantheria_tree.rda")

# Body mass / length / gestation / longevity span 5+ orders of magnitude
# in mammals (shrew to whale) and are conventionally log-transformed
# (Felsenstein 1985 Am Nat 125:1; Garland et al. 1993 Syst Biol 42:265).
LOG_TRAITS <- c("body_mass_g", "head_body_length_mm",
                "gestation_d",  "max_longevity_m")

result <- benchmark_dataset(
  traits_df    = pantheria_traits,
  tree         = pantheria_tree,
  dataset_name = "pantheria",
  log_traits   = LOG_TRAITS,
  subset_n     = 2000L,
  nitt = 20000, thin = 15, burnin = 4000,
  runs = 5, n_final = 10,
  max_attempts = 2, n_cores = 4L,
  verbose = TRUE
)
