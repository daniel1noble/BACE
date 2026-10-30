## =============================================================================
## Bird social behaviour example: shared settings, sourced by the run scripts
## =============================================================================

# Three analyses use the same data (vignettes/bird_social/data/):
#   full    impute the social behaviour of the species without a description in
#           BIRDBASE, using all known values (nothing hidden); these are the
#           imputations to use in downstream analyses
#   validation  check the imputations: as full, but with the known behaviour
#           of the held-out species (holdout.rds) also hidden, so that their
#           imputations can be scored (ordinal variable only)
#   mixed   check joint imputation: as validation, and also hide the hand-wing
#           index of 20% of species (hwi_mask.rds), so that BACE imputes an
#           ordinal and a continuous variable that predict each other
analyses <- c("full", "validation", "mixed")

# Imputation models. In the validation analysis (run first) the formula is a
# single string, so BACE also builds (and fits) models for the complete predictors;
# they have nothing to impute. The full analysis gives the same model as a
# one-element list, which fits only the social model and saves time; the social
# model itself is identical.
validation_formula <- "social_max ~ log_mass + hwi"
full_formulas  <- list("social_max ~ log_mass + hwi")
mixed_formulas <- list("social_max ~ log_mass + hwi", "hwi ~ log_mass + social_max")

# MCMC settings. The chained phase only needs point imputations to settle; the
# final imputations need long chains for the phylogenetic variance and the
# cutpoints of the ordinal (threshold) model. In the mixed analysis the second
# element of each setting applies to the hand-wing-index model, which is
# Gaussian and mixes quickly.
chain_settings <- list(full   = list(runs = 5, nitt = 50000, burnin = 10000, thin = 40),
                       validation = list(runs = 5, nitt = 50000, burnin = 10000, thin = 40),
                       mixed  = list(runs = 25, nitt = list(50000, 20000),
                                     burnin = list(10000, 5000), thin = list(40, 15)))
final_settings <- list(full   = list(nitt = 500000, burnin = 50000, thin = 450),
                       validation = list(nitt = 500000, burnin = 50000, thin = 450),
                       mixed  = list(nitt = list(500000, 100000), burnin = list(50000, 10000),
                                     thin = list(450, 90)))
n_final <- c(full = 50, validation = 20, mixed = 20)

# Seeds: chained phase, and final imputation i uses final_seed + i
chain_seed <- c(full = 3001, validation = 5001, mixed = 8001)
final_seed <- c(full = 3100, validation = 5100, mixed = 8100)

## -----------------------------------------------------------------------------
## Data for one analysis
## -----------------------------------------------------------------------------
# Returns the BACE input (species, social_max, log_mass, hwi), the tree, and the
# rows whose social behaviour or hand-wing index is missing in this analysis.
make_analysis_data <- function(analysis) {
  social_data <- readRDS("vignettes/bird_social/data/social_data.rds")
  social_tree <- readRDS("vignettes/bird_social/data/social_tree.rds")
  holdout     <- readRDS("vignettes/bird_social/data/holdout.rds")
  hwi_mask    <- readRDS("vignettes/bird_social/data/hwi_mask.rds")

  d <- social_data[, c("species", "social_max", "log_mass", "hwi")]
  if (analysis != "full") d$social_max[d$species %in% holdout$species] <- NA
  if (analysis == "mixed") d$hwi[d$species %in% hwi_mask] <- NA

  list(data = d, tree = social_tree,
       social_rows = which(is.na(d$social_max)),                 # imputed by BACE
       scored_rows = if (analysis == "full") integer(0) else match(holdout$species, d$species),
       hwi_rows    = which(is.na(d$hwi)),
       formulas    = switch(analysis, full = full_formulas, validation = validation_formula, mixed = mixed_formulas))
}
