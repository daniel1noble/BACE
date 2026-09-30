## =============================================================================
## Amphibian CTmax example, step 4: pre-flight check of phylogenetic signal
## =============================================================================
##
## Fits one phylogenetic mixed model per continuous variable, as recommended in
## the BACE tutorial before a long imputation run. Uses the experimental rows
## only, because the standardised rows carry no observed values. Run from the
## BACE repository root:
##
##   Rscript vignettes/amphibian_ctmax/04_phylo_signal.R
## =============================================================================

suppressMessages({
  pkgload::load_all(".", quiet = TRUE)   # the BACE source in this repository
  library(ape)
})
source("vignettes/amphibian_ctmax/00_setup.R")

ctmax_data <- ctmax_fill_medium(readRDS("vignettes/amphibian_ctmax/data/ctmax_data.rds"))
ctmax_tree <- readRDS("vignettes/amphibian_ctmax/data/ctmax_tree.rds")

exp_data <- ctmax_data[ctmax_data$row_type == "experimental", ctmax_vars]
exp_tree <- keep.tip(ctmax_tree, unique(exp_data$species))

set.seed(2026)
signal <- phylo_signal_summary(data = exp_data, tree = exp_tree, species_col = "species",
                               variables = c("CTmax", "ln_sd_CTmax", "ln_acclimation_time", "ramping"),
                               species = TRUE, nitt = 105000, burnin = 5000, thin = 50,
                               keep_models = FALSE, verbose = TRUE)
print(signal)
saveRDS(signal, "vignettes/amphibian_ctmax/results/phylo_signal.rds")
