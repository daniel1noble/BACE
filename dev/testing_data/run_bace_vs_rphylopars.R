
## Compare BACE (bace_imp) vs Rphylopars on AVONET (2000 spp) with artificial masking
## Inputs (prepared by ChatGPT):
##  - avonet_2000_masked.csv : masked trait matrix (species + 8 continuous traits)
##  - avonet_2000_truth.csv  : held-out true values (species, trait, true_value)
##  - Hackett_tree_2000.tre  : pruned Hackett phylogeny with exactly those 2000 tips
##
## Outputs:
##  - bace_last.csv, rphylopars_recon.csv
##  - compare_metrics_by_trait.csv, compare_predictions_long.csv

## ---- 0) Packages ----
pkgs <- c("ape", "dplyr", "tidyr", "readr", "MCMCglmm", "Rphylopars")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) {
  stop("Missing packages: ", paste(to_install, collapse=", "),
       "\nInstall them first, e.g. install.packages(c(", paste0('"', to_install, '"', collapse=", "), ")).")
}

library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(MCMCglmm)
library(Rphylopars)
library(here)
library(BACE)

## ---- 1) Source BACE functions ----
## Adjust these paths if you have the package installed; this is for running from the provided scripts.
# source("/mnt/data/prep_functions.R")
# source("/mnt/data/build_functions.R")
# source("/mnt/data/model_functions.R")
# source("/mnt/data/bace_imp.R")

## ---- 2) Load prepared data + tree ----
dat <- read_csv(here::here("dev", "testing_data", "avonet_2000_masked.csv"), show_col_types = FALSE)
truth <- read_csv(here::here("dev", "testing_data", "avonet_2000_truth.csv"), show_col_types = FALSE)
tree <- read.tree(here::here("dev", "testing_data", "Hackett_tree_2000.tre"))

stopifnot(all(dat$species %in% tree$tip.label), length(tree$tip.label) == nrow(dat))

## Ensure the order matches tree tips (helps some methods)
dat <- dat %>% slice(match(tree$tip.label, species))
dat <- as.data.frame(dat)

traits <- setdiff(names(dat), "species")

## ---- 3) Run BACE imputation ----
set.seed(1)

## Use any variable as the LHS; bace_imp will generate chained formulas for all vars in the formula
fixformula <- paste0("Mass ~ ", paste(setdiff(traits, "Mass"), collapse = " + "))

fit_bace <- bace_imp(
  fixformula     = fixformula,
  ran_phylo_form = "~ 1 | species",
  phylo          = tree,
  data           = dat,
  runs           = 5,       
  nitt           = 4000,    
  burnin         = 1000,
  thin           = 10,
  verbose        = TRUE
)

bace_last <- fit_bace$data[[length(fit_bace$data)]]
bace_last <- as.data.frame(bace_last)
write_csv(bace_last, "bace_last.csv")

## ---- 4) Run Rphylopars imputation ----
## Rphylopars expects rownames = species and a trait matrix/data.frame of numeric traits
trait_mat <- dat %>%
  select(-species) %>%
  as.data.frame()
rownames(trait_mat) <- dat$species

rp <- phylopars(trait_data = trait_mat, tree = tree, model = "BM", pheno_error = FALSE)

## Extract reconstructed tip traits (includes imputations)
## Different versions expose slightly different slot names; handle both.
if (!is.null(rp$tip_recon)) {
  rp_recon <- rp$tip_recon
} else if (!is.null(rp$trait_recon)) {
  rp_recon <- rp$trait_recon
} else {
  stop("Could not find reconstructed tip traits in phylopars output. Try str(rp) to locate the matrix.")
}

rp_recon <- as.data.frame(rp_recon)
rp_recon$species <- rownames(rp_recon)
rp_recon <- rp_recon %>% select(species, all_of(traits))
write_csv(rp_recon, "rphylopars_recon.csv")

## ---- 5) Compare on held-out (masked) entries ----
## Long-form predictions for held-out cells
bace_long <- bace_last %>%
  select(species, all_of(traits)) %>%
  pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "bace_pred")

rp_long <- rp_recon %>%
  pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "rphy_pred")

preds <- truth %>%
  left_join(bace_long, by = c("species_tip" = "species", "trait" = "trait")) %>%
  left_join(rp_long,   by = c("species_tip" = "species", "trait" = "trait")) %>%
  rename(species = species_tip)

## Metrics by trait
metrics <- preds %>%
  group_by(trait) %>%
  summarise(
    n = n(),
    ## BACE
    cor_bace  = cor(true_value, bace_pred, use = "complete.obs"),
    rmse_bace = sqrt(mean((bace_pred - true_value)^2, na.rm = TRUE)),
    mae_bace  = mean(abs(bace_pred - true_value), na.rm = TRUE),
    ## Rphylopars
    cor_rphy  = cor(true_value, rphy_pred, use = "complete.obs"),
    rmse_rphy = sqrt(mean((rphy_pred - true_value)^2, na.rm = TRUE)),
    mae_rphy  = mean(abs(rphy_pred - true_value), na.rm = TRUE),
    ## Direct agreement between methods on held-out cells
    cor_methods = cor(bace_pred, rphy_pred, use = "complete.obs"),
    .groups = "drop"
  ) %>% arrange(trait)

write_csv(metrics, "compare_metrics_by_trait.csv")
write_csv(preds, "compare_predictions_long.csv")

print(metrics)
