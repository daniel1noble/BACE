## Compare BACE (bace_imp) vs Rphylopars on AVONET (2000 spp) with artificial masking
## Includes: point metrics + 95% interval coverage + interval width
##
## Inputs (expected in dev/testing_data/):
##  - avonet_2000_masked.csv : masked trait matrix (species + continuous traits)
##  - avonet_2000_truth.csv  : held-out true values (species_tip, trait, true_value)
##  - Hackett_tree_2000.tre  : pruned phylogeny with exactly those 2000 tips
##
## Assumptions:
##  - bace_imp() returns fit_bace$pred_list_last_run[[trait]] with columns:
##      post_mean, post_sd, ci_lower, ci_upper
##  - dat is reordered to tree tips, so per-observation prediction rows match dat$species

## ---- 0) Packages ----
pkgs <- c("ape", "dplyr", "tidyr", "readr", "MCMCglmm", "Rphylopars", "here")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) {
  stop("Missing packages: ", paste(to_install, collapse = ", "),
       "\nInstall them first, e.g. install.packages(c(",
       paste0('"', to_install, '"', collapse = ", "), ")).")
}

library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(MCMCglmm)
library(Rphylopars)
library(here)

## ---- 1) Source BACE functions (edit paths as needed) ----
## If you are running from scripts in your repo:
# source("/mnt/data/prep_functions.R")
# source("/mnt/data/build_functions.R")
# source("/mnt/data/model_functions.R")
# source("/mnt/data/bace_imp.R")

## If you have BACE installed:
# library(BACE)

## ---- 2) Load prepared data + tree ----
dat   <- read_csv(here::here("dev", "testing_data", "avonet_2000_masked.csv"),
                  show_col_types = FALSE)
truth <- read_csv(here::here("dev", "testing_data", "avonet_2000_truth.csv"),
                  show_col_types = FALSE)
tree  <- read.tree(here::here("dev", "testing_data", "Hackett_tree_2000.tre"))

stopifnot(all(dat$species %in% tree$tip.label),
          length(tree$tip.label) == nrow(dat))

## Ensure order matches tree tips (helps both methods, and crucial for BACE CI join)
dat <- dat %>% slice(match(tree$tip.label, species))
dat <- as.data.frame(dat)

traits <- setdiff(names(dat), "species")

## Optional: ensure branch lengths exist / are compatible
tree <- ape::compute.brlen(tree, method = "Grafen")

## ---- 3) Run BACE imputation ----
set.seed(1)

## Use any variable as the LHS; bace_imp will generate chained formulas for all vars in the formula
fixformula <- paste0("Mass ~ ", paste(setdiff(traits, "Mass"), collapse = " + "))

fit_bace <- bace_imp(
  fixformula     = fixformula,
  ran_phylo_form = "~ 1 | species",
  phylo          = tree,
  data           = dat,
  runs           = 10,
  nitt           = 4000 * 100,
  burnin         = 1000 * 100,
  thin           = 10 * 50,
  verbose        = TRUE
)

## Last completed dataset (point imputations)
bace_last <- fit_bace$data[[length(fit_bace$data)]]
bace_last <- as.data.frame(bace_last)

## ---- 4) Run Rphylopars imputation ----
## Rphylopars expects a data.frame with a species column + numeric traits
trait_df <- dat %>%
  dplyr::select(species, dplyr::everything()) %>%
  as.data.frame()

rp <- Rphylopars::phylopars(
  trait_data  = trait_df,
  tree        = tree,
  model       = "BM",
  pheno_error = FALSE
)

## Robustly extract reconstructed TIP traits from a phylopars fit
get_tip_recon_phylopars <- function(rp, what = c("mean", "var")) {
  what <- match.arg(what)
  
  if (is.null(rp$tree) || is.null(rp$tree$tip.label)) {
    stop("rp$tree$tip.label not found; cannot identify tips.")
  }
  tip_labels <- rp$tree$tip.label
  
  mat <- switch(
    what,
    mean = rp$anc_recon,
    var  = rp$anc_var
  )
  if (is.null(mat)) {
    stop(sprintf("Could not find rp$anc_%s in phylopars output.", what))
  }
  if (is.null(rownames(mat))) {
    stop("Reconstruction matrix has no rownames; cannot match to tip labels.")
  }
  
  idx <- match(tip_labels, rownames(mat))
  if (anyNA(idx)) {
    missing <- tip_labels[is.na(idx)]
    stop("Some tip labels were not found in reconstruction matrix rownames. Examples: ",
         paste(utils::head(missing, 10), collapse = ", "))
  }
  
  mat[idx, , drop = FALSE]
}

rp_tip_recon <- get_tip_recon_phylopars(rp, what = "mean")  # 2000 x p
rp_tip_var   <- get_tip_recon_phylopars(rp, what = "var")   # 2000 x p

rp_recon <- as.data.frame(rp_tip_recon)
rp_recon$species <- rownames(rp_tip_recon)
rp_recon <- rp_recon %>% dplyr::select(species, all_of(traits))

## ---- 4b) Rphylopars SE + CI (normal approx) ----
rp_var <- as.data.frame(rp_tip_var)
rp_var$species <- rownames(rp_tip_var)
rp_var <- rp_var %>% dplyr::select(species, all_of(traits))

rp_long <- rp_recon %>%
  tidyr::pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "rphy_mean")

rpv_long <- rp_var %>%
  tidyr::pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "rphy_var")

rp_long <- rp_long %>%
  dplyr::left_join(rpv_long, by = c("species", "trait")) %>%
  dplyr::mutate(
    rphy_se  = sqrt(pmax(rphy_var, 0)),
    rphy_lwr = rphy_mean - 1.96 * rphy_se,
    rphy_upr = rphy_mean + 1.96 * rphy_se
  )

## ---- 5) Compare on held-out (masked) entries + COVERAGE ----

## BACE: extract per-trait fitted summaries (mean/sd/CI) from last run
bace_sum_long <- do.call(rbind, lapply(traits, function(tr) {
  s <- fit_bace$pred_list_last_run[[tr]]
  if (is.null(s)) stop("Missing fit_bace$pred_list_last_run[[", tr, "]]")
  
  data.frame(
    species   = dat$species,
    trait     = tr,
    bace_mean = s$post_mean,
    bace_sd   = s$post_sd,
    bace_lwr  = s$ci_lower,
    bace_upr  = s$ci_upper
  )
}))

## BACE point-imputed values (optional, for point accuracy metrics)
bace_point_long <- bace_last %>%
  dplyr::select(species, all_of(traits)) %>%
  tidyr::pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "bace_pred")

preds <- truth %>%
  dplyr::rename(species = species_tip) %>%
  dplyr::left_join(bace_point_long, by = c("species", "trait")) %>%
  dplyr::left_join(bace_sum_long,   by = c("species", "trait")) %>%
  dplyr::left_join(rp_long,         by = c("species", "trait"))

metrics <- preds %>%
  dplyr::group_by(trait) %>%
  dplyr::summarise(
    n = dplyr::n(),
    
    ## Accuracy (point)
    cor_bace  = cor(true_value, bace_pred, use = "complete.obs"),
    rmse_bace = sqrt(mean((bace_pred - true_value)^2, na.rm = TRUE)),
    mae_bace  = mean(abs(bace_pred - true_value), na.rm = TRUE),
    
    cor_rphy  = cor(true_value, rphy_mean, use = "complete.obs"),
    rmse_rphy = sqrt(mean((rphy_mean - true_value)^2, na.rm = TRUE)),
    mae_rphy  = mean(abs(rphy_mean - true_value), na.rm = TRUE),
    
    ## Coverage (95% intervals)
    cov_bace = mean(true_value >= bace_lwr & true_value <= bace_upr, na.rm = TRUE),
    cov_rphy = mean(true_value >= rphy_lwr & true_value <= rphy_upr, na.rm = TRUE),
    
    ## Interval width (sharpness)
    width_bace = mean(bace_upr - bace_lwr, na.rm = TRUE),
    width_rphy = mean(rphy_upr - rphy_lwr, na.rm = TRUE),
    
    ## Agreement between methods on held-out cells
    cor_methods = cor(bace_pred, rphy_mean, use = "complete.obs"),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(trait)

print(metrics)

## Optional: save outputs
# readr::write_csv(metrics, "compare_metrics_by_trait_with_coverage.csv")
# readr::write_csv(preds,   "compare_predictions_long_with_intervals.csv")