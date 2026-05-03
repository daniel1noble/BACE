## Compare BACE (bace_imp) vs Rphylopars on AVONET (2000 spp) with artificial masking
## Includes: point metrics + 95% interval coverage + interval width
##
## Inputs:
##  - dev/testing_data/data/avonet_traits.rda : 9993 spp x 16 cols, built
##                                              by data-raw/make_avonet.R
##  - dev/testing_data/data/avonet_tree.rda   : 9993-tip Hackett phylogeny
## We sample 2000 species at runtime (set.seed-controlled), keep the
## 8 continuous AVONET morphometric / range / centroid traits, then
## hide a fixed 10% of cells per column. Held-out cells are stored in
## `truth` for downstream scoring.
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

## ---- 2) Load bundled data + tree, sample 2000 spp, mask 10% per col ----
load(here::here("dev", "testing_data", "data", "avonet_traits.rda"))
load(here::here("dev", "testing_data", "data", "avonet_tree.rda"))

# AVONET-native names so the rest of the script reads naturally.
col_map <- c(
  mass_g                = "Mass",
  wing_length_mm        = "Wing.Length",
  beak_length_culmen_mm = "Beak.Length_Culmen",
  tarsus_length_mm      = "Tarsus.Length",
  tail_length_mm        = "Tail.Length",
  range_size_km2        = "Range.Size",
  centroid_lat          = "Centroid.Latitude",
  centroid_lon          = "Centroid.Longitude"
)

set.seed(1)
N_SUBSET  <- 2000
MISS_RATE <- 0.10

keep_sp <- sample(rownames(avonet_traits), N_SUBSET)
dat     <- avonet_traits[keep_sp, names(col_map), drop = FALSE]
colnames(dat) <- unname(col_map)
dat     <- data.frame(species = rownames(dat), dat,
                       stringsAsFactors = FALSE,
                       check.names      = FALSE)
tree    <- ape::keep.tip(avonet_tree, keep_sp)

# Hide MISS_RATE per continuous column, build a long-format truth frame.
traits <- unname(col_map)
truth_list <- list()
for (v in traits) {
  vals <- dat[[v]]
  obs  <- which(!is.na(vals))
  miss <- sample(obs, floor(length(obs) * MISS_RATE))
  truth_list[[v]] <- data.frame(
    species_tip = dat$species[miss],
    trait       = v,
    true_value  = vals[miss],
    stringsAsFactors = FALSE
  )
  dat[[v]][miss] <- NA
}
truth <- do.call(rbind, truth_list)

stopifnot(all(dat$species %in% tree$tip.label),
          length(tree$tip.label) == nrow(dat))

## Reorder to tree tips (crucial for BACE CI join in section 5)
dat <- dat %>% slice(match(tree$tip.label, species)) %>% as.data.frame()

## Optional: ensure branch lengths exist / are compatible
tree <- ape::compute.brlen(tree, method = "Grafen")

## ---- 2.5) Optional: log-transform traits BEFORE imputation ----

# Choose which traits to log-transform:
# For AVONET continuous traits, often: all positive traits (e.g., Mass, Wing, Tarsus, etc.)
log_traits <- traits  # or specify: c("Mass","Wing.Length","Tarsus.Length", ...)

# Helper: build a safe log transform (handles zeros/negatives via per-trait shift)
make_log_transform <- function(df, vars, eps = 0) {
  # eps can be 0 for strictly positive variables; use small eps if you want log(x + eps)
  shifts <- setNames(numeric(length(vars)), vars)
  
  for (v in vars) {
    x <- df[[v]]
    x_obs <- x[is.finite(x)]
    if (!length(x_obs)) next
    
    min_x <- min(x_obs, na.rm = TRUE)
    
    # If any observed values <= 0, shift up so min becomes (1 + eps), then log
    # (Using 1 rather than tiny values makes back-transforming more stable.)
    if (min_x <= 0) {
      shifts[v] <- (1 + eps) - min_x
    } else {
      shifts[v] <- eps
    }
  }
  
  forward <- function(df_in) {
    df_out <- df_in
    for (v in vars) df_out[[v]] <- log(df_out[[v]] + shifts[v])
    df_out
  }
  
  inverse <- function(df_in) {
    df_out <- df_in
    for (v in vars) df_out[[v]] <- exp(df_out[[v]]) - shifts[v]
    df_out
  }
  
  list(shifts = shifts, forward = forward, inverse = inverse)
}

log_tf <- make_log_transform(dat, log_traits, eps = 0)

# Apply forward transform to the data used for imputation
dat_imp <- log_tf$forward(dat)

# IMPORTANT: if your truth values are on the original scale (they are),
# keep them as-is for evaluation. We'll back-transform predictions later.

## Use any variable as the LHS; bace_imp will generate chained formulas for all vars in the formula
fixformula <- paste0("Mass ~ ", paste(setdiff(traits, "Mass"), collapse = " + "))

fit_bace <- bace_imp(
  fixformula     = fixformula,
  ran_phylo_form = "~ 1 | species",
  phylo          = tree,
  data           = dat_imp,
  runs           = 5,       
  nitt           = 4000*5,    
  burnin         = 1000*5,
  thin           = 10*2,
  verbose        = TRUE
)

bace_last <- fit_bace$data[[length(fit_bace$data)]]
bace_last <- as.data.frame(bace_last)
#write_csv(bace_last, "bace_last.csv")

## ---- 4) Run Rphylopars imputation (on log scale) ----
trait_df <- dat_imp %>%  # <- changed
  dplyr::select(species, dplyr::everything()) %>%
  as.data.frame()
  
  system.time(
    rp <- Rphylopars::phylopars(
      trait_data  = trait_df,
      tree        = tree,
      model       = "BM",
      pheno_error = FALSE
    ))


## Extract reconstructed tip traits (includes imputations)
## Different versions expose slightly different slot names; handle both.
## Robustly extract reconstructed TIP traits from a phylopars fit
get_tip_recon_phylopars <- function(rp, what = c("mean", "var")) {
  what <- match.arg(what)
  
  if (is.null(rp$tree) || is.null(rp$tree$tip.label)) {
    stop("rp$tree$tip.label not found; cannot identify tips.")
  }
  
  tip_labels <- rp$tree$tip.label
  
  # Choose source matrix
  mat <- switch(
    what,
    mean = rp$anc_recon,
    var  = rp$anc_var
  )
  
  if (is.null(mat)) {
    stop(sprintf("Could not find rp$anc_%s in phylopars output.", what))
  }
  
  # Need rownames to match tips; anc_recon in your object *does* have them
  if (is.null(rownames(mat))) {
    stop("Reconstruction matrix has no rownames; cannot match to tip labels.")
  }
  
  # Subset/reorder to tips
  idx <- match(tip_labels, rownames(mat))
  
  if (anyNA(idx)) {
    missing <- tip_labels[is.na(idx)]
    stop(
      "Some tip labels were not found in reconstruction matrix rownames. ",
      "Examples: ", paste(utils::head(missing, 10), collapse = ", ")
    )
  }
  
  mat[idx, , drop = FALSE]
}

## Usage:
rp_tip_recon <- get_tip_recon_phylopars(rp, what = "mean")
rp_tip_var   <- get_tip_recon_phylopars(rp, what = "var")
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