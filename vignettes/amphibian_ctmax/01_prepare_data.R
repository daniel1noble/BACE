## =============================================================================
## Amphibian CTmax example, step 1: prepare the data
## =============================================================================
##
## Builds the compact inputs used by the amphibian worked example in the
## vignette (vignettes/amphibian_ctmax/_amphibian_ctmax.qmd). Run once, from the
## BACE repository root:
##
##   Rscript vignettes/amphibian_ctmax/01_prepare_data.R
##
## Source: Pottier et al. (2025) Nature, doi:10.1038/s41586-025-08665-0; data and code at
## https://github.com/p-pottier/Vulnerability_amphibians_global_warming (GPL-3).
## Files are read from a local clone when `src_local` exists, otherwise they are
## downloaded from GitHub (Git LFS) at the pinned commit below.
##
## Outputs (vignettes/amphibian_ctmax/data/):
##   ctmax_data.rds      18,270 rows: 2,661 experimental CTmax estimates for 524
##                       species, plus 15,609 rows (5,203 species x 3
##                       acclimation temperatures) at standardised assay
##                       conditions whose CTmax is to be imputed
##   ctmax_tree.rds      phylogeny pruned to the 5,203 species, made ultrametric
##   cv_folds.rds        the five cross-validation folds of the original study
## =============================================================================

suppressMessages({
  library(ape)
  library(phytools)
})

commit    <- "f47fb6b4935ac4265c0302e2a1c7e51d7b49fdf9"
src_url   <- paste0("https://media.githubusercontent.com/media/p-pottier/",
                    "Vulnerability_amphibians_global_warming/", commit, "/")
src_local <- Sys.getenv("AMPHIBIAN_REPO", unset = "")   # optional local clone
out_dir   <- "vignettes/amphibian_ctmax/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_src <- function(path) {
  local_file <- file.path(src_local, path)
  if (nzchar(src_local) && file.exists(local_file)) return(readRDS(local_file))
  tmp <- tempfile(fileext = ".rds")
  utils::download.file(paste0(src_url, path), tmp, mode = "wb", quiet = TRUE)
  readRDS(tmp)
}

## -----------------------------------------------------------------------------
## 1. Data: keep the variables used by the original imputation models
## -----------------------------------------------------------------------------
raw <- read_src("RData/General_data/data_for_imputation_with_temp.rds")
# `row_n` is unique for experimental rows only; the three standardised rows of a
# species share one value. `row_id` gives every row its own key.
stopifnot(nrow(raw) == 18270, !anyDuplicated(raw$row_n[raw$imputed == "no"]))

ctmax_data <- data.frame(
  row_id              = seq_len(nrow(raw)),
  row_n               = as.integer(raw$row_n),
  species             = raw$tip.label,
  row_type            = factor(ifelse(raw$imputed == "yes", "standardised", "experimental"),
                               levels = c("experimental", "standardised")),
  CTmax               = raw$mean_UTL,
  acclimation_temp    = raw$acclimation_temp,
  ln_acclimation_time = log(raw$acclimation_time),
  ramping             = raw$ramping,
  ln_sd_CTmax         = log(raw$sd_UTL),
  medium              = factor(raw$medium_test_temp),
  endpoint            = factor(raw$endpoint,
                               levels = c("LRR", "OS", "LOE", "prodding", "other", "death"),
                               labels = c("LRR", "OS", "LRR", "other", "other", "other")),
  acclimated          = factor(raw$acclimated),
  life_stage          = factor(raw$life_stage_tested,
                               levels = c("adult", "adults", "larvae"),
                               labels = c("adults", "adults", "larvae")),
  ecotype             = factor(raw$ecotype),
  stringsAsFactors = FALSE)

# The three standardised rows per species sit at the 5th, 50th and 95th
# percentiles of the operative body temperatures across its range.
std <- ctmax_data$row_type == "standardised"
ctmax_data$temp_level <- NA_character_
ctmax_data$temp_level[std] <- ave(ctmax_data$acclimation_temp[std],
                                  ctmax_data$species[std],
                                  FUN = function(x) c("p05", "p50", "p95")[rank(x, ties.method = "first")])
ctmax_data$temp_level <- factor(ctmax_data$temp_level, levels = c("p05", "p50", "p95"))

stopifnot(
  sum(!std) == 2661, length(unique(ctmax_data$species[!std])) == 524,
  sum(std) == 15609, all(table(ctmax_data$species[std]) == 3),
  !any(is.infinite(ctmax_data$ln_sd_CTmax)), !any(is.infinite(ctmax_data$ln_acclimation_time)))

## -----------------------------------------------------------------------------
## 2. Phylogeny: prune to the data and remove rounding-level non-ultrametricity,
##    as in the original study
## -----------------------------------------------------------------------------
tree <- read_src("RData/General_data/tree_for_imputation.rds")
ctmax_tree <- keep.tip(tree, unique(ctmax_data$species))
ctmax_tree <- suppressMessages(force.ultrametric(ctmax_tree, method = "extend"))
stopifnot(Ntip(ctmax_tree) == 5203, setequal(ctmax_tree$tip.label, ctmax_data$species))

## -----------------------------------------------------------------------------
## 3. Cross-validation folds of the original study
## -----------------------------------------------------------------------------
# Each fold masked the experimental CTmax of 16 species tested under the
# standardised conditions (adults, onset of spasms, 1 C/min ramping) and
# dropped 234 species with no such data. We store, per fold, the experimental
# rows kept and masked, keyed by their unique row_n.
fold_names <- c("1st", "2nd", "3rd", "4th", "5th")
cv_folds <- lapply(seq_along(fold_names), function(k) {
  fold   <- read_src(paste0("RData/Imputation/data/Data_crossV_", fold_names[k], "_set.rds"))
  masked <- as.integer(fold$row_n[fold$dat_to_validate %in% "yes"])
  stopifnot(all(masked %in% ctmax_data$row_n[!std]))
  list(fold        = k,
       rows_kept   = as.integer(fold$row_n[fold$imputed == "no"]),
       rows_masked = masked)
})
names(cv_folds) <- paste0("fold", 1:5)

# Check the reconstruction against the held-out data reported by Pottier et
# al. (2025): experimental CTmax 36.19 +/- 2.67 (mean +/- s.d.) across the
# folds, with rows masked in two folds counted once.
held_out <- unique(unlist(lapply(cv_folds, `[[`, "rows_masked")))
chk <- ctmax_data[!std & ctmax_data$row_n %in% held_out, ]
cat(sprintf("Held-out rows: n = %d from %d species; experimental CTmax %.2f +/- %.2f
",
            nrow(chk), length(unique(chk$species)), mean(chk$CTmax), sd(chk$CTmax)))
stopifnot(nrow(chk) == 375, round(mean(chk$CTmax), 2) == 36.19, round(sd(chk$CTmax), 2) == 2.67)

saveRDS(ctmax_data, file.path(out_dir, "ctmax_data.rds"))
saveRDS(ctmax_tree, file.path(out_dir, "ctmax_tree.rds"))
saveRDS(cv_folds,   file.path(out_dir, "cv_folds.rds"))
cat("Wrote", paste(list.files(out_dir), collapse = ", "), "to", out_dir, "
")
