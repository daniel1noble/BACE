## =============================================================================
## Bird migration example, step 1: prepare the data
## =============================================================================
##
## Builds the compact inputs used by the bird migration worked example in the
## vignette (vignettes/bird_migration/_bird_migration.qmd). Run once, from the
## BACE repository root:
##
##   Rscript vignettes/bird_migration/01_prepare_data.R
##
## Sources
##   Traits: AVONET (Tobias et al. 2022, Ecology Letters, doi:10.1111/ele.13898),
##     eBird-taxonomy table (TraitData/AVONET2_eBird.xlsx in ELEData); figshare doi:10.6084/m9.figshare.16586228.v7,
##     licence CC BY 4.0. Read from a local copy of the ELEData folder when the
##     environment variable AVONET_ELEDATA points to it, otherwise downloaded
##     (ELEData.zip, figshare file 38429873) and checked against its MD5.
##   Phylogeny: Clements/eBird-taxonomy tree of the clootl package (version 1.6,
##     2025 taxonomy; McTavish et al. 2025, PNAS, doi:10.1073/pnas.2409658122),
##     bundled with clootl (GPL-3); the underlying tree data (AvesData,
##     https://github.com/McTavishLab/AvesData) are licensed CC BY 4.0.
##
## Licence and attribution of the files written here: they are derived from the
## two CC BY 4.0 sources above and are redistributed under CC BY 4.0
## (https://creativecommons.org/licenses/by/4.0/). Changes made: species names
## matched to the 2025 taxonomy, species dropped (lumped or unmatched), traits
## subset and transformed (log body mass), tree pruned to the species and made
## ultrametric. Cite Tobias et al. (2022) and McTavish et al. (2025) when reusing.
##
## Outputs (vignettes/bird_migration/data/):
##   migration_data.rds  one row per species: migratory behaviour (ordered),
##                       log body mass, hand-wing index, taxonomic family and order
##   migration_tree.rds  clootl tree pruned to these species, made ultrametric
##   cv_folds.rds        five cross-validation folds of species with known
##                       migratory behaviour
##   hwi_mask.rds        species whose hand-wing index is masked in the
##                       mixed-variable example
## =============================================================================

suppressMessages({
  library(readxl)
  library(ape)
  library(phytools)
  library(clootl)   # attach: extractTree() looks up its bundled tree data
})

out_dir <- "vignettes/bird_migration/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------------------------------------------------------
## 1. AVONET traits (eBird taxonomy)
## -----------------------------------------------------------------------------
avonet_dir <- Sys.getenv("AVONET_ELEDATA", unset = "")
if (!nzchar(avonet_dir)) {
  zip_file <- tempfile(fileext = ".zip")
  utils::download.file("https://ndownloader.figshare.com/files/38429873", zip_file, mode = "wb", quiet = TRUE)
  stopifnot(unname(tools::md5sum(zip_file)) == "ab2db72c48e0200a697c88d95e362380")
  avonet_dir <- file.path(tempdir(), "ELEData")
  utils::unzip(zip_file, exdir = tempdir())
}
# The published table is an Excel sheet; read it and give the columns R-safe names
avonet <- as.data.frame(readxl::read_excel(file.path(avonet_dir, "TraitData", "AVONET2_eBird.xlsx"),
                                           sheet = "AVONET2_eBird", na = c("", "NA")))
names(avonet) <- make.names(names(avonet))      # e.g. "Hand-Wing.Index" -> "Hand.Wing.Index"
avonet$Migration <- as.integer(avonet$Migration)
stopifnot(nrow(avonet) == 10661)

## -----------------------------------------------------------------------------
## 2. Match AVONET names (eBird/Clements 2021) to the 2025 taxonomy of the tree
## -----------------------------------------------------------------------------
# Most names are unchanged. For the rest, eBird species codes link the 2021 and
# 2025 taxonomies (e.g. Accipiter badius -> Tachyspiza badia).
tax_2021 <- clootl_data$taxonomies$year2021
tax_2025 <- clootl_data$taxonomies$year2025
tax_2025 <- tax_2025[tax_2025$CATEGORY == "species", ]

species_code   <- tax_2021$SPECIES_CODE[match(avonet$Species2, tax_2021$SCI_NAME)]
name_by_code   <- tax_2025$SCI_NAME[match(species_code, tax_2025$SPECIES_CODE)]
avonet$name_2025 <- ifelse(avonet$Species2 %in% tax_2025$SCI_NAME, avonet$Species2, name_by_code)

# Species lumped since 2021 map to the same 2025 name; their traits are
# ambiguous, so we drop them, together with names that cannot be matched.
lumped <- avonet$name_2025 %in% avonet$name_2025[duplicated(avonet$name_2025) & !is.na(avonet$name_2025)]
cat(sprintf("AVONET species: %d; matched directly: %d; via species code: %d; dropped (lumped): %d; unmatched: %d\n",
            nrow(avonet), sum(avonet$Species2 %in% tax_2025$SCI_NAME),
            sum(!(avonet$Species2 %in% tax_2025$SCI_NAME) & !is.na(name_by_code)),
            sum(lumped), sum(is.na(avonet$name_2025))))
avonet <- avonet[!is.na(avonet$name_2025) & !lumped, ]

## -----------------------------------------------------------------------------
## 3. Phylogeny: clootl v1.6, 2025 taxonomy
## -----------------------------------------------------------------------------
migration_tree <- suppressMessages(extractTree(species = avonet$name_2025, label_type = "scientific",
                                               taxonomy_year = 2025, version = "1.6", force = TRUE))
migration_tree$node.label <- NULL
migration_tree$tip.label  <- gsub(" ", "_", migration_tree$tip.label, fixed = TRUE)
# The tree is ultrametric up to rounding (~1e-6); MCMCglmm::inverseA() needs it exact.
migration_tree <- suppressMessages(force.ultrametric(migration_tree, method = "extend"))

## -----------------------------------------------------------------------------
## 4. Analysis data: one row per species, in tree order
## -----------------------------------------------------------------------------
avonet$species <- gsub(" ", "_", avonet$name_2025, fixed = TRUE)
avonet <- avonet[match(migration_tree$tip.label, avonet$species), ]
stopifnot(identical(avonet$species, migration_tree$tip.label))

# Migration (AVONET, after Tobias & Pigot 2019): 1 = sedentary; 2 = partially
# migratory (a minority of the population migrates long distances, or most of
# it moves short distances, nomadically or altitudinally); 3 = migratory (most
# of the population migrates long distances). The order reflects increasing
# migratory tendency, so we treat it as an ordinal trait.
migration_data <- data.frame(
  species   = avonet$species,
  migration = factor(avonet$Migration, levels = 1:3,
                     labels = c("sedentary", "partial", "migratory"), ordered = TRUE),
  log_mass  = log(avonet$Mass),
  hwi       = avonet$Hand.Wing.Index,
  # "family" would clash with MCMCglmm's family argument, hence the prefix
  taxon_family = avonet$Family2,
  taxon_order  = avonet$Order2,
  # TRUE when AVONET estimated the hand-wing index from relatives
  hwi_inferred = grepl("Hand-Wing Index", avonet$Traits.inferred, fixed = TRUE),
  stringsAsFactors = FALSE)

stopifnot(!anyNA(migration_data$log_mass), !anyNA(migration_data$hwi))
print(table(migration_data$migration, useNA = "ifany"))

## -----------------------------------------------------------------------------
## 5. Cross-validation folds and the hand-wing-index mask
## -----------------------------------------------------------------------------
# Five folds of the species with known migratory behaviour; each fold's
# species are masked once and their behaviour predicted from the rest.
set.seed(2026)
known <- migration_data$species[!is.na(migration_data$migration)]
fold_id  <- sample(rep(1:5, length.out = length(known)))
cv_folds <- lapply(1:5, function(k) known[fold_id == k])
names(cv_folds) <- paste0("fold", 1:5)

# Mixed-variable example: mask the hand-wing index of 20% of species, drawn
# from species whose index was measured (not inferred by AVONET), so that every
# masked value can be compared with a measurement.
set.seed(2027)
measured_hwi <- migration_data$species[!migration_data$hwi_inferred]
hwi_mask <- sample(measured_hwi, round(0.2 * nrow(migration_data)))

saveRDS(migration_data, file.path(out_dir, "migration_data.rds"))
saveRDS(migration_tree, file.path(out_dir, "migration_tree.rds"))
saveRDS(cv_folds,       file.path(out_dir, "cv_folds.rds"))
saveRDS(hwi_mask,       file.path(out_dir, "hwi_mask.rds"))
cat("Wrote", paste(list.files(out_dir), collapse = ", "), "to", out_dir, "\n")
