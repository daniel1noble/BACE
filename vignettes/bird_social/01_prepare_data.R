## =============================================================================
## Bird social behaviour example, step 1: prepare the data
## =============================================================================
##
## Builds the compact inputs used by the bird social behaviour worked example in
## the vignette (vignettes/bird_social/_bird_social.qmd). Run once, from the
## BACE repository root:
##
##   Rscript vignettes/bird_social/01_prepare_data.R
##
## Sources
##   Social behaviour: BIRDBASE v2025.1 (Sekercioglu et al. 2025, Scientific Data,
##     doi:10.1038/s41597-025-05615-3); figshare doi:10.6084/m9.figshare.27051040.v1,
##     licence CC BY 4.0. Read from a local copy when the environment variable
##     BIRDBASE_XLSX points to it, otherwise downloaded (figshare file 55634729)
##     and checked against its MD5.
##   Body mass, hand-wing index, habitat, family and order: AVONET (Tobias et al.
##     2022, Ecology Letters, doi:10.1111/ele.13898), eBird-taxonomy table
##     (TraitData/AVONET2_eBird.xlsx in ELEData); figshare
##     doi:10.6084/m9.figshare.16586228.v7, licence CC BY 4.0. Read from a local
##     copy of the ELEData folder when the environment variable AVONET_ELEDATA
##     points to it, otherwise downloaded (ELEData.zip, figshare file 38429873)
##     and checked against its MD5.
##   Phylogeny: Clements/eBird-taxonomy tree of the clootl package (version 1.6,
##     2025 taxonomy; McTavish et al. 2025, PNAS, doi:10.1073/pnas.2409658122),
##     bundled with clootl (GPL-3); the underlying tree data (AvesData,
##     https://github.com/McTavishLab/AvesData) are licensed CC BY 4.0.
##
## Licence and attribution of the files written here: they are derived from the
## three CC BY 4.0 sources above and are redistributed under CC BY 4.0
## (https://creativecommons.org/licenses/by/4.0/). Changes made: species names
## matched to the 2025 taxonomy, species dropped (lumped or unmatched), traits
## subset, recoded (social behaviour) and transformed (log body mass), tree
## pruned to the species and made ultrametric. Cite Sekercioglu et al. (2025),
## Tobias et al. (2022) and McTavish et al. (2025) when reusing.
##
## Outputs (vignettes/bird_social/data/):
##   social_data.rds  one row per species: social behaviour (ordered; all known
##                    values, including the held-out ones), the six BIRDBASE
##                    social columns, log body mass, hand-wing index, habitat,
##                    taxonomic family and order
##   social_tree.rds  clootl tree pruned to these species, made ultrametric
##   holdout.rds      species whose known social behaviour is hidden and scored
##                    (whole families and random species)
##   hwi_mask.rds     species whose hand-wing index is hidden in the
##                    mixed-variable example
## =============================================================================

suppressMessages({
  library(readxl)
  library(ape)
  library(phytools)
  library(clootl)   # attach: extractTree() looks up its bundled tree data
})

out_dir <- "vignettes/bird_social/data"
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
avonet <- avonet[!is.na(avonet$name_2025) & !lumped, ]

## -----------------------------------------------------------------------------
## 3. Phylogeny: clootl v1.6, 2025 taxonomy
## -----------------------------------------------------------------------------
social_tree <- suppressMessages(extractTree(species = avonet$name_2025, label_type = "scientific",
                                            taxonomy_year = 2025, version = "1.6", force = TRUE))
social_tree$node.label <- NULL
social_tree$tip.label  <- gsub(" ", "_", social_tree$tip.label, fixed = TRUE)
# The tree is ultrametric up to rounding (~1e-6); MCMCglmm::inverseA() needs it exact.
social_tree <- suppressMessages(force.ultrametric(social_tree, method = "extend"))

avonet$species <- gsub(" ", "_", avonet$name_2025, fixed = TRUE)
avonet <- avonet[match(social_tree$tip.label, avonet$species), ]
stopifnot(identical(avonet$species, social_tree$tip.label))

## -----------------------------------------------------------------------------
## 4. BIRDBASE social behaviour, matched to the tips of the tree
## -----------------------------------------------------------------------------
birdbase_file <- Sys.getenv("BIRDBASE_XLSX", unset = "")
if (!nzchar(birdbase_file)) {
  birdbase_file <- tempfile(fileext = ".xlsx")
  utils::download.file("https://ndownloader.figshare.com/files/55634729", birdbase_file, mode = "wb", quiet = TRUE)
}
stopifnot(unname(tools::md5sum(birdbase_file)) == "36d53d59af1e449db1fc4c30062d9d1f")
birdbase <- as.data.frame(read_excel(birdbase_file, sheet = "Data", skip = 1,
                                     .name_repair = "unique", guess_max = 20000))

# Clements v2024b names match most tips. A few genera were renamed in the 2025
# taxonomy of the tree, so for the remaining tips we try the other name columns
# in turn, using a name only where it is unique and its row is not already used.
name_cols <- c("eBird/Clements (V2024b)", "AviList v1 2025", "IOC World Bird List (v15.1)",
               "Latin (BirdLife > IOC > Clements>AviList)", "HBW/BirdLife International (v9.1)")
row_for_tip <- rep(NA_integer_, length(social_tree$tip.label))
for (col in name_cols) {
  tip_names <- gsub(" ", "_", trimws(birdbase[[col]]), fixed = TRUE)
  ok <- !is.na(tip_names) & !(tip_names %in% tip_names[duplicated(tip_names)]) &
        !(seq_len(nrow(birdbase)) %in% row_for_tip)
  todo <- is.na(row_for_tip)
  row_for_tip[todo] <- match(social_tree$tip.label[todo], ifelse(ok, tip_names, NA))
}
stopifnot(!anyDuplicated(na.omit(row_for_tip)))
cat("Tips matched to BIRDBASE:", sum(!is.na(row_for_tip)), "of", length(row_for_tip), "\n")
bb <- birdbase[row_for_tip, ]

# BIRDBASE social columns: Social_1 colonial; Social_2 social (large numbers,
# mixed-species or seasonal flocks); Social_3 pairs and family groups; Social_4
# singly and pairs; Social_5 solitary; Social_6 lekking. A 0 means that the
# species was not described as behaving that way (not that it never does), and
# a species can have several 1s. Species without any 1 have no description.
social <- as.data.frame(lapply(bb[paste0("Social_", 1:6)], as.numeric))
names(social) <- c("colonial", "flocks", "pairs_family", "singly_pairs", "solitary", "lek")
described <- !is.na(row_for_tip) & rowSums(social == 1, na.rm = TRUE) > 0
social[!described, ] <- NA

# Ordinal social behaviour: the largest grouping in which the species is
# described, solitary < singly and pairs < pairs and family groups < flocks or
# colonies. The columns mix breeding and non-breeding contexts (some species
# are described as both solitary and flocking), so using the largest grouping
# is an analytic choice. Colonial breeding aggregates many birds, so it joins
# the top level. Leks are display aggregations of males and are left out:
# species described only as lekking are missing here.
social_max <- ifelse(social$flocks == 1 | social$colonial == 1, 4,
              ifelse(social$pairs_family == 1, 3,
              ifelse(social$singly_pairs == 1, 2,
              ifelse(social$solitary == 1, 1, NA))))

## -----------------------------------------------------------------------------
## 5. Analysis data: one row per species, in tree order
## -----------------------------------------------------------------------------
social_data <- data.frame(
  species    = avonet$species,
  social_max = factor(social_max, levels = 1:4,
                      labels = c("solitary", "singly_pairs", "pairs_family", "flocks_colonies"),
                      ordered = TRUE),
  log_mass   = log(avonet$Mass),
  hwi        = avonet$Hand.Wing.Index,
  # "family" would clash with MCMCglmm's family argument, hence the prefix
  taxon_family = avonet$Family2,
  taxon_order  = avonet$Order2,
  habitat      = avonet$Habitat,
  # TRUE when AVONET estimated the hand-wing index from relatives
  hwi_inferred = grepl("Hand-Wing Index", avonet$Traits.inferred, fixed = TRUE),
  in_birdbase  = !is.na(row_for_tip),
  described    = described,
  social,
  stringsAsFactors = FALSE)

stopifnot(!anyNA(social_data$log_mass), !anyNA(social_data$hwi))
print(table(social_data$social_max, useNA = "ifany"))

## -----------------------------------------------------------------------------
## 6. Held-out species and the hand-wing-index mask
## -----------------------------------------------------------------------------
# Two disjoint sets of species with known social behaviour, whose values are
# hidden from BACE and used only to score the imputations:
#   family block: whole families (drawn at random among families with at most
#                 100 species of known behaviour) until about 10% of the known
#                 species, mimicking the family-wide gaps in BIRDBASE;
#   random:       10% of the remaining species with known behaviour.
known   <- !is.na(social_data$social_max)
n_known <- sum(known)
set.seed(20261027)
fam_n <- table(social_data$taxon_family[known])
candidates <- sample(names(fam_n)[fam_n <= 100])
block_families <- character(0)
for (f in candidates) {
  if (sum(fam_n[block_families]) >= 0.10 * n_known) break
  block_families <- c(block_families, f)
}
in_block  <- known & social_data$taxon_family %in% block_families
rest      <- which(known & !in_block)
in_random <- seq_len(nrow(social_data)) %in% sample(rest, round(0.10 * n_known))
holdout <- data.frame(species = social_data$species[in_block | in_random],
                      set = ifelse(in_block, "family_block", "random")[in_block | in_random],
                      stringsAsFactors = FALSE)
cat("Held out:", sum(in_block), "species in", length(block_families), "whole families and",
    sum(in_random), "random species\n")

# Mixed-variable example: hide the hand-wing index of 20% of species, drawn
# from species whose index was measured (not inferred by AVONET), so that every
# hidden value can be compared with a measurement.
set.seed(2027)
measured_hwi <- social_data$species[!social_data$hwi_inferred]
hwi_mask <- sample(measured_hwi, round(0.2 * nrow(social_data)))

saveRDS(social_data, file.path(out_dir, "social_data.rds"))
saveRDS(social_tree, file.path(out_dir, "social_tree.rds"))
saveRDS(holdout,     file.path(out_dir, "holdout.rds"))
saveRDS(hwi_mask,    file.path(out_dir, "hwi_mask.rds"))
cat("Wrote", paste(list.files(out_dir), collapse = ", "), "to", out_dir, "\n")
