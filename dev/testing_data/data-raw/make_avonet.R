#!/usr/bin/env Rscript
# data-raw/make_avonet.R
#
# Build avonet_traits + avonet_tree from the bundled AVONET trait
# database (Tobias et al. 2022, Ecol Lett 25:581) and matching
# Hackett bird phylogeny (Hackett et al. 2008, Science 320:1763)
# at dev/testing_data/AVONET.csv + Hackett_tree.tre.
#
# Trait selection mirrors dev/00_benchmark_AVONET.R but adds two
# ordinals (habitat_density, migration) and one extra categorical
# (habitat) so the dataset spans the same trait-type mix as the
# other bundled benchmarks: 8 continuous, 2 ordinal, 3 categorical,
# plus higher taxonomy. Continuous morphometrics + range size are
# kept on the raw scale here; users typically log-transform them at
# fit time (see 00_benchmark_AVONET.R header for the rationale).

suppressPackageStartupMessages({
  library(ape)
})

cache_dir <- "/Users/noble/Library/CloudStorage/Dropbox/1_Research/1_Manuscripts/In_Preparation/BACE/dev/testing_data"
out_dir   <- "/Users/noble/Library/CloudStorage/Dropbox/1_Research/1_Manuscripts/In_Preparation/BACE/dev/testing_data/data"

csv_path  <- file.path(cache_dir, "AVONET.csv")
tree_path <- file.path(cache_dir, "Hackett_tree.tre")
stopifnot(file.exists(csv_path), file.exists(tree_path))

av   <- utils::read.csv(csv_path, stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8-BOM")
tree <- ape::read.tree(tree_path)

av$species_key <- gsub("[^A-Za-z0-9_]", "_", trimws(av$Species3))
tree$tip.label <- gsub("[^A-Za-z0-9_]", "_", tree$tip.label)

overlap <- intersect(tree$tip.label, av$species_key)
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, overlap))
av   <- av[av$species_key %in% overlap, ]
av   <- av[!duplicated(av$species_key), ]
av   <- av[match(tree$tip.label, av$species_key), ]
rownames(av) <- av$species_key

df <- data.frame(row.names = av$species_key)
df$mass_g                <- as.numeric(av$Mass)
df$wing_length_mm        <- as.numeric(av$Wing.Length)
df$beak_length_culmen_mm <- as.numeric(av$Beak.Length_Culmen)
df$tarsus_length_mm      <- as.numeric(av$Tarsus.Length)
df$tail_length_mm        <- as.numeric(av$Tail.Length)
df$range_size_km2        <- as.numeric(av$Range.Size)
df$centroid_lat          <- as.numeric(av$Centroid.Latitude)
df$centroid_lon          <- as.numeric(av$Centroid.Longitude)
for (v in c("mass_g", "wing_length_mm", "beak_length_culmen_mm",
            "tarsus_length_mm", "tail_length_mm", "range_size_km2")) {
  x <- df[[v]]
  x[x <= 0 | !is.finite(x)] <- NA
  df[[v]] <- x
}

# Ordinals. Habitat.Density: 1 = dense, 2 = semi-open, 3 = open.
# Migration: 1 = sedentary, 2 = partial, 3 = full migrant.
df$habitat_density <- ordered(as.integer(av$Habitat.Density),
                               levels = 1:3)
df$migration       <- ordered(as.integer(av$Migration),
                               levels = 1:3)

# Categoricals. Trim trailing whitespace in Habitat (e.g. "Shrubland ").
df$trophic_level     <- factor(trimws(av$Trophic.Level))
df$primary_lifestyle <- factor(trimws(av$Primary.Lifestyle))
df$habitat           <- factor(trimws(av$Habitat))

# Higher taxonomy (related data, useful for taxonomic random effects)
df$order_name  <- factor(av$Order3)
df$family_name <- factor(av$Family3)
df$genus_name  <- factor(sub("_.*$", "", av$species_key))

# Drop species with no observed traits at all
non_missing <- rowSums(!is.na(df)) > 0L
df   <- df[non_missing, , drop = FALSE]
tree <- ape::keep.tip(tree, rownames(df))
df   <- df[tree$tip.label, , drop = FALSE]

avonet_traits <- df
avonet_tree   <- tree

cat(sprintf("avonet: %d species, %d traits, %d tips\n",
             nrow(avonet_traits), ncol(avonet_traits),
             length(avonet_tree$tip.label)))

save(avonet_traits, file = file.path(out_dir, "avonet_traits.rda"),
     compress = "xz")
save(avonet_tree,   file = file.path(out_dir, "avonet_tree.rda"),
     compress = "xz")
