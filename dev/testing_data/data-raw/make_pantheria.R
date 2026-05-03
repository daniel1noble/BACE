#!/usr/bin/env Rscript
# data-raw/make_pantheria.R
#
# Build pantheria_traits + pantheria_tree from the cached PanTHERIA
# raw text + tip-resolved mammal tree at
# ~/Dropbox/Github Local/pigauto/script/data-cache/.
#
# Trait selection mirrors script/bench_pantheria_full.R and
# script/bench_phase_c_cross_dataset.R: 8 mixed-type traits (4 log-cont,
# 1 count, 2 ordinal, 1 categorical).

suppressPackageStartupMessages({
  library(ape)
})

cache_dir <- "/Users/z3437171/Dropbox/Github Local/pigauto/script/data-cache"
out_dir   <- "/Users/z3437171/Dropbox/Github Local/phyloTraitData/data"

stopifnot(file.exists(file.path(cache_dir, "pantheria.txt")),
          file.exists(file.path(cache_dir, "mammal_tree.tre")))

pan <- utils::read.table(file.path(cache_dir, "pantheria.txt"),
                          header = TRUE, sep = "\t",
                          na.strings = "-999", quote = "",
                          stringsAsFactors = FALSE, comment.char = "")
pan$species_key <- gsub("[^A-Za-z0-9_]", "",
                        paste(pan$MSW93_Genus, pan$MSW93_Species, sep = "_"))

tree <- ape::read.tree(file.path(cache_dir, "mammal_tree.tre"))
tree$tip.label <- gsub("[^A-Za-z0-9_]", "", tree$tip.label)

overlap <- intersect(tree$tip.label, pan$species_key)
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, overlap))
pan  <- pan[pan$species_key %in% overlap, ]
pan  <- pan[match(tree$tip.label, pan$species_key), ]
rownames(pan) <- pan$species_key

df <- data.frame(row.names = pan$species_key)
df$body_mass_g          <- as.numeric(pan$X5.1_AdultBodyMass_g)
df$head_body_length_mm  <- as.numeric(pan$X13.1_AdultHeadBodyLen_mm)
df$gestation_d          <- as.numeric(pan$X9.1_GestationLen_d)
df$max_longevity_m      <- as.numeric(pan$X17.1_MaxLongevity_m)
df$litter_size          <- as.integer(round(as.numeric(pan$X15.1_LitterSize)))
df$diet_breadth         <- ordered(as.integer(pan$X6.1_DietBreadth),
                                    levels = 1:5)
df$habitat_breadth      <- ordered(as.integer(pan$X12.1_HabitatBreadth),
                                    levels = 1:3)
df$terrestriality       <- factor(as.integer(pan$X12.2_Terrestriality),
                                  levels = 1:2)
# Higher taxonomy (related data, useful for taxonomic random effects)
df$order_name  <- factor(pan$MSW93_Order)
df$family_name <- factor(pan$MSW93_Family)
df$genus_name  <- factor(pan$MSW93_Genus)
for (v in c("body_mass_g", "head_body_length_mm", "gestation_d",
            "max_longevity_m")) {
  x <- df[[v]]
  x[x <= 0 | !is.finite(x)] <- NA
  df[[v]] <- x
}

# Drop species with no observed traits at all
non_missing <- rowSums(!is.na(df)) > 0L
df   <- df[non_missing, , drop = FALSE]
tree <- ape::keep.tip(tree, rownames(df))
df   <- df[tree$tip.label, , drop = FALSE]

# MCMCglmm::inverseA needs the tree rooted and every tip + node label
# distinct. The Bininda-Emonds 2007 supertree ships unrooted in this
# cache, and most internal nodes carry empty-string labels that all
# collide. Root if needed, then make.unique node labels, preserving
# informative names where present.
if (!ape::is.rooted(tree)) {
  tree <- ape::root.phylo(tree, outgroup = 1L, resolve.root = TRUE)
}
if (!is.null(tree$node.label)) {
  nl <- ifelse(nzchar(tree$node.label) & !is.na(tree$node.label),
               tree$node.label, "N")
  tree$node.label <- make.unique(nl)
}

pantheria_traits <- df
pantheria_tree   <- tree

cat(sprintf("pantheria: %d species, %d traits, %d tips\n",
             nrow(pantheria_traits), ncol(pantheria_traits),
             length(pantheria_tree$tip.label)))

save(pantheria_traits, file = file.path(out_dir, "pantheria_traits.rda"),
     compress = "xz")
save(pantheria_tree,   file = file.path(out_dir, "pantheria_tree.rda"),
     compress = "xz")
