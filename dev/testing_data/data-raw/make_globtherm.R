#!/usr/bin/env Rscript
# data-raw/make_globtherm.R
#
# Build globtherm_traits + globtherm_tree from globTherm.rda at
# ~/Dropbox/Github Local/pigauto/script/data-cache/external/.  Tree
# is built taxonomically (Class / Order / Family / Genus / Species)
# with Grafen rank-based branch lengths because GlobTherm does not
# ship a phylogeny.

suppressPackageStartupMessages({
  library(ape)
})

cache_dir <- "/Users/z3437171/Dropbox/Github Local/pigauto/script/data-cache"
out_dir   <- "/Users/z3437171/Dropbox/Github Local/phyloTraitData/data"

gt_path <- file.path(cache_dir, "external", "globTherm.rda")
stopifnot(file.exists(gt_path))

e <- new.env(parent = emptyenv())
load(gt_path, envir = e)
gt <- as.data.frame(e[[ls(e)[1]]], stringsAsFactors = FALSE)

# Need a clean species key + complete higher-tax for the tree
ok_tax <- !is.na(gt$Class) & !is.na(gt$Order) & !is.na(gt$Family) &
           !is.na(gt$Genus) & !is.na(gt$scientificNameStd) &
           nzchar(gt$Class) & nzchar(gt$Order) &
           nzchar(gt$Family) & nzchar(gt$Genus) &
           nzchar(gt$scientificNameStd)
gt <- gt[ok_tax, , drop = FALSE]
gt$species_key <- gsub(" ", "_", gt$scientificNameStd)
gt <- gt[!duplicated(gt$species_key), , drop = FALSE]

df <- data.frame(row.names = gt$species_key)
df$Tmax            <- as.numeric(gt$Tmax)
df$Tmin            <- suppressWarnings(as.numeric(gt$tmin))
df$lat_max         <- as.numeric(gt$lat_max)
df$long_max        <- as.numeric(gt$long_max)
df$elevation_max   <- as.numeric(gt$elevation_max)
df$class_name      <- factor(gt$Class)
df$order_name      <- factor(gt$Order)
df$family_name     <- factor(gt$Family)
df$genus_name      <- factor(gt$Genus)

# Drop rows with no thermal trait (Tmax / Tmin both NA)
keep <- !is.na(df$Tmax) | !is.na(df$Tmin)
df <- df[keep, , drop = FALSE]
gt <- gt[keep, , drop = FALSE]

# Build taxonomic tree from gt rows aligned to df
tax_df <- gt[, c("Class", "Order", "Family", "Genus", "species_key"),
             drop = FALSE]
tax_df[] <- lapply(tax_df, factor)
tree <- as.phylo(~Class/Order/Family/Genus/species_key, data = tax_df,
                  collapse = FALSE)
tree <- ape::collapse.singles(tree)
if (!ape::is.rooted(tree)) {
  tree <- ape::root.phylo(tree, outgroup = 1L, resolve.root = TRUE)
}
set.seed(2026L)
tree <- ape::multi2di(tree, random = TRUE)
tree <- ape::compute.brlen(tree, method = "Grafen")
tree$edge.length[tree$edge.length <= 0] <- 1e-8

stopifnot(setequal(tree$tip.label, rownames(df)))
df <- df[tree$tip.label, , drop = FALSE]

globtherm_traits <- df
globtherm_tree   <- tree

cat(sprintf("globtherm: %d species, %d traits/covariates, %d tips\n",
             nrow(globtherm_traits), ncol(globtherm_traits),
             length(globtherm_tree$tip.label)))

save(globtherm_traits, file = file.path(out_dir, "globtherm_traits.rda"),
     compress = "xz")
save(globtherm_tree,   file = file.path(out_dir, "globtherm_tree.rda"),
     compress = "xz")
