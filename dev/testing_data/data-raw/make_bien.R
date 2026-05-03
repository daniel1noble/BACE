#!/usr/bin/env Rscript
# data-raw/make_bien.R
#
# Build bien_traits + bien_tree from cached BIEN species-level trait
# means + V.PhyloMaker2 tree at
# ~/Dropbox/Github Local/pigauto/script/data-cache/.
#
# bien_trait_means.rds is a list of one-trait data.frames (long format,
# one row per species per trait); we widen it to a single per-species
# data.frame keyed by row.names and intersect with the supplied tree.
# All five traits are continuous and conventionally log-transformed
# downstream (height_m, leaf_area, sla, seed_mass, wood_density).

suppressPackageStartupMessages({
  library(ape)
})

cache_dir <- "/Users/z3437171/Dropbox/Github Local/pigauto/script/data-cache"
out_dir   <- "/Users/z3437171/Dropbox/Github Local/phyloTraitData/data"

means_path <- file.path(cache_dir, "bien_trait_means.rds")
tree_path  <- file.path(cache_dir, "bien_tree.rds")
stopifnot(file.exists(means_path), file.exists(tree_path))

trait_means <- readRDS(means_path)
tree_obj    <- readRDS(tree_path)

# V.PhyloMaker2 returns a list (scenario.1/2/3) or a phylo.  Take the
# scenario.3 phylo if it is a list, otherwise the object itself.
if (inherits(tree_obj, "multiPhylo") || is.list(tree_obj) &&
    !inherits(tree_obj, "phylo")) {
  if (!is.null(tree_obj$scenario.3)) {
    tree <- tree_obj$scenario.3
  } else {
    tree <- tree_obj[[length(tree_obj)]]
  }
} else {
  tree <- tree_obj
}
stopifnot(inherits(tree, "phylo"))

# Wide format: union of all species across the 5 trait dfs
species_pool <- unique(unlist(lapply(trait_means, function(d) d$species)))
species_pool <- species_pool[!is.na(species_pool) & nzchar(species_pool)]

df <- data.frame(row.names = species_pool)
for (tn in names(trait_means)) {
  d <- trait_means[[tn]]
  v <- rep(NA_real_, length(species_pool))
  v[match(d$species, species_pool)] <- d$mean_value
  df[[tn]] <- v
}

# Tree tips use underscore separators: "Genus_species".
df_keys   <- gsub(" ", "_", rownames(df))
tip_keys  <- gsub(" ", "_", tree$tip.label)
overlap   <- intersect(df_keys, tip_keys)

if (length(overlap) < 200) {
  stop("BIEN/tree intersection too small: ", length(overlap),
       " species. Check naming convention.", call. = FALSE)
}

df <- df[match(overlap, df_keys), , drop = FALSE]
rownames(df) <- overlap

tree$tip.label <- tip_keys
tree <- ape::keep.tip(tree, overlap)
df   <- df[tree$tip.label, , drop = FALSE]

# Drop species with all-NA traits
non_missing <- rowSums(!is.na(df)) > 0L
df   <- df[non_missing, , drop = FALSE]
tree <- ape::keep.tip(tree, rownames(df))
df   <- df[tree$tip.label, , drop = FALSE]

bien_traits <- df
bien_tree   <- tree

cat(sprintf("bien: %d species, %d traits, %d tips\n",
             nrow(bien_traits), ncol(bien_traits),
             length(bien_tree$tip.label)))

save(bien_traits, file = file.path(out_dir, "bien_traits.rda"),
     compress = "xz")
save(bien_tree,   file = file.path(out_dir, "bien_tree.rda"),
     compress = "xz")
