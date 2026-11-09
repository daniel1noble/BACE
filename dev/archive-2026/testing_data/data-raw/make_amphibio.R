#!/usr/bin/env Rscript
# data-raw/make_amphibio.R
#
# Build amphibio_traits + amphibio_tree from the cached AmphiBIO csv +
# tree.rda at ~/Dropbox/Github Local/pigauto/script/data-cache/.
#
# Trait selection mirrors script/bench_phase_c_cross_dataset.R (Phase C
# v2): 6 mixed-type columns -- 2 log-cont, 1 ordinal K=5, 1 categorical
# K=4, plus binary Diu / Noc included here for users who want the raw
# presence-only encoding (drop them at use-site if you want a clean
# binary signal).

suppressPackageStartupMessages({
  library(ape)
})

cache_dir <- "/Users/z3437171/Dropbox/Github Local/pigauto/script/data-cache"
out_dir   <- "/Users/z3437171/Dropbox/Github Local/phyloTraitData/data"

csv_local  <- file.path(cache_dir, "AmphiBIO_v1.csv")
tree_local <- file.path(cache_dir, "amphibio_tree.rda")
stopifnot(file.exists(csv_local), file.exists(tree_local))

amph <- utils::read.csv(csv_local, stringsAsFactors = FALSE)
e <- new.env(parent = emptyenv())
load(tree_local, envir = e)
tree <- e[[ls(e)[1]]]

amph$species_key <- gsub("[^A-Za-z0-9_]", "_", amph$Species)
tree$tip.label <- gsub("[^A-Za-z0-9_]", "_", tree$tip.label)

overlap <- intersect(tree$tip.label, amph$species_key)
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, overlap))
amph <- amph[amph$species_key %in% overlap, ]
amph <- amph[match(tree$tip.label, amph$species_key), ]
rownames(amph) <- amph$species_key

df <- data.frame(row.names = amph$species_key)
df$body_size_mm <- {
  x <- as.numeric(amph$Body_size_mm); x[x <= 0 | !is.finite(x)] <- NA; x
}
df$body_mass_g  <- {
  x <- as.numeric(amph$Body_mass_g);  x[x <= 0 | !is.finite(x)] <- NA; x
}
# Diu / Noc are presence-only in AmphiBIO (1 = recorded as diurnal /
# nocturnal; NA = no record).  Keep as factor so the user can decide
# whether to treat NA as 0 or as truly missing.
to_bin <- function(x) {
  raw <- as.character(x)
  raw[raw %in% c("1", "TRUE", "T")] <- "yes"
  raw[raw %in% c("0", "FALSE", "F")] <- "no"
  raw[!raw %in% c("yes", "no")] <- NA
  factor(raw, levels = c("no", "yes"))
}
df$diurnal   <- to_bin(amph$Diu)
df$nocturnal <- to_bin(amph$Noc)

# Diet breadth from 6 indicator columns
diet_cols <- intersect(c("Leaves", "Flowers", "Seeds", "Fruits",
                         "Arthro", "Vert"), names(amph))
diet_mat <- sapply(diet_cols, function(cn) {
  as.integer(!is.na(amph[[cn]]) & amph[[cn]] == 1L)
})
any_diet_data <- rowSums(!is.na(amph[, diet_cols, drop = FALSE])) > 0L
breadth <- rowSums(diet_mat, na.rm = TRUE)
breadth[!any_diet_data] <- NA_integer_
breadth[breadth == 0L]  <- NA_integer_
breadth <- pmin(breadth, 5L)
df$diet_breadth <- ordered(as.integer(breadth), levels = 1:5)

# Habitat -- pick first-marked (Fos / Ter / Aqu / Arb)
hab_cols <- intersect(c("Fos", "Ter", "Aqu", "Arb"), names(amph))
hab_assignment <- apply(amph[, hab_cols, drop = FALSE], 1, function(row) {
  idx <- which(!is.na(row) & row == 1L)
  if (length(idx) == 0L) NA_character_ else hab_cols[idx[1]]
})
df$habitat <- factor(hab_assignment, levels = hab_cols)
# Higher taxonomy (related data, useful for taxonomic random effects)
if ("Order"  %in% names(amph)) df$order_name  <- factor(amph$Order)
if ("Family" %in% names(amph)) df$family_name <- factor(amph$Family)
if ("Genus"  %in% names(amph)) df$genus_name  <- factor(amph$Genus)

non_missing <- rowSums(!is.na(df)) > 0L
df   <- df[non_missing, , drop = FALSE]
tree <- ape::keep.tip(tree, rownames(df))
df   <- df[tree$tip.label, , drop = FALSE]

# MCMCglmm::inverseA needs every tip + node label distinct. Source
# amphibio tree carries collapsed-singleton internals labelled "" that
# all duplicate. make.unique-ing preserves informative Family names
# where present and kills the empty-string duplicates.
if (!is.null(tree$node.label)) {
  nl <- ifelse(nzchar(tree$node.label) & !is.na(tree$node.label),
               tree$node.label, "N")
  tree$node.label <- make.unique(nl)
}

amphibio_traits <- df
amphibio_tree   <- tree

cat(sprintf("amphibio: %d species, %d traits, %d tips\n",
             nrow(amphibio_traits), ncol(amphibio_traits),
             length(amphibio_tree$tip.label)))

save(amphibio_traits, file = file.path(out_dir, "amphibio_traits.rda"),
     compress = "xz")
save(amphibio_tree,   file = file.path(out_dir, "amphibio_tree.rda"),
     compress = "xz")
