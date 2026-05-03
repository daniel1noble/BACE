#!/usr/bin/env Rscript
# data-raw/make_leptraits.R
#
# Build leptraits_traits + leptraits_tree from LepTraits/consensus/
# consensus.csv at
# ~/Dropbox/Github Local/pigauto/script/data-cache/external/.  Tree
# is built taxonomically (Family / Genus / Species) with Grafen
# rank-based branch lengths because LepTraits does not ship a
# phylogeny.
#
# Trait selection mirrors script/bench_leptraits_covariates.R:
#   wingspan_lower (mm; log)
#   forewing_length_lower (mm; log)
#   flight_duration (months / yr)
#   n_hostplant_families (count; log1p)
# Plus 12 monthly flight indicators (Jan..Dec) as 0/1 covariates.

suppressPackageStartupMessages({
  library(ape)
})

cache_dir <- "/Users/z3437171/Dropbox/Github Local/pigauto/script/data-cache"
out_dir   <- "/Users/z3437171/Dropbox/Github Local/phyloTraitData/data"

csv_path <- file.path(cache_dir, "external", "LepTraits",
                      "consensus", "consensus.csv")
stopifnot(file.exists(csv_path))

lt <- utils::read.csv(csv_path, stringsAsFactors = FALSE)

trait_cols <- c(wingspan_lower        = "WS_L",
                forewing_length_lower = "FW_L",
                flight_duration       = "FlightDuration",
                n_hostplant_families  = "NumberOfHostplantFamilies")
month_cols <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
trait_cols <- trait_cols[trait_cols %in% colnames(lt)]
month_cols <- intersect(month_cols, colnames(lt))
stopifnot(length(trait_cols) >= 2L, length(month_cols) == 12L)

# Tax filter
tax_ok <- !is.na(lt$Family) & !is.na(lt$Genus) & !is.na(lt$Species) &
           nzchar(lt$Family) & nzchar(lt$Genus) & nzchar(lt$Species)
lt <- lt[tax_ok, , drop = FALSE]

lt$species_key <- gsub(" ", "_", trimws(lt$Species))
lt <- lt[!duplicated(lt$species_key), , drop = FALSE]
rownames(lt) <- lt$species_key

# Coerce traits + months to numeric; raw scale (no log here)
df <- data.frame(row.names = lt$species_key)
for (nm in names(trait_cols)) {
  v <- suppressWarnings(as.numeric(lt[[trait_cols[[nm]]]]))
  v[!is.finite(v)] <- NA_real_
  df[[nm]] <- v
}
for (m in month_cols) {
  v <- suppressWarnings(as.numeric(lt[[m]]))
  v[!is.finite(v)] <- NA_real_
  df[[m]] <- as.integer(v > 0)
}
# Higher taxonomy (related data, useful for taxonomic random effects)
df$family_name <- factor(lt$Family)
df$genus_name  <- factor(lt$Genus)

# Drop species with no observed trait
non_missing <- rowSums(!is.na(df[, names(trait_cols), drop = FALSE])) > 0L
df <- df[non_missing, , drop = FALSE]
lt <- lt[rownames(df), , drop = FALSE]

# Build taxonomic tree
tax_df <- lt[, c("Family", "Genus", "species_key"), drop = FALSE]
tax_df[] <- lapply(tax_df, factor)
tree <- as.phylo(~Family/Genus/species_key, data = tax_df, collapse = FALSE)
tree <- ape::collapse.singles(tree)
if (!ape::is.rooted(tree))
  tree <- ape::root.phylo(tree, outgroup = 1L, resolve.root = TRUE)
set.seed(2026L)
tree <- ape::multi2di(tree, random = TRUE)
tree <- ape::compute.brlen(tree, method = "Grafen")
tree$edge.length[tree$edge.length <= 0] <- 1e-8

stopifnot(setequal(tree$tip.label, rownames(df)))
df <- df[tree$tip.label, , drop = FALSE]

leptraits_traits <- df
leptraits_tree   <- tree

cat(sprintf("leptraits: %d species, %d cols, %d tips\n",
             nrow(leptraits_traits), ncol(leptraits_traits),
             length(leptraits_tree$tip.label)))

save(leptraits_traits, file = file.path(out_dir, "leptraits_traits.rda"),
     compress = "xz")
save(leptraits_tree,   file = file.path(out_dir, "leptraits_tree.rda"),
     compress = "xz")
