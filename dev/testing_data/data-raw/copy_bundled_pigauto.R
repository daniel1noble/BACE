#!/usr/bin/env Rscript
# data-raw/copy_bundled_pigauto.R
#
# Copy AVONET-independent phylogenetic datasets shipped inside pigauto
# (real Delhey bird-plumage data, simulated CTmax multi-obs data, plus
# the 300-tip mammal tree and its posterior sample) into this package
# so they ride along.  The .rda files use the original variable names
# (delhey5809, tree_delhey, ctmax_sim, tree300, trees300) so that
# code written against pigauto continues to work after `library(phyloTraitData)`.

src_dir <- "/Users/z3437171/Dropbox/Github Local/pigauto/data"
out_dir <- "/Users/z3437171/Dropbox/Github Local/phyloTraitData/data"

bundled <- c("delhey5809.rda", "tree_delhey.rda",
             "ctmax_sim.rda",  "tree300.rda",
             "trees300.rda")

for (f in bundled) {
  src <- file.path(src_dir, f)
  dst <- file.path(out_dir, f)
  stopifnot(file.exists(src))
  file.copy(src, dst, overwrite = TRUE)
  cat(sprintf("  copied %s (%.1f KB)\n", f, file.size(dst) / 1024))
}
