## =============================================================================
## Amphibian CTmax example, step 2: cross-validation
## =============================================================================
##
## Re-runs the cross-validation of Pottier et al. (2025) with BACE, using the
## same five sets of held-out estimates. Each set masks the experimental CTmax
## of 16 species tested as adults at 1 C/min with CTmax recorded at the onset of
## spasms; BACE then predicts the masked values, which we score against the
## experimental estimates.
##
## Run from the BACE repository root, one fold per process (folds are
## independent, so they can run in parallel):
##
##   Rscript vignettes/amphibian_ctmax/02_run_cv.R 1
##   ...
##   Rscript vignettes/amphibian_ctmax/02_run_cv.R 5
##
## All 2,661 experimental rows enter each run; only the standardised rows,
## which have no CTmax, are left out. Pottier et al. (2025) also dropped 234
## species without experimental data per fold, to keep the proportion of
## missing data constant; rows without CTmax do not inform the CTmax model, so
## that step is not needed here. Dropping species without data from the tree is
## exact under Brownian motion: pruning unobserved tips leaves the covariances
## among the remaining species unchanged.
## =============================================================================

fold_id <- as.integer(commandArgs(trailingOnly = TRUE)[1])
stopifnot(fold_id %in% 1:5)

suppressMessages({
  pkgload::load_all(".", quiet = TRUE)   # the BACE source in this repository
  library(ape)
})
source("vignettes/amphibian_ctmax/00_setup.R")

ctmax_data <- ctmax_fill_medium(readRDS("vignettes/amphibian_ctmax/data/ctmax_data.rds"))
ctmax_tree <- readRDS("vignettes/amphibian_ctmax/data/ctmax_tree.rds")
fold       <- readRDS("vignettes/amphibian_ctmax/data/cv_folds.rds")[[fold_id]]

d <- ctmax_data[ctmax_data$row_type == "experimental", ]
masked <- d$row_n %in% fold$rows_masked
truth  <- d$CTmax[masked]
d$CTmax[masked] <- NA
tree <- keep.tip(ctmax_tree, unique(d$species))

cat(sprintf("Fold %d: %d rows, %d species, %d masked rows from %d species\n",
            fold_id, nrow(d), Ntip(tree), sum(masked), length(unique(d$species[masked]))))

chain <- bace_chain_step(d, tree, nitt = 15000, burnin = 3000, thin = 24,
                         runs = 10, seed = 2026 + fold_id)
final <- bace_final_step(chain, tree, n_final = 50, chunk = 10,
                         seed = 3026 + 100 * fold_id, verbose = TRUE)

out <- list(fold        = fold_id,
            row_n       = d$row_n[masked],
            species     = d$species[masked],
            truth       = truth,
            draws       = final$ctmax_draws[masked, , drop = FALSE],
            chain       = chain$ctmax_chain[masked, , drop = FALSE],
            convergence = chain$convergence,
            mcmc_diag   = chain$mcmc_diag,
            settings    = c(chain$settings, n_final = 50),
            minutes     = c(chained = chain$minutes, final = final$minutes))

dir.create("vignettes/amphibian_ctmax/results", showWarnings = FALSE)
saveRDS(out, sprintf("vignettes/amphibian_ctmax/results/cv_fold%d.rds", fold_id))
cat(sprintf("Fold %d done in %.1f min\n", fold_id, sum(out$minutes)))
