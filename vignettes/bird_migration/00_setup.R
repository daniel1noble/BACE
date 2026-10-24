## =============================================================================
## Bird migration example: shared settings, sourced by the run scripts
## =============================================================================

# Three analyses use the same data (vignettes/bird_migration/data/):
#   full    impute the migratory behaviour of the 17 species for which AVONET
#           has no record (ordinal variable only)
#   cv1-cv5 five-fold cross-validation: additionally mask the behaviour of one
#           fold of species with known behaviour, and score the imputations
#   mixed   mixed-variable imputation: mask fold 1 of migratory behaviour AND the
#           hand-wing index of 20% of species, so that BACE imputes an ordinal
#           and a continuous variable that predict each other
analyses <- c("full", paste0("cv", 1:5), "mixed")

# One imputation model per variable with missing data
migration_formula <- "migration ~ log_mass + hwi"
hwi_formula       <- "hwi ~ log_mass + migration"

# MCMC settings. The chained phase only needs point imputations to settle; the
# final imputations need long chains because the phylogenetic variance of the
# ordinal (threshold) model mixes slowly (pilot: effective sample size of about
# 130 per 500,000 iterations for 10,563 species). For the mixed analysis the
# second element of each setting applies to the hand-wing-index model, which is
# Gaussian and mixes quickly.
chain_settings <- list(nitt = list(50000, 20000), burnin = list(10000, 5000),
                       thin = list(40, 15), runs_ordinal = 5, runs_mixed = 10)
final_settings <- list(nitt = list(500000, 100000), burnin = list(50000, 10000),
                       thin = list(450, 90))
n_final <- c(full = 50, cv1 = 20, cv2 = 20, cv3 = 20, cv4 = 20, cv5 = 20, mixed = 50)

## -----------------------------------------------------------------------------
## Data for one analysis
## -----------------------------------------------------------------------------
# Returns the BACE input (species, migration, log_mass, hwi), the tree, and the
# rows whose migration or hand-wing index is missing in this analysis.
make_analysis_data <- function(analysis) {
  migration_data <- readRDS("vignettes/bird_migration/data/migration_data.rds")
  migration_tree <- readRDS("vignettes/bird_migration/data/migration_tree.rds")
  cv_folds       <- readRDS("vignettes/bird_migration/data/cv_folds.rds")
  hwi_mask       <- readRDS("vignettes/bird_migration/data/hwi_mask.rds")

  d <- migration_data[, c("species", "migration", "log_mass", "hwi")]
  scored_species <- character(0)
  if (analysis %in% paste0("cv", 1:5)) {
    scored_species <- cv_folds[[sub("cv", "fold", analysis)]]
  }
  if (analysis == "mixed") {
    scored_species <- cv_folds$fold1
    d$hwi[d$species %in% hwi_mask] <- NA
  }
  d$migration[d$species %in% scored_species] <- NA

  list(data = d, tree = migration_tree,
       migration_rows = which(is.na(d$migration)),              # imputed by BACE
       scored_rows    = which(d$species %in% scored_species),   # masked known values
       hwi_rows       = which(is.na(d$hwi)),
       formulas = if (analysis == "mixed") list(migration_formula, hwi_formula) else list(migration_formula))
}

## -----------------------------------------------------------------------------
## Category probabilities from one fitted threshold model
## -----------------------------------------------------------------------------
# For each requested row, the probability of each category averaged over the
# posterior draws: P(y = k) = mean over draws of
#   pnorm(c_k - eta) - pnorm(c_(k-1) - eta),  with eta = X b + u_species,
# cutpoints (0, c_2) and the residual variance fixed at 1. (BACE's own
# prob_preds instead plug the posterior means of eta and c_2 into pnorm(),
# which is sharper than the posterior predictive distribution.)
category_probabilities <- function(fit, rows, species) {
  X    <- as.matrix(fit$X)[rows, c("(Intercept)", "log_mass", "hwi"), drop = FALSE]
  beta <- as.matrix(fit$Sol[, c("(Intercept)", "log_mass", "hwi")])
  u    <- as.matrix(fit$Sol[, paste0("species.", species), drop = FALSE])
  eta  <- beta %*% t(X) + u                                   # draws x rows
  c2   <- as.numeric(fit$CP[, 1])
  cbind(sedentary = colMeans(pnorm(-eta)),
        partial   = colMeans(pnorm(c2 - eta) - pnorm(-eta)),
        migratory = colMeans(1 - pnorm(c2 - eta)))
}

## -----------------------------------------------------------------------------
## Chained phase
## -----------------------------------------------------------------------------
bace_chain_step <- function(a, runs, seed) {
  set.seed(seed)
  t0  <- Sys.time()
  fit <- bace_imp(fixformula = a$formulas, ran_phylo_form = "~ 1 | species",
                  phylo = a$tree, data = a$data,
                  runs = runs, nitt = chain_settings$nitt[seq_along(a$formulas)],
                  burnin = chain_settings$burnin[seq_along(a$formulas)],
                  thin = chain_settings$thin[seq_along(a$formulas)],
                  species = FALSE, verbose = TRUE)
  conv <- assess_convergence(fit, method = "summary")
  # Drop the heavy fitted models of the last chained iteration
  fit$models_last_run <- NULL
  fit$pred_list_last_run <- NULL
  list(fit = fit, convergence = conv, runs = runs, seed = seed,
       minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

## -----------------------------------------------------------------------------
## One final (posterior-predictive) imputation, reduced to what we keep
## -----------------------------------------------------------------------------
bace_final_one <- function(a, chain, seed) {
  set.seed(seed)
  t0  <- Sys.time()
  k   <- seq_along(a$formulas)
  fin <- bace_final_imp(chain$fit, fixformula = a$formulas, ran_phylo_form = "~ 1 | species",
                        phylo = a$tree, nitt = final_settings$nitt[k],
                        burnin = final_settings$burnin[k], thin = final_settings$thin[k],
                        n_final = 1, species = FALSE, verbose = FALSE, n_cores = 1)
  completed <- fin$all_datasets[[1]]
  fit_mig   <- fin$all_models[[1]][["migration"]]
  out <- list(
    imputed_migration = completed$migration[a$migration_rows],
    probs = category_probabilities(fit_mig, a$migration_rows, a$data$species[a$migration_rows]),
    plugin_probs = as.matrix(fin$all_prob_preds[[1]][["migration"]]),
    sol = fit_mig$Sol[, seq_len(fit_mig$Fixed$nfl), drop = FALSE],
    vcv = fit_mig$VCV, cp = fit_mig$CP,
    fit_mb = as.numeric(utils::object.size(fit_mig)) / 2^20)
  if ("hwi" %in% names(fin$all_models[[1]])) {
    fit_hwi <- fin$all_models[[1]][["hwi"]]
    out$imputed_hwi <- completed$hwi[a$hwi_rows]
    out$hwi_sol <- fit_hwi$Sol[, seq_len(fit_hwi$Fixed$nfl), drop = FALSE]
    out$hwi_vcv <- fit_hwi$VCV
    # The (scaled) hand-wing index that the migration model of this imputation
    # conditioned on, at the rows where it was imputed
    out$hwi_used_by_migration <- as.matrix(fit_mig$X)[a$hwi_rows, "hwi"]
  }
  out$minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  out
}
