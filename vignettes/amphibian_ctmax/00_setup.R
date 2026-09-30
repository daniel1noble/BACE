## =============================================================================
## Amphibian CTmax example: shared settings, sourced by the run scripts
## =============================================================================

# One imputation model per variable with missing data, mirroring the models of
# Pottier et al. (2025) with the random slopes on acclimation temperature
# removed. The standard deviation of CTmax (ln_sd_CTmax) is an auxiliary
# variable: it carries information about assay methods, which helps impute the
# methodological covariates.
ctmax_formulas <- list(
  "ln_sd_CTmax ~ life_stage + acclimation_temp + endpoint + acclimated + ln_acclimation_time + medium + ramping + CTmax",
  "ln_acclimation_time ~ life_stage + ln_sd_CTmax + CTmax",
  "ramping ~ life_stage + ln_sd_CTmax + CTmax",
  "CTmax ~ acclimation_temp + ln_acclimation_time + ramping + medium + endpoint + acclimated + life_stage + ecotype"
)

# The test medium (ambient vs body/water temperature) is unrecorded for 36 of
# 2,661 experimental rows. Only 51 of the 2,625 recorded assays used ambient
# temperature, too few to estimate the phylogenetic variance of a binary
# threshold model: in a pilot, its variance components reached an effective
# sample size of 12-14 of 500 draws under BACE's default prior and under a
# chi-square(1) prior. Both priors predicted body/water for all 36 rows, with
# probability >= 0.97 (0.99 on average), so we set them to body/water and keep
# `medium` as a complete covariate rather than impute it.
ctmax_fill_medium <- function(data) {
  data$medium[is.na(data$medium)] <- "body_water"
  data
}

ctmax_vars <- c("species", "CTmax", "acclimation_temp", "ln_acclimation_time", "ramping",
                "ln_sd_CTmax", "medium", "endpoint", "acclimated", "life_stage", "ecotype")

# The two helpers below run the BACE workflow step by step. Together they are
# statistically the same as bace(..., skip_conv = TRUE): the chained phase
# (bace_imp) runs first, then n_final independent posterior-predictive
# imputations start from its last iteration (bace_final_imp).
#
# The final phase runs in chunks because keeping all n_final x 4 MCMCglmm fits,
# with their random effects and latent variables, needs more memory than a
# desktop has for the full 18,270-row dataset. From each chunk we keep the
# imputed CTmax values and the fixed-effect and variance-component draws of the
# CTmax model, which is what pool_posteriors() would stack for that model.

# Chained phase: point imputations until the imputed values settle.
bace_chain_step <- function(data, tree, nitt, burnin, thin, runs, seed, verbose = FALSE) {
  set.seed(seed)
  t0  <- Sys.time()
  fit <- bace_imp(fixformula = ctmax_formulas, ran_phylo_form = "~ 1 | species",
                  phylo = tree, data = data[, ctmax_vars],
                  runs = runs, nitt = nitt, burnin = burnin, thin = thin,
                  species = TRUE, verbose = verbose)
  conv <- assess_convergence(fit, method = "summary")

  # MCMC diagnostics of the last chained iteration, then drop the heavy models
  mcmc_diag <- do.call(rbind, lapply(names(fit$models_last_run), function(v) {
    m <- fit$models_last_run[[v]]
    data.frame(variable = v,
               min_ess_fixed = min(coda::effectiveSize(m$Sol[, seq_len(m$Fixed$nfl), drop = FALSE])),
               min_ess_vcv = min(coda::effectiveSize(m$VCV)),
               n_samples = nrow(m$VCV))
  }))
  fit$models_last_run <- NULL
  fit$pred_list_last_run <- NULL

  list(fit = fit, convergence = conv, mcmc_diag = mcmc_diag,
       ctmax_chain = sapply(fit$data[-1], function(d) d$CTmax),   # rows x runs
       settings = list(nitt = nitt, burnin = burnin, thin = thin, runs = runs, seed = seed),
       minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

# Final phase: n_final posterior-predictive imputations, in chunks.
bace_final_step <- function(chain, tree, n_final, chunk, seed, verbose = FALSE) {
  s  <- chain$settings
  t0 <- Sys.time()
  draws <- list(); sol <- list(); vcv <- list(); mem_mb <- numeric(0)
  for (k in seq_len(ceiling(n_final / chunk))) {
    set.seed(seed + k)
    n_k <- min(chunk, n_final - (k - 1) * chunk)
    fin <- bace_final_imp(chain$fit, fixformula = ctmax_formulas, ran_phylo_form = "~ 1 | species",
                          phylo = tree, nitt = s$nitt, burnin = s$burnin, thin = s$thin,
                          n_final = n_k, species = TRUE, verbose = FALSE)
    mem_mb[k] <- as.numeric(utils::object.size(fin)) / 2^20 / n_k   # per imputation
    if (verbose) cat(sprintf("chunk %d: %.0f MB per imputation\n", k, mem_mb[k]))
    draws[[k]] <- sapply(fin$all_datasets, function(d) d$CTmax)
    m   <- lapply(fin$all_models, `[[`, "CTmax")
    sol <- c(sol, lapply(m, function(x) x$Sol[, seq_len(x$Fixed$nfl), drop = FALSE]))
    vcv <- c(vcv, lapply(m, function(x) x$VCV))
    rm(fin, m); invisible(gc())
    if (verbose) cat("final imputations done:", min(k * chunk, n_final), "/", n_final, "\n")
  }
  list(ctmax_draws = do.call(cbind, draws),   # rows x n_final, raw scale
       sol = sol, vcv = vcv,                  # CTmax model, one element per imputation
       mem_mb_per_imputation = mean(mem_mb),  # all four fitted models of one imputation
       minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
