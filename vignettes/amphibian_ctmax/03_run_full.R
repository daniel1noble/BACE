## =============================================================================
## Amphibian CTmax example, step 3: impute CTmax for all 5,203 species
## =============================================================================
##
## Runs BACE on the full dataset (18,270 rows). The final imputations are
## independent given the chained phase, so they can be split across processes.
## Run from the BACE repository root:
##
##   Rscript vignettes/amphibian_ctmax/03_run_full.R chain        # chained phase
##   Rscript vignettes/amphibian_ctmax/03_run_full.R final 1 6    # final part 1 of 6
##   ...                                                          # (parts in parallel;
##   Rscript vignettes/amphibian_ctmax/03_run_full.R final 6 6    #  each needs ~5 GB RAM)
##   Rscript vignettes/amphibian_ctmax/03_run_full.R combine      # summarise
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
stopifnot(mode %in% c("chain", "final", "combine"))

suppressMessages({
  pkgload::load_all(".", quiet = TRUE)   # the BACE source in this repository
  library(ape)
})
source("vignettes/amphibian_ctmax/00_setup.R")

ctmax_data <- ctmax_fill_medium(readRDS("vignettes/amphibian_ctmax/data/ctmax_data.rds"))
ctmax_tree <- readRDS("vignettes/amphibian_ctmax/data/ctmax_tree.rds")
tmp_dir <- "vignettes/amphibian_ctmax/results/tmp"   # large intermediates, not tracked
res_dir <- "vignettes/amphibian_ctmax/results"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

n_final <- 50

if (mode == "chain") {
  chain <- bace_chain_step(ctmax_data, ctmax_tree, nitt = 15000, burnin = 3000, thin = 24,
                           runs = 10, seed = 2026, verbose = TRUE)
  saveRDS(chain, file.path(tmp_dir, "full_chain.rds"))
  cat(sprintf("Chained phase done in %.1f min\n", chain$minutes))
}

if (mode == "final") {
  part <- as.integer(args[2]); n_parts <- as.integer(args[3])
  part_sizes <- diff(round(seq(0, n_final, length.out = n_parts + 1)))   # e.g. 8, 9, 8, ...
  chain <- readRDS(file.path(tmp_dir, "full_chain.rds"))
  final <- bace_final_step(chain, ctmax_tree, n_final = part_sizes[part], chunk = 1,
                           seed = 5000 + 100 * part, verbose = TRUE)
  saveRDS(final, file.path(tmp_dir, sprintf("full_final_part%d.rds", part)))
  cat(sprintf("Final part %d done in %.1f min\n", part, final$minutes))
}

if (mode == "combine") {
  chain <- readRDS(file.path(tmp_dir, "full_chain.rds"))
  parts <- lapply(sort(list.files(tmp_dir, "^full_final_part", full.names = TRUE)), readRDS)
  draws <- do.call(cbind, lapply(parts, `[[`, "ctmax_draws"))
  stopifnot(ncol(draws) == n_final)
  std <- ctmax_data$row_type == "standardised"

  # Per-row summaries of the n_final imputations for the standardised rows
  imputed <- data.frame(
    row_id     = ctmax_data$row_id[std],
    species    = ctmax_data$species[std],
    temp_level = ctmax_data$temp_level[std],
    acclimation_temp = ctmax_data$acclimation_temp[std],
    CTmax_mean = rowMeans(draws[std, ]),
    CTmax_sd   = apply(draws[std, ], 1, sd),
    lower95    = apply(draws[std, ], 1, quantile, 0.025),
    upper95    = apply(draws[std, ], 1, quantile, 0.975))

  # Pooled posterior of the CTmax model: stack the per-imputation draws, as
  # pool_posteriors() does. Coefficients are on BACE's internal z-scale; the
  # scale factors let the vignette back-transform them.
  pooled <- list(Sol = do.call(rbind, lapply(parts, function(p) do.call(rbind, p$sol))),
                 VCV = do.call(rbind, lapply(parts, function(p) do.call(rbind, p$vcv))),
                 sd_CTmax = sd(ctmax_data$CTmax, na.rm = TRUE),
                 sd_acclimation_temp = sd(ctmax_data$acclimation_temp))

  # One extra final imputation, run with the same settings. The chunked run
  # kept only the CTmax draws, so this fit provides (a) complete chains of all
  # four imputation models, for trace plots, and (b) the MCMCglmm structure of
  # the pooled model below.
  s <- chain$settings
  set.seed(1)
  extra <- bace_final_imp(chain$fit, fixformula = ctmax_formulas,
                          ran_phylo_form = "~ 1 | species", phylo = ctmax_tree,
                          nitt = s$nitt, burnin = s$burnin, thin = s$thin, n_final = 1,
                          species = TRUE, verbose = FALSE)$all_models[[1]]
  extra_chains <- lapply(extra, function(m)
    list(Sol = m$Sol[, seq_len(m$Fixed$nfl), drop = FALSE], VCV = m$VCV))

  # The pooled posterior as an object of the class get_pooled_model() returns,
  # so that summary() and plot() work as in the rest of the tutorial. As in
  # pool_posteriors(), one fitted model supplies the structure; every posterior
  # draw in the object comes from the n_final final imputations.
  pooled_model <- extra$CTmax
  rm(extra)
  # Keep only what summary() and plot() need. Formulas inside the fit carry the
  # environment of BACE's internal functions (and with it the full data), which
  # would inflate the saved file to hundreds of MB, so reset them.
  pooled_model[c("Liab", "Deviance", "DIC", "X", "Z", "ZR", "XL", "ginverse", "y.additional")] <- NULL
  strip_env <- function(x) {
    if (inherits(x, "formula")) environment(x) <- globalenv()
    else if (is.list(x)) x[] <- lapply(x, strip_env)
    x
  }
  pooled_model <- strip_env(pooled_model)
  pooled_model$Sol <- coda::as.mcmc(pooled$Sol)
  pooled_model$VCV <- coda::as.mcmc(pooled$VCV)
  pooled_model$BACE_pooling <- list(n_imputations = n_final,
                                    n_samples_per_imputation = nrow(pooled$Sol) / n_final,
                                    original_samples_per_imputation = nrow(pooled$Sol) / n_final,
                                    total_samples = nrow(pooled$Sol), variable = "CTmax",
                                    pooled = TRUE, sampled = FALSE)
  class(pooled_model) <- c("bace_pooled_MCMCglmm", "MCMCglmm")
  pooled <- pooled[c("sd_CTmax", "sd_acclimation_temp")]   # draws now live in pooled_model

  # Each imputation's CTmax model is an independent MCMC chain. On the full
  # data every chain mixes slowly, because MCMCglmm re-samples the 15,609
  # missing responses at each iteration, so check that the chains agree
  # (Gelman-Rubin R-hat) and how many effective draws they give together.
  sol_list <- unlist(lapply(parts, `[[`, "sol"), recursive = FALSE)
  vcv_list <- unlist(lapply(parts, `[[`, "vcv"), recursive = FALSE)
  chains <- coda::mcmc.list(lapply(seq_along(sol_list), function(i)
    coda::mcmc(cbind(as.matrix(sol_list[[i]]), as.matrix(vcv_list[[i]])))))
  rhat <- coda::gelman.diag(chains, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  ess_each <- sapply(chains, coda::effectiveSize)     # parameters x chains
  chain_diag <- data.frame(parameter = names(rhat), rhat = unname(rhat),
                           ess_per_chain = apply(ess_each, 1, median),
                           ess_total = rowSums(ess_each), n_draws = coda::niter(chains[[1]]) * length(chains))

  # Burn-in check. If burn-in were too short, every chain would still be drifting
  # away from similar starting values, so its early draws would differ from its
  # late draws in the same direction across chains. For each parameter, test the
  # shift between the first 10% and the second half of the kept draws across the
  # independent chains (in posterior SDs).
  n_it <- coda::niter(chains[[1]])
  shift <- t(sapply(chain_diag$parameter, function(p) {
    x <- sapply(chains, function(m) as.numeric(m[, p]))          # draws x chains
    d <- (colMeans(x[seq_len(n_it %/% 10), , drop = FALSE]) -
            colMeans(x[(n_it %/% 2 + 1):n_it, , drop = FALSE])) / stats::sd(as.vector(x))
    c(mean(d), stats::t.test(d)$p.value)
  }))
  chain_diag$early_shift_sd <- shift[, 1]
  chain_diag$early_shift_p  <- shift[, 2]

  out <- list(imputed = imputed,
              pooled = pooled,
              pooled_model = pooled_model,
              types = chain$fit$types,
              mcmc_diagnostics = chain$fit$diagnostics,   # BACE's summary of the last chained iteration
              extra_imputation = extra_chains,   # complete chains of all four models, for trace plots
              convergence = chain$convergence,
              mcmc_diag = chain$mcmc_diag,
              chain_diag = chain_diag,
              settings = c(chain$settings, n_final = n_final),
              mem_mb_per_imputation = mean(vapply(parts, `[[`, numeric(1), "mem_mb_per_imputation")),
              minutes = c(chained = chain$minutes,
                          final = sum(vapply(parts, `[[`, numeric(1), "minutes"))))
  saveRDS(out, file.path(res_dir, "full_imputation.rds"))
  cat("Saved", file.path(res_dir, "full_imputation.rds"), "\n")
}
