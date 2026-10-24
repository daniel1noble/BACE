## =============================================================================
## Bird migration example, step 3: summarise the BACE runs
## =============================================================================
##
## Combines the chained phase and final imputations of each analysis (from
## 02_run_bace.R, in results/tmp/) into the compact results read by the
## vignette. Run from the BACE repository root:
##
##   Rscript vignettes/bird_migration/03_summarise.R
##
## Outputs (vignettes/bird_migration/results/):
##   full_imputation.rds  imputed behaviour of the 17 species without records,
##                        pooled posterior and MCMC diagnostics
##   cv_results.rds       five-fold cross-validation: scores of BACE and of
##                        three simple baselines, calibration, diagnostics
##   mixed_results.rds    mixed-variable example: scores for migration and the
##                        hand-wing index, diagnostics
## =============================================================================

suppressMessages({
  library(ape)
  library(coda)
  library(MASS)
})
source("vignettes/bird_migration/00_setup.R")

tmp_dir <- "vignettes/bird_migration/results/tmp"
res_dir <- "vignettes/bird_migration/results"
log_dir <- Sys.getenv("BACE_LOG_DIR", unset = "")   # optional: /usr/bin/time -v logs
categories <- c("sedentary", "partial", "migratory")

migration_data <- readRDS("vignettes/bird_migration/data/migration_data.rds")
migration_tree <- readRDS("vignettes/bird_migration/data/migration_tree.rds")

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
# Ranked probability score of one prediction (lower is better)
rps <- function(probs, truth) {
  observed_cdf <- outer(truth, 1:3, `<=`) * 1
  rowSums((t(apply(probs, 1, cumsum))[, 1:2, drop = FALSE] - observed_cdf[, 1:2, drop = FALSE])^2) / 2
}

read_finals <- function(analysis) {
  files <- sort(list.files(tmp_dir, sprintf("^%s_final[0-9]+\\.rds$", analysis), full.names = TRUE))
  stopifnot(length(files) == n_final[[analysis]])
  lapply(files, readRDS)
}

# MCMC diagnostics across the final imputations. In this example every final
# migration model is fitted to the same data (the response's missing values are
# re-sampled inside MCMCglmm; predictors are complete or fixed at the chained
# phase), so the final fits are independent chains for one posterior and the
# Gelman-Rubin R-hat across them is meaningful.
final_diagnostics <- function(finals) {
  chains <- mcmc.list(lapply(finals, function(f) {
    Vp <- as.numeric(f$vcv[, "species"]); c2 <- as.numeric(f$cp[, 1]); s <- sqrt(Vp + 1)
    mcmc(cbind(as.matrix(f$sol), phylo_var = Vp, cutpoint2 = c2,
               h2_latent = Vp / (Vp + 1), std_cutpoint2 = c2 / s,
               std_log_mass = f$sol[, "log_mass"] / s, std_hwi = f$sol[, "hwi"] / s))
  }))
  ess_each <- sapply(chains, effectiveSize)
  rhat <- gelman.diag(chains, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  # Burn-in check (as in the amphibian example): shift between the first 10% and
  # the second half of each chain, in posterior SDs, tested across chains
  n_it <- niter(chains[[1]])
  shift <- t(sapply(colnames(chains[[1]]), function(p) {
    x <- sapply(chains, function(m) as.numeric(m[, p]))
    d <- (colMeans(x[seq_len(n_it %/% 10), , drop = FALSE]) - colMeans(x[(n_it %/% 2 + 1):n_it, , drop = FALSE])) /
      sd(as.vector(x))
    c(mean(d), t.test(d)$p.value)
  }))
  data.frame(parameter = names(rhat), rhat = unname(rhat),
             ess_per_chain = apply(ess_each, 1, median), ess_min_chain = apply(ess_each, 1, min),
             ess_total = rowSums(ess_each), n_draws = n_it * length(chains),
             early_shift_sd = shift[, 1], early_shift_p = shift[, 2], row.names = NULL)
}

pooled_draws <- function(finals) {
  list(sol = do.call(rbind, lapply(finals, function(f) as.matrix(f$sol))),
       vcv = do.call(rbind, lapply(finals, function(f) as.matrix(f$vcv))),
       cp  = do.call(rbind, lapply(finals, function(f) as.matrix(f$cp))))
}

# Run time and peak memory per process from /usr/bin/time -v logs, if available
resource_use <- function(analysis) {
  if (!nzchar(log_dir)) return(NULL)
  read_one <- function(file) {
    x <- readLines(file)
    wall <- sub(".*: ", "", grep("Elapsed \\(wall clock\\)", x, value = TRUE))
    parts <- as.numeric(rev(strsplit(wall, ":")[[1]]))
    c(minutes = sum(parts * c(1 / 60, 1, 60)[seq_along(parts)]),
      max_rss_gb = as.numeric(sub(".*: ", "", grep("Maximum resident set size", x, value = TRUE))) / 2^20)
  }
  chain  <- read_one(file.path(log_dir, sprintf("%s_chain.log", analysis)))
  finals <- t(sapply(sort(list.files(log_dir, sprintf("^%s_final", analysis), full.names = TRUE)), read_one))
  list(chain = chain,
       final = c(minutes_mean = mean(finals[, "minutes"]), core_hours = sum(finals[, "minutes"]) / 60,
                 max_rss_gb = max(finals[, "max_rss_gb"])))
}

# Baselines for the masked species, fitted to the species with known behaviour
baseline_probs <- function(d, scored_rows) {
  train <- which(!is.na(d$migration))
  # (a) observed category frequencies
  freq <- as.numeric(table(d$migration[train])) / length(train)
  marginal <- matrix(freq, length(scored_rows), 3, byrow = TRUE)
  # (b) ordinal probit regression on the same predictors, without phylogeny
  polr_fit <- polr(migration ~ log_mass + hwi, data = d[train, ], method = "probit")
  polr_p   <- predict(polr_fit, newdata = d[scored_rows, ], type = "probs")
  # (c) category frequencies among the five closest relatives with known
  #     behaviour (patristic distance), plus 0.5 per category
  dist <- cophenetic(migration_tree)[d$species[scored_rows], d$species[train]]
  nearest <- t(apply(dist, 1, function(x) {
    counts <- table(factor(d$migration[train][order(x)[1:5]], levels = categories)) + 0.5
    as.numeric(counts) / sum(counts)
  }))
  list(marginal = marginal, polr = unname(polr_p), nearest_relatives = nearest)
}

score_table <- function(prob_list, truth) {
  do.call(rbind, lapply(names(prob_list), function(m) {
    p <- prob_list[[m]]; modal <- max.col(p, ties.method = "first")
    # proportion of species of each observed category whose most likely category was correct
    recall <- sapply(1:3, function(k) mean(modal[truth == k] == k))
    data.frame(method = m, rps = mean(rps(p, truth)), accuracy = mean(modal == truth),
               recall_sedentary = recall[1], recall_partial = recall[2], recall_migratory = recall[3],
               balanced_accuracy = mean(recall),
               mean_abs_category_error = mean(abs(modal - truth)),
               mean_max_probability = mean(apply(p, 1, max)))   # how confident the predictions are
  }))
}

# Paired bootstrap of a mean difference, resampling whole taxonomic families.
# Species within a family are not independent (shared ancestry), so resampling
# species would understate the uncertainty; resampling families accounts for
# dependence within families (though not for deeper phylogenetic dependence).
family_bootstrap <- function(differences, families, n_boot = 2000) {
  by_family <- split(differences, families)
  boot <- replicate(n_boot, {
    picked <- sample(names(by_family), replace = TRUE)
    mean(unlist(by_family[picked], use.names = FALSE))
  })
  c(lower95 = unname(quantile(boot, 0.025)), upper95 = unname(quantile(boot, 0.975)))
}

## -----------------------------------------------------------------------------
## 1. Full imputation: the 17 species without records
## -----------------------------------------------------------------------------
a      <- make_analysis_data("full")
chain  <- readRDS(file.path(tmp_dir, "full_chain.rds"))
finals <- read_finals("full")
probs  <- Reduce(`+`, lapply(finals, `[[`, "probs")) / length(finals)
imputed_counts <- t(apply(sapply(finals, function(f) as.character(f$imputed_migration)), 1,
                          function(x) table(factor(x, levels = categories))))
colnames(imputed_counts) <- paste0("n_", categories)   # imputed datasets per category
full <- list(
  imputed = data.frame(species = a$data$species[a$migration_rows],
                       taxon_family = migration_data$taxon_family[a$migration_rows],
                       taxon_order  = migration_data$taxon_order[a$migration_rows],
                       probs, imputed_counts, check.names = FALSE, row.names = NULL),
  pooled = pooled_draws(finals),
  final_diagnostics = final_diagnostics(finals),
  chain_convergence = chain$convergence,
  bace_diagnostics = chain$fit$diagnostics,
  settings = list(chain = chain_settings, final = final_settings, n_final = n_final[["full"]]),
  minutes = c(chained = chain$minutes, final_mean = mean(sapply(finals, `[[`, "minutes"))),
  fit_mb = mean(sapply(finals, `[[`, "fit_mb")),   # one fitted migration model (MB)
  resources = resource_use("full"))
saveRDS(full, file.path(res_dir, "full_imputation.rds"))

## -----------------------------------------------------------------------------
## 2. Cross-validation
## -----------------------------------------------------------------------------
cv <- lapply(1:5, function(k) {
  analysis <- paste0("cv", k)
  a      <- make_analysis_data(analysis)
  chain  <- readRDS(file.path(tmp_dir, sprintf("%s_chain.rds", analysis)))
  finals <- read_finals(analysis)
  pos    <- match(a$scored_rows, a$migration_rows)          # scored rows within the imputed rows
  truth  <- as.integer(migration_data$migration[a$scored_rows])
  bace_p   <- Reduce(`+`, lapply(finals, function(f) f$probs[pos, ])) / length(finals)
  plugin_p <- Reduce(`+`, lapply(finals, function(f) f$plugin_probs[pos, ])) / length(finals)
  imputed  <- sapply(finals, function(f) as.integer(f$imputed_migration[pos]))   # rows x n_final
  imputed_freq <- t(apply(imputed, 1, function(x) tabulate(x, 3) / length(x)))
  base_p <- baseline_probs(a$data, a$scored_rows)
  list(fold = k, species = a$data$species[a$scored_rows], truth = truth,
       probs = c(list(bace = bace_p, bace_plugin = plugin_p, bace_imputed_frequency = imputed_freq), base_p),
       final_diagnostics = final_diagnostics(finals),
       chain_convergence = chain$convergence,
       minutes = c(chained = chain$minutes, final_mean = mean(sapply(finals, `[[`, "minutes"))),
       resources = resource_use(analysis))
})
truth_all <- unlist(lapply(cv, `[[`, "truth"))
methods   <- names(cv[[1]]$probs)
probs_all <- lapply(setNames(methods, methods), function(m) do.call(rbind, lapply(cv, function(f) f$probs[[m]])))
scores    <- score_table(probs_all, truth_all)
scores_by_fold <- do.call(rbind, lapply(cv, function(f) cbind(fold = f$fold, score_table(f$probs, f$truth))))

# Paired bootstrap of the RPS difference between each baseline and BACE:
# resampling species (descriptive; ignores shared ancestry) and resampling
# families (accounts for dependence within families)
set.seed(99)
rps_bace <- rps(probs_all$bace, truth_all)
families_all <- migration_data$taxon_family[match(unlist(lapply(cv, `[[`, "species")), migration_data$species)]
rps_gain <- do.call(rbind, lapply(c("marginal", "polr", "nearest_relatives"), function(m) {
  rps_m <- rps(probs_all[[m]], truth_all)
  boot  <- replicate(2000, { i <- sample(length(truth_all), replace = TRUE); mean(rps_m[i]) - mean(rps_bace[i]) })
  fam   <- family_bootstrap(rps_m - rps_bace, families_all)
  data.frame(baseline = m, rps_gain = mean(rps_m) - mean(rps_bace),
             lower95_species = quantile(boot, 0.025), upper95_species = quantile(boot, 0.975),
             lower95_family = fam[["lower95"]], upper95_family = fam[["upper95"]], row.names = NULL)
}))

# Confusion matrices and calibration of the BACE probabilities
confusion <- lapply(probs_all[c("bace", "polr", "nearest_relatives")], function(p)
  table(observed = factor(categories[truth_all], levels = categories),
        predicted = factor(categories[max.col(p, ties.method = "first")], levels = categories)))
calibration <- do.call(rbind, lapply(1:3, function(k) {
  p <- probs_all$bace[, k]; bin <- cut(p, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  data.frame(category = categories[k], bin = levels(bin),
             mean_predicted = as.numeric(tapply(p, bin, mean)),
             observed = as.numeric(tapply(truth_all == k, bin, mean)),
             n = as.numeric(table(bin)), row.names = NULL)
}))

cv_results <- list(folds = lapply(cv, function(f) f[c("fold", "species", "truth", "probs", "minutes", "resources")]),
                   scores = scores, scores_by_fold = scores_by_fold, rps_gain = rps_gain,
                   confusion = confusion, calibration = calibration,
                   final_diagnostics = lapply(cv, `[[`, "final_diagnostics"),
                   chain_convergence = lapply(cv, `[[`, "chain_convergence"),
                   settings = list(chain = chain_settings, final = final_settings, n_final = n_final[["cv1"]]))
saveRDS(cv_results, file.path(res_dir, "cv_results.rds"))

## -----------------------------------------------------------------------------
## 3. Mixed-variable imputation
## -----------------------------------------------------------------------------
a      <- make_analysis_data("mixed")
chain  <- readRDS(file.path(tmp_dir, "mixed_chain.rds"))
finals <- read_finals("mixed")
pos    <- match(a$scored_rows, a$migration_rows)
truth  <- as.integer(migration_data$migration[a$scored_rows])
mig_p  <- Reduce(`+`, lapply(finals, function(f) f$probs[pos, ])) / length(finals)
hwi_masked <- a$data$species[a$scored_rows] %in% a$data$species[a$hwi_rows]

# Hand-wing index: imputed values (rows x n_final) against the measured values
hwi_true <- migration_data$hwi[a$hwi_rows]
hwi_imp  <- sapply(finals, `[[`, "imputed_hwi")
train    <- setdiff(seq_len(nrow(a$data)), a$hwi_rows)
lm_fit   <- lm(hwi ~ log_mass, data = migration_data[train, ])
dist     <- cophenetic(migration_tree)[a$data$species[a$hwi_rows], a$data$species[train]]
hwi_pred <- list(bace = rowMeans(hwi_imp),
                 observed_mean = rep(mean(migration_data$hwi[train]), length(a$hwi_rows)),
                 lm_log_mass = unname(predict(lm_fit, newdata = migration_data[a$hwi_rows, ])),
                 nearest_relatives = apply(dist, 1, function(x) mean(migration_data$hwi[train][order(x)[1:5]])))
hwi_scores <- data.frame(method = names(hwi_pred),
                         rmse = sapply(hwi_pred, function(p) sqrt(mean((p - hwi_true)^2))),
                         correlation = sapply(hwi_pred, function(p) suppressWarnings(cor(p, hwi_true))),
                         row.names = NULL)
hwi_scores$coverage90 <- c(mean(hwi_true >= apply(hwi_imp, 1, quantile, 0.05) &
                                  hwi_true <= apply(hwi_imp, 1, quantile, 0.95)), NA, NA, NA)

# Compare with the ordinal-only analysis of the same masked species (fold 1)
cv1_p <- cv[[1]]$probs$bace
mixed_migration_scores <- rbind(
  data.frame(analysis = "mixed", species_set = c("all", "hwi_also_missing", "hwi_observed"),
             n = c(length(truth), sum(hwi_masked), sum(!hwi_masked)),
             rps = c(mean(rps(mig_p, truth)), mean(rps(mig_p[hwi_masked, ], truth[hwi_masked])),
                     mean(rps(mig_p[!hwi_masked, ], truth[!hwi_masked])))),
  data.frame(analysis = "ordinal_only (cv fold 1)", species_set = c("all", "hwi_also_missing", "hwi_observed"),
             n = c(length(truth), sum(hwi_masked), sum(!hwi_masked)),
             rps = c(mean(rps(cv1_p, truth)), mean(rps(cv1_p[hwi_masked, ], truth[hwi_masked])),
                     mean(rps(cv1_p[!hwi_masked, ], truth[!hwi_masked])))))

# Paired bootstrap of the RPS difference, mixed minus ordinal-only, resampling
# families (see family_bootstrap above)
set.seed(98)
rps_mixed <- rps(mig_p, truth); rps_ordinal <- rps(cv1_p, truth)
fold1_families <- migration_data$taxon_family[a$scored_rows]
mixed_vs_ordinal <- do.call(rbind, lapply(list(all = rep(TRUE, length(truth)), hwi_also_missing = hwi_masked,
                                               hwi_observed = !hwi_masked), function(keep) {
  d <- rps_mixed[keep] - rps_ordinal[keep]
  fam <- family_bootstrap(d, fold1_families[keep])
  data.frame(n = sum(keep), rps_difference = mean(d), lower95 = fam[["lower95"]], upper95 = fam[["upper95"]])
}))
mixed_vs_ordinal$species_set <- c("all", "hwi_also_missing", "hwi_observed")

# The hand-wing index that the migration models conditioned on, per final
# imputation: zero spread across imputations means that uncertainty in the
# imputed index was not passed on to the migration imputations
hwi_used_sd <- apply(sapply(finals, `[[`, "hwi_used_by_migration"), 1, sd)

mixed_results <- list(
  migration_scores = mixed_migration_scores, mixed_vs_ordinal = mixed_vs_ordinal,
  migration_probs = mig_p, migration_truth = truth, hwi_masked = hwi_masked,
  hwi_scores = hwi_scores, hwi_true = hwi_true, hwi_imputed_mean = rowMeans(hwi_imp),
  hwi_imputed_sd = apply(hwi_imp, 1, sd),
  hwi_used_by_migration_max_sd = max(hwi_used_sd),
  final_diagnostics = final_diagnostics(finals),
  hwi_model_draws = list(sol = do.call(rbind, lapply(finals, function(f) as.matrix(f$hwi_sol))),
                         vcv = do.call(rbind, lapply(finals, function(f) as.matrix(f$hwi_vcv)))),
  chain_convergence = chain$convergence,
  settings = list(chain = chain_settings, final = final_settings, n_final = n_final[["mixed"]]),
  minutes = c(chained = chain$minutes, final_mean = mean(sapply(finals, `[[`, "minutes"))),
  fit_mb = mean(sapply(finals, `[[`, "fit_mb")),
  resources = resource_use("mixed"))
saveRDS(mixed_results, file.path(res_dir, "mixed_results.rds"))

# Provenance of the runs (package version, R and MCMCglmm versions, machine)
writeLines(c(sprintf("BACE source: %s (git %s)", as.character(utils::packageVersion("BACE")),
                     tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) "unknown")),
             sprintf("MCMCglmm: %s", as.character(utils::packageVersion("MCMCglmm"))),
             sprintf("Summarised on: %s, %s", Sys.info()[["nodename"]], format(Sys.time(), "%Y-%m-%d")),
             "", capture.output(sessionInfo())),
           file.path(res_dir, "run_info.txt"))

print(scores); print(rps_gain); print(mixed_migration_scores); print(mixed_vs_ordinal); print(hwi_scores)
cat("Saved", paste(list.files(res_dir, "\\.rds$"), collapse = ", "), "to", res_dir, "\n")
