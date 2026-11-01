## =============================================================================
## Bird social behaviour example, step 3: summarise the BACE runs
## =============================================================================
##
## Combines the chained phase and final imputations of each analysis (from
## 02_run_bace.R, in results/tmp/) into the compact results read by the
## vignette. Run from the BACE repository root:
##
##   Rscript vignettes/bird_social/03_summarise.R
##
## Outputs (vignettes/bird_social/results/):
##   data_summary.rds        who is missing in BIRDBASE: missingness by habitat
##                           and order, deviance explained by family versus traits
##   full_results.rds        full analysis (all known values): imputed behaviour of
##                           the species without a description, pooled posterior,
##                           MCMC diagnostics, downstream analysis with Rubin's rules
##   validation_results.rds  validation analysis: scores of BACE and three simple
##                           baselines on the held-out species, MCMC diagnostics
##   mixed_results.rds       mixed-variable validation: scores for social behaviour
##                           and the hand-wing index, convergence diagnostics
## =============================================================================

suppressMessages({
  pkgload::load_all(".", quiet = TRUE)   # assess_convergence(), with_imputations(), pool_mi()
  library(ape)
  library(coda)
  library(MASS)
  library(phylolm)
})
source("vignettes/bird_social/00_setup.R")

tmp_dir <- "vignettes/bird_social/results/tmp"
res_dir <- "vignettes/bird_social/results"
social_data <- readRDS("vignettes/bird_social/data/social_data.rds")
social_tree <- readRDS("vignettes/bird_social/data/social_tree.rds")
holdout     <- readRDS("vignettes/bird_social/data/holdout.rds")
categories  <- levels(social_data$social_max)
K <- length(categories)

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
# Ranked probability score of each prediction, scaled to 0-1 by K - 1 (lower is better)
rps <- function(probs, truth) {
  observed_cdf <- outer(truth, 1:K, `<=`) * 1
  cum_probs <- t(apply(probs, 1, cumsum))
  rowSums((cum_probs[, 1:(K - 1), drop = FALSE] - observed_cdf[, 1:(K - 1), drop = FALSE])^2) / (K - 1)
}

read_finals <- function(analysis) {
  files <- sort(list.files(tmp_dir, sprintf("^%s_final[0-9]+\\.rds$", analysis), full.names = TRUE))
  stopifnot(length(files) == n_final[[analysis]])
  lapply(files, readRDS)
}

score_table <- function(prob_list, truth) {
  do.call(rbind, lapply(names(prob_list), function(m) {
    p <- prob_list[[m]]; modal <- max.col(p, ties.method = "first")
    # proportion of species of each observed category whose most likely category was correct
    recall <- sapply(1:K, function(k) mean(modal[truth == k] == k))
    data.frame(method = m, n = length(truth), rps = mean(rps(p, truth)),
               accuracy = mean(modal == truth),
               recall_solitary = recall[1], recall_singly_pairs = recall[2],
               recall_pairs_family = recall[3], recall_flocks_colonies = recall[4],
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

# Posterior draws of one final imputation as an mcmc object: fixed effects,
# phylogenetic variance, the two estimated cutpoints and the latent H2 of the
# social model, plus the variances of the hand-wing-index model if present
final_draws <- function(f) {
  Vp <- as.numeric(f$vcv[, "species"])
  m <- cbind(as.matrix(f$sol), phylo_var = Vp, cutpoint2 = as.numeric(f$cp[, 1]),
             cutpoint3 = as.numeric(f$cp[, 2]), h2_latent = Vp / (Vp + 1))
  if (!is.null(f$hwi_vcv)) m <- cbind(m, hwi_phylo_var = as.numeric(f$hwi_vcv[, "species"]),
                                      hwi_residual_var = as.numeric(f$hwi_vcv[, "units"]))
  mcmc(m)
}

# Closest relatives (patristic distance) of one species among the training
# species, as weights. In an ultrametric tree, all members of a sister clade
# are at the same distance, so several relatives are often tied at the k-th
# closest distance. Relatives strictly closer than the k-th distance get weight
# 1, and those tied at it share the remaining weight equally, so the weights
# always sum to k ("weighted ties"; the average over random tie-breaking). With
# ties = "all", every tied relative gets weight 1 instead. Distances that differ
# by less than tie_tol (in the units of the tree, million years) count as tied:
# in this tree, differences from floating-point rounding are below 1e-12 and
# real differences above 1e-6.
tie_tol <- 1e-8
relative_weights <- function(x, k, ties = c("weighted", "all")) {
  ties <- match.arg(ties)
  d_k    <- sort(x)[k]
  closer <- x < d_k - tie_tol
  tied   <- abs(x - d_k) <= tie_tol
  w <- numeric(length(x))
  w[closer] <- 1
  w[tied] <- if (ties == "weighted") (k - sum(closer)) / sum(tied) else 1
  w
}
# Category probabilities from the weighted relatives, plus a pseudo-count per
# category (default 0.5 whatever the total weight; the sensitivity analysis
# also scales it with k, 0.5 * k / 5, so that its relative strength is constant)
relatives_probs <- function(dist, train_cat, k, ties = "weighted", pseudo = 0.5) {
  t(apply(dist, 1, function(x) {
    w <- relative_weights(x, k, ties)
    counts <- tapply(w, factor(train_cat, levels = categories), sum)
    counts[is.na(counts)] <- 0
    as.numeric(counts + pseudo) / sum(counts + pseudo)
  }))
}

# MCMC diagnostics across the final imputations (one chain each):
#   ESS per chain; Geweke z per chain (first 10% vs last 50%); split R-hat
#   within each chain (first versus second half); early shift (mean of the
#   first 10% minus mean of the second half, in posterior SDs, averaged over
#   chains: the size of any drift behind Geweke flags); a placebo Geweke test
#   with the same window sizes as the real one but away from the start of the
#   chain (draws 401-500 vs 501-1000, against 1-100 vs 501-1000): its flag
#   rate shows how often the test flags these chains when no burn-in can be
#   involved; the number of flagged chains with positive z (drift direction); the change in the pooled mean when
#   the first 10% of every chain is dropped; and R-hat across chains.
#   R-hat compares B, the variance of the chain means, with W, the mean
#   within-chain variance (R-hat is about sqrt(1 + B/W)).
#   In bace_final_imp() the missing values of the response are set back to NA
#   in every final model and imputed by MCMCglmm, and the predictors are kept
#   at their imputed values from the chained phase or the previous sweep. When only social behaviour is missing, every final model is fitted
#   to the same data, so the chains target one posterior and the R-hat across
#   them is a convergence diagnostic. In the mixed analysis each final model
#   conditions on different imputed hand-wing indices (final-phase sweeps), so
#   the R-hat across them mixes imputation variability with convergence.
final_diagnostics <- function(finals) {
  chains <- mcmc.list(lapply(finals, final_draws))
  n_it <- niter(chains[[1]])
  ess_each <- sapply(chains, effectiveSize)
  geweke_z <- sapply(chains, function(m) geweke.diag(m)$z)
  split_rhat <- sapply(chains, function(m) {
    h <- n_it %/% 2
    gelman.diag(mcmc.list(mcmc(m[1:h, ]), mcmc(m[(h + 1):(2 * h), ])),
                autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  })
  early_shift <- sapply(colnames(chains[[1]]), function(p) {
    x <- sapply(chains, function(m) as.numeric(m[, p]))
    mean((colMeans(x[seq_len(n_it %/% 10), , drop = FALSE]) - colMeans(x[(n_it %/% 2 + 1):n_it, , drop = FALSE])) /
           apply(x, 2, sd))
  })
  placebo_z <- sapply(chains, function(m) {
    later <- mcmc(as.matrix(m)[(4 * n_it %/% 10 + 1):n_it, , drop = FALSE])   # draws 401-1000
    geweke.diag(later, frac1 = 1 / 6, frac2 = 5 / 6)$z
  })
  pooled  <- do.call(rbind, lapply(chains, as.matrix))
  trimmed <- do.call(rbind, lapply(chains, function(m) as.matrix(m)[-seq_len(n_it %/% 10), , drop = FALSE]))
  chain_means <- sapply(chains, function(m) colMeans(m))
  chain_sds   <- sapply(chains, function(m) apply(m, 2, sd))
  B <- apply(chain_means, 1, var)
  W <- rowMeans(chain_sds^2)
  data.frame(parameter = colnames(chains[[1]]),
             rhat = gelman.diag(chains, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1],
             ess_per_chain = apply(ess_each, 1, median), ess_min_chain = apply(ess_each, 1, min),
             ess_total = rowSums(ess_each),
             n_chains_ess_below400 = rowSums(ess_each < 400),
             geweke_flag_rate = rowMeans(abs(geweke_z) > qnorm(0.975)),
             geweke_flag_p = sapply(rowSums(abs(geweke_z) > qnorm(0.975)),
                                    function(k) binom.test(k, length(chains), 0.05)$p.value),
             geweke_flag_positive = rowSums(geweke_z > qnorm(0.975)),
             geweke_flag_negative = rowSums(geweke_z < -qnorm(0.975)),
             geweke_placebo_rate = rowMeans(abs(placebo_z) > qnorm(0.975)),
             drop10_shift_sd = (colMeans(trimmed) - colMeans(pooled)) / apply(pooled, 2, sd),
             split_rhat_max = apply(split_rhat, 1, max),
             early_shift_sd = early_shift,
             between_within = apply(chain_means, 1, sd) / rowMeans(chain_sds),
             B_over_W = B / W,
             # B expected from Monte Carlo error alone (variance of a chain mean, W / ESS)
             B_mcmc_over_W = 1 / rowMeans(ess_each),
             n_draws = n_it * length(chains), row.names = NULL)
}

# Chained phase: imputed category shares and the share of species whose
# imputed category changed from one iteration to the next
chain_series <- function(chain) {
  rows <- chain$fit$miss_dat$row[chain$fit$miss_dat$colname == "social_max"]
  imp  <- sapply(chain$fit$data[-1], function(d) as.integer(d$social_max[rows]))
  shares <- t(apply(imp, 2, function(x) tabulate(x, K) / length(x)))
  colnames(shares) <- categories
  data.frame(iteration = seq_len(ncol(imp)), shares,
             changed = c(NA, colMeans(imp[, -1, drop = FALSE] != imp[, -ncol(imp), drop = FALSE])))
}

## -----------------------------------------------------------------------------
## 1. Who is missing? (all species, before holding any out)
## -----------------------------------------------------------------------------
missing_social <- is.na(social_data$social_max)
by_habitat <- aggregate(missing ~ habitat, FUN = mean,
                        data = data.frame(missing = missing_social, habitat = social_data$habitat))
by_habitat$n <- as.numeric(table(social_data$habitat)[by_habitat$habitat])
by_order <- aggregate(missing ~ taxon_order, FUN = mean,
                      data = data.frame(missing = missing_social, taxon_order = social_data$taxon_order))
by_order$n <- as.numeric(table(social_data$taxon_order)[by_order$taxon_order])
# Share of deviance in "missing or not" explained by traits and by family
d_miss <- data.frame(missing = missing_social, social_data[, c("log_mass", "hwi", "habitat", "taxon_family")])
d_miss <- d_miss[!is.na(d_miss$habitat), ]
null_dev <- glm(missing ~ 1, family = binomial, data = d_miss)$deviance
aic_missing <- c(null = AIC(glm(missing ~ 1, family = binomial, data = d_miss)),
                 traits = AIC(glm(missing ~ log_mass + hwi + habitat, family = binomial, data = d_miss)),
                 family = AIC(glm(missing ~ taxon_family, family = binomial, data = d_miss)))
dev_explained <- c(
  traits = 1 - glm(missing ~ log_mass + hwi + habitat, family = binomial, data = d_miss)$deviance / null_dev,
  family = 1 - glm(missing ~ taxon_family, family = binomial, data = d_miss)$deviance / null_dev)
# Species without records whose family has no species with a record: these
# resemble the whole-family held-out set rather than the random one
known_in_family <- tapply(!missing_social, social_data$taxon_family, sum)
no_relative_in_family <- known_in_family[social_data$taxon_family] == 0
# and, more generally, the share of species with a record in each missing
# species' family
family_known_share <- tapply(!missing_social, social_data$taxon_family, mean)[social_data$taxon_family]
missing_by_family_coverage <- table(cut(family_known_share[missing_social], c(-0.01, 0, 0.25, 0.5, 1),
                                        labels = c("0", "0-25%", "25-50%", "50-100%")))

data_summary <- list(
  counts = table(social_data$social_max, useNA = "ifany"),
  n_species = nrow(social_data),
  n_not_in_birdbase = sum(!social_data$in_birdbase),
  n_not_described = sum(social_data$in_birdbase & !social_data$described),
  n_lek_only = sum(social_data$described & missing_social),
  n_solitary_and_flocks = sum(social_data$solitary == 1 & social_data$flocks == 1, na.rm = TRUE),
  by_habitat = by_habitat[order(by_habitat$missing), ],
  by_order = by_order[order(-by_order$missing), ],
  dev_explained = dev_explained, aic_missing = aic_missing,
  n_families = length(unique(d_miss$taxon_family)),
  n_hwi_inferred = sum(social_data$hwi_inferred),
  n_missing_no_family_record = sum(missing_social & no_relative_in_family),
  n_families_no_record = sum(known_in_family == 0),
  missing_by_family_coverage = missing_by_family_coverage,
  holdout = table(holdout$set),
  n_block_families = length(unique(social_data$taxon_family[social_data$species %in%
                                                            holdout$species[holdout$set == "family_block"]])))
saveRDS(data_summary, file.path(res_dir, "data_summary.rds"))

## -----------------------------------------------------------------------------
## 2. Full analysis: imputations of the species without a description
## -----------------------------------------------------------------------------
a      <- make_analysis_data("full")
chain  <- readRDS(file.path(tmp_dir, "full_chain.rds"))
finals <- read_finals("full")
probs  <- Reduce(`+`, lapply(finals, `[[`, "probs")) / length(finals)
stopifnot(identical(rownames(probs), a$data$species[a$social_rows]))
imputed_cat <- sapply(finals, function(f) as.integer(f$imputed_social))   # species x imputations
imputed_counts <- t(apply(imputed_cat, 1, function(x) tabulate(x, K)))
colnames(imputed_counts) <- paste0("n_", categories)   # imputed datasets per category
imputed <- data.frame(species = a$data$species[a$social_rows],
                      taxon_family = social_data$taxon_family[a$social_rows],
                      taxon_order  = social_data$taxon_order[a$social_rows],
                      no_relative_in_family = no_relative_in_family[a$social_rows],
                      probs[, categories], imputed_counts, check.names = FALSE, row.names = NULL)

# Category shares among all species, per completed dataset: their spread
# across datasets shows the imputation uncertainty in a simple summary
shares_all <- t(apply(imputed_cat, 2, function(x) {
  s <- as.integer(social_data$social_max); s[a$social_rows] <- x
  tabulate(s, K) / length(s)
}))
colnames(shares_all) <- categories

# Downstream analysis: does wing shape differ with social behaviour? A
# phylogenetic regression of the hand-wing index on social behaviour (as an
# unordered factor, contrasts against "solitary") and log body mass, fitted to
# each completed dataset and pooled with Rubin's rules, against the same model
# fitted to the species with records only (complete-case analysis).
completed <- lapply(seq_along(finals), function(i) {
  d <- a$data
  d$social_max[a$social_rows] <- finals[[i]]$imputed_social
  d$social <- factor(d$social_max, ordered = FALSE)
  rownames(d) <- d$species
  d
})
downstream_fun <- function(d, tree) phylolm(hwi ~ social + log_mass, data = d, phy = tree, model = "BM")
# df_fun gives the complete-data residual df, so that pool_mi() uses the
# Barnard-Rubin small-sample degrees of freedom
downstream_mi <- pool_mi(with_imputations(completed, downstream_fun, tree = a$tree),
                         df_fun = function(f) f$n - f$d)
cc <- completed[[1]][-a$social_rows, ]
cc_tree <- drop.tip(a$tree, setdiff(a$tree$tip.label, cc$species))
cc_fit <- downstream_fun(cc, cc_tree)
downstream_cc <- data.frame(term = names(coef(cc_fit)), estimate = unname(coef(cc_fit)),
                            std.error = sqrt(diag(vcov(cc_fit))), row.names = NULL)

full_results <- list(
  imputed = imputed, shares_all = shares_all,
  downstream_mi = as.data.frame(downstream_mi), downstream_cc = downstream_cc,
  n_complete_case = nrow(cc),
  pooled = list(sol = do.call(rbind, lapply(finals, function(f) as.matrix(f$sol))),
                vcv = do.call(rbind, lapply(finals, function(f) as.matrix(f$vcv))),
                cp  = do.call(rbind, lapply(finals, function(f) as.matrix(f$cp)))),
  final_diagnostics = final_diagnostics(finals),
  chain_convergence = chain$convergence, chain_series = chain_series(chain),
  chain_model_diagnostics = chain$fit$diagnostics$all_diagnostics,
  settings = list(chain = chain_settings$full, final = final_settings$full, n_final = n_final[["full"]]),
  minutes = c(chained = chain$minutes, final_mean = mean(sapply(finals, `[[`, "minutes"))))
saveRDS(full_results, file.path(res_dir, "full_results.rds"))

## -----------------------------------------------------------------------------
## 3. Validation analysis: scores on the held-out species
## -----------------------------------------------------------------------------
a      <- make_analysis_data("validation")
chain  <- readRDS(file.path(tmp_dir, "validation_chain.rds"))
finals <- read_finals("validation")
probs  <- Reduce(`+`, lapply(finals, `[[`, "probs")) / length(finals)
stopifnot(identical(rownames(probs), a$data$species[a$social_rows]))

truth    <- as.integer(social_data$social_max[a$scored_rows])
families <- social_data$taxon_family[a$scored_rows]
pos      <- match(a$scored_rows, a$social_rows)
bace_p   <- unname(probs[pos, categories])
imputed_freq <- t(apply(sapply(finals, function(f) as.integer(f$imputed_social[pos])), 1,
                        function(x) tabulate(x, K) / length(x)))
# Baselines fitted to the species with known behaviour in the training data
train <- a$data[!is.na(a$data$social_max), ]
test  <- a$data[a$scored_rows, ]
#   observed category frequencies
freq <- as.numeric(table(train$social_max)) / nrow(train)
marginal <- matrix(freq, nrow(test), K, byrow = TRUE)
#   ordinal probit regression on the same predictors, without phylogeny
polr_fit <- polr(social_max ~ log_mass + hwi, data = train, method = "probit")
polr_p   <- unname(predict(polr_fit, newdata = test, type = "probs"))
#   category frequencies among the closest relatives with known behaviour
#   (k = 5, weighted ties; see relative_weights()), plus 0.5 per category.
#   k = 5 was carried over from an earlier example, not chosen from these results.
tree_sub <- drop.tip(social_tree, setdiff(social_tree$tip.label, c(train$species, test$species)))
dist <- cophenetic(tree_sub)[test$species, train$species]
nearest <- relatives_probs(dist, train$social_max, k = 5)
# number of relatives with non-zero weight (k = 5): can exceed k with ties
n_relatives <- apply(dist, 1, function(x) sum(relative_weights(x, 5) > 0))
probs_all <- list(bace = bace_p, bace_imputed_frequency = imputed_freq, marginal = marginal,
                  polr = polr_p, nearest_relatives = unname(nearest))

set.seed(28)
scores <- list(); rps_gain <- list(); confusion <- list(); prob_observed <- list()
for (s in c("random", "family_block", "all")) {
  keep <- if (s == "all") rep(TRUE, nrow(holdout)) else holdout$set == s
  p_s  <- lapply(probs_all, function(p) p[keep, , drop = FALSE])
  scores[[s]] <- cbind(set = s, score_table(p_s, truth[keep]))
  rps_bace <- rps(p_s$bace, truth[keep])
  # RPS of each baseline minus RPS of BACE (positive = BACE better)
  rps_gain[[s]] <- do.call(rbind, lapply(c("marginal", "polr", "nearest_relatives"), function(m) {
    diffs <- rps(p_s[[m]], truth[keep]) - rps_bace
    ci <- family_bootstrap(diffs, families[keep])
    data.frame(set = s, baseline = m, rps_gain = mean(diffs), lower95 = ci[["lower95"]], upper95 = ci[["upper95"]])
  }))
  confusion[[s]] <- table(observed = factor(categories[truth[keep]], levels = categories),
                          predicted = factor(categories[max.col(p_s$bace, ties.method = "first")], levels = categories))
  # Mean probability given to the observed category, by observed category
  prob_observed[[s]] <- do.call(rbind, lapply(c("bace", "nearest_relatives", "marginal"), function(m) {
    p_obs <- p_s[[m]][cbind(seq_len(sum(keep)), truth[keep])]
    data.frame(set = s, method = m, category = categories,
               n = as.numeric(table(factor(truth[keep], levels = 1:K))),
               mean_prob_observed = as.numeric(tapply(p_obs, factor(truth[keep], levels = 1:K), mean)))
  }))
}
scores <- do.call(rbind, scores); rps_gain <- do.call(rbind, rps_gain)
prob_observed <- do.call(rbind, prob_observed)

# Sensitivity of the comparison to the number of relatives (k) and to the tie
# rule; RPS of the relatives rule minus RPS of BACE, with family bootstrap
set.seed(32)
relatives_sensitivity <- do.call(rbind, lapply(c("weighted", "all"), function(ties) {
  do.call(rbind, lapply(c(1, 3, 5, 10, 20), function(k) {
    p_k <- relatives_probs(dist, train$social_max, k, ties)
    n_k <- apply(dist, 1, function(x) sum(relative_weights(x, k, ties) > 0))
    do.call(rbind, lapply(c("random", "family_block"), function(s) {
      keep <- holdout$set == s
      sc_k <- score_table(list(relatives = p_k[keep, , drop = FALSE]), truth[keep])
      diffs <- rps(p_k[keep, , drop = FALSE], truth[keep]) - rps(bace_p[keep, , drop = FALSE], truth[keep])
      ci <- family_bootstrap(diffs, families[keep])
      data.frame(ties = ties, smoothing = "fixed", k = k, set = s, rps = sc_k$rps,
                 balanced_accuracy = sc_k$balanced_accuracy,
                 recall_solitary = sc_k$recall_solitary, rps_gain_bace = mean(diffs),
                 lower95 = ci[["lower95"]], upper95 = ci[["upper95"]],
                 relatives_median = median(n_k[keep]), relatives_max = max(n_k[keep]))
    }))
  }))
}))
# The same with weighted ties and a pseudo-count proportional to k
# (0.5 * k / 5 per category; identical to the main baseline at k = 5)
set.seed(33)
relatives_sensitivity <- rbind(relatives_sensitivity, do.call(rbind, lapply(c(1, 3, 5, 10, 20), function(k) {
  p_k <- relatives_probs(dist, train$social_max, k, "weighted", pseudo = 0.5 * k / 5)
  n_k <- apply(dist, 1, function(x) sum(relative_weights(x, k, "weighted") > 0))
  do.call(rbind, lapply(c("random", "family_block"), function(s) {
    keep <- holdout$set == s
    sc_k <- score_table(list(relatives = p_k[keep, , drop = FALSE]), truth[keep])
    diffs <- rps(p_k[keep, , drop = FALSE], truth[keep]) - rps(bace_p[keep, , drop = FALSE], truth[keep])
    ci <- family_bootstrap(diffs, families[keep])
    data.frame(ties = "weighted", smoothing = "scaled", k = k, set = s, rps = sc_k$rps,
               balanced_accuracy = sc_k$balanced_accuracy,
               recall_solitary = sc_k$recall_solitary, rps_gain_bace = mean(diffs),
               lower95 = ci[["lower95"]], upper95 = ci[["upper95"]],
               relatives_median = median(n_k[keep]), relatives_max = max(n_k[keep]))
  }))
})))

# Calibration of the BACE probabilities (all held-out species)
calibration <- do.call(rbind, lapply(1:K, function(k) {
  p <- bace_p[, k]; bin <- cut(p, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  data.frame(category = categories[k], bin = levels(bin),
             mean_predicted = as.numeric(tapply(p, bin, mean)),
             observed = as.numeric(tapply(truth == k, bin, mean)),
             n = as.numeric(table(bin)), row.names = NULL)
}))

validation_results <- list(
  scores = scores, rps_gain = rps_gain, confusion = confusion, calibration = calibration,
  prob_observed = prob_observed, holdout_probs = probs_all, holdout_truth = truth,
  n_relatives = n_relatives, relatives_sensitivity = relatives_sensitivity, tie_tol = tie_tol,
  train_shares = freq,
  final_diagnostics = final_diagnostics(finals),
  chain_convergence = chain$convergence, chain_series = chain_series(chain),
  chain_model_diagnostics = chain$fit$diagnostics$all_diagnostics,
  settings = list(chain = chain_settings$validation, final = final_settings$validation,
                  n_final = n_final[["validation"]]),
  minutes = c(chained = chain$minutes, final_mean = mean(sapply(finals, `[[`, "minutes"))))
saveRDS(validation_results, file.path(res_dir, "validation_results.rds"))

## -----------------------------------------------------------------------------
## 4. Mixed-variable validation
## -----------------------------------------------------------------------------
a      <- make_analysis_data("mixed")
chain  <- readRDS(file.path(tmp_dir, "mixed_chain.rds"))
mixed_finals <- read_finals("mixed")

# Chained phase. assess_convergence() tests the share of the first category
# among the imputed values over all iterations; we rerun the same function on
# the iterations from 5 (and from 10) onwards, after the initial transient.
chain_late <- lapply(c(from5 = 5, from10 = 10), function(k) {
  fit_late <- chain$fit
  fit_late$data <- chain$fit$data[c(1, (k + 1):length(chain$fit$data))]   # keep the initial data first
  conv <- assess_convergence(fit_late, method = "summary")
  c(converged = conv$converged,
    unlist(conv$method_results$summary$stationarity_tests$social_max[c("acf", "trend_converged",
                                                                       "mean_pct_change")]))
})

mixed_p <- Reduce(`+`, lapply(mixed_finals, `[[`, "probs")) / length(mixed_finals)
mixed_p <- unname(mixed_p[a$data$species[a$scored_rows], categories])
hwi_masked <- a$scored_rows %in% a$hwi_rows

# Social behaviour: mixed analysis minus validation analysis, same held-out species
set.seed(29)
rps_mixed  <- rps(mixed_p, truth)
rps_single <- rps(bace_p, truth)
mixed_vs_validation <- do.call(rbind, lapply(c("random", "family_block"), function(s) {
  do.call(rbind, lapply(c("all", "hwi_also_missing", "hwi_observed"), function(h) {
    keep <- holdout$set == s & switch(h, all = TRUE, hwi_also_missing = hwi_masked, hwi_observed = !hwi_masked)
    d  <- rps_mixed[keep] - rps_single[keep]
    ci <- family_bootstrap(d, families[keep])
    data.frame(set = s, species_set = h, n = sum(keep), rps_mixed = mean(rps_mixed[keep]),
               rps_validation = mean(rps_single[keep]), rps_difference = mean(d),
               lower95 = ci[["lower95"]], upper95 = ci[["upper95"]])
  }))
}))

# Hand-wing index: imputed values (rows x n_final) against the measured values
hwi_true <- social_data$hwi[a$hwi_rows]
hwi_imp  <- sapply(mixed_finals, `[[`, "imputed_hwi")
hwi_train <- setdiff(seq_len(nrow(a$data)), a$hwi_rows)
lm_fit   <- lm(hwi ~ log_mass, data = social_data[hwi_train, ])
dist     <- cophenetic(social_tree)[a$data$species[a$hwi_rows], a$data$species[hwi_train]]
hwi_pred <- list(bace = rowMeans(hwi_imp),
                 observed_mean = rep(mean(social_data$hwi[hwi_train]), length(a$hwi_rows)),
                 lm_log_mass = unname(predict(lm_fit, newdata = social_data[a$hwi_rows, ])),
                 nearest_relatives = apply(dist, 1, function(x) sum(relative_weights(x, 5) * social_data$hwi[hwi_train]) / 5))
# the same with other k and with all ties at full weight
hwi_relatives_sensitivity <- do.call(rbind, lapply(c("weighted", "all"), function(ties) {
  do.call(rbind, lapply(c(1, 3, 5, 10, 20), function(k) {
    pred <- apply(dist, 1, function(x) { w <- relative_weights(x, k, ties); sum(w * social_data$hwi[hwi_train]) / sum(w) })
    data.frame(ties = ties, k = k, rmse = sqrt(mean((pred - hwi_true)^2)))
  }))
}))
hwi_scores <- data.frame(method = names(hwi_pred),
                         rmse = sapply(hwi_pred, function(p) sqrt(mean((p - hwi_true)^2))),
                         correlation = sapply(hwi_pred, function(p) suppressWarnings(cor(p, hwi_true))),
                         row.names = NULL)
hwi_scores$coverage90 <- c(mean(hwi_true >= apply(hwi_imp, 1, quantile, 0.05) &
                                  hwi_true <= apply(hwi_imp, 1, quantile, 0.95)), NA, NA, NA)
# Squared error of the closest-relatives mean minus that of BACE, resampling families
set.seed(30)
hwi_gain_d  <- (hwi_pred$nearest_relatives - hwi_true)^2 - (hwi_pred$bace - hwi_true)^2
hwi_gain_ci <- family_bootstrap(hwi_gain_d, social_data$taxon_family[a$hwi_rows])
social_missing_too <- is.na(a$data$social_max[a$hwi_rows])
hwi_by_social <- data.frame(social_missing = c(TRUE, FALSE),
                            n = c(sum(social_missing_too), sum(!social_missing_too)),
                            rmse_bace = c(sqrt(mean((hwi_pred$bace - hwi_true)[social_missing_too]^2)),
                                          sqrt(mean((hwi_pred$bace - hwi_true)[!social_missing_too]^2))),
                            rmse_nearest = c(sqrt(mean((hwi_pred$nearest_relatives - hwi_true)[social_missing_too]^2)),
                                             sqrt(mean((hwi_pred$nearest_relatives - hwi_true)[!social_missing_too]^2))))

# The hand-wing index that the social models conditioned on, per final
# imputation: a spread across imputations shows that uncertainty in the imputed
# index is passed on to the social imputations (BACE's final-phase sweeps)
hwi_used_sd <- apply(sapply(mixed_finals, `[[`, "hwi_used_by_social"), 1, sd)

# Spread of the hand-wing-index effect across final imputations, in the first
# and second halves of the chains: if the chains were still moving towards a
# common value, B/W would shrink from the first to the second half
hwi_effect_halves <- sapply(list(first = 1:500, second = 501:1000), function(rows) {
  x <- sapply(mixed_finals, function(f) f$sol[rows, "hwi"])
  var(colMeans(x)) / mean(apply(x, 2, var))
})

# Traces of the slowest-mixing parameters in four final imputations
trace <- do.call(rbind, lapply(1:4, function(k) {
  m <- final_draws(mixed_finals[[k]])
  data.frame(chain = k, draw = seq_len(nrow(m)),
             parameter = rep(c("phylo_var", "cutpoint2", "cutpoint3"), each = nrow(m)),
             value = c(m[, "phylo_var"], m[, "cutpoint2"], m[, "cutpoint3"]))
}))

mixed_results <- list(
  hwi_relatives_sensitivity = hwi_relatives_sensitivity,
  mixed_vs_validation = mixed_vs_validation,
  hwi_scores = hwi_scores, hwi_gain = c(mse_gain = mean(hwi_gain_d), hwi_gain_ci),
  hwi_by_social = hwi_by_social, hwi_true = hwi_true,
  hwi_imputed_mean = rowMeans(hwi_imp), hwi_imputed_sd = apply(hwi_imp, 1, sd),
  hwi_used_by_social_median_sd = median(hwi_used_sd),
  hwi_used_by_social_max_sd = max(hwi_used_sd),
  n_hwi_and_social_missing = sum(social_missing_too),
  sweeps = mixed_finals[[1]]$sweeps,   # stochastic passes per final imputation
  final_diagnostics = final_diagnostics(mixed_finals), trace = trace,
  hwi_effect_halves = hwi_effect_halves,
  chain_convergence = chain$convergence, chain_series = chain_series(chain), chain_late = chain_late,
  chain_hwi_mean = sapply(chain$fit$data[-1], function(d) mean(d$hwi[a$hwi_rows])),
  chain_model_diagnostics = chain$fit$diagnostics$all_diagnostics,
  settings = list(chain = chain_settings$mixed, final = final_settings$mixed, n_final = n_final[["mixed"]]),
  minutes = c(chained = chain$minutes, final_mean = mean(sapply(mixed_finals, `[[`, "minutes"))))
saveRDS(mixed_results, file.path(res_dir, "mixed_results.rds"))

# Provenance of the runs
writeLines(c("Runs on totoro (Ubuntu 24.04, R 4.6.1, BACE git dc79be1, MCMCglmm 2.36), one process per final imputation.",
             "validation, mixed: 2026-10-27 to 2026-10-28, run with scripts in an investigation workspace using the",
             "  same data, calls, settings and seeds as 02_run_bace.R (input data checked identical with identical());",
             "  their outputs were copied to results/tmp/ and renamed (social_* -> validation_*; field 'imputed' ->",
             "  'imputed_social' in the validation finals).",
             "full: 2026-10-30, run with 02_run_bace.R (chained phase 22 min; its log ends with a parse error",
             "  because the script file was replaced (header comments edited) while R was still reading it;",
             "  full_chain.rds had already been saved).",
             sprintf("Summarised on: %s, %s, R %s", Sys.info()[["nodename"]], format(Sys.time(), "%Y-%m-%d"),
                     getRversion())),
           file.path(res_dir, "run_info.txt"))

options(width = 200)
print(data_summary[c("counts", "dev_explained", "n_missing_no_family_record", "n_families_no_record")])
print(full_results$final_diagnostics, digits = 3); print(full_results$downstream_mi); print(downstream_cc)
print(scores, digits = 3, row.names = FALSE); print(rps_gain, digits = 3, row.names = FALSE)
print(prob_observed, digits = 3, row.names = FALSE)
print(validation_results$final_diagnostics, digits = 3)
print(mixed_results$final_diagnostics, digits = 3); print(chain_late)
print(mixed_vs_validation, digits = 3, row.names = FALSE); print(hwi_scores, digits = 3)
cat("Saved", paste(list.files(res_dir, "\\.rds$"), collapse = ", "), "to", res_dir, "\n")
