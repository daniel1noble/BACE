# code/rubin_stack_posthoc.R
#
# Post-hoc arm for the completed 2026-10 campaign: score BACE's NATIVE combiner
# (posterior stacking, as in pool_posteriors()) on the downstream slope, next to
# the Rubin's-rules arm already in the rds files.
#
# For every dataset with stored completed datasets ($imputations$bace; saved at
# n <= 300 by the launcher), run BACE's own decoupled-analysis workflow:
#   BACE::with_imputations()  fits the Bayesian ANALYSIS model
#                             MCMCglmm(c2 ~ c1, random = ~Species, BM phylo ginverse)
#                             to each of the M = 20 completed datasets;
#   BACE::pool_posteriors()   stacks the M posteriors (the fits are wrapped in the
#                             bace_final shape it expects);
#   BACE::get_pooled_model()  reads the pooled model; its stacked Sol[, "c1"] gives
#                             the posterior mean and central 95% interval,
# scored against the known truth (rho).
#
#   Rscript code/rubin_stack_posthoc.R [RESULTS_DIR] [OUT_DIR] [CORES]
#
# Resumable: one small rds per dataset in OUT_DIR/stack_cache/; existing files
# are skipped. Final outputs: stack_down.csv (dataset level), agg_stack.csv
# (scenario level), printed sanity lines.

suppressMessages({ library(MCMCglmm); library(ape); library(BACE) })
RNGkind("L'Ecuyer-CMRG")   # rubin_cell.R:31 — required to reproduce make_cell()'s datasets exactly

a <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(a) >= 1) a[1] else "results_2026-10"
out_dir <- if (length(a) >= 2) a[2] else res_dir
cores   <- if (length(a) >= 3) as.integer(a[3]) else 7L
cache   <- file.path(out_dir, "stack_cache"); dir.create(cache, showWarnings = FALSE, recursive = TRUE)

here <- "code"
source(file.path(here, "campaign_gnn_off_lib.R"))   # make_cell()

NITT <- 13000; BURNIN <- 3000; THIN <- 10           # 1000 draws per fit x M = 20000 stacked
prior <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),  # par-expanded
              R = list(V = 1, nu = 0.002))

fs <- list.files(res_dir, "^rubin_types_mixed_.*_n(100|300)_.*[.]rds$", full.names = TRUE)
cat(sprintf("candidate files (n<=300): %d\n", length(fs)))

one <- function(f) {
  key <- sub("[.]rds$", "", basename(f))
  cf <- file.path(cache, paste0(key, ".rds"))
  if (file.exists(cf)) return(readRDS(cf))
  x <- readRDS(f)
  if (is.null(x$imputations$bace)) { saveRDS(NULL, cf); return(NULL) }   # failed dataset: no arm
  cell <- make_cell("types_mixed", x$n, x$seed, miss_frac = x$miss_frac, miss = x$miss,
                    lambda = x$lambda, rho = x$rho, thresholds = x$thresholds, driver = x$driver)
  # tripwire: the rebuilt cell must reproduce the stored dataset exactly
  stopifnot(isTRUE(all.equal(cell$truth$c1, x$truth$c1)),
            isTRUE(all.equal(cell$truth$c2, x$truth$c2)),
            identical(rownames(cell$truth), rownames(x$truth)))
  tree_b <- cell$tree
  if (any(tree_b$edge.length == 0)) tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8  # as in fit_bace_mi
  Ainv <- inverseA(tree_b, nodes = "TIPS")$Ainv
  rn_fallback <- rownames(x$truth)

  # 1) BACE's model-agnostic MI layer fits the analysis model to each completed dataset
  set.seed(x$seed * 1000L)   # fits run sequentially, so one seed makes the whole set reproducible
  fits <- suppressWarnings(with_imputations(x$imputations$bace, function(d) {
    rn <- rownames(d)
    if (is.null(rn) || all(rn == as.character(seq_len(nrow(d))))) rn <- rn_fallback
    dd <- data.frame(c1 = d$c1, c2 = d$c2, Species = factor(rn))
    MCMCglmm(c2 ~ c1, random = ~Species, ginverse = list(Species = Ainv),
             data = dd, prior = prior, nitt = NITT, burnin = BURNIN, thin = THIN,
             verbose = FALSE, pr = FALSE)
  }, .progress = FALSE))
  n_fail <- attr(fits, "n_failed")
  ok <- Filter(function(f) inherits(f, "MCMCglmm"), fits)
  if (!length(ok)) { saveRDS(NULL, cf); return(NULL) }

  # 2) BACE's native combiner stacks the per-imputation posteriors
  bf <- structure(list(all_models = lapply(ok, function(f) list(c2 = f))), class = "bace_final")
  pooled <- pool_posteriors(bf)
  pm <- get_pooled_model(pooled, variable = "c2")
  draws <- as.numeric(pm$Sol[, "c1"])
  sols <- ok   # length used below as m_ok
  tr <- x$estimands$truth[x$estimands$arm == "bace" & x$estimands$estimand == "slope"][1]
  if (!length(tr) || is.na(tr)) tr <- x$rho
  lo <- unname(quantile(draws, 0.025)); hi <- unname(quantile(draws, 0.975))
  row <- data.frame(n = x$n, lambda = x$lambda, rho = x$rho, seed = x$seed,
                    estimand = "slope", method = "stack",
                    estimate = mean(draws), lower = lo, upper = hi, truth = tr,
                    covered = tr >= lo & tr <= hi, ciw = hi - lo,
                    m_ok = length(sols), m_fail = n_fail)
  saveRDS(row, cf); row
}

rows <- parallel::mclapply(fs, function(f) tryCatch(one(f), error = function(e) {
  cat("ERROR", basename(f), ":", conditionMessage(e), "\n"); NULL }), mc.cores = cores)
sd_tab <- do.call(rbind, Filter(Negate(is.null), rows)); rownames(sd_tab) <- NULL
write.csv(sd_tab, file.path(out_dir, "stack_down.csv"), row.names = FALSE)

se <- function(v) sd(v) / sqrt(length(v))
agg <- do.call(rbind, lapply(split(sd_tab, sd_tab[c("n", "lambda", "rho")], drop = TRUE),
  function(g) data.frame(n = g$n[1], lambda = g$lambda[1], rho = g$rho[1],
    coverage = mean(g$covered), coverage_se = se(as.numeric(g$covered)),
    bias = mean(g$estimate - g$truth), bias_se = se(g$estimate - g$truth),
    rmse = sqrt(mean((g$estimate - g$truth)^2)), ci_width = mean(g$ciw),
    datasets = nrow(g), fit_failures = sum(g$m_fail), row.names = NULL)))
write.csv(agg, file.path(out_dir, "agg_stack.csv"), row.names = FALSE)
cat("datasets scored:", nrow(sd_tab), " total MCMCglmm fit failures:", sum(sd_tab$m_fail), "\n")
print(agg, digits = 3, row.names = FALSE)
