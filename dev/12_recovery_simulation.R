# =============================================================================
# 12_recovery_simulation.R
#
# Does BACE recover the true VALUES and PARAMETERS when data are missing —
# under MCAR *and* MAR? Multiple imputation is valid (unbiased, nominal
# coverage) under MAR when the imputation model is correct (Rubin 1976, 1987);
# this script demonstrates that empirically for BACE and contrasts it with a
# complete-case analysis.
#
# Design (per replicate):
#   tree ~ pure-birth; x  ~ Brownian motion on the tree (fully observed)
#   y = b0 + b1*x + u + e,  u ~ BM phylo effect,  e ~ N(0, sigma_e^2)
#   Missingness applied to y under one of:
#     MCAR : P(miss) constant
#     MAR  : P(miss) = plogis(a + c*x)   [depends on OBSERVED x -> MAR]
#
# Analyses (all target the slope b1 and the marginal mean of y):
#   oracle        : fit on the complete data          (benchmark)
#   complete_case : fit on rows with observed y        (the naive baseline)
#   bace          : bace_imp -> bace_final_imp ->
#                   with_imputations(gls) -> pool_mi()  (Track D pathway)
#
# Metrics aggregated across replicates, per mechanism x method:
#   slope_bias      mean(est_b1) - b1_true
#   slope_cov95     fraction of 95% CIs containing b1_true   (target 0.95)
#   slope_ci_width  mean 95% CI width
#   mean_bias       mean(estimated E[y]) - true E[y]         (the MAR signature)
#   cell_cor        cor(imputed y, truth) on hidden cells    (BACE only)
#
# Expectations a reviewer should check:
#   - Under BOTH MCAR and MAR, BACE slope_bias ~ 0 and slope_cov95 ~ 0.95.
#   - Under MAR, complete_case mean_bias is NON-zero (the observed sample is a
#     biased sample of y) while BACE mean_bias ~ 0.
#   - cell_cor is high when phylogenetic signal + x-signal are strong.
#
# Scale via env vars (smoke defaults are tiny; set for production):
#   RECOV_REPS RECOV_NSPP RECOV_NITT RECOV_BURNIN RECOV_THIN RECOV_RUNS RECOV_NFINAL
# =============================================================================

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape); library(MASS); library(nlme)
})

.envint <- function(nm, d) {
  v <- suppressWarnings(as.integer(Sys.getenv(nm, unset = NA_character_)))
  if (is.na(v)) d else v
}
CFG <- list(
  reps    = .envint("RECOV_REPS",   4L),
  nspp    = .envint("RECOV_NSPP",   60L),
  nitt    = .envint("RECOV_NITT",   3000L),
  burnin  = .envint("RECOV_BURNIN", 600L),
  thin    = .envint("RECOV_THIN",   3L),
  runs    = .envint("RECOV_RUNS",   4L),
  n_final = .envint("RECOV_NFINAL", 15L),
  b0 = 0, b1 = 0.8, sigma_e = 0.5, miss_rate = 0.30, mar_slope = 2.0
)

# ---- one simulated dataset --------------------------------------------------
simulate_one <- function(nspp, cfg, seed) {
  set.seed(seed)
  tr <- ape::rtree(nspp)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  R  <- ape::vcv(tr, corr = TRUE)
  x  <- as.numeric(MASS::mvrnorm(1, rep(0, nspp), R))          # phylo predictor
  u  <- as.numeric(MASS::mvrnorm(1, rep(0, nspp), R))          # phylo effect on y
  y  <- cfg$b0 + cfg$b1 * x + u + rnorm(nspp, 0, cfg$sigma_e)
  data.frame(y = y, x = x, Species = tr$tip.label,
             stringsAsFactors = FALSE) -> d
  list(tree = tr, data = d, y_true = y)
}

apply_missing <- function(d, mechanism, cfg, seed) {
  set.seed(seed)
  n <- nrow(d)
  p <- switch(mechanism,
    MCAR = rep(cfg$miss_rate, n),
    MAR  = {
      lin <- cfg$mar_slope * as.numeric(scale(d$x))
      # calibrate intercept so the marginal missing rate ~ miss_rate
      a <- stats::uniroot(function(a) mean(plogis(a + lin)) - cfg$miss_rate,
                          interval = c(-10, 10))$root
      plogis(a + lin)
    })
  miss <- stats::rbinom(n, 1, p) == 1
  miss[all(miss)] <- FALSE
  dm <- d; dm$y[miss] <- NA
  list(data = dm, miss = miss)
}

# ---- analysis helpers -------------------------------------------------------
fit_gls <- function(dat, tree) {
  dat$Species <- factor(dat$Species, levels = tree$tip.label)
  tryCatch(
    nlme::gls(y ~ x, data = dat,
              correlation = ape::corBrownian(1, tree, form = ~ Species),
              method = "ML"),
    error = function(e) stats::lm(y ~ x, data = dat))  # graceful fallback
}
# Residual df robust to gls (df.residual.gls returns NULL) and lm.
resid_df <- function(fit)
  as.numeric(stats::nobs(fit)) - length(stats::coef(fit))

slope_ci <- function(fit) {
  ct  <- summary(fit)$tTable          # gls: Value/Std.Error ; lm falls back below
  if (is.null(ct)) ct <- summary(fit)$coefficients
  est <- ct["x", 1]; se <- ct["x", 2]
  q   <- stats::qt(0.975, resid_df(fit))
  c(est = est, lo = est - q * se, hi = est + q * se)
}

# ---- one replicate ----------------------------------------------------------
run_rep <- function(mechanism, rep_id, cfg) {
  seed <- 1000L * rep_id + switch(mechanism, MCAR = 1L, MAR = 2L)
  sim  <- simulate_one(cfg$nspp, cfg, seed)
  mm   <- apply_missing(sim$data, mechanism, cfg, seed + 5L)
  dm   <- mm$data; miss <- mm$miss

  # oracle + complete-case (frequentist gls)
  or   <- slope_ci(fit_gls(sim$data, sim$tree))
  cc   <- slope_ci(fit_gls(dm[!miss, , drop = FALSE], sim$tree))

  # BACE imputation + Track D pooling
  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
    phylo = sim$tree, data = dm, runs = cfg$runs,
    nitt = cfg$nitt, thin = cfg$thin, burnin = cfg$burnin, verbose = FALSE)))
  fin  <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = init, fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
    phylo = sim$tree, n_final = cfg$n_final,
    nitt = cfg$nitt, thin = cfg$thin, burnin = cfg$burnin, verbose = FALSE)))
  fits <- suppressWarnings(suppressMessages(
    with_imputations(fin, function(d, tree) fit_gls(d, tree),
                     tree = sim$tree, .progress = FALSE)))
  pooled <- suppressWarnings(suppressMessages(
    pool_mi(fits, df_fun = function(f) tryCatch(resid_df(f),
                                                error = function(e) NA_real_))))
  brow <- pooled[pooled$term == "x", ]
  bace <- c(est = brow$estimate, lo = brow$conf.low, hi = brow$conf.high)

  # marginal mean of y: truth, complete-case (observed only), BACE (imputed)
  mean_true <- mean(sim$y_true)
  mean_cc   <- mean(dm$y, na.rm = TRUE)
  imp_means <- vapply(fin$all_datasets, function(z) mean(z$y), numeric(1))
  mean_bace <- mean(imp_means)

  # cell recovery on hidden cells (average imputed vs truth)
  cell_cor <- NA_real_
  if (sum(miss) >= 3L) {
    imp_hidden <- rowMeans(vapply(fin$all_datasets,
                                  function(z) z$y[miss], numeric(sum(miss))))
    cell_cor <- suppressWarnings(stats::cor(imp_hidden, sim$y_true[miss]))
  }

  data.frame(
    mechanism = mechanism, rep_id = rep_id,
    method = c("oracle", "complete_case", "bace"),
    slope_est = c(or["est"], cc["est"], bace["est"]),
    covered   = c(or["lo"] <= cfg$b1 & cfg$b1 <= or["hi"],
                  cc["lo"] <= cfg$b1 & cfg$b1 <= cc["hi"],
                  bace["lo"] <= cfg$b1 & cfg$b1 <= bace["hi"]),
    ci_width  = c(or["hi"] - or["lo"], cc["hi"] - cc["lo"],
                  bace["hi"] - bace["lo"]),
    mean_bias = c(0, mean_cc - mean_true, mean_bace - mean_true),
    cell_cor  = c(NA, NA, cell_cor),
    row.names = NULL, stringsAsFactors = FALSE)
}

# ---- main -------------------------------------------------------------------
cores <- .envint("RECOV_CORES", 6L)
cat(sprintf("Recovery sim: reps=%d nspp=%d nitt=%d n_final=%d b1=%.2f cores=%d\n",
            CFG$reps, CFG$nspp, CFG$nitt, CFG$n_final, CFG$b1, cores))
wl <- expand.grid(mech = c("MCAR", "MAR"), rep = seq_len(CFG$reps),
                  stringsAsFactors = FALSE)
rows <- parallel::mclapply(seq_len(nrow(wl)), function(i) {
  tryCatch(run_rep(wl$mech[i], wl$rep[i], CFG),
    error = function(e) { message(wl$mech[i], " rep ", wl$rep[i], " failed: ",
      e$message); NULL })
}, mc.cores = cores)
res <- do.call(rbind, rows)

agg <- do.call(rbind, lapply(split(res, list(res$mechanism, res$method)), function(g) {
  data.frame(
    mechanism      = g$mechanism[1], method = g$method[1], n = nrow(g),
    slope_bias     = round(mean(g$slope_est) - CFG$b1, 3),
    slope_cov95    = round(mean(g$covered), 3),
    slope_ci_width = round(mean(g$ci_width), 3),
    mean_bias      = round(mean(g$mean_bias), 3),
    cell_cor       = round(mean(g$cell_cor, na.rm = TRUE), 3),
    row.names = NULL)
}))
agg <- agg[order(agg$mechanism, agg$method), ]
cat("\n=== Recovery summary (aggregated across reps) ===\n")
print(agg, row.names = FALSE)
cat("\nTargets: slope_bias ~ 0, slope_cov95 ~ 0.95, BACE mean_bias ~ 0.\n")
cat("Under MAR, complete_case mean_bias should be non-zero; BACE should correct it.\n")

out_dir <- file.path("dev", "simulation_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
saveRDS(list(per_rep = res, summary = agg, cfg = CFG),
        file.path(out_dir, "recovery_simulation.rds"))
