# =============================================================================
# 15_sweeps_calibration.R
#
# Calibrate the final-phase `sweeps` default (2026-10 proper-MI fix) in the
# regime the fix targets: TWO variables missing at once. dev/12 and dev/13
# only ever made y missing - there a single pass is exact and sweeps are
# auto-collapsed - so this script is the first simulation that exercises the
# multi-variable final phase at all.
#
# Design (per replicate; paired across sweeps arms):
#   tree ~ pure-birth; x ~ BM on tree; y = b0 + b1*x + u + e (u ~ BM, e iid)
#   30% missingness on BOTH x and y:
#     x : MCAR
#     y : MCAR, or MAR via P(miss) = plogis(a + 2*z(x)) computed from x's
#         OBSERVED values (x-missing rows get the marginal rate, keeping the
#         mechanism MAR, never conditioned on unobserved x)
#   ONE chained phase per replicate (it does not depend on sweeps), then one
#   bace_final_imp per sweeps arm k in K_GRID, same seed -> paired arms.
#
# Metrics per mechanism x k (targets in brackets):
#   slope_cov95 via with_imputations(gls) + pool_mi  - Rubin's rules
#                 DELIBERATELY: understated between-imputation variance B
#                 shows up directly in T = W + (1+1/M)B        [~0.95]
#   slope_ci_width                                            [grows with k,
#                 plateaus at the converged-FCS reference k = max(K_GRID)]
#   B_m         between-imputation variance of per-dataset gls slopes
#                 [rises with k, plateaus]
#   sd_cell_y/x mean between-run SD of imputed values at hidden cells
#                 [k = 1 understates for y (first-visited); plateaus]
#   joint_cor   mean over rows hidden in BOTH x and y of cor across the
#                 n_final runs between the row's y and x draws
#                 [~0 at k = 1 (y never sees x's draws); non-zero after]
#   acc_cor_y   cor(mean imputed, truth) at hidden y cells  [unchanged by k]
#
# The k = 1 column doubles as the published quantification of the pre-fix
# behaviour. Acceptance for the shipped default k*: cov95(k*) within MC error
# of cov95(k_max) and >= 0.93; B_m(k*)/B_m(k_max) >= 0.85;
# sd_cell(k*)/sd_cell(k_max) >= 0.90; profiles monotone in k.
#
# Scale via env vars (smoke defaults tiny; production used 2026-10):
#   SWEEP_REPS=100 SWEEP_NSPP=60 SWEEP_NITT=3000 SWEEP_BURNIN=600 SWEEP_THIN=3
#   SWEEP_RUNS=4 SWEEP_NFINAL=15 SWEEP_CORES=7 SWEEP_KGRID="1,2,3,5,8"
# =============================================================================

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape); library(MASS); library(nlme)
})

.envint <- function(nm, d) {
  v <- suppressWarnings(as.integer(Sys.getenv(nm, unset = NA_character_)))
  if (is.na(v)) d else v
}
.envgrid <- function(nm, d) {
  v <- Sys.getenv(nm, unset = NA_character_)
  if (is.na(v)) d else as.integer(strsplit(v, ",")[[1]])
}
CFG <- list(
  reps    = .envint("SWEEP_REPS",   2L),
  nspp    = .envint("SWEEP_NSPP",   60L),
  nitt    = .envint("SWEEP_NITT",   3000L),
  burnin  = .envint("SWEEP_BURNIN", 600L),
  thin    = .envint("SWEEP_THIN",   3L),
  runs    = .envint("SWEEP_RUNS",   4L),
  n_final = .envint("SWEEP_NFINAL", 15L),
  k_grid  = .envgrid("SWEEP_KGRID", c(1L, 3L, 8L)),
  b0 = 0, b1 = 0.8, sigma_e = 0.5, miss_rate = 0.30, mar_slope = 2.0
)

simulate_one <- function(nspp, cfg, seed) {
  set.seed(seed)
  tr <- ape::rtree(nspp)
  tr <- ape::compute.brlen(tr, method = "Grafen")
  tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  R  <- ape::vcv(tr, corr = TRUE)
  x  <- as.numeric(MASS::mvrnorm(1, rep(0, nspp), R))
  u  <- as.numeric(MASS::mvrnorm(1, rep(0, nspp), R))
  y  <- cfg$b0 + cfg$b1 * x + u + rnorm(nspp, 0, cfg$sigma_e)
  list(tree = tr,
       data = data.frame(y = y, x = x, Species = tr$tip.label,
                         stringsAsFactors = FALSE),
       y_true = y, x_true = x)
}

# x MCAR everywhere; y MCAR or MAR-on-observed-x (marginal rate where x hidden)
apply_missing2 <- function(d, mechanism, cfg, seed) {
  set.seed(seed)
  n <- nrow(d)
  miss_x <- stats::rbinom(n, 1, cfg$miss_rate) == 1
  if (all(miss_x)) miss_x[1] <- FALSE

  p_y <- rep(cfg$miss_rate, n)
  if (mechanism == "MAR") {
    obs_x <- d$x[!miss_x]
    lin <- cfg$mar_slope * (d$x - mean(obs_x)) / stats::sd(obs_x)
    a <- stats::uniroot(function(a) mean(plogis(a + lin[!miss_x])) - cfg$miss_rate,
                        interval = c(-10, 10))$root
    p_y <- plogis(a + lin)
    p_y[miss_x] <- cfg$miss_rate       # never condition on unobserved x
  }
  miss_y <- stats::rbinom(n, 1, p_y) == 1
  if (all(miss_y)) miss_y[1] <- FALSE

  dm <- d
  dm$x[miss_x] <- NA
  dm$y[miss_y] <- NA
  list(data = dm, miss_x = miss_x, miss_y = miss_y)
}

fit_gls <- function(dat, tree) {
  dat$Species <- factor(dat$Species, levels = tree$tip.label)
  tryCatch(
    nlme::gls(y ~ x, data = dat,
              correlation = ape::corBrownian(1, tree, form = ~ Species),
              method = "ML"),
    error = function(e) stats::lm(y ~ x, data = dat))
}
resid_df <- function(fit)
  as.numeric(stats::nobs(fit)) - length(stats::coef(fit))

run_rep <- function(mechanism, rep_id, cfg) {
  seed <- 1000L * rep_id + switch(mechanism, MCAR = 1L, MAR = 2L)
  sim  <- simulate_one(cfg$nspp, cfg, seed)
  mm   <- apply_missing2(sim$data, mechanism, cfg, seed + 5L)
  dm   <- mm$data

  # One chained phase per replicate; the sweeps arms share it (paired design).
  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x", ran_phylo_form = "~ 1 | Species",
    phylo = sim$tree, data = dm, runs = cfg$runs,
    nitt = cfg$nitt, thin = cfg$thin, burnin = cfg$burnin, verbose = FALSE)))

  out <- vector("list", length(cfg$k_grid))
  for (ki in seq_along(cfg$k_grid)) {
    k <- cfg$k_grid[ki]
    t0 <- Sys.time()
    fin <- suppressWarnings(suppressMessages(bace_final_imp(
      bace_object = init, fixformula = "y ~ x",
      ran_phylo_form = "~ 1 | Species", phylo = sim$tree,
      n_final = cfg$n_final, sweeps = k, seed = seed + 9L,
      nitt = cfg$nitt, thin = cfg$thin, burnin = cfg$burnin, verbose = FALSE)))
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    fits <- suppressWarnings(suppressMessages(
      with_imputations(fin, function(d, tree) fit_gls(d, tree),
                       tree = sim$tree, .progress = FALSE)))
    pooled <- suppressWarnings(suppressMessages(
      pool_mi(fits, df_fun = function(f) tryCatch(resid_df(f),
                                                  error = function(e) NA_real_))))
    brow <- pooled[pooled$term == "x", ]

    # with_imputations returns the list of fits itself (class bace_mi_fits)
    slopes <- vapply(fits, function(f) stats::coef(f)[["x"]], numeric(1))

    y_mat <- vapply(fin$all_datasets, function(z) z$y[mm$miss_y],
                    numeric(sum(mm$miss_y)))
    x_mat <- vapply(fin$all_datasets, function(z) z$x[mm$miss_x],
                    numeric(sum(mm$miss_x)))
    sd_cell_y <- mean(apply(y_mat, 1, stats::sd))
    sd_cell_x <- mean(apply(x_mat, 1, stats::sd))

    both <- which(mm$miss_y & mm$miss_x)
    joint_cor <- NA_real_
    if (length(both) >= 3L) {
      yb <- vapply(fin$all_datasets, function(z) z$y[both], numeric(length(both)))
      xb <- vapply(fin$all_datasets, function(z) z$x[both], numeric(length(both)))
      rr <- vapply(seq_along(both), function(r)
        suppressWarnings(stats::cor(yb[r, ], xb[r, ])), numeric(1))
      joint_cor <- mean(rr, na.rm = TRUE)
    }

    acc_cor_y <- suppressWarnings(stats::cor(rowMeans(y_mat),
                                             sim$y_true[mm$miss_y]))

    out[[ki]] <- data.frame(
      mechanism = mechanism, rep_id = rep_id, sweeps = k,
      slope_est = brow$estimate,
      covered   = brow$conf.low <= cfg$b1 & cfg$b1 <= brow$conf.high,
      ci_width  = brow$conf.high - brow$conf.low,
      B_m       = stats::var(slopes),
      sd_cell_y = sd_cell_y, sd_cell_x = sd_cell_x,
      joint_cor = joint_cor, acc_cor_y = acc_cor_y,
      wall_sec  = wall,
      row.names = NULL, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# ---- main -------------------------------------------------------------------
cores <- .envint("SWEEP_CORES", 2L)
cat(sprintf(
  "Sweeps calibration: reps=%d nspp=%d nitt=%d n_final=%d k={%s} cores=%d\n",
  CFG$reps, CFG$nspp, CFG$nitt, CFG$n_final,
  paste(CFG$k_grid, collapse = ","), cores))

wl <- expand.grid(mech = c("MCAR", "MAR"), rep = seq_len(CFG$reps),
                  stringsAsFactors = FALSE)
t_all <- Sys.time()
rows <- parallel::mclapply(seq_len(nrow(wl)), function(i) {
  tryCatch(run_rep(wl$mech[i], wl$rep[i], CFG),
    error = function(e) { message(wl$mech[i], " rep ", wl$rep[i], " failed: ",
      e$message); NULL })
}, mc.cores = cores)
res <- do.call(rbind, rows)
cat(sprintf("Total wall time: %.1f min; %d/%d replicate rows OK\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins")),
            length(unique(paste(res$mechanism, res$rep_id))), nrow(wl)))

agg <- do.call(rbind, lapply(split(res, list(res$mechanism, res$sweeps)),
  function(g) data.frame(
    mechanism   = g$mechanism[1], sweeps = g$sweeps[1], n = nrow(g),
    slope_bias  = round(mean(g$slope_est) - CFG$b1, 3),
    slope_cov95 = round(mean(g$covered), 3),
    ci_width    = round(mean(g$ci_width), 3),
    B_m         = signif(mean(g$B_m), 3),
    sd_cell_y   = round(mean(g$sd_cell_y), 3),
    sd_cell_x   = round(mean(g$sd_cell_x), 3),
    joint_cor   = round(mean(g$joint_cor, na.rm = TRUE), 3),
    acc_cor_y   = round(mean(g$acc_cor_y, na.rm = TRUE), 3),
    wall_sec    = round(mean(g$wall_sec), 1),
    row.names = NULL)))
agg <- agg[order(agg$mechanism, agg$sweeps), ]
cat("\n=== Sweeps calibration summary ===\n")
print(agg, row.names = FALSE)
cat("\nRead: k = 1 is the pre-fix behaviour; the shipped default must sit on",
    "\nthe k = max plateau for cov95, B_m and sd_cell, at acceptable wall_sec.\n")

out_dir <- file.path("dev", "simulation_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
saveRDS(list(per_rep = res, summary = agg, cfg = CFG,
             pkg_commit = tryCatch(system("git rev-parse --short HEAD",
                                          intern = TRUE),
                                   error = function(e) NA_character_),
             date = Sys.time()),
        file.path(out_dir, "sweeps_calibration.rds"))
cat("\nSaved: dev/simulation_results/sweeps_calibration.rds\n")
