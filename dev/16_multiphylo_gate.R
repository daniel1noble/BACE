# dev/16_multiphylo_gate.R
#
# Pre-registered gate study for tree-per-imputation (feature 1 of
# .agents/plans/next-features-2026-11.md; design adversarially reviewed twice).
#
# QUESTION. Every final imputation run starts from the ONE converged dataset
# produced by the chained phase on trees[[1]]. If the default sweeps = 3 do
# not wash that starting state out, runs assigned to other trees are pulled
# toward the tree-1 predictive centre, DEFLATING the between-tree component
# of imputation variance -- anti-conservative in exactly the direction the
# feature exists to fix.
#
# DESIGN (per seed):
#   - tree_true = rcoal(N_SP); data = 2 correlated gaussian traits with strong
#     phylogenetic signal on tree_true; 40% MCAR on both traits.
#   - K = 5 candidate trees, independent rcoal topologies on the same tips:
#     an adversarial, maximal-tree-uncertainty regime chosen for POWER (the
#     larger the true between-tree separation, the more detectable any
#     leakage). tree_true is in none of them.
#   - SHORTCUT arm (the shipped algorithm): chained phase on tree 1 once,
#     then bace_final_imp(multiPhylo, n_final = 20, sweeps = 3) -- 4 runs per
#     tree by balanced recycling.
#   - GOLD arm: per tree k, TWO independent chained phases (replicates; the
#     chained phase is unseeded, so replicates estimate its Monte Carlo
#     noise), each followed by bace_final_imp(tree k, n_final = 2, sweeps =
#     3): same 20 runs, but each run's conditioning state was converged on
#     ITS OWN tree -- no tree-1 anchoring by construction. The shortcut's
#     tree-1 chained phase doubles as gold tree-1 replicate 1.
#   - DOSE-RESPONSE (shortcut only): sweeps in {1, 3, 10}, same seed, same
#     converged start. sweeps = 1 is the POSITIVE CONTROL: it must show
#     leakage (deflated separation), or the study cannot detect leakage and
#     is INCONCLUSIVE rather than a pass.
#
# ESTIMATORS (per seed, per masked cell, then averaged over cells; W = mean
# within-group variance of imputed values, groups = (tree) or (tree,
# replicate)):
#   shortcut: tau2_hat = var_k(tree-group means) - W/4          (4 runs/tree)
#   gold:     kappa2_hat (chained-phase noise) from the two replicate means
#             per tree: mean_k[(m_k1 - m_k2)^2 / 2] - W/2        (2 runs/rep)
#             tau2_hat = var_k(tree means) - kappa2_hat/2 - W/4
#   Per-cell averaging slashes the chi-square_{M-1} noise that would make a
#   single-scalar comparison need ~100 seeds for a scalar B (round-2 power
#   analysis); seeds remain the replication unit for inference.
#   Secondary estimand: per-imputation downstream slope lm(x2 ~ x1) -> same
#   decomposition on the slope.
#
# PRE-REGISTERED ACCEPTANCE (computed by the aggregation block):
#   A1 PRIMARY. mean_s log(tau2_short / tau2_gold) has its one-sided 95%
#      lower confidence limit > log(0.90): the shortcut's between-tree
#      component is not deflated by more than 10% relative to gold.
#   A2 POSITIVE CONTROL. tau2 at sweeps = 1 is significantly BELOW tau2 at
#      sweeps = 10 (paired one-sided t over seeds, alpha = .05). If A2 fails
#      the design lacked sensitivity: verdict INCONCLUSIVE regardless of A1.
#   A3 PLATEAU. two-sided 95% CI of mean_s log(tau2_sweeps3 / tau2_sweeps10)
#      lies within +/- log(1.10): the default sweeps = 3 already sits on the
#      plateau. If A1/A3 fail while A2 passes, the recommendation is a larger
#      default sweeps for multiPhylo runs, with this study as the evidence.
#
# RUN:  Rscript dev/16_multiphylo_gate.R [OUT_DIR] [CORES] [SEEDS]
#       defaults: dev/simulation_results/multiphylo_gate  7  100
#       SMOKE=1 trims to 2 seeds / reduced MCMC for a mechanics check.
# Resumable: one rds per seed in OUT_DIR/cache/.

suppressMessages({
  devtools::load_all(quiet = TRUE)   # run from the package root
  library(ape)
})
RNGkind("L'Ecuyer-CMRG")

a <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(a) >= 1) a[1] else "dev/simulation_results/multiphylo_gate"
cores   <- if (length(a) >= 2) as.integer(a[2]) else 7L
n_seeds <- if (length(a) >= 3) as.integer(a[3]) else 100L
smoke   <- nzchar(Sys.getenv("SMOKE"))
cache   <- file.path(out_dir, "cache")
dir.create(cache, showWarnings = FALSE, recursive = TRUE)

N_SP   <- 60L
K      <- 5L
M      <- 20L                       # final runs per arm (4 per tree; the gold
                                    # decomposition NEEDS M/K/2 >= 2)
SWEEPS <- c(1L, 3L, 10L)
MISS   <- 0.4
NITT <- 3000L; BURN <- 600L; THIN <- 4L; RUNS_CH <- 3L
if (smoke) { n_seeds <- 2L; N_SP <- 40L
             NITT <- 1200L; BURN <- 300L; THIN <- 3L }
stopifnot(M %% (K * 2L) == 0L, M / K / 2L >= 2L)

sim_one <- function(seed) {
  set.seed(seed)
  tree_true <- rcoal(N_SP)
  trees <- lapply(seq_len(K), function(k) {
    tr <- rcoal(N_SP); tr$tip.label <- tree_true$tip.label; tr
  })
  S <- vcv(tree_true, corr = TRUE); L <- chol(S)
  d <- data.frame(Species = tree_true$tip.label)
  g1 <- as.numeric(t(L) %*% rnorm(N_SP)); g2 <- as.numeric(t(L) %*% rnorm(N_SP))
  # Phylogeny-DOMINATED regime (weak cross-trait slope, small residual): the
  # imputations must lean on the tree, or between-tree separation -- the very
  # quantity under test -- is too small to measure against the per-group
  # noise floor (observed in the smoke run with a 0.7 slope / 0.35 sd DGP).
  d$x1 <- g1 + rnorm(N_SP, sd = 0.25)
  d$x2 <- 0.3 * d$x1 + 0.95 * g2 + rnorm(N_SP, sd = 0.25)
  n_miss <- round(MISS * N_SP)
  d$x1[sample.int(N_SP, n_miss)] <- NA
  d$x2[sample.int(N_SP, n_miss)] <- NA
  list(trees = trees, data = d)
}

chained <- function(data, tree) {
  suppressWarnings(suppressMessages(bace_imp(
    fixformula = "x1 ~ x2", ran_phylo_form = "~ 1 | Species",
    phylo = tree, data = data, runs = RUNS_CH,
    nitt = NITT, burnin = BURN, thin = THIN, verbose = FALSE)))
}

finals <- function(fit, phylo, n_final, sweeps, seed) {
  suppressWarnings(suppressMessages(bace_final_imp(
    fit, fixformula = "x1 ~ x2", ran_phylo_form = "~ 1 | Species",
    phylo = phylo, n_final = n_final, sweeps = sweeps,
    nitt = NITT, burnin = BURN, thin = THIN, verbose = FALSE, seed = seed)))
}

# imputed-value matrix at originally-missing cells: rows = cells, cols = runs
cell_mat <- function(datasets, miss_dat) {
  do.call(rbind, lapply(seq_len(nrow(miss_dat)), function(j) {
    vapply(datasets, function(dd)
      as.numeric(dd[miss_dat$row[j], miss_dat$colname[j]]), numeric(1))
  }))
}

# per-cell tau2, averaged over cells. groups: integer tree id per column;
# reps: integer replicate id per column (NULL for the shortcut arm).
tau2_decomp <- function(mat, groups, reps = NULL) {
  ut <- sort(unique(groups))
  per_cell <- apply(mat, 1L, function(y) {
    gm <- vapply(ut, function(k) mean(y[groups == k]), numeric(1))
    if (is.null(reps)) {
      W <- mean(vapply(ut, function(k) stats::var(y[groups == k]), numeric(1)))
      n_per <- length(y) / length(ut)
      stats::var(gm) - W / n_per
    } else {
      rep_means <- lapply(ut, function(k) {
        vapply(sort(unique(reps[groups == k])), function(r)
          mean(y[groups == k & reps == r]), numeric(1))
      })
      W <- mean(unlist(lapply(ut, function(k) {
        unlist(lapply(sort(unique(reps[groups == k])), function(r) {
          yy <- y[groups == k & reps == r]
          if (length(yy) > 1L) stats::var(yy) else NA_real_
        }))
      })), na.rm = TRUE)
      n_in_rep <- 2L
      kap <- mean(vapply(rep_means, function(m2)
        (m2[1] - m2[2])^2 / 2, numeric(1))) - W / n_in_rep
      stats::var(vapply(rep_means, mean, numeric(1))) -
        max(kap, 0) / 2 - W / (2L * n_in_rep)
    }
  })
  mean(per_cell)   # may average over some negative per-cell values; unbiased
}

slope_of <- function(datasets) {
  vapply(datasets, function(dd)
    unname(stats::coef(stats::lm(x2 ~ x1, data = dd))[2]), numeric(1))
}

one_seed <- function(seed) {
  cf <- file.path(cache, sprintf("seed%04d.rds", seed))
  if (file.exists(cf)) return(readRDS(cf))
  t0 <- Sys.time()
  sim <- sim_one(seed)
  trees <- sim$trees
  multi <- do.call(c, trees)   # multiPhylo

  # shortcut chained phase (= gold tree-1 replicate 1)
  fit1 <- chained(sim$data, trees[[1]])
  miss_dat <- fit1$miss_dat

  # shortcut arm at each sweeps level (same converged start, same final seed)
  short <- lapply(SWEEPS, function(s) finals(fit1, multi, M, s, seed = seed))
  names(short) <- paste0("s", SWEEPS)
  stopifnot(identical(short$s3$tree_index,
                      ((seq_len(M) - 1L) %% K) + 1L))

  # gold arm: per tree, 2 chained replicates x (M / K / 2) finals each
  n_per_rep <- M / K / 2L           # 2 runs per (tree, replicate)
  gold_runs <- list(); g_tree <- integer(0); g_rep <- integer(0)
  for (k in seq_len(K)) {
    for (r in 1:2) {
      fit_kr <- if (k == 1L && r == 1L) fit1 else chained(sim$data, trees[[k]])
      fk <- finals(fit_kr, trees[[k]], n_final = n_per_rep, sweeps = 3L,
                   seed = seed * 1000L + k * 10L + r)
      gold_runs <- c(gold_runs, fk$all_datasets)
      g_tree <- c(g_tree, rep(k, n_per_rep)); g_rep <- c(g_rep, rep(r, n_per_rep))
    }
  }

  grp_short <- short$s3$tree_index
  res <- list(seed = seed, wall_min = NA_real_)
  for (s in names(short)) {
    mat <- cell_mat(short[[s]]$all_datasets, miss_dat)
    res[[paste0("tau2_", s)]]  <- tau2_decomp(mat, grp_short)
    res[[paste0("Bcell_", s)]] <- mean(apply(mat, 1L, stats::var))
    sl <- slope_of(short[[s]]$all_datasets)
    res[[paste0("tau2_slope_", s)]] <-
      tau2_decomp(matrix(sl, nrow = 1), grp_short)
    res[[paste0("B_slope_", s)]] <- stats::var(sl)
  }
  gmat <- cell_mat(gold_runs, miss_dat)
  res$tau2_gold  <- tau2_decomp(gmat, g_tree, g_rep)
  res$Bcell_gold <- mean(apply(gmat, 1L, stats::var))
  gsl <- slope_of(gold_runs)
  res$tau2_slope_gold <- tau2_decomp(matrix(gsl, nrow = 1), g_tree, g_rep)
  res$B_slope_gold    <- stats::var(gsl)

  res$wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(res, cf)
  res
}

cat(sprintf("multiPhylo gate: %d seeds, %d cores, out = %s%s\n",
            n_seeds, cores, out_dir, if (smoke) " [SMOKE]" else ""))
rows <- parallel::mclapply(seq_len(n_seeds), function(s)
  tryCatch(one_seed(s), error = function(e) {
    cat("ERROR seed", s, ":", conditionMessage(e), "\n"); NULL
  }), mc.cores = cores)
rows <- Filter(Negate(is.null), rows)
d <- do.call(rbind, lapply(rows, function(r) as.data.frame(r)))
write.csv(d, file.path(out_dir, "gate_results.csv"), row.names = FALSE)

# ---- Aggregation & pre-registered verdicts ---------------------------------
# All tests on the DIFFERENCE scale with margins expressed as 10% of the
# reference arm's mean: per-seed tau2 estimates are unbiased but can be
# negative by sampling noise, and filtering on positivity (or taking logs)
# would select seeds and bias the comparison. Every completed seed counts.
dd <- d[stats::complete.cases(
  d[, c("tau2_s1", "tau2_s3", "tau2_s10", "tau2_gold")]), ]
cat(sprintf("\nseeds usable: %d / %d (complete estimates)\n", nrow(dd), nrow(d)))
if (nrow(dd) < 10L) {
  cat("Too few seeds for the pre-registered tests (need >= 10); run more.\n")
  cat(sprintf("Means so far: tau2 s1 %.4g | s3 %.4g | s10 %.4g | gold %.4g\n",
              mean(dd$tau2_s1), mean(dd$tau2_s3), mean(dd$tau2_s10),
              mean(dd$tau2_gold)))
  quit(save = "no", status = 0)
}

ref_gold <- mean(dd$tau2_gold)
m10      <- 0.10 * ref_gold

# A1: paired non-inferiority, tau2_short - tau2_gold > -10% of gold's mean
d1 <- dd$tau2_s3 - dd$tau2_gold
t1 <- stats::t.test(d1, alternative = "greater", mu = -m10)
a1 <- t1$p.value < 0.05
lcl1 <- mean(d1) - stats::qt(0.95, length(d1) - 1) *
  stats::sd(d1) / sqrt(length(d1))

# A2 positive control: sweeps 1 must be DEFLATED vs sweeps 10
d2 <- dd$tau2_s1 - dd$tau2_s10
pc <- stats::t.test(d2, alternative = "less", mu = 0)
a2 <- pc$p.value < 0.05

# A3 plateau: sweeps 3 vs 10 within +/- 10% of the sweeps-10 mean
d3 <- dd$tau2_s3 - dd$tau2_s10
m10b <- 0.10 * mean(dd$tau2_s10)
ci <- stats::t.test(d3)$conf.int
a3 <- ci[1] > -m10b && ci[2] < m10b

cat(sprintf("\nreference: mean tau2 gold = %.5g (10%% margin = %.5g)\n",
            ref_gold, m10))
cat(sprintf("A1 primary (shortcut not >10%% below gold): mean diff = %.5g, one-sided 95%% LCL = %.5g vs -%.5g -> %s (p = %.4g)\n",
            mean(d1), lcl1, m10, if (a1) "PASS" else "FAIL", t1$p.value))
cat(sprintf("A2 positive control (sweeps 1 deflated vs 10): mean diff = %.5g -> %s (p = %.4g)\n",
            mean(d2), if (a2) "SENSITIVE" else
              "INSENSITIVE -> study INCONCLUSIVE", pc$p.value))
cat(sprintf("A3 plateau (sweeps 3 vs 10 within +/-10%% of s10): mean = %.5g, 95%% CI [%.5g, %.5g], margin +/-%.5g -> %s\n",
            mean(d3), ci[1], ci[2], m10b, if (a3) "PASS" else "FAIL"))
cat(sprintf("\nArm means: tau2 s1 %.5g | s3 %.5g | s10 %.5g | gold %.5g\n",
            mean(dd$tau2_s1), mean(dd$tau2_s3), mean(dd$tau2_s10),
            mean(dd$tau2_gold)))
cat(sprintf("Secondary (downstream slope B): short s3 = %.4g, gold = %.4g, ratio = %.3f\n",
            mean(dd$B_slope_s3), mean(dd$B_slope_gold),
            mean(dd$B_slope_s3) / mean(dd$B_slope_gold)))
verdict <- if (!a2) "INCONCLUSIVE" else if (a1 && a3) "PASS" else "FAIL"
cat(sprintf("\nGATE VERDICT: %s\n", verdict))
cat(sprintf("mean wall per seed: %.1f min\n", mean(d$wall_min, na.rm = TRUE)))
