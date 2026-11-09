# =============================================================================
# 13_response_type_simulation.R  (v3, 2026-08)
#
# v3 changes (post bug-fix re-validation):
#  * Metrics now separate MODEL calibration from FINITE-DRAW artifacts.
#    The v2 continuous "cover95" (empirical 2.5-97.5% quantiles of n_final
#    draws) has a hard ceiling of ~0.845 at n_final=15 even for a PERFECTLY
#    calibrated imputer (order statistics: the [min,max] of m draws covers
#    (m-1)/(m+1)); the v2 discrete "95% set" degenerates to the full draw
#    support at 1/15 resolution. v3 keeps both v2 metrics (comparability with
#    the 2026-07 run) and adds:
#      cover95_t     (cont): mean +/- t_{m-1,.975} * sd(draws) * sqrt(1+1/m)
#                    — a Rubin-style interval NOT capped by finite m.
#      recovery_spearman, recovery_log (cont): rank / log1p-scale recovery,
#                    since Pearson-on-raw-counts is tail-dominated for poisson.
#      *_prob (disc): metrics from the models' AVERAGED class-probability
#                    matrices (fin$all_prob_preds) instead of draw frequencies:
#                    point = argmax, 95% set = smallest cumulative-prob set.
#      bayes_bal_acc (categorical): per-rep Bayes-classifier balanced accuracy
#                    computed from the TRUE simulated class probabilities
#                    (sim_bace$response_probs) — the DGP ceiling to report
#                    balanced accuracy against.
#  * Categorical DGP de-degenerated: explicit strong per-liability betas via
#    beta_resp_liab (the v2 run silently used runif(-0.5,0.5) shared slopes,
#    whose Bayes ceiling was ~0.50 balanced accuracy for 3 classes).
#
# Comprehensive per-response-type recovery study. For EVERY response type
# (gaussian, poisson/count, binary, categorical/multinomial, ordinal/threshold)
# and both MCAR and MAR, over many independent replicates, we report:
#
#   value recovery : correlation (cont/count) or accuracy + balanced accuracy
#                    (discrete), on hidden cells.
#   BIAS           : cont/count  = mean(point-estimate - truth)      [signed]
#                    binary       = P_hat(positive) - P(positive)     [signed]
#                    ordinal      = mean(E[rank]) - mean(true rank)    [signed]
#                    categorical  = total-variation of class props     [>= 0]
#   INTERVAL/SET CALIBRATION (95%):
#                    cont/count  = P(truth in [2.5%, 97.5%] of the n_final draws)
#                    discrete    = P(true class in the smallest class set whose
#                                  imputation-frequency reaches 0.95)
#   timing         : BACE wall-clock seconds.
#
# Data-generating model (sim_bace): birth-death tree; two gaussian predictors
# x1, x2; response of the target type from a KNOWN linear predictor
# (beta_resp) plus a phylogenetic random effect. Phylogenetic signal HIGH
# (Pagel lambda = 0.90) on all variables. Only y is made missing (30%).
# MCAR = constant; MAR = plogis(a + 2*z(x1)) (depends on observed x1).
#
# n = 1000 replicates per (type x mechanism) cell by default. The run is
# RESUMABLE: each cell's results are written to chunks/<type>_<mech>.rds as
# soon as it finishes; a restart skips completed cells.
#
# Env vars: RESP_REPS RESP_NSPP RESP_NITT RESP_BURNIN RESP_THIN RESP_RUNS
#           RESP_NFINAL RESP_CORES RESP_PHYLO
# =============================================================================

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape); library(MASS); library(parallel)
})

.envint <- function(nm, d) { v <- suppressWarnings(as.integer(Sys.getenv(nm, unset = NA))); if (is.na(v)) d else v }
.envnum <- function(nm, d) { v <- suppressWarnings(as.numeric(Sys.getenv(nm, unset = NA))); if (is.na(v)) d else v }
CFG <- list(
  reps    = .envint("RESP_REPS",   1000L),
  nspp    = .envint("RESP_NSPP",   80L),
  nitt    = .envint("RESP_NITT",   4000L),
  burnin  = .envint("RESP_BURNIN", 800L),
  thin    = .envint("RESP_THIN",   4L),
  runs    = .envint("RESP_RUNS",   3L),
  n_final = .envint("RESP_NFINAL", 15L),
  cores   = .envint("RESP_CORES",  7L),
  phylo   = .envnum("RESP_PHYLO",  0.90),
  rate = 0.30, mar_slope = 2.0,
  beta = list(x1 = 0.9, x2 = -0.6)      # KNOWN predictor effects (strong signal)
)

TYPES <- list(
  gaussian    = list(resp = "gaussian",     kind = "cont"),
  poisson     = list(resp = "poisson",      kind = "cont"),
  binary      = list(resp = "binary",       kind = "disc"),
  ordinal     = list(resp = "threshold3",   kind = "disc"),
  categorical = list(resp = "multinomial3", kind = "disc"))
MECHS <- c("MCAR", "MAR")

OUT   <- file.path("dev", "simulation_results")
CHUNK <- file.path(OUT, "response_type_chunks")
if (!dir.exists(CHUNK)) dir.create(CHUNK, recursive = TRUE)

# ---- simulate + missingness -------------------------------------------------
post_process <- function(d, type) {
  if (type == "poisson")     d$y <- as.integer(round(d$y))
  if (type == "binary")      d$y <- factor(d$y, levels = sort(unique(d$y)))
  if (type == "categorical") d$y <- factor(as.character(d$y))
  if (type == "ordinal")     d$y <- factor(d$y, levels = sort(unique(d$y)), ordered = TRUE)
  d
}
simulate_one <- function(type, cfg, seed) {
  set.seed(seed)
  # KNOWN predictor effects. Poisson uses a gentler linear predictor so exp()
  # does not produce runaway counts. Categorical uses KNOWN strong
  # per-liability betas (distinct per liability, so the multinomial DGP is
  # non-degenerate; see sim_bace beta_resp_liab).
  is_cat <- TYPES[[type]]$resp == "multinomial3"
  br <- if (type == "poisson") list(x1 = 0.35, x2 = -0.25)
        else list(x1 = 0.9, x2 = -0.6)
  bl <- if (is_cat) list(x1 = c(0.9, -0.4), x2 = c(-0.6, 0.5)) else NULL
  sim <- suppressMessages(sim_bace(
    response_type   = TYPES[[type]]$resp,
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    beta_resp       = br, beta_resp_liab = bl, beta_sparsity = 0.0,
    phylo_signal    = rep(cfg$phylo, 3),
    n_cases = cfg$nspp, n_species = cfg$nspp, missingness = c(0, 0, 0)))
  d <- sim$complete_data
  names(d)[names(d) == "species"] <- "Species"
  list(data = post_process(d, type), tree = sim$tree,
       true_probs = sim$response_probs)
}
apply_missing <- function(d, mechanism, cfg, seed) {
  set.seed(seed); n <- nrow(d)
  p <- if (mechanism == "MCAR") rep(cfg$rate, n) else {
    lin <- cfg$mar_slope * as.numeric(scale(d$x1))
    a <- stats::uniroot(function(a) mean(plogis(a + lin)) - cfg$rate, c(-12, 12))$root
    plogis(a + lin)
  }
  miss <- stats::rbinom(n, 1, p) == 1
  if (is.factor(d$y)) for (lv in levels(d$y)) {          # keep >=2 observed/level
    idx <- which(d$y == lv); obs <- idx[!miss[idx]]
    if (length(obs) < 2L && length(idx) >= 2L) miss[sample(idx, 2L)] <- FALSE
  }
  if (all(miss)) miss[sample(n, 1L)] <- FALSE
  miss
}

# ---- metrics ----------------------------------------------------------------
# Bias is reported as STANDARDISED bias = (mean predicted - mean true) / SD(true)
# on a numeric encoding, so it is comparable across gaussian/poisson/binary/
# ordinal. Categorical (unordered) has no signed scalar bias, so it is reported
# as the total-variation distance of the marginal class proportions (>= 0).
.sdiv <- function(num, s) if (!is.finite(s) || s == 0) NA_real_ else num / s
metrics_cont <- function(true_v, imp_mat) {
  point <- rowMeans(imp_mat)
  m <- ncol(imp_mat)
  lo <- apply(imp_mat, 1, stats::quantile, 0.025, names = FALSE)
  hi <- apply(imp_mat, 1, stats::quantile, 0.975, names = FALSE)
  # t-interval: not capped by the finite-draw quantile ceiling (m-1)/(m+1)
  sdd  <- apply(imp_mat, 1, stats::sd)
  half <- stats::qt(0.975, df = m - 1) * sdd * sqrt(1 + 1 / m)
  list(recovery = suppressWarnings(stats::cor(point, true_v)),
       recovery_spearman = suppressWarnings(
         stats::cor(point, true_v, method = "spearman")),
       recovery_log = suppressWarnings(
         stats::cor(log1p(pmax(point, 0)), log1p(pmax(true_v, 0)))),
       bias = .sdiv(mean(point - true_v), stats::sd(true_v)),
       cover95 = mean(true_v >= lo & true_v <= hi),
       cover95_t = mean(true_v >= point - half & true_v <= point + half))
}
.greedy_setcov <- function(P, true_v, levs) {
  mean(vapply(seq_len(nrow(P)), function(i) {
    o <- order(P[i, ], decreasing = TRUE); k <- which(cumsum(P[i, o]) >= 0.95)[1]
    true_v[i] %in% levs[o[seq_len(k)]] }, logical(1)))
}
.bal_acc <- function(point, true_v, levs) {
  rec <- vapply(levs, function(l) { i <- true_v == l
    if (!any(i)) NA_real_ else mean(point[i] == l) }, numeric(1))
  mean(rec, na.rm = TRUE)
}
metrics_disc <- function(true_v, imp_mat, levs, kind, P_model = NULL) {
  P <- t(apply(imp_mat, 1, function(r)
    as.numeric(table(factor(r, levels = levs))) / length(r)))
  colnames(P) <- levs
  point <- levs[max.col(P, ties.method = "first")]
  acc <- mean(point == true_v)
  bal <- .bal_acc(point, true_v, levs)
  setcov <- .greedy_setcov(P, true_v, levs)
  # Probability-based variants from the models' averaged class-probability
  # matrices (not limited by 1/n_final draw-frequency resolution).
  acc_p <- bal_p <- setcov_p <- NA_real_
  if (!is.null(P_model)) {
    point_p  <- levs[max.col(P_model, ties.method = "first")]
    acc_p    <- mean(point_p == true_v)
    bal_p    <- .bal_acc(point_p, true_v, levs)
    setcov_p <- .greedy_setcov(P_model, true_v, levs)
  }
  if (kind == "categorical") {
    phat <- colMeans(P)
    ptrue <- as.numeric(table(factor(true_v, levels = levs))) / length(true_v)
    bias <- 0.5 * sum(abs(phat - ptrue))                  # total variation
  } else if (kind == "binary") {
    pos <- levs[length(levs)]
    tnum <- as.numeric(true_v == pos); pnum <- P[, pos]   # predicted P(positive)
    bias <- .sdiv(mean(pnum) - mean(tnum), stats::sd(tnum))
  } else {                                                # ordinal: expected rank
    ranks <- seq_along(levs); tnum <- match(true_v, levs)
    pnum <- as.numeric(P %*% ranks)
    bias <- .sdiv(mean(pnum) - mean(tnum), stats::sd(tnum))
  }
  list(recovery = acc, balanced_accuracy = bal, bias = bias, cover95 = setcov,
       recovery_prob = acc_p, balanced_accuracy_prob = bal_p,
       cover95_prob = setcov_p)
}

# ---- one replicate ----------------------------------------------------------
run_rep <- function(type, mechanism, rep_id, cfg) {
  seed <- 5000L * rep_id + switch(mechanism, MCAR = 1L, MAR = 2L) +
    100000L * match(type, names(TYPES))
  sim <- simulate_one(type, cfg, seed)
  d <- sim$data; tree <- sim$tree
  miss <- apply_missing(d, mechanism, cfg, seed + 7L)
  dm <- d; dm$y[miss] <- NA

  t0 <- Sys.time()
  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x1 + x2", ran_phylo_form = "~ 1 | Species", phylo = tree,
    data = dm, runs = cfg$runs, nitt = cfg$nitt, thin = cfg$thin,
    burnin = cfg$burnin, verbose = FALSE)))
  fin <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = init, fixformula = "y ~ x1 + x2",
    ran_phylo_form = "~ 1 | Species", phylo = tree, n_final = cfg$n_final,
    nitt = cfg$nitt, thin = cfg$thin, burnin = cfg$burnin, verbose = FALSE)))
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  hidden <- which(miss)
  na <- NA_real_
  if (TYPES[[type]]$kind == "cont") {
    true_v <- as.numeric(d$y[hidden])
    imp_mat <- vapply(fin$all_datasets, function(z) as.numeric(z$y[hidden]),
                      numeric(length(hidden)))
    m <- metrics_cont(true_v, imp_mat)
    out <- list(recovery = m$recovery, balanced_accuracy = na,
      bias = m$bias, cover95 = m$cover95, cover95_t = m$cover95_t,
      recovery_spearman = m$recovery_spearman, recovery_log = m$recovery_log,
      recovery_prob = na, balanced_accuracy_prob = na, cover95_prob = na,
      bayes_bal_acc = na)
  } else {
    levs <- levels(d$y)
    true_v <- as.character(d$y[hidden])
    imp_mat <- vapply(fin$all_datasets, function(z) as.character(z$y[hidden]),
                      character(length(hidden)))
    # Average class-probability matrix over the n_final runs (rows = hidden
    # cells in miss_dat row order = ascending, matching which(miss)).
    P_model <- tryCatch({
      mats <- lapply(fin$all_prob_preds, function(pp) {
        pm <- pp[["y"]]
        if (is.null(pm) || nrow(pm) != length(hidden) ||
            !all(levs %in% colnames(pm))) return(NULL)
        as.matrix(pm[, levs, drop = FALSE])
      })
      mats <- Filter(Negate(is.null), mats)
      if (length(mats) == 0) NULL else Reduce(`+`, mats) / length(mats)
    }, error = function(e) NULL)
    m <- metrics_disc(true_v, imp_mat, levs, kind = type, P_model = P_model)
    # DGP ceiling (categorical only): Bayes classifier from the TRUE
    # simulated class probabilities, scored on the same hidden cells.
    bayes <- na
    if (type == "categorical" && !is.null(sim$true_probs)) {
      tp <- sim$true_probs[hidden, , drop = FALSE]
      bayes_point <- colnames(tp)[max.col(tp, ties.method = "first")]
      bayes <- .bal_acc(bayes_point, true_v, levs)
    }
    out <- list(recovery = m$recovery, balanced_accuracy = m$balanced_accuracy,
      bias = m$bias, cover95 = m$cover95, cover95_t = na,
      recovery_spearman = na, recovery_log = na,
      recovery_prob = m$recovery_prob,
      balanced_accuracy_prob = m$balanced_accuracy_prob,
      cover95_prob = m$cover95_prob, bayes_bal_acc = bayes)
  }
  data.frame(type = type, mechanism = mechanism, rep_id = rep_id, status = "ok",
    runtime_sec = runtime, recovery = out$recovery,
    balanced_accuracy = out$balanced_accuracy,
    bias = out$bias, cover95 = out$cover95, cover95_t = out$cover95_t,
    recovery_spearman = out$recovery_spearman, recovery_log = out$recovery_log,
    recovery_prob = out$recovery_prob,
    balanced_accuracy_prob = out$balanced_accuracy_prob,
    cover95_prob = out$cover95_prob, bayes_bal_acc = out$bayes_bal_acc,
    n_missing = length(hidden),
    stringsAsFactors = FALSE, row.names = NULL)
}

err_row <- function(type, mech, rep_id, msg) data.frame(type = type,
  mechanism = mech, rep_id = rep_id, status = paste("error:", msg),
  runtime_sec = NA, recovery = NA, balanced_accuracy = NA, bias = NA,
  cover95 = NA, cover95_t = NA, recovery_spearman = NA, recovery_log = NA,
  recovery_prob = NA, balanced_accuracy_prob = NA, cover95_prob = NA,
  bayes_bal_acc = NA, n_missing = NA, stringsAsFactors = FALSE)

# ---- main: run each cell, save its chunk, resume if present -----------------
cat(sprintf("Response-type sim v2: %d types x %d mech x %d reps | nspp=%d phylo=%.2f nitt=%d n_final=%d cores=%d\n",
            length(TYPES), length(MECHS), CFG$reps, CFG$nspp, CFG$phylo,
            CFG$nitt, CFG$n_final, CFG$cores))
cat("Started:", format(Sys.time()), "\n")
t_all <- Sys.time()

for (type in names(TYPES)) for (mech in MECHS) {
  cf <- file.path(CHUNK, sprintf("%s_%s.rds", type, mech))
  if (file.exists(cf)) { cat(sprintf("[skip] %s %s (chunk exists)\n", type, mech)); next }
  tc <- Sys.time()
  rows <- mclapply(seq_len(CFG$reps), function(r)
    tryCatch(run_rep(type, mech, r, CFG),
             error = function(e) err_row(type, mech, r, conditionMessage(e))),
    mc.cores = CFG$cores)
  cell <- do.call(rbind, rows)
  saveRDS(cell, cf)
  ok <- sum(cell$status == "ok")
  cat(sprintf("[%s] %-11s %-4s : %d/%d ok, %.1f min\n", format(Sys.time(), "%H:%M"),
      type, mech, ok, CFG$reps, as.numeric(difftime(Sys.time(), tc, units = "mins"))))
}

# ---- aggregate all chunks ---------------------------------------------------
res <- do.call(rbind, lapply(list.files(CHUNK, pattern = "\\.rds$", full.names = TRUE), readRDS))
saveRDS(list(per_rep = res, cfg = CFG), file.path(OUT, "response_type_simulation.rds"))
utils::write.csv(res, file.path(OUT, "response_type_per_rep.csv"), row.names = FALSE)

ok <- res[res$status == "ok", ]
se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
mn <- function(x) round(mean(x, na.rm = TRUE), 3)
agg <- do.call(rbind, lapply(split(ok, list(ok$type, ok$mechanism), drop = TRUE), function(g)
  data.frame(type = g$type[1], mechanism = g$mechanism[1], n = nrow(g),
    runtime_s = round(mean(g$runtime_sec), 1),
    recovery = mn(g$recovery),
    recovery_spearman = mn(g$recovery_spearman),
    recovery_log = mn(g$recovery_log),
    balanced_accuracy = mn(g$balanced_accuracy),
    balanced_accuracy_prob = mn(g$balanced_accuracy_prob),
    bayes_bal_acc = mn(g$bayes_bal_acc),
    bias = round(median(g$bias[is.finite(g$bias)]), 3),   # robust: poisson tail
    bias_mean = mn(g$bias),
    bias_mcse = round(se(g$bias), 4),
    cover95 = mn(g$cover95),
    cover95_t = mn(g$cover95_t),
    cover95_prob = mn(g$cover95_prob),
    cover95_mcse = round(se(g$cover95), 4), row.names = NULL)))
agg <- agg[order(match(agg$type, names(TYPES)), agg$mechanism), ]
utils::write.csv(agg, file.path(OUT, "response_type_summary.csv"), row.names = FALSE)
cat("\n=== Per-response-type recovery (mean over replicates) ===\n"); print(agg, row.names = FALSE)
cat(sprintf("\n%d ok / %d total. Wallclock: %.1f min.\n", nrow(ok), nrow(res),
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
if (any(res$status != "ok")) { cat("Failures:\n"); print(table(res$status[res$status != "ok"])) }
