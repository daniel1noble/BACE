# =============================================================================
# 13_response_type_simulation.R
#
# Comprehensive simulation study: how well does BACE recover MISSING VALUES and
# the RESPONSE's marginal structure across response distributions (gaussian,
# count/poisson, binary, categorical/multinomial, ordinal/threshold), under
# MCAR and MAR, over many replicates?
#
# Data-generating model (per replicate), via the package engine sim_bace():
#   tree ~ birth-death; two fully-observed gaussian predictors x1, x2 with
#   phylogenetic signal; a response y of the target type generated from
#   x1, x2 (auto-generated effect sizes) + a phylogenetic random effect
#   (Pagel lambda = 0.7 on all variables). Only y is made missing.
#
# Missingness on y:
#   MCAR : P(miss) = rate (constant)
#   MAR  : P(miss) = plogis(a + s * z(x1))   [depends on OBSERVED x1 -> MAR]
#
# BACE model fit (per replicate):
#   bace_imp(y ~ x1 + x2, ~1|Species)  ->  bace_final_imp(n_final)
#   (chained equations; x1, x2 complete so only y is imputed)
#
# Metrics per (type, mechanism), aggregated across replicates:
#   value recovery (type-appropriate, via evaluate_imputation_ensemble):
#     gaussian/count : correlation, nrmse, coverage95 (of 95% PI on hidden cells)
#     binary/categ.  : accuracy, balanced_accuracy, brier
#     ordinal        : accuracy, brier, ordinal_mae
#   marginal recovery (the MAR signature):
#     marg_bias_bace : E[imputed y summary] - true summary
#     marg_bias_cc   : E[complete-case y summary] - true summary
#     (summary = mean for gaussian/count, P(level2) for binary, mean-rank for
#      ordinal; categorical uses accuracy only.)
#   timing: bace wall-clock seconds per replicate.
#
# Scale via env vars (defaults are a moderate run):
#   SIMTYPE_REPS SIMTYPE_NSPP SIMTYPE_NITT SIMTYPE_BURNIN SIMTYPE_THIN
#   SIMTYPE_RUNS SIMTYPE_NFINAL SIMTYPE_CORES
# =============================================================================

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape); library(MASS); library(parallel)
})

# Reuse the shared cell-metric evaluator (single source of truth with dev/10/11).
source_funcs <- function(path, fn_names) {
  src <- readLines(path)
  for (fn in fn_names) {
    s <- grep(paste0("^", fn, " <- function"), src)[1]
    if (is.na(s)) stop("fn not found: ", fn)
    depth <- 0; e <- NA
    for (i in s:length(src)) {
      depth <- depth + sum(gregexpr("\\{", src[i])[[1]] > 0) -
        sum(gregexpr("\\}", src[i])[[1]] > 0)
      if (i > s && depth == 0) { e <- i; break }
    }
    eval(parse(text = src[s:e]), envir = globalenv())
  }
}
source_funcs(file.path("dev", "02_benchmark_simulated_full.R"),
             c("evaluate_imputation_ensemble"))

.envint <- function(nm, d) {
  v <- suppressWarnings(as.integer(Sys.getenv(nm, unset = NA_character_)))
  if (is.na(v)) d else v
}
CFG <- list(
  reps    = .envint("SIMTYPE_REPS",   25L),
  nspp    = .envint("SIMTYPE_NSPP",   80L),
  nitt    = .envint("SIMTYPE_NITT",   5000L),
  burnin  = .envint("SIMTYPE_BURNIN", 1000L),
  thin    = .envint("SIMTYPE_THIN",   5L),
  runs    = .envint("SIMTYPE_RUNS",   4L),
  n_final = .envint("SIMTYPE_NFINAL", 20L),
  cores   = .envint("SIMTYPE_CORES",  6L),
  rate = 0.30, mar_slope = 2.0, phylo = 0.7
)

TYPES <- list(
  gaussian    = list(resp = "gaussian",     eval = "gaussian"),
  poisson     = list(resp = "poisson",      eval = "count"),
  binary      = list(resp = "binary",       eval = "categorical"),
  categorical = list(resp = "multinomial3", eval = "categorical"),
  ordinal     = list(resp = "threshold3",   eval = "ordered")
)
MECHS <- c("MCAR", "MAR")

# ---- helpers ----------------------------------------------------------------
post_process <- function(d, type) {
  if (type == "poisson")     d$y <- as.integer(round(d$y))
  if (type == "binary")      d$y <- factor(d$y, levels = sort(unique(d$y)))
  if (type == "categorical") d$y <- factor(as.character(d$y))
  if (type == "ordinal")     d$y <- factor(d$y,
                                           levels = sort(unique(d$y)), ordered = TRUE)
  d
}

simulate_one <- function(type, cfg, seed) {
  set.seed(seed)
  sim <- suppressMessages(sim_bace(
    response_type   = TYPES[[type]]$resp,
    predictor_types = c("gaussian", "gaussian"),
    var_names       = c("y", "x1", "x2"),
    phylo_signal    = rep(cfg$phylo, 3),
    n_cases         = cfg$nspp, n_species = cfg$nspp,
    beta_sparsity   = 0.0,                 # predictors always affect the response
    missingness     = c(0, 0, 0)))
  d <- sim$complete_data
  names(d)[names(d) == "species"] <- "Species"
  list(data = post_process(d, type), tree = sim$tree)
}

apply_missing <- function(d, mechanism, cfg, seed) {
  set.seed(seed); n <- nrow(d)
  p <- if (mechanism == "MCAR") rep(cfg$rate, n) else {
    lin <- cfg$mar_slope * as.numeric(scale(d$x1))
    a <- stats::uniroot(function(a) mean(plogis(a + lin)) - cfg$rate,
                        c(-12, 12))$root
    plogis(a + lin)
  }
  miss <- stats::rbinom(n, 1, p) == 1
  # Guard: discrete responses need >= 2 observed per level for MCMCglmm.
  if (is.factor(d$y)) {
    for (lv in levels(d$y)) {
      idx <- which(d$y == lv)
      obs <- idx[!miss[idx]]
      if (length(obs) < 2L && length(idx) >= 2L) {
        miss[sample(idx, 2L)] <- FALSE
      }
    }
  }
  if (all(miss)) miss[sample(n, 1L)] <- FALSE
  miss
}

# numeric summary of the response for marginal-bias tracking
resp_summary <- function(y, type) {
  if (type %in% c("gaussian", "poisson")) return(mean(as.numeric(y)))
  if (type == "binary") return(mean(as.numeric(y) == max(as.numeric(y))))
  if (type == "ordinal") return(mean(as.integer(y)))
  NA_real_                                # categorical: no scalar summary
}

# ---- one replicate ----------------------------------------------------------
run_rep <- function(type, mechanism, rep_id, cfg) {
  seed <- 5000L * rep_id + switch(mechanism, MCAR = 1L, MAR = 2L) +
    100000L * match(type, names(TYPES))
  sim  <- simulate_one(type, cfg, seed)
  d    <- sim$data; tree <- sim$tree
  miss <- apply_missing(d, mechanism, cfg, seed + 7L)

  dm <- d; dm$y[miss] <- NA
  miss_mask <- data.frame(y = miss, x1 = FALSE, x2 = FALSE)
  var_types <- c(y = TYPES[[type]]$eval, x1 = "gaussian", x2 = "gaussian")

  t0 <- Sys.time()
  init <- suppressWarnings(suppressMessages(bace_imp(
    fixformula = "y ~ x1 + x2", ran_phylo_form = "~ 1 | Species",
    phylo = tree, data = dm, runs = cfg$runs,
    nitt = cfg$nitt, thin = cfg$thin, burnin = cfg$burnin, verbose = FALSE)))
  fin <- suppressWarnings(suppressMessages(bace_final_imp(
    bace_object = init, fixformula = "y ~ x1 + x2",
    ran_phylo_form = "~ 1 | Species", phylo = tree, n_final = cfg$n_final,
    nitt = cfg$nitt, thin = cfg$thin, burnin = cfg$burnin, verbose = FALSE)))
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # value recovery on hidden y cells
  ev <- evaluate_imputation_ensemble(true_data = d, imp_list = fin$all_datasets,
                                     miss_mask = miss_mask, var_types = var_types)
  ev <- ev[ev$variable == "y", ]
  getm <- function(m) { v <- ev$value[ev$metric == m]; if (length(v)) v[1] else NA_real_ }

  # marginal recovery
  s_true <- resp_summary(d$y, type)
  s_cc   <- resp_summary(d$y[!miss], type)
  s_bace <- mean(vapply(fin$all_datasets,
                        function(z) resp_summary(z$y, type), numeric(1)))

  data.frame(
    type = type, mechanism = mechanism, rep_id = rep_id, status = "ok",
    runtime_sec       = runtime,
    correlation       = getm("correlation"),
    nrmse             = getm("nrmse"),
    coverage95        = getm("coverage95"),
    accuracy          = getm("accuracy"),
    balanced_accuracy = getm("balanced_accuracy"),
    brier             = getm("brier"),
    ordinal_mae       = getm("ordinal_mae"),
    marg_bias_bace    = s_bace - s_true,
    marg_bias_cc      = s_cc   - s_true,
    n_missing         = sum(miss),
    stringsAsFactors  = FALSE, row.names = NULL)
}

# ---- main -------------------------------------------------------------------
worklist <- expand.grid(type = names(TYPES), mechanism = MECHS,
                        rep_id = seq_len(CFG$reps), stringsAsFactors = FALSE)
cat(sprintf("Response-type sim: %d types x %d mechs x %d reps = %d runs | nspp=%d nitt=%d n_final=%d | cores=%d\n",
            length(TYPES), length(MECHS), CFG$reps, nrow(worklist),
            CFG$nspp, CFG$nitt, CFG$n_final, CFG$cores))
cat("Started:", format(Sys.time()), "\n")
t_all <- Sys.time()

rows <- mclapply(seq_len(nrow(worklist)), function(i) {
  w <- worklist[i, ]
  out <- tryCatch(run_rep(w$type, w$mechanism, w$rep_id, CFG),
    error = function(e) data.frame(type = w$type, mechanism = w$mechanism,
      rep_id = w$rep_id, status = paste("error:", conditionMessage(e)),
      runtime_sec = NA, correlation = NA, nrmse = NA, coverage95 = NA,
      accuracy = NA, balanced_accuracy = NA, brier = NA, ordinal_mae = NA,
      marg_bias_bace = NA, marg_bias_cc = NA, n_missing = NA,
      stringsAsFactors = FALSE))
  rt <- if (is.na(out$runtime_sec)) 0 else out$runtime_sec
  cat(sprintf("[%s] %-11s %-4s rep %02d -> %s (%.0fs)\n",
              format(Sys.time(), "%H:%M:%S"), w$type, w$mechanism, w$rep_id,
              substr(out$status, 1, 20), rt))
  out
}, mc.cores = CFG$cores)

res <- do.call(rbind, rows)
ok  <- res[res$status == "ok", ]

se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
agg <- do.call(rbind, lapply(split(ok, list(ok$type, ok$mechanism), drop = TRUE),
  function(g) data.frame(
    type = g$type[1], mechanism = g$mechanism[1], n = nrow(g),
    runtime_s   = round(mean(g$runtime_sec, na.rm = TRUE), 1),
    correlation = round(mean(g$correlation, na.rm = TRUE), 3),
    nrmse       = round(mean(g$nrmse, na.rm = TRUE), 3),
    coverage95  = round(mean(g$coverage95, na.rm = TRUE), 3),
    accuracy    = round(mean(g$accuracy, na.rm = TRUE), 3),
    bal_acc     = round(mean(g$balanced_accuracy, na.rm = TRUE), 3),
    brier       = round(mean(g$brier, na.rm = TRUE), 3),
    ord_mae     = round(mean(g$ordinal_mae, na.rm = TRUE), 3),
    marg_bias_bace = round(mean(g$marg_bias_bace, na.rm = TRUE), 3),
    marg_bias_cc   = round(mean(g$marg_bias_cc, na.rm = TRUE), 3),
    row.names = NULL)))
type_order <- names(TYPES)
agg <- agg[order(match(agg$type, type_order), agg$mechanism), ]

cat("\n=== Per-response-type recovery (mean over replicates) ===\n")
print(agg, row.names = FALSE)
cat(sprintf("\nRuns: %d ok / %d total. Wallclock: %.1f min.\n",
            nrow(ok), nrow(res),
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
if (any(res$status != "ok")) {
  cat("Failures:\n"); print(table(res$status[res$status != "ok"]))
}

out_dir <- file.path("dev", "simulation_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
saveRDS(list(per_rep = res, summary = agg, cfg = CFG),
        file.path(out_dir, "response_type_simulation.rds"))
utils::write.csv(agg, file.path(out_dir, "response_type_summary.csv"),
                 row.names = FALSE)
cat("Saved:", file.path(out_dir, "response_type_simulation.rds"), "\n")
