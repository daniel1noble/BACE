# code/coverage_experiment.R
#
# Why do the discrete 95% category sets undercover at n = 100, lambda = 0.3
# (ord 0.886-0.909, cat3 0.906-0.914), when binary and all continuous traits are
# nominal in the same runs? (Dan, 2026-11-05.) Ruled out already:
#   - the set rule itself: a perfectly calibrated imputer at M = 20 scores >= 0.95
#     under the campaign's exact rule (scratch set_rule_calibration.R demo);
#   - generic chain length: the well-covered lambda = 1 cells mix WORSE
#     (ess_med 224 vs 730) than the undercovering cells.
# Two remaining hypotheses, one arm each, paired to the campaign datasets:
#   multi   ovr_categorical = FALSE (true multinomial, the new package default),
#           campaign MCMC settings. Tests the one-vs-rest construction for cat3.
#   long3x  ovr_categorical = TRUE, 3x chain (nitt 150k, burnin 30k, thin 75;
#           same 1,600 retained samples). Tests cutpoint/threshold mixing for ord.
#
# Cells: types_mixed, n = 100, lambda = 0.3, rho in {0, 0.5}, seeds 1-25 -- the
# SAME datasets as the campaign (make_cell rebuild + truth tripwire), so each
# variant is compared to the campaign's stored per-dataset bace scores, paired
# by seed. Discrete scoring only (score_discrete; continuous was already nominal).
#
#   Rscript code/coverage_experiment.R [RESULTS_DIR] [OUT_DIR] [CORES]
#
# Resumable: one rds per dataset x variant in OUT_DIR/coverage_exp_cache/.
# Output: coverage_exp.csv (per dataset x variant x trait) and
# coverage_exp_summary.csv (per cell x variant x trait, with paired diffs).

suppressMessages({ library(MCMCglmm); library(ape); library(BACE) })
RNGkind("L'Ecuyer-CMRG")   # required to reproduce make_cell()'s datasets exactly

a <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(a) >= 1) a[1] else "results_2026-10"
out_dir <- if (length(a) >= 2) a[2] else res_dir
cores   <- if (length(a) >= 3) as.integer(a[3]) else 7L
cache   <- file.path(out_dir, "coverage_exp_cache"); dir.create(cache, showWarnings = FALSE, recursive = TRUE)

source(file.path("code", "campaign_gnn_off_lib.R"))
source(file.path("code", "rubin_lib.R"))
source(file.path("code", "rubin_bace.R"))
source(file.path("code", "rubin_discrete.R"))

SEEDS <- 1:25; RHOS <- c("r0", "r0.5")
if (nzchar(Sys.getenv("SMOKE"))) { SEEDS <- 1L; RHOS <- "r0.5" }   # one dataset, multi arm only
VARIANTS <- list(
  multi  = list(ovr = FALSE, nitt = 50000,  burnin = 10000, thin = 25, voff = 1L),
  long3x = list(ovr = TRUE,  nitt = 150000, burnin = 30000, thin = 75, voff = 2L))
if (nzchar(Sys.getenv("SMOKE"))) VARIANTS <- VARIANTS["multi"]

fs <- unlist(lapply(RHOS, function(r) file.path(res_dir,
  sprintf("rubin_types_mixed_BM_l0.3_%s_mcar0.3_n100_M20_s%d.rds", r, SEEDS))))
fs <- fs[file.exists(fs)]
grid <- expand.grid(f = fs, variant = names(VARIANTS), stringsAsFactors = FALSE)
cat(sprintf("datasets: %d x variants: %d = %d runs\n", length(fs), length(VARIANTS), nrow(grid)))

one <- function(f, vn) {
  v <- VARIANTS[[vn]]
  key <- paste0(sub("[.]rds$", "", basename(f)), "_", vn)
  cf <- file.path(cache, paste0(key, ".rds"))
  if (file.exists(cf)) return(readRDS(cf))
  x <- readRDS(f)
  cell <- make_cell("types_mixed", x$n, x$seed, miss_frac = x$miss_frac, miss = x$miss,
                    lambda = x$lambda, rho = x$rho, thresholds = x$thresholds, driver = x$driver)
  stopifnot(isTRUE(all.equal(cell$truth$c1, x$truth$c1)),
            identical(rownames(cell$truth), rownames(x$truth)),
            all(is.na(cell$df_miss[["ord"]]) == x$mask[, "ord"]))
  set.seed(x$seed * 100L + v$voff)
  fb <- fit_bace_mi(cell, M = 20L, nitt = v$nitt, burnin = v$burnin, thin = v$thin,
                    runs = 5L, ovr_categorical = v$ovr)
  sets <- mi_bace_shipped(fb$outb, cell$df_miss)
  fd <- fill_degenerate(sets, cell$df_miss, disc_traits_of(x$truth))
  dc <- score_discrete(vn, fd$sets, x$truth, x$mask)
  dc$n <- x$n; dc$lambda <- x$lambda; dc$rho <- x$rho; dc$seed <- x$seed
  dc$variant <- vn; dc$wall_s <- fb$wall_s; dc$ess_med <- fb$ess_med
  # campaign baseline for the same dataset, for the paired comparison
  # ($disc_cells is one data.frame with an `arm` column, not a list by arm)
  base <- x$disc_cells[x$disc_cells$arm == "bace", ]
  dc$base_set_coverage <- base$set_coverage[match(dc$trait, base$trait)]
  dc$base_accuracy    <- base$accuracy[match(dc$trait, base$trait)]
  saveRDS(dc, cf); dc
}

rows <- parallel::mclapply(seq_len(nrow(grid)), function(i)
  tryCatch(one(grid$f[i], grid$variant[i]), error = function(e) {
    cat("ERROR", basename(grid$f[i]), grid$variant[i], ":", conditionMessage(e), "\n"); NULL
  }), mc.cores = cores)
d <- do.call(rbind, Filter(Negate(is.null), rows)); rownames(d) <- NULL
write.csv(d, file.path(out_dir, "coverage_exp.csv"), row.names = FALSE)

se <- function(v) sd(v) / sqrt(length(v))
summ <- do.call(rbind, lapply(split(d, d[c("variant", "rho", "trait")], drop = TRUE), function(g) {
  dd <- g$set_coverage - g$base_set_coverage
  data.frame(variant = g$variant[1], rho = g$rho[1], trait = g$trait[1], datasets = nrow(g),
             set_coverage = mean(g$set_coverage), base_set_coverage = mean(g$base_set_coverage),
             paired_diff = mean(dd), paired_diff_se = se(dd),
             accuracy = mean(g$accuracy), base_accuracy = mean(g$base_accuracy),
             mean_wall_min = mean(g$wall_s) / 60, row.names = NULL)
}))
write.csv(summ, file.path(out_dir, "coverage_exp_summary.csv"), row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)
