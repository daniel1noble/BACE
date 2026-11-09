# =============================================================================
# RETIRED 2026-11-09, part-run (17 of 40). This pilot asked which chain
# settings make BACE's convergence flag pass. That flag turned out to be
# defective (it could never pass for mean-centred traits; replaced in BACE
# 1.0.1 - see .agents/investigation-convergence-2026-11.md), so the question
# dissolves. What the 17 finished runs DID show, and what retired it:
#   - coverage 17/17 and bias ~0 at EVERY setting tested
#   - runs 5 -> 15 and final_nitt 50k -> 150k changed nothing but cost
# So the campaign settings (runs 5, nitt 50k) already deliver nominal
# coverage; longer chains buy nothing measurable. Kept for the record.
# =============================================================================
# =================================================================================================
# pilot_settings.R -- how long should the chains be?
#
#   Rscript pilot_settings.R --out pilot --reps 8
#
# WHY THIS EXISTS. The campaign is three weeks of compute, and its MCMC settings were inherited
# from the 2026-09/10 campaign without ever being checked. Two measurements say they may be light:
#
#   - the chained loop converged in only 28% of n = 100 datasets (and 51% at n = 300, 82% at
#     n = 1000), so at the small sizes this design is built around, it mostly does not settle;
#   - 23% of the fixed effects in the FINAL phase -- the models that produce the published
#     numbers -- had an effective sample size below 100.
#
# The two real case studies in the package both use far heavier settings, and both split the
# phases, which is what `final_nitt` was added for:
#
#   2026-09/10 campaign   chained nitt 50000, runs  5    final: same chain,      n_final 20
#   amphibian_ctmax       chained nitt 15000, runs 10    final_nitt 105000,      n_final 50
#   bird_social           chained nitt 50000, runs 25    final_nitt 500000,      n_final 50
#
# We cannot simply adopt those. Cost goes as runs x nitt for the chained phase and
# n_final x final_nitt for the final one, so bird_social's settings are roughly 21x the
# campaign's -- 1.5 million core-hours. The question is where between the two the honest answer
# sits, and this pilot answers it with about 20 core-hours against a 72,000 core-hour campaign.
#
# WHAT IT MEASURES, per setting:
#   converged      how often the chained loop settles          -- fixes the 28% problem or not
#   ess_min        worst effective sample size in the final models  -- fixes the mixing problem
#   covered        does the 95% interval catch the true slope  -- the thing that actually matters
#   wall_s         what it costs
#
# HOW TO READ IT: take the CHEAPEST setting whose coverage is ~0.95 and whose convergence and ESS
# are respectable. More chain than that is wasted money; less is indefensible.
# =================================================================================================

suppressPackageStartupMessages({ library(ape); library(nlme); library(phytools) })

a <- commandArgs(trailingOnly = TRUE)
arg <- function(f, d = NULL) { i <- match(paste0("--", f), a); if (is.na(i)) d else a[i + 1L] }
out_dir <- arg("out", "pilot")
reps    <- as.integer(arg("reps", 8L))
id      <- arg("id", "n100_ps0.7_typical_mcar")
## --only <n> runs just the nth setting, so the five can go on five cores at once:
##   for k in 1 2 3 4 5; do Rscript pilot_settings.R --only $k --out pilot & done
## Each (setting, seed) writes its own file and is skipped if it exists, so this is resumable and
## the copies never duplicate each other's work.
only    <- if (is.null(arg("only"))) NULL else as.integer(arg("only"))

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
for (f in sprintf("R/0%d_%s.R", 1:4, c("scenarios","simulate","fit","score"))) source(file.path(here, f))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- the settings to compare -------------------------------------------------------------------
## `runs` drives whether the chained loop settles. `final_nitt` drives how well the final models
## mix. Each is varied on its own so we can see which one buys what, rather than changing both
## and not knowing which helped.
settings <- rbind(
  data.frame(label = "campaign (current)",  runs =  5, final_nitt =  50000),
  data.frame(label = "more chained runs",   runs = 10, final_nitt =  50000),
  data.frame(label = "many chained runs",   runs = 15, final_nitt =  50000),
  data.frame(label = "longer final chain",  runs =  5, final_nitt = 150000),
  data.frame(label = "both",                runs = 10, final_nitt = 150000)
)

scenario <- make_scenarios()[make_scenarios()$id == id, ]
stopifnot(nrow(scenario) == 1L)
message("pilot on ", id, ", ", reps, " datasets x ", nrow(settings), " settings")

## ---- run -----------------------------------------------------------------------------------
rows <- list()
for (k in if (is.null(only)) seq_len(nrow(settings)) else only) {
  st <- settings[k, ]
  for (s in seq_len(reps)) {
    f <- file.path(out_dir, sprintf("%s_runs%d_fn%d_s%d.rds", id, st$runs, st$final_nitt, s))
    if (file.exists(f)) { rows[[length(rows)+1L]] <- readRDS(f); next }   # resumable

    d <- simulate_dataset(scenario, s)
    ## final_thin is scaled with the chain, so a longer chain stores the same number of draws
    ## rather than 3x as many -- every retained draw is kept for all n_final runs, so not scaling
    ## it would multiply memory as well as time.
    ft <- as.integer(25 * st$final_nitt / 50000)
    fit <- tryCatch(fit_bace(d, n_final = 20L, nitt = 50000L, burnin = 10000L, thin = 25L,
                             runs = st$runs, seed = s,
                             final_nitt = st$final_nitt, final_thin = ft),
                    error = function(e) NULL)
    if (is.null(fit)) next

    sl <- score_slope_bace(fit$imputed, d$tree)
    row <- data.frame(label = st$label, runs = st$runs, final_nitt = st$final_nitt, seed = s,
                      converged = fit$converged, wall_s = fit$wall_s,
                      covered = if (is.null(sl)) NA else
                        sl$lower <= scenario$trait_beta & scenario$trait_beta <= sl$upper,
                      bias = if (is.null(sl)) NA_real_ else sl$estimate - scenario$trait_beta)
    saveRDS(row, f); rows[[length(rows)+1L]] <- row
    message(sprintf("  %-20s seed %d: %.0f min, converged=%s", st$label, s, fit$wall_s/60, fit$converged))
  }
}

## ---- report --------------------------------------------------------------------------------
## With --only, read every setting's files back so the summary covers the whole pilot however it
## was split up.
if (!is.null(only)) rows <- lapply(list.files(out_dir, pattern = "_runs.*[.]rds$", full.names = TRUE), readRDS)
res <- do.call(rbind, rows)
if (is.null(res) || !nrow(res)) { message("no results yet"); quit(save = "no") }
utils::write.csv(res, file.path(out_dir, "pilot_results.csv"), row.names = FALSE)

summ <- stats::aggregate(cbind(converged, covered, bias, wall_s) ~ label + runs + final_nitt,
                         res, function(z) mean(z, na.rm = TRUE))
summ$cost_vs_current <- round(summ$wall_s / summ$wall_s[summ$runs == 5 & summ$final_nitt == 50000], 2)
print(summ, digits = 3)
message("\nTake the cheapest row whose coverage is ~0.95 and whose convergence is respectable.")
