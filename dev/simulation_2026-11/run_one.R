# =================================================================================================
# run_one.R -- ONE scenario, ONE seed, end to end. This is the unit of work.
#
#   Rscript run_one.R --id n100_ps0.7_typical_clade --seed 1 --out results/
#
# Everything the campaign does is this script, run 54,000 times with different arguments. Run it
# once by hand and read the output to see exactly what each dataset produces.
#
# Resumable: if the output file already exists, it stops. So the campaign can be interrupted and
# relaunched freely, and replicates can be added later without redoing earlier ones.
# =================================================================================================

suppressPackageStartupMessages({ library(ape); library(nlme); library(phytools) })

## ---- which BACE? Checked before anything else. ------------------------------------------------
## A campaign measures one version of BACE. An older build sitting in the R library produces
## results that look entirely normal and mean nothing -- that happened on 2026-11-09, with the
## repository at 1.0.0 and the library still at 0.1.0.9000.
BACE_MIN_VERSION <- "1.0.0"
bace_version <- as.character(utils::packageVersion("BACE"))
if (utils::compareVersion(bace_version, BACE_MIN_VERSION) < 0)
  stop("BACE ", bace_version, " is installed; this campaign needs >= ", BACE_MIN_VERSION,
       ". Run  R CMD INSTALL .  from the BACE repo root.", call. = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

## ---- arguments --------------------------------------------------------------------------------
a <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) { i <- match(paste0("--", flag), a)
                                        if (is.na(i)) default else a[i + 1L] }
id      <- arg("id",   stop("--id is required"))
seed    <- as.integer(arg("seed", 1L))
out_dir <- arg("out",  "results")
M       <- as.integer(arg("M", 20L))
nitt    <- as.integer(arg("nitt", 50000L))
burnin  <- as.integer(arg("burnin", 10000L))
thin    <- as.integer(arg("thin", 25L))
runs    <- as.integer(arg("runs", 5L))
## The final imputation runs may have their own, usually longer, chains. NULL = reuse the
## chained-phase settings. The two real case studies both split them: amphibian_ctmax runs
## nitt 15000 chained with final_nitt 105000, bird_social nitt 50000 with final_nitt 500000.
int_or_null <- function(x) if (is.null(x)) NULL else as.integer(x)
final_nitt   <- int_or_null(arg("final_nitt"))
final_burnin <- int_or_null(arg("final_burnin"))
final_thin   <- int_or_null(arg("final_thin"))
## Checkpointing is ON by default so an interrupted dataset resumes instead of restarting; it is
## deleted as soon as the dataset succeeds. --no_checkpoint turns it off (e.g. on a disk too small
## to hold one per concurrent job: budget a few hundred MB per core).
use_ckpt <- !("--no_checkpoint" %in% a)

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
for (f in sprintf("R/0%d_%s.R", 1:4, c("scenarios", "simulate", "fit", "score")))
  source(file.path(here, f))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, paste0(id, "_s", seed, ".rds"))
if (file.exists(out_file)) { message("already done: ", basename(out_file)); quit(save = "no") }

## ---- find the scenario ------------------------------------------------------------------------
scenarios <- make_scenarios()
scenario  <- scenarios[scenarios$id == id, ]
if (nrow(scenario) != 1L) stop("no such scenario: ", id, call. = FALSE)

message(id, " seed ", seed)

## ---- 0. a dataset may fail, and that must not kill the job --------------------------------
## A campaign of 54,000 datasets will hit failures, and two are known in BACE 1.0.0, both
## reproduced on this grid:
##   sim_bace()  "NA in probability vector" -- the softmax in mnom_liab2cat() overflows when the
##               liabilities are large (exp(710) is Inf, and Inf/Inf is NaN). Hit 3 of 54
##               scenarios, all at n = 500 with high signal and related traits.   [issue #19]
##   bace()      "zero-variance series" -- coda::geweke.diag() is called unguarded at
##               R/model_functions.R:1635, where the same call is wrapped in tryCatch at
##               R/phylo_signal_summary.R:850. Hit 2 of 54 scenarios.             [issue #20]
##
## A crashed job writes nothing, so it leaves a silent gap and the failure rate is invisible.
## Instead every stage is wrapped: a failure is RECORDED and the dataset still writes its file,
## so the failure rate is a number we report rather than a hole nobody notices.
##
## Note these are DETERMINISTIC failures: the same scenario and seed fails the same way every
## time, so re-running or resuming cannot help. Only fixing them in BACE will.
fail <- function(stage, e) {
  message("FAILED at ", stage, ": ", conditionMessage(e))
  saveRDS(list(id = id, seed = seed, scenario = scenario, failed_at = stage,
               error = conditionMessage(e), bace_version = bace_version, time = Sys.time()),
          out_file)
  quit(save = "no", status = 0)
}

## ---- 1. simulate ------------------------------------------------------------------------------
dataset <- tryCatch(simulate_dataset(scenario, seed), error = function(e) fail("simulate", e))

## ---- 2. run the three arms --------------------------------------------------------------------
## The true slope of x2 on x1 is the regression slope the data were generated with. The oracle arm
## is the check on that: if the oracle does not recover it with about 95% coverage, the assumed
## truth is wrong and every other arm's coverage is meaningless.
truth_slope <- scenario$trait_beta

oracle <- score_slope_single(dataset$truth, dataset$tree)
cc     <- fit_complete_case(dataset)
cc_slope <- if (cc$status == "scored") score_slope_single(cc$data, dataset$tree) else NULL

ckpt <- if (use_ckpt) file.path(out_dir, "checkpoints", paste0(id, "_s", seed)) else NULL
bace <- tryCatch(fit_bace(dataset, n_final = M, nitt = nitt, burnin = burnin,
                          thin = thin, runs = runs, seed = seed, checkpoint_dir = ckpt,
                          final_nitt = final_nitt, final_burnin = final_burnin,
                          final_thin = final_thin),
                 error = function(e) fail("bace", e))
bace_slope <- tryCatch(score_slope_bace(bace$imputed, dataset$tree), error = function(e) NULL)

## ---- 3. score ---------------------------------------------------------------------------------
recovery <- tryCatch(score_recovery(dataset$truth, dataset$mask, bace$imputed),
                     error = function(e) NULL)
coverage <- tryCatch(score_cell_coverage(bace$intervals, dataset$truth), error = function(e) NULL)

## One tidy row per arm for the slope, with coverage of the TRUE value.
slope_row <- function(arm, s) {
  if (is.null(s)) return(NULL)
  data.frame(arm = arm, s, truth = truth_slope,
             bias = s$estimate - truth_slope,
             covered = s$lower <= truth_slope & truth_slope <= s$upper)
}
slopes <- rbind(slope_row("oracle", oracle),
                slope_row("complete_case", cc_slope),
                slope_row("bace", bace_slope))

## ---- 4. write a SMALL result ------------------------------------------------------------------
## Never the fit: a fitted bace object is about 35 MB, and 54,000 of them is not storable.
saveRDS(list(
  id = id, seed = seed, scenario = scenario, failed_at = NA_character_,
  slopes = slopes, recovery = recovery, cell_coverage = coverage,
  realised = dataset$realised, realised_frac = dataset$realised_frac,
  complete_case = cc[c("n_complete", "frac_kept", "status")],
  bace_converged = bace$converged, bace_conv_detail = bace$conv_detail,
  settings = c(M = M, nitt = nitt, burnin = burnin, thin = thin, runs = runs,
               final_nitt = final_nitt %||% nitt, final_thin = final_thin %||% thin),
  wall_s = bace$wall_s,
  bace_version = bace_version, time = Sys.time()
), out_file)

message("wrote ", basename(out_file), "  (", round(bace$wall_s), "s)")
print(slopes, digits = 3)
print(recovery, digits = 3)
