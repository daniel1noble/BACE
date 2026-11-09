#!/bin/bash
# =================================================================================================
# rubin_campaign_local.sh -- build the task list for the 2026-11 campaign and run it.
#
# WHAT THIS SCRIPT DOES: writes one line per dataset into a plain text file, then feeds that file
# to xargs, which runs CORES of them at a time. Nothing clever happens here. If you want to know
# exactly what the campaign will run, run this script and read $OUT/tasks.txt -- every line is a
# complete, copy-pasteable rubin_cell.R command you can run on its own.
#
# PLAN: .agents/plans/simulation-plan-2026-11.md  (the design and the reasoning behind it)
#
# THE DESIGN, in words:
#
#   2 tree sizes          n = 100 and 500        the sizes comparative biologists actually have
#   x 3 signal levels     lambda = 0.3, 0.7, 0.9 weak, moderate and strong phylogenetic signal
#   x 3 information       floor / typical / good covariate   -- see below, this is NOT a factorial
#   x 3 gap mechanisms    mcar, mar, clade       random gaps, driven gaps, whole clades missing
#   x REPS datasets each
#
#   = 2 x 3 x 3 x 3 = 54 scenarios, x REPS datasets.
#
# THE THREE INFORMATION REGIMES. This is one axis with three levels, not two crossed axes. Each
# level says how much a missing trait can borrow from somewhere other than the phylogeny:
#
#   floor           traits uncorrelated (rho 0), no covariate   -- nothing to borrow but the tree
#   typical         traits correlated 0.5, no covariate         -- the ordinary situation
#   good covariate  traits correlated 0.5, PLUS one fully       -- you have a body-size-like
#                   observed variable correlated 0.8 with c1       variable measured for everything
#
# Claim 2 of the manuscript is read off the last two: at lambda = 0.3, where the phylogeny cannot
# help, how much does the good covariate buy you?
#
# RESUMABLE, ALWAYS. rubin_cell.R skips any dataset whose .rds already exists. So you can stop
# this script at any moment, rerun it, and it picks up where it left off; and you can RAISE the
# replicate count later and only the new seeds will run.
#
# USAGE
#   CORES=7 bash code/rubin_campaign_local.sh              # the full design
#   CORES=7 REPS=2 bash code/rubin_campaign_local.sh       # a tiny rehearsal of the same thing
#   CORES=7 N_LIST=100 bash code/rubin_campaign_local.sh   # just the small trees
#   DRY_RUN=1 bash code/rubin_campaign_local.sh            # write tasks.txt and stop
# =================================================================================================
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"

# ---- knobs, all overridable from the environment ------------------------------------------------
CORES="${CORES:-7}"                       # how many datasets to run at once
OUT="${OUT:-$HERE/../results_2026-11}"    # where the .rds files go
REPS="${REPS:-1000}"                      # datasets per scenario
N_LIST="${N_LIST:-100 500}"               # tree sizes
LAMBDAS="${LAMBDAS:-0.3 0.7 0.9}"         # phylogenetic signal
MECHS="${MECHS:-mcar mar clade}"          # how the gaps are made
ARMS="${ARMS:-bace,complete_case}"        # BACE, and the drop-the-gappy-species competitor
DRIVER_AS_PRED="${DRIVER_AS_PRED:-0}"     # 1 = do not waste a model on the fully observed driver

# MCMC settings, identical to the 2026-10 campaign so cost comparisons stay meaningful
NITT="${NITT:-50000}"; BURNIN="${BURNIN:-10000}"; THIN="${THIN:-25}"; RUNS="${RUNS:-5}"; M="${M:-20}"

mkdir -p "$OUT"

# ---- refuse to launch against the wrong BACE ----------------------------------------------------
# rubin_cell.R checks this too, per dataset. It is checked here as well so a wrong version stops
# the campaign in one second rather than 54,000 times over.
BACE_MIN_VERSION="${BACE_MIN_VERSION:-1.0.0}"

# NB: the check is its OWN command, not part of a pipeline. In a pipeline the exit status is the
# LAST command's (here `tee`), so `Rscript ... | tee ... || exit 1` would never fire -- the guard
# would print its error and the campaign would launch anyway. Capture first, then decide.
VERSIONS=$(Rscript -e "v <- as.character(utils::packageVersion('BACE'))
  if (utils::compareVersion(v, '$BACE_MIN_VERSION') < 0)
    stop('installed BACE ', v, ' is older than the required $BACE_MIN_VERSION -- run  R CMD INSTALL .  from the BACE repo root', call. = FALSE)
  cat('BACE', v, '| MCMCglmm', as.character(utils::packageVersion('MCMCglmm')), '| R', as.character(getRversion()), '\n')" 2>&1)
if [ $? -ne 0 ]; then
  echo "$VERSIONS" >&2
  echo "REFUSING TO LAUNCH: wrong BACE version." >&2
  exit 1
fi
echo "$VERSIONS" | tee "$OUT/versions.txt"

# ---- build the task list ------------------------------------------------------------------------
TASKS="$OUT/tasks.txt"; : > "$TASKS"

DRIVER_FLAG=""; [ "$DRIVER_AS_PRED" = "1" ] && DRIVER_FLAG="--driver_as_predictor"

for n in $N_LIST; do

  # Results are split into one directory per tree size. 54,000 files in a single flat folder is
  # slow to list and awkward for the aggregation scripts; this costs nothing and avoids that.
  N_OUT="$OUT/n$n"
  mkdir -p "$N_OUT"

  for lambda in $LAMBDAS; do

    # The three information regimes, one per line so you can read each one whole. A regime is
    # fully identified downstream by (rho, n_cov), both of which appear in the result filename,
    # so the label here is for the human reading this file.
    # Format: "<label> <rho between traits> <number of covariates> <covariate correlation>"
    for regime in "floor 0 0 0" "typical 0.5 0 0" "good_covariate 0.5 1 0.8"; do
      set -- $regime
      rho="$2"; n_cov="$3"; rho_cov="$4"

      # Only pass covariate flags when there is a covariate. With --n_cov 0 the filename carries no
      # covariate token, which is what keeps the pre-2026-11 results readable.
      COV_FLAGS=""
      [ "$n_cov" != "0" ] && COV_FLAGS="--n_cov $n_cov --rho_cov $rho_cov --cov_target c1"

      for miss in $MECHS; do
        for seed in $(seq 1 "$REPS"); do
          echo "--n $n --seed $seed --lambda $lambda --rho $rho $COV_FLAGS --miss $miss --frac 0.30 --M $M --arms $ARMS --discrete --save_imp $DRIVER_FLAG --bace_nitt $NITT --bace_burnin $BURNIN --bace_thin $THIN --bace_runs $RUNS --out $N_OUT" >> "$TASKS"
        done
      done

    done
  done
done

# ---- report what is about to happen, then do it -------------------------------------------------
TOTAL=$(wc -l < "$TASKS" | tr -d " ")
DONE=$(find "$OUT" -name '*.rds' 2>/dev/null | wc -l | tr -d " ")
echo "---------------------------------------------------------------"
echo "scenarios : $(( $(echo $N_LIST | wc -w) * $(echo $LAMBDAS | wc -w) * 3 * $(echo $MECHS | wc -w) ))"
echo "datasets  : $TOTAL queued  ($REPS per scenario)"
echo "already   : $DONE .rds present -- those lines will be skipped, not redone"
echo "cores     : $CORES"
echo "output    : $OUT/n<size>/"
echo "task list : $TASKS   <- read this to see every command that will run"
echo "---------------------------------------------------------------"
[ "${DRY_RUN:-0}" = "1" ] && { echo "DRY_RUN=1: task list written, nothing launched."; exit 0; }

echo "progress: find $OUT -name '*.rds' | wc -l   |   stop any time: kill this script"
< "$TASKS" xargs -P "$CORES" -L 1 Rscript "$HERE/rubin_cell.R" >> "$OUT/campaign.log" 2>&1
echo "campaign pass finished: $(find "$OUT" -name '*.rds' | wc -l | tr -d " ") rds present"
