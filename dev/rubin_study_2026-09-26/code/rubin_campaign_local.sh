#!/bin/bash
# Local (single-machine) launcher for the 2026-10 campaign re-run on the CORRECTED BACE package.
# 2026-10-27 scope (Dan): BACE ONLY -- the goal is accuracy, low bias and nominal coverage for
# Current BACE across the full scenario grid (n x lambda x rho), scored against the per-dataset
# complete-data reference; competitor arms are off by default (ARMS env restores them). Fully
# resumable: rubin_cell.R skips any cell whose rds already exists, so this script can be stopped
# and relaunched at will, and seed counts can be RAISED later (new seeds append, old ones skip).
#
# Usage (defaults: BACE only, 50 seeds at n=100/300, 25 at n=1000,
# nitt 50000 / burnin 10000 / thin 25 / runs 5, M = 20, discrete scoring at n <= 300):
#   CORES=7 bash code/rubin_campaign_local.sh
# Tiered runs, e.g. the n=100 grid only:
#   CORES=7 N300_SEEDS=0 N1000_SEEDS=0 bash code/rubin_campaign_local.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
CORES="${CORES:-7}"
OUT="${OUT:-$HERE/../results_2026-10}"
ARMS="${ARMS:-bace}"
NITT="${NITT:-50000}"; BURNIN="${BURNIN:-10000}"; THIN="${THIN:-25}"; RUNS="${RUNS:-5}"; M="${M:-20}"
N100_SEEDS="${N100_SEEDS:-50}"; N300_SEEDS="${N300_SEEDS:-50}"; N1000_SEEDS="${N1000_SEEDS:-25}"

mkdir -p "$OUT"
Rscript -e 'cat("BACE", as.character(utils::packageVersion("BACE")), "| MCMCglmm", as.character(utils::packageVersion("MCMCglmm")), "| Rphylopars", as.character(utils::packageVersion("Rphylopars")), "| castor", as.character(utils::packageVersion("castor")), "\n")' | tee "$OUT/versions.txt"

TASKS="$OUT/tasks.txt"; : > "$TASKS"
for n in 100 300 1000; do
  case $n in
    100)  S=$N100_SEEDS ;;
    300)  S=$N300_SEEDS ;;
    1000) S=$N1000_SEEDS ;;
  esac
  [ "$S" -eq 0 ] && continue
  DISC=""; [ "$n" -le 300 ] && DISC="--discrete --save_imp"   # 2026-09: discrete at n = 100, 300 only
  for lambda in 0.3 0.7 1; do
    for rho in 0 0.5; do
      for seed in $(seq 1 "$S"); do
        echo "--n $n --seed $seed --lambda $lambda --rho $rho --M $M --arms $ARMS $DISC --bace_nitt $NITT --bace_burnin $BURNIN --bace_thin $THIN --bace_runs $RUNS --out $OUT" >> "$TASKS"
      done
    done
  done
done

TOTAL=$(wc -l < "$TASKS" | tr -d " ")
echo "$TOTAL cells queued -> $OUT  (cores = $CORES; resume = on: existing rds are skipped)"
echo "progress: grep -c 'wrote' $OUT/campaign.log   |   stop any time: kill this script"
< "$TASKS" xargs -P "$CORES" -L 1 Rscript "$HERE/rubin_cell.R" >> "$OUT/campaign.log" 2>&1
echo "campaign pass finished: $(ls "$OUT"/*.rds 2>/dev/null | wc -l | tr -d " ") rds present"
