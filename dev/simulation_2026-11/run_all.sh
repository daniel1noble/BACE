#!/bin/bash
# =================================================================================================
# run_all.sh -- run the whole campaign by calling run_one.R once per dataset.
#
#   bash run_all.sh                    # the full design
#   REPS=2 bash run_all.sh             # a tiny rehearsal of exactly the same thing
#   CORES=20 bash run_all.sh           # 20 at a time
#   DRY_RUN=1 bash run_all.sh          # write the task list and stop -- read it
#
# It writes every command to tasks.txt before running anything, so you can see precisely what the
# campaign will do. Each line is a complete, copy-pasteable run_one.R call.
#
# Stop it whenever you like. run_one.R skips any dataset whose .rds already exists, so relaunching
# picks up where it left off, and raising REPS later only runs the new seeds.
# =================================================================================================
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"

CORES="${CORES:-7}"
REPS="${REPS:-1000}"
OUT="${OUT:-$HERE/results}"
M="${M:-20}"; NITT="${NITT:-50000}"; BURNIN="${BURNIN:-10000}"; THIN="${THIN:-25}"; RUNS="${RUNS:-5}"

mkdir -p "$OUT"

# Refuse to start on the wrong BACE. Checked here as its own command (not in a pipeline, where the
# exit status would be tee's) so a wrong version costs one second rather than 54,000 failures.
VERS=$(Rscript -e "v <- as.character(utils::packageVersion('BACE'))
  if (utils::compareVersion(v, '1.0.0') < 0) stop('BACE ', v, ' is too old; run  R CMD INSTALL .  from the repo root', call. = FALSE)
  cat('BACE', v, '| R', as.character(getRversion()), '\n')" 2>&1)
if [ $? -ne 0 ]; then echo "$VERS" >&2; echo "REFUSING TO LAUNCH." >&2; exit 1; fi
echo "$VERS" | tee "$OUT/versions.txt"

# The scenario list comes from R/01_scenarios.R -- one definition, used here and by run_one.R.
IDS=$(Rscript -e "source('$HERE/R/01_scenarios.R'); cat(make_scenarios()\$id, sep = '\n')")

TASKS="$OUT/tasks.txt"; : > "$TASKS"
for id in $IDS; do
  for seed in $(seq 1 "$REPS"); do
    echo "--id $id --seed $seed --out $OUT --M $M --nitt $NITT --burnin $BURNIN --thin $THIN --runs $RUNS" >> "$TASKS"
  done
done

TOTAL=$(wc -l < "$TASKS" | tr -d ' ')
DONE=$(ls "$OUT"/*.rds 2>/dev/null | wc -l | tr -d ' ')
echo "-------------------------------------------------------"
echo "scenarios : $(echo "$IDS" | wc -w | tr -d ' ')"
echo "datasets  : $TOTAL queued ($REPS per scenario)"
echo "already   : $DONE done -- those will be skipped, not redone"
echo "cores     : $CORES"
echo "task list : $TASKS   <- every command that will run"
echo "-------------------------------------------------------"
[ "${DRY_RUN:-0}" = "1" ] && { echo "DRY_RUN=1: nothing launched."; exit 0; }

< "$TASKS" xargs -P "$CORES" -L 1 Rscript "$HERE/run_one.R" >> "$OUT/campaign.log" 2>&1
echo "finished: $(ls "$OUT"/*.rds 2>/dev/null | wc -l | tr -d ' ') results present"
