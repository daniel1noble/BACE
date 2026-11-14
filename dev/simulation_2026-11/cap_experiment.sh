#!/usr/bin/env bash
# =================================================================================================
# cap_experiment.sh -- does BACE's count rate cap bias the imputed counts, and by how much?
#
#   CORES=10 bash cap_experiment.sh results/        # pick the 10 campaign datasets the cap hit most
#   CORES=10 bash cap_experiment.sh                 # no campaign results to hand: a fixed list
#   R_LIBS=/path/to/lib105 CORES=10 bash cap_experiment.sh results/
#
# BACE caps an imputed count's rate at 3 * max(observed) + 5 to stop runaway draws. The cap is set
# from the counts it can SEE, so it binds -- and pulls the imputations down -- exactly when the
# missing counts are the large ones. Which datasets that happens in is recorded in every campaign
# result file, so with a campaign results directory this reruns the ten datasets where the cap
# fired most (cap_pick.R); without one it uses signal-0.3 datasets, where one measurement showed
# 39 cap events (at signal 0.9 the cap did not fire at all: the observed counts already include
# giants). Each dataset runs twice: the cap as in the campaign (_capon) and switched off
# (_capoff). Same seed, same data, same MCMC settings; only the cap differs. Then:
#
#   Rscript cap_compare.R results_cap/
#
# Needs BACE >= 1.0.5 (the cap became an option there). The campaign keeps running on 1.0.4: do
# NOT reinstall over a library that live workers are using -- install 1.0.5 into its own library
# and point R_LIBS at it, as in the third example above.
# =================================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:-$HERE/results_cap}"; CORES="${CORES:-4}"
mkdir -p "$OUT"
Rscript -e "v <- as.character(utils::packageVersion('BACE')); if (utils::compareVersion(v, '1.0.5') < 0) stop('BACE ', v, ' is too old for --rate_cap; this needs 1.0.5', call. = FALSE); cat('BACE', v, '| R', as.character(getRversion()), '| library', .libPaths()[1], '\n')" | tee "$OUT/versions.txt"
CAMPAIGN="${1:-}"
if [ -n "$CAMPAIGN" ]; then
  DATASETS=$(Rscript "$HERE/cap_pick.R" "$CAMPAIGN" 10)
else
  DATASETS=$(for seed in 1 2 3 4 5; do echo "n100_ps0.3_typical_clade $seed"; echo "n100_ps0.3_typical_mar $seed"; done)
fi
TASKS="$OUT/tasks.txt"; : > "$TASKS"
echo "$DATASETS" | while read -r id seed; do
  [ -z "$id" ] && continue
  echo "Rscript '$HERE/run_one.R' --id $id --seed $seed --out '$OUT' --no_checkpoint --tag _capon" >> "$TASKS"
  echo "Rscript '$HERE/run_one.R' --id $id --seed $seed --out '$OUT' --no_checkpoint --rate_cap off --tag _capoff" >> "$TASKS"
done
echo "$(wc -l < "$TASKS") runs (each dataset cap on and cap off), $CORES at a time -> $OUT"
xargs -P "$CORES" -I CMD bash -c CMD < "$TASKS"
echo "done. Now: Rscript '$HERE/cap_compare.R' '$OUT'"
