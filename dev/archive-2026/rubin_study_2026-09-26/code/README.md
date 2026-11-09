# Simulation code — read it in this order

Everything here does one job: generate datasets with known truth and known gaps, impute them,
and score how well the imputation did. The plan behind the current design is
[`.agents/plans/simulation-plan-2026-11.md`](../../../.agents/plans/simulation-plan-2026-11.md).

**If you only read one thing, read `rubin_cell.R`.** It takes one dataset from start to finish
and everything else is either something it calls or something that runs it many times.

## The chain, in order

| # | file | what it does |
|---|---|---|
| 1 | `campaign_gnn_off_lib.R` | **Makes the data.** `make_cell()` builds a tree, 7 traits with known correlations and known phylogenetic signal, the auxiliary covariate, and then punches the holes. The truth is kept, so every imputation can be scored against what the value really was. |
| 2 | `rubin_covariate.R` | **The auxiliary covariate** — one fully observed variable (think body size) that helps impute the others. Also the check that it actually reached the model, because when it did not, nothing complained. |
| 3 | `rubin_bace.R` | **The BACE arm.** Builds one imputation model per trait and hands them to `BACE::bace()`. |
| 4 | `rubin_completecase.R` | **The competitor.** Throw away every species with a gap, analyse the rest. This is what the manuscript's first claim is measured against. |
| 5 | `rubin_lib.R` | **The estimators and the maths.** PGLS slope, phylogenetic correlation, Rubin's pooling rules, and `cell_tag()`, which names every result file. |
| 6 | `rubin_discrete.R` | **Scoring for binary, ordinal and categorical traits** (accuracy with credit for ties, and category sets instead of intervals). |
| 7 | `rubin_cell.R` | **One dataset, end to end.** Builds it, runs every arm on it, scores them, writes one `.rds`. Everything above is called from here. |
| 8 | `rubin_campaign_local.sh` | **Runs it many times.** Writes one `rubin_cell.R` command per dataset into `tasks.txt`, then runs them N at a time. Read `tasks.txt` to see exactly what will run. |
| 9 | `rubin_missing.R` | **What has not run yet.** Compares `tasks.txt` against the files present and writes the leftovers back out, ready to relaunch. |

## Running it

```bash
# install the package being measured FIRST -- the scripts refuse to start on an older one
cd <BACE repo root> && R CMD INSTALL .

# one dataset, to see what happens (about a minute at n = 60 with --smoke)
Rscript code/rubin_cell.R --n 100 --seed 1 --lambda 0.7 --rho 0.5 \
  --n_cov 1 --rho_cov 0.8 --miss clade --arms bace,complete_case --discrete --out /tmp/try

# the whole campaign (see the header of the launcher for every knob)
DRY_RUN=1 bash code/rubin_campaign_local.sh      # write the task list and stop -- read it
CORES=7  bash code/rubin_campaign_local.sh       # actually run it; stop and restart freely

# what is left to do
Rscript code/rubin_missing.R results_2026-11

# the tests
Rscript -e 'testthat::test_dir("code/tests-rubin")'
```

## Supporting files

| file | what it does |
|---|---|
| `rubin_freq.R` | the frequentist comparison arms (Rphylopars). **Off** in the 2026-11 campaign; loaded only if `--arms` asks for `freqA`/`freqB`. |
| `rubin_campaign_aggregate.R`, `rubin_discrete_aggregate.R`, `rubin_full_aggregate.R` | turn thousands of `.rds` files into the summary tables the report and manuscript read |
| `oracle_floor.R` | the "best anyone could do" and "no better than the average" reference lines for the figures |
| `rubin_discrete_truth.R` | Monte Carlo estimate of the true discrete downstream slope, which is not simply rho |
| `rubin_totoro.sh`, `rubin_campaign_nibi.sh` | cluster launchers (one core per dataset; Nibi needs Shinichi's sign-off per its own header) |
| `tests-rubin/` | the test suite. Run it before launching anything long. |
| `archive-2026-09/` | finished tooling from the 2026-09 bundle — the HTML report builders, the pre-run stage, the one-off experiments. Kept for reference, used by nothing here. |

## Two things that bite, both silent

1. **The wrong BACE.** An older build sitting in the R library produces results that look entirely
   normal and mean nothing. This happened on 2026-11-09. Both `rubin_cell.R` and the launcher now
   refuse to start below `BACE_MIN_VERSION`, and every result file records the version that
   actually ran.
2. **Two scenarios sharing a filename.** The runner skips any dataset whose `.rds` already exists —
   that is what makes a three-week run restartable. It also means a filename collision causes the
   second scenario to be skipped silently and never run. `cell_tag()` in `rubin_lib.R` is the one
   place filenames are built, and `tests-rubin/test-campaign-grid.R` checks every scenario in the
   design gets its own name.
