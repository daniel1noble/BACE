# The simulation — final file set

**Agreed 2026-11-09.** These are the files that make up the BACE simulation
campaign. Anything not listed here is not part of it.

Plan and design: [`.agents/plans/simulation-plan-2026-11.md`](../../.agents/plans/simulation-plan-2026-11.md)
How to read the code: [`code/README.md`](code/README.md)
To see it work: `Rscript code/walkthrough.R`

## Run the campaign — these nine files, in this order

| # | file | lines | what it does |
|---|---|---|---|
| 1 | `code/campaign_gnn_off_lib.R` | 838 | Makes the data: tree, seven traits with known correlations and signal, the covariate, then punches the holes. Keeps the truth. |
| 2 | `code/rubin_covariate.R` | 110 | The auxiliary covariate (think body size), and the check that it reached the model. |
| 3 | `code/rubin_bace.R` | 156 | The BACE arm: one imputation model per trait, handed to `BACE::bace()`. |
| 4 | `code/rubin_completecase.R` | 150 | The competitor: drop every species with a gap, analyse the rest. |
| 5 | `code/rubin_lib.R` | 295 | The estimators and the maths: PGLS, phylogenetic correlation, Rubin pooling, `cell_tag()`. |
| 6 | `code/rubin_discrete.R` | 252 | Scoring for binary, ordinal and categorical traits. |
| 7 | `code/rubin_cell.R` | 365 | **One dataset, end to end.** Everything above is called from here. Read this one first. |
| 8 | `code/rubin_campaign_local.sh` | 133 | Runs it many times. Writes every command to `tasks.txt` first. |
| 9 | `code/rubin_missing.R` | 89 | What has not run yet; writes the leftovers back out ready to relaunch. |

## Aggregate and report

| file | what it does |
|---|---|
| `code/rubin_campaign_aggregate.R` | continuous-trait results -> summary tables |
| `code/rubin_discrete_aggregate.R` | discrete-trait results -> summary tables |
| `code/rubin_full_aggregate.R` | the combined table the manuscript reads |
| `code/oracle_floor.R` | the "best anyone could do" and "no better than the average" reference lines |
| `code/rubin_discrete_truth.R` | Monte Carlo estimate of the true discrete downstream slope (it is not simply rho) |

## Run it on a cluster

| file | what it does |
|---|---|
| `code/rubin_totoro.sh` | Totoro: no queue, one dataset per core, capped at 150 cores (D-143) |
| `code/rubin_campaign_nibi.sh` | Nibi: SLURM array. **Needs Shinichi's sign-off per its own header.** |

## Tests — run before launching anything long

`Rscript -e 'testthat::test_dir("code/tests-rubin")'` — 16 test files, **455 assertions, 0 failures**.

The three that guard this campaign's new work:

| file | guards against |
|---|---|
| `code/tests-rubin/test-covariate.R` | the covariate not existing, being the wrong strength, or never reaching the model |
| `code/tests-rubin/test-completecase.R` | the survivor rule, the refusal floor, and attrition by mechanism |
| `code/tests-rubin/test-campaign-grid.R` | two scenarios writing the same filename, so one is silently skipped |

## Not part of the simulation

| path | why it is here |
|---|---|
| `code/rubin_freq.R` | the frequentist comparison arms. **Off** in this campaign; loaded only if `--arms` asks for them. Kept for the competitor benchmark. |
| `code/archive-2026-09/` | finished 2026-09 tooling: HTML report builders, the pre-run stage, one-off experiments. Used by nothing here. |
| `results_2026-10/` | the previous campaign's 750 result files. Superseded, kept for the cost measurements. |
| `report/`, `pages/`, `records/` | the 2026-09 write-up and its records. |

## Before you launch

1. `R CMD INSTALL .` from the BACE repo root. The scripts refuse to start on an
   older BACE than `BACE_MIN_VERSION` — on 2026-11-09 the repo was at 1.0.0
   while the installed library was still 0.1.0.9000, and nothing warned anybody.
2. `DRY_RUN=1 bash code/rubin_campaign_local.sh` and read `tasks.txt`.
3. Run the tests.
4. The four decisions in §8 of the plan, two of which are cheap checks worth
   doing first.
