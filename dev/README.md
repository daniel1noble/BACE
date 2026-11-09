# dev/

Development code. Not built into the package.

## → `simulation_2026-11/` — the live simulation

**This is the one.** The simulation behind the manuscript: it generates comparative datasets with
known truth, hides values in three different ways, imputes them with `BACE::bace()`, and scores
how well they came back.

Plan and design: [`.agents/plans/simulation-plan-2026-11.md`](../.agents/plans/simulation-plan-2026-11.md)

Four steps, one file each, in order:

| file | step |
|---|---|
| `R/01_scenarios.R` | the list of scenarios — print it to see what will run |
| `R/02_simulate.R` | one simulated dataset from one scenario's settings |
| `R/03_fit.R` | run BACE on it, and the complete-case comparison |
| `R/04_score.R` | pull out the statistics |

`ms/ms.qmd` reads its numbers from `simulation_2026-11/results/`.

## `archive-2026/` — finished work

Everything that came before, moved here in November 2026 so the live simulation is easy to find.
Nothing here is used to build the package, the manuscript or the vignettes. Nothing was deleted,
and git keeps the history.

| what | why it is here |
|---|---|
| `rubin_study_2026-09-26/` | the 2026-09/10 campaign (BACE vs Rphylopars vs castor). Superseded by `simulation_2026-11/`; its aggregated tables are kept because the manuscript quoted them until the new campaign runs. |
| `simulation_results/` | the 2026-08 simulation report (`SIMULATION_REPORT.qmd`) and its figures |
| `00_`–`17_*.R` | the numbered benchmark and validation scripts, each a finished one-off |
| `benchmark_*`, `testing_data/` | benchmark engine, results and reference datasets |
| `sim_accuracy_*` | the accuracy sweep the vignette's results table is built from |
| `notes/`, `figure_prompts/` | working notes |
