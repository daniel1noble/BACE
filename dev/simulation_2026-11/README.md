# The simulation

Generates comparative datasets with known truth, hides values three different ways, imputes them
with `BACE::bace()`, and scores what came back.

Design and reasoning: [`.agents/plans/simulation-plan-2026-11.md`](../../.agents/plans/simulation-plan-2026-11.md)

## Four steps, one file each

| file | step |
|---|---|
| `R/01_scenarios.R` | **the scenarios.** `make_scenarios()` returns a table — print it to see what will run. |
| `R/02_simulate.R` | **one dataset.** `sim_bace()` makes the traits and tree; the three missingness mechanisms are ours, because `sim_bace()` only hides values at random. |
| `R/03_fit.R` | **run the methods.** `bace()` for imputation, plus the complete-case comparison. |
| `R/04_score.R` | **the statistics.** Recovery of the hidden values, and the pooled PGLS slope. |

Then three scripts that use them:

| file | what |
|---|---|
| `run_one.R` | **one dataset, end to end.** The unit of work — run it by hand to see exactly what happens. |
| `run_all.sh` | runs `run_one.R` once per dataset. Writes `tasks.txt` first, so you can read every command before anything starts. |
| `collect.R` | gathers the results into five CSVs. |

## Running it

```bash
# 1. install the package being measured. The scripts refuse to start on an older one.
cd <BACE repo root> && R CMD INSTALL .

# 2. one dataset, to watch it work
Rscript run_one.R --id n100_ps0.7_good_covariate_clade --seed 1 --out results

# 3. see what the campaign would do, without doing it
DRY_RUN=1 bash run_all.sh        # writes results/tasks.txt

# 4. run it. Stop and restart whenever you like.
CORES=8 bash run_all.sh

# 5. turn the results into tables
Rscript collect.R results
```

## The three arms

| arm | what it is | what it is for |
|---|---|---|
| `oracle` | the truth, before any holes | the best anyone could do — and the **calibration check**: its bias must be ~0 and its coverage ~0.95, or the assumed true slope is wrong and nothing else means anything |
| `complete_case` | drop every species with a gap | what people do now; the comparison claim 1 rests on |
| `bace` | impute with `bace()` | the method the paper is about |

## Things that will bite, and what is done about them

**The wrong BACE.** An older build in the R library gives results that look entirely normal and
mean nothing — it happened on 2026-11-09, repo at 1.0.0 and library at 0.1.0.9000. Both
`run_one.R` and `run_all.sh` refuse to start below 1.0.0, and every result records the version
that produced it.

**Two known BACE crashes**, both reproduced on this grid and both deterministic:

| where | error | hit |
|---|---|---|
| `sim_bace()` | `NA in probability vector` — the softmax in `mnom_liab2cat()` overflows (`exp(710)` is `Inf`, `Inf/Inf` is `NaN`) | 3 of 54 scenarios, all n = 500 with high signal |
| `bace()` | `zero-variance series` — `coda::geweke.diag()` is called unguarded at `R/model_functions.R:1635`, while the same call is wrapped in `tryCatch` at `R/phylo_signal_summary.R:850` | 2 of 54 scenarios |

A crashed job writes nothing, which leaves a silent gap. So every stage is wrapped: a failure is
**recorded** and the dataset still writes its file. `collect.R` counts them, writes
`failures.csv`, and reports the rate. Being deterministic, these do not resume — only fixing them
in BACE will.

**Interruptions are different.** A wall-time limit or a dead node can strike any dataset, and at
n = 500 that is 2.2 hours lost. `bace()`'s `checkpoint_dir` is on by default: a resumed run takes
about 5% of the time and gives an identical result. Checkpoints are large (51.6 MB even for a tiny
fit), so the directory is **deleted as soon as the dataset succeeds** — only interrupted datasets
leave one behind, which is exactly what should be resumed. `--no_checkpoint` turns it off.

**Resumability at the dataset level** is separate and always on: `run_one.R` skips any dataset
whose `.rds` already exists, so the campaign can be stopped and relaunched freely, and replicates
can be raised later without redoing earlier seeds.

## What the numbers mean

`phylo_signal` is phylogenetic heritability on the **latent** scale, which is the scale BACE's
threshold models work on. It is delivered as asked and equally for every trait type. Measured on
the observed values a binary trait looks far less phylogenetic than a continuous one at the same
setting — that is what dichotomising a measurement does, not a flaw — so every dataset records its
realised observed-scale lambda.

`trait_beta` and `aux_beta` are **regression slopes, not correlations**. A slope buys less
correlation than its own value: slope 0.5 gave r = 0.45 and slope 0.8 gave r = 0.62 over 25
datasets at n = 200. Every dataset records its own realised correlation.
