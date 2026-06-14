# Simulated reference benchmark — running on your own hardware

The four-step simulated reference benchmark (generate → evaluate → aggregate)
is designed to run unattended on a multi-core machine or HPC node. Off-CI you
do **not** need the per-chunk sharding the GitHub Actions workflow uses — there
is no 6h job ceiling, so you just give `mclapply` all your cores and run the
whole sweep in one go.

All paths below are relative to the repo root. Run from there after checking
out the branch with these scripts.

## 0. Prerequisites (once)

```r
install.packages(c(
  "devtools", "ape", "MASS", "coda", "MCMCglmm",
  "dplyr", "magrittr", "stringi", "methods",
  "phytools", "caper"          # Suggests; used by sim + phylo-signal paths
))
```

`MCMCglmm` compiles from source — allow a few minutes the first time.

## 1. Generate the four reference datasets (~5 min)

```bash
Rscript dev/09_generate_reference_datasets.R
```

Writes `dev/simulation_results/reference_datasets/{sim_ideal,sim_typical,
sim_heterogeneous,sim_hard}/rep_NN.rds` plus a `meta.rds` per dataset.
180 reps total (30 / 50 / 50 / 50). This directory is gitignored, so it is
regenerated locally and never committed.

## 2. Evaluate BACE on every replicate

```bash
# Use as many cores as you have RAM for (see the memory note below).
BACE_EVAL_CORES=32 Rscript dev/10_evaluate_reference_datasets.R
```

With no positional arguments this evaluates **all reps of all datasets**.
`BACE_EVAL_CORES` controls how many reps run concurrently via `mclapply`
(default 6). Each `bace()` call itself runs single-threaded, so total load is
roughly `BACE_EVAL_CORES` busy cores.

Production MCMC budget is baked into the script defaults and matches the CI
production run:

| field        | default | env override   |
|:-------------|--------:|:---------------|
| nitt         |  50000  | `BACE_NITT`    |
| burnin       |  10000  | `BACE_BURNIN`  |
| thin         |     25  | `BACE_THIN`    |
| runs         |     10  | `BACE_RUNS`    |
| n_final      |     50  | `BACE_NFINAL`  |

Example quick pre-flight (tiny budget, 1 rep of one dataset, in-process):

```bash
BACE_NITT=2000 BACE_BURNIN=500 BACE_THIN=5 BACE_RUNS=2 BACE_NFINAL=4 \
  Rscript dev/10_evaluate_reference_datasets.R sim_ideal 1
```

Other CLI forms (all intersect with the reps actually present on disk):

```bash
Rscript dev/10_evaluate_reference_datasets.R sim_hard          # all reps of one dataset
Rscript dev/10_evaluate_reference_datasets.R sim_hard 5 12     # reps 5..12 of one dataset
```

**Resumable.** Each rep writes `eval_rep_NN.rds`; a rep already saved with
`status == "ok"` is skipped on re-run. If the machine is interrupted, just
re-run the same command — it picks up where it left off.

**Memory note.** Each concurrent rep holds pooled MCMC chains
(`n_final` × stored draws across 5 equations), which can run to a few GB at the
production budget. A safe ceiling is roughly
`BACE_EVAL_CORES ≈ min(physical_cores, total_RAM_GB / 3)`. If you see workers
killed, lower `BACE_EVAL_CORES` or `n_final`.

## 3. Aggregate into a cross-dataset summary

```bash
Rscript dev/11_aggregate_reference_evaluations.R
```

Writes to `dev/simulation_results/evaluation_results/`:

- `summary_report.md` — headline numbers (status, convergence rate, per-rep
  wallclock, cell-level accuracy, beta coverage vs oracle and dialled truth)
- `per_rep_long.csv` / `.rds` — one row per (dataset, rep, variable, metric)
- `beta_compare_long.csv` — per-coefficient BACE vs oracle vs dialled
- `timings.csv` — per-rep BACE / oracle runtimes

## One-liner for the full sweep

```bash
Rscript dev/09_generate_reference_datasets.R && \
BACE_EVAL_CORES=32 Rscript dev/10_evaluate_reference_datasets.R && \
Rscript dev/11_aggregate_reference_evaluations.R
```

Budget guide: the full sweep is on the order of 80–100 CPU-hours at the
production budget, so wall time ≈ that divided by `BACE_EVAL_CORES`
(e.g. ~3 h on 32 cores). Because step 2 is resumable, it is safe to run under
`nohup`/`tmux` and reconnect.
