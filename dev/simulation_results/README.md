# BACE simulation results — file guide

All artifacts here are **tracked in git** and regenerable. Large raw bundles
(`reference_datasets/`, `evaluation_results/`, per-run RDS) are gitignored.

## The two studies

| Study | Script | What it does |
|---|---|---|
| **A — parameter recovery** | `dev/12_recovery_simulation.R` | gaussian response, known slope b1=0.8, 40 reps × {MCAR, MAR}, n=80. Oracle / complete-case / BACE. Slope bias + 95% coverage, marginal-mean bias, cell correlation. |
| **B — response-type recovery** | `dev/13_response_type_simulation.R` | 5 response types × {MCAR, MAR} × 30 reps = 300 datasets, n=80. Per-type value recovery, marginal bias, timing. |

Figures are built from the RDS by `dev/14_make_simulation_report.R`.

## Files

| File | Contents |
|---|---|
| `SIMULATION_REPORT.qmd` | The report source (Quarto): aims, questions, measures, model + formulas, results with mean±MCSE plots. |
| `SIMULATION_REPORT.html` | Rendered, **fully self-contained** HTML (open in any browser). Rebuild with `quarto render SIMULATION_REPORT.qmd`. |
| `recovery_simulation.rds` | Study A: `list(per_rep, summary, cfg)`. |
| `recovery_per_rep.csv` | Study A, one row per (mechanism, rep, method). |
| `recovery_summary.csv` | Study A aggregate. |
| `response_type_simulation.rds` | Study B: `list(per_rep, summary, cfg)`. |
| `response_type_per_rep.csv` | Study B, one row per (type, mechanism, rep). |
| `response_type_summary.csv` | Study B aggregate. |
| `figures/*.png` | Report figures (Okabe-Ito colourblind-safe palette). |
| `logs/*.log` | Full stdout of each run (per-replicate status + timing). |

## Column dictionaries

**Study A per-rep** (`recovery_per_rep.csv`): `mechanism` (MCAR/MAR), `rep_id`,
`method` (oracle/complete_case/bace), `slope_est`, `covered` (true b1=0.8 inside
95% CI?), `ci_width`, `mean_bias` (est. mean of y − true mean; 0 for oracle),
`cell_cor` (imputed vs truth on hidden cells; bace only).

**Study B per-rep** (`response_type_per_rep.csv`): `type`, `mechanism`, `rep_id`,
`status`, `runtime_sec`, then type-appropriate metrics (`correlation`, `nrmse`,
`coverage95` for gaussian/count; `accuracy`, `balanced_accuracy`, `brier` for
binary/categorical; `+ordinal_mae` for ordinal), `marg_bias_bace`,
`marg_bias_cc`, `n_missing`. Metrics not applicable to a type are `NA`.

## Inspect in R

```r
A <- readRDS("dev/simulation_results/recovery_simulation.rds")
B <- readRDS("dev/simulation_results/response_type_simulation.rds")
A$summary            # Study A aggregate
B$per_rep            # every replicate's metrics (Study B)
subset(B$per_rep, type == "categorical")   # audit the weak categorical case

# Regenerate everything (env vars scale reps/MCMC/cores; see script headers):
# Rscript dev/12_recovery_simulation.R
# Rscript dev/13_response_type_simulation.R
# Rscript dev/14_make_simulation_report.R   # figures + CSVs from the RDS
```

## Not yet captured (future, if deeper auditing is needed)

The per-rep CSVs hold aggregate **metrics** per replicate, not per-**cell**
truth-vs-imputed values. To audit individual imputations against truth, the
scripts would need to additionally save each replicate's complete data,
missingness mask, and imputed datasets — a re-run with per-rep bundle saving
(the pattern already used by `dev/09`/`dev/10`).
