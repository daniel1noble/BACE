# Rubin MI simulation study: BACE vs frequentist imputation (continuous and discrete traits)

Bundle date 2026-09-26. Source: `arc/rubin-freq-bace`.
(The original bundle README said "Private; not for circulation" — pushing/publishing this folder
is held until Shinichi confirms; the pristine imported bundle is git commit `eaca7e3`.)

## 2026-10 re-run — this directory is now the live campaign

The 2026-09 results below were produced with the PRE-FIX BACE package (its "BACE as shipped" arm
measured the improper-MI defect this study itself identified, independently reported in PR #15
and fixed in BACE 0.1.0 on 2026-10-26). The campaign is being re-run on the corrected package:

- Arms: `freqA`, `freqB` (Rphylopars on continuous traits, castor on discrete — display names
  match the ms software table), `bace` (**"Current BACE"**: each final dataset from an
  independent run with 3 internal Gibbs sweeps — the released algorithm) and `bace_chain`
  (the 2026-09 chained variant, pinned to `sweeps = 1` per link, kept for comparison).
  The `bace_resid` negative control is retired with the defect it probed.
- Local launcher: `code/rubin_campaign_local.sh` (tiered, resumable; see its header). The cluster
  launchers still work with the corrected BACE installed cluster-side.
- New results land in `results_2026-10/`; the report rebuilds from `report/` once the campaign's
  aggregation (`code/rubin_campaign_aggregate.R`, `code/rubin_study/build_data.R`) is re-run.
  Until then, `report/study.html` and `report/data/` describe the 2026-09 (pre-fix) run.

## Start here

`report/study.html` is the full write-up (methods, results, supplement), self-contained; open it in a browser.
`report/study.qmd` is its source. It rebuilds from the tables in `report/data/` with
`quarto render study.qmd --to html` run inside `report/` (needs R and knitr only; no refitting).

## Contents

| folder | what |
|---|---|
| `report/` | `study.qmd`, `study.html`, and `data/`, the per-fit and aggregated result tables (the file table in section 5 of the report describes each) |
| `code/` | every script of the study: data generation (`campaign_gnn_off_lib.R`, `make_cell()`), one dataset through all methods (`rubin_cell.R`), methods (`rubin_freq.R`, `rubin_bace.R`, `rubin_discrete.R`, `rubin_lib.R`), cluster launchers (`rubin_campaign_nibi.sh`, `rubin_totoro.sh`, `rubin_prerun_*.sh`), pooling of result files (`rubin_pool.sh`, `rubin_disc_pool.sh`), aggregation (`rubin_campaign_aggregate.R`, `rubin_discrete_aggregate.R`, `rubin_study/build_data.R`), Monte Carlo target (`rubin_discrete_truth.R`), rescoring (`rubin_discrete_rescore.R`), page build (`rubin_report_build.sh`, `rubin_disc_pages_data.py`, `*_template.html`) and tests (`tests-rubin/`) |
| `pages/` | the three campaign pages (report, summary, accuracy) and the pre-run results page, as built HTML |
| `records/` | plans, pre-run plan, review, after-task reports for each stage |

## Reproducing

The scripts source each other by their repository paths (`script/...`), so run them from a checkout of the
pigauto repository on branch `arc/rubin-freq-bace`, where `code/` here is `script/`. One dataset, all methods, with
the campaign's settings:
`Rscript script/rubin_cell.R --n 100 --lambda 0.7 --rho 0.5 --seed 1 --out results/ --discrete --save_imp --bace_runs 5`
(flags in the header of `rubin_cell.R`; the script's own default is 10 BACE runs, the campaign used 5).
Tests: `Rscript -e 'testthat::test_dir("script/tests-rubin")'`.
BACE 0.1.0, MCMCglmm 2.36, Rphylopars 0.3.10, castor 1.8.7, ape 5.8.1, R 4.5.

## Not included

The raw per-fit result files (`.rds`, one per dataset and method set, including all 20 imputations of every
method in the discrete part; 1.4 GB) are in a separate archive, `rubin_raw_pools_2026-09-26.tar.gz`.
