# AGENTS.md — BACE

Canonical brief for anyone (human or AI agent) working in this repo. Tool-agnostic;
Claude Code also reads `CLAUDE.md`, which imports this file. Detailed shared notes
live in [`.agents/`](.agents/).

## What BACE is

**B**ayesian **A**ugmentation by **C**hained **E**quations: an R package for multiple
imputation of missing data in **phylogenetic comparative** datasets. Every per-variable
imputation model is a phylogenetic mixed model fit with `MCMCglmm`, so related taxa
inform each other's imputations. Handles mixed types (continuous/gaussian, count/poisson,
binary & ordinal on the threshold scale, unordered multinomial). Feeds a manuscript
(`ms/ms.qmd`). Companion package: **pigauto** (Nakagawa) — see [.agents/pigauto.md](.agents/pigauto.md).

## Working norms (all contributors)

BACE feeds a manuscript; every methodological choice must be **defensible to peer
reviewers**. Be a critical reviewer / statistician, not a default-agreeable assistant:

- Evaluate proposals critically before implementing; raise better alternatives or known
  weaknesses. Constructive pushback = **(1) the issue, (2) the consequence, (3) a
  demonstration** (a number, a 3-line derivation, a short simulation — not an assertion).
- Back non-trivial statistical claims with peer-reviewed references (author, year, journal,
  specific result). Never invent citations.
- **Verify, don't assert.** When a claim is checkable (is it a bug in X or Y? does the
  imputation overshoot?), reproduce it and show the numbers. This repo's findings (e.g.
  the poisson count bug) were confirmed by worked examples, not reasoning alone.

## Build / test / dev

```r
devtools::load_all()          # interactive dev
devtools::document()          # regenerate NAMESPACE + man/ after roxygen edits
devtools::test()              # run all tests   (currently 1300+ pass / 0 fail)
devtools::test(filter = "pool_mi")   # one test file
devtools::check()             # full R CMD check (warnings = errors)
```
From the shell: `Rscript -e "devtools::test()"`. Tests skip on CRAN via
`skip_on_cran()`, so run with `NOT_CRAN=true` (devtools::test sets this).

## Architecture — the load-bearing invariants

The pipeline has **two phases**; conflating them breaks the method:

```
bace_imp()           convergence chain: point estimates (sample = FALSE).
                     NEVER pass sample = TRUE here — the chained-equations
                     convergence diagnostic needs the imputed values to settle,
                     which only deterministic point estimates do.
assess_convergence() chained-equations convergence check (single-chain stabilisation:
                     %-change + trend + ACF; NOT literally Gelman-Rubin R-hat).
bace_final_imp()     n_final independent posterior-PREDICTIVE draws (sample = TRUE).
                     This is what makes BACE multiple imputation; without it, between-
                     imputation variance = 0 and intervals ignore missing-data uncertainty.
pool_posteriors()    stacks per-imputation MCMC chains -> a Monte Carlo approximation to
                     the MARGINAL posterior. This is NOT Rubin's rules (Zhou & Reiter 2010;
                     BDA3 Ch.18). Differs from Rubin by O(1/M) at small M; default n_final=50.
```
`bace()` runs all four with auto-retry on non-convergence.

**Type detection** (`R/prep_functions.R::.detect_type`): numeric→gaussian, integer→poisson,
factor(2)→binary(threshold), factor(≥3)→categorical(multinomial probit), ordered→threshold.
`ovr_categorical=TRUE` (default) fits J one-vs-rest binary threshold models per categorical.

**Prediction core** `.predict_bace()` (`R/model_functions.R`): `sample=FALSE` → posterior
mean/argmax; `sample=TRUE` → one posterior-predictive draw (K=1; median-of-3 was dropped —
see `dev/benchmark_report_2026-04-19.md`).

## Repo layout

| Path | What |
|---|---|
| `R/`, `tests/`, `man/`, `NAMESPACE`, `DESCRIPTION` | the package (tracked, shipped) |
| `ms/`, `bib/` | manuscript + bibliography |
| `dev/` | dev scripts, benchmarks, simulations (tracked; NOT built into the package) |
| `dev/simulation_results/` | simulation report (`SIMULATION_REPORT.qmd`), figures, data |
| `.agents/` | **shared agent/onboarding knowledge** — start at `.agents/README.md` |
| `CLAUDE.md` | Claude-Code entry (imports this file) |

## Current state (2026-07-12)

- Package healthy: `devtools::test()` = **1300+ pass / 0 fail / 0 warn**. R 4.4.2, MCMCglmm 2.36.
- **Done this cycle:** Track A numerical-correctness hardening; Track D Rubin's-rules
  pathway (`with_imputations()` + `pool_mi()`, golden-tested vs `mice::pool`); two
  simulation studies validating recovery + calibration under MCAR/MAR (n up to 1000/cell).
- **Open / priority:** see [.agents/roadmap.md](.agents/roadmap.md). Notably a **poisson
  imputation bug** (log-link extrapolation on MAR-hidden high-leverage cells → catastrophic
  counts; fix = PMM or latent-scale clipping) and continuous-interval under-coverage.
- Full simulation write-up: [.agents/simulation-report.html](.agents/simulation-report.html)
  (self-contained) / source `dev/simulation_results/SIMULATION_REPORT.qmd`.

## Where to look

- **What to do next** → `.agents/roadmap.md`
- **Simulation results** → `.agents/simulation-report.html`
- **pigauto relationship** → `.agents/pigauto.md`
- **Rubin's-rules pathway design** → `.agents/plans/pool_mi_rubin.md`
