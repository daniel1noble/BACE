# .agents/ — shared knowledge for BACE

Onboarding for collaborators (Dan, Shinichi, Szymek) and AI agents, so nobody
starts from scratch. Read `../AGENTS.md` first (what BACE is, the invariants,
build/test). Then:

| File | What it is |
|---|---|
| [`roadmap.md`](roadmap.md) | Everything remaining to reach a submittable manuscript + release. Tracks A–D, priorities, known bugs. **Start here for "what next".** |
| [`project-state.md`](project-state.md) | Snapshot of where things stand: what's done, test status, open items. |
| [`investigation-convergence-2026-11.md`](investigation-convergence-2026-11.md) | **The convergence-check repair (2026-11)**: the chained-equations convergence verdict measured each variable's scale, not its stability (never passed centred traits, passed drifting ones, called a collapsed sampler "converged"). Replaced in 1.0.1 by a calibrated, scale-free check (drift / oscillation / stuck); bace() retries only on drift. Full evidence, literature and calibration tables. |
| [`investigation-audit-2026-11.md`](investigation-audit-2026-11.md) | **The 2026-11 numerical audits**: three rounds of adversarial review for numbers that are *wrong* rather than code that crashes; every finding reproduced with a printed number before it counted. The recurring defect class: a derived quantity treated as a measurement. |
| [`investigation-simulator-reproducibility-2026-11.md`](investigation-simulator-reproducibility-2026-11.md) | **Same seed, different data on different machines (fixed 2026-11-13, 1.0.3)**: `MASS::mvrnorm` builds its matrix square root from eigenvectors, which have no fixed sign; three maths libraries made three valid choices. Now a Cholesky factor, which is unique. Golden test pins it. |
| [`investigation-ordinal-generator-2026-11.md`](investigation-ordinal-generator-2026-11.md) | **The simulated ordinal trait collapsed at low signal (fixed 2026-11-14, 1.0.4)**: fixed cutpoints on a latent scale that grows with the signal. A first fix (sample quantiles) was rejected by review for coupling hidden to observed cells; shipped: population thresholds, `sim_bace(ordinal_thresholds=)`, level counts recorded per dataset. Every trait type pinned by a golden test. |
| [`investigation-2026-08.md`](investigation-2026-08.md) | **The 2026-08 performance investigation**: six root causes of the Study B coverage/accuracy failures (improper final-phase MI, finite-draw metric ceiling, poisson point-estimate/clip, mis-scaled priors, degenerate multinomial DGP, AIC-race type detection) and the fixes. |
| [`pigauto.md`](pigauto.md) | The companion package (Nakagawa): how it works, and how BACE differs. Manuscript positioning. |
| [`plans/pool_mi_rubin.md`](plans/pool_mi_rubin.md) | Design of the Rubin's-rules pooling pathway (`with_imputations()` + `pool_mi()`). |
| [`simulation-report.html`](simulation-report.html) | Self-contained simulation study (aims, model + formulas, results with mean±MCSE). Open in a browser. Source: `../dev/archive-2026/simulation_results/SIMULATION_REPORT.qmd`. |

## Conventions

- **This directory is the durable, shared record.** When you finish a substantive
  piece of work or learn something non-obvious, update `project-state.md` and
  `roadmap.md` so the next person (or agent) inherits it.
- `dev/` holds working scripts/data (tracked, not built into the package). The
  curated, human-facing summaries live here in `.agents/`.
- Simulation data files (`*.rds`) and rendered `*.html` under `dev/` are gitignored
  (regenerable / large); the self-contained `simulation-report.html` is kept here so
  it ships even though loose `*.html` is ignored.

## Reproduce the simulation study

```bash
Rscript dev/archive-2026/12_recovery_simulation.R        # parameter recovery (gaussian slope)
Rscript dev/archive-2026/13_response_type_simulation.R   # all response types, resumable (per-cell chunks)
Rscript dev/archive-2026/14_make_simulation_report.R     # per-rep CSVs + figures
quarto render dev/archive-2026/simulation_results/SIMULATION_REPORT.qmd   # -> self-contained HTML
```
