# .agents/ — shared knowledge for BACE

Onboarding for collaborators (Dan, Shinichi, Szymek) and AI agents, so nobody
starts from scratch. Read `../AGENTS.md` first (what BACE is, the invariants,
build/test). Then:

| File | What it is |
|---|---|
| [`roadmap.md`](roadmap.md) | Everything remaining to reach a submittable manuscript + release. Tracks A–D, priorities, known bugs. **Start here for "what next".** |
| [`project-state.md`](project-state.md) | Snapshot of where things stand: what's done, test status, open items. |
| [`pigauto.md`](pigauto.md) | The companion package (Nakagawa): how it works, and how BACE differs. Manuscript positioning. |
| [`plans/pool_mi_rubin.md`](plans/pool_mi_rubin.md) | Design of the Rubin's-rules pooling pathway (`with_imputations()` + `pool_mi()`). |
| [`simulation-report.html`](simulation-report.html) | Self-contained simulation study (aims, model + formulas, results with mean±MCSE). Open in a browser. Source: `../dev/simulation_results/SIMULATION_REPORT.qmd`. |

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
Rscript dev/12_recovery_simulation.R        # parameter recovery (gaussian slope)
Rscript dev/13_response_type_simulation.R   # all response types, resumable (per-cell chunks)
Rscript dev/14_make_simulation_report.R     # per-rep CSVs + figures
quarto render dev/simulation_results/SIMULATION_REPORT.qmd   # -> self-contained HTML
```
