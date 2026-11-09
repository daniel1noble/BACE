# Archive: 2026-08 simulation outputs (pre 2026-10 fixes)

Snapshot taken 2026-10-27, before re-running the validation studies on the
package with the five 2026-10 fixes (Ayumi's PR #15 review; see
`.agents/investigation-2026-10.md`).

- `recovery_*` — Study A (dev/12) outputs from the 2026-08 run. Re-run as a
  no-regression guard (y-only missing, so the sweeps change auto-collapses;
  numbers should match within Monte Carlo error).
- `categorical_{MCAR,MAR}_chunk.rds` — the two Study B (dev/13) categorical
  cells, **stale by construction**: they were produced while categorical
  prediction silently dropped all phylogenetic effects (problem 5). Moved out
  of `../response_type_chunks/` so dev/13's resume logic regenerates exactly
  these two cells.
- `response_type_simulation.rds` + csvs — the 2026-08 full-grid aggregate,
  archived because the re-aggregation after the categorical re-run overwrites
  it. Note the binary/ordinal `*_prob` metrics in the old aggregate were
  computed from plug-in probabilities (problem 2); their imputed values are
  unaffected. Full-grid regeneration stays on the manuscript roadmap.
