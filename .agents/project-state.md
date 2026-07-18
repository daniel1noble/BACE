# Project state (updated 2026-07-12)

Keep this current — it's the shared "where are we" so nobody re-derives it.

## Health
- `devtools::test()` = **1300+ pass / 0 fail / 0 skip / 0 warn**. R 4.4.2, MCMCglmm 2.36.
- Coverage-calibration fix merged to main earlier: R-structure prior `nu = 2`,
  `par_expand = TRUE`, single posterior-predictive draw (K=1). Rationale +
  evidence in `dev/benchmark_report_2026-04-19.md`.

## Done this cycle (2026-07)
- **Track A — numerical-correctness hardening:** stale `K=3` comments fixed in all
  5 `.predict_bace`/`.fit_predict_ovr` sample branches (K=1 machinery kept as an
  annotated no-op hook, not refactored); explicit `sample = FALSE` at
  `bace_imp.R:352`; oracle-comparator prior in `dev/10` aligned to BACE's; new
  `tests/testthat/test-track-a-hardening.R` (predict determinism + stacking-variance
  identity).
- **Track D — Rubin's-rules pathway:** `R/with_imputations.R` (fit any model per
  imputed dataset; class `bace_mi_fits`) + `R/pool_mi.R` (Rubin's rules, accepts
  MCMCglmm via posterior mean/var, Barnard-Rubin df, reports `fmi`/`riv`).
  **Golden test matches `mice::pool` to ~1e-15** (`tests/testthat/test-pool_mi.R`).
  D3 (generalising `pool_posteriors`) deferred. Plan: `plans/pool_mi_rubin.md`.
- **Simulation studies** (`dev/12`, `dev/13`, `dev/14`; report in
  `simulation-report.html`):
  - Study A (parameter recovery, gaussian, 40 reps): BACE slope 95% coverage nominal
    (0.925 MAR / 1.00 MCAR), unbiased vs the oracle, and **corrects the MAR
    complete-case mean bias** (-0.26 → +0.005).
  - Study B (all 5 response types, MCAR/MAR, **1000 reps/cell = 10,000 datasets**,
    phylo signal λ=0.90): gaussian/binary/ordinal **unbiased** (median bias ≈ 0),
    recovery 0.86-0.88; binary/ordinal calibration **nominal** (0.945/0.968).
    Poisson recovery ~0.6; categorical the hardest (balanced acc 0.43-0.46 vs
    chance 0.33) and slowest (~30 s/dataset).

## Open / priority (details in roadmap.md)
- **P1 — poisson imputation bug.** ~3.4% of poisson replicates produce catastrophic
  imputed counts (e.g. true = 2, imputed ≈ 3000; max standardised bias 217,838).
  Cause verified: log-link extrapolation for high-leverage cells, which MAR
  preferentially hides; the `rate` clip at 1e6 in `.predict_bace` is far too loose.
  Fix: predictive mean matching (PMM) for counts, or latent-scale clipping. **Not
  a `sim_bace` artefact** — the true data is well-behaved.
- **P2 — continuous/count 95% PI under-coverage** (gaussian ≈ 0.76): predictive
  intervals too narrow. Residual-variance posterior too tight in this regime.
- **Track B** — run the production reference benchmark + a competitor arm
  (missForest / Rphylopars). **Track C** — manuscript (restore missing bib assets;
  draft Methods/Results). **Release** — NEWS.md, trim exported internals.
