# BACE completion roadmap

Consolidated view of everything remaining to take BACE from "usable research code"
to "submittable manuscript + release-ready package". Created 2026-07-11.
Priorities weight the user's stated emphasis: **numerical accuracy + careful testing**.

Current health: `devtools::test()` = 1268 pass / 0 fail / 0 skip / 0 warn (R 4.4.2,
MCMCglmm 2.36). Coverage-calibration fix (R nu=2, par_expand, K=1) merged to main.
Reference datasets generated (180 reps); production evaluation NOT yet run.
Manuscript `ms/ms.qmd` is a 98-line skeleton whose bib/template assets are missing.

---

## Robustness bugs surfaced by the n=1000 simulation (2026-07-12)

- [ ] **P1 (priority). Poisson imputation emits catastrophic counts.** 3.4% of
  poisson replicates produce astronomically large imputed counts (max standardised
  bias 217,838). Cause: `.predict_bace()` poisson branch clips the RATE at 1e6
  (`rate_k <- pmin(rate_k, 1e6)`) so a large `Liab` draw gives `rpois(1e6)` ≈ a
  million-count. Fix: clip the imputed COUNT to a data-adaptive ceiling (small
  multiple of max observed count), not the rate to 1e6. Add a regression test.
  See `dev/simulation_results/SIMULATION_REPORT.qmd` §5.
- [ ] **P2. Continuous/count 95% PI under-coverage** (gaussian cell-PI 0.76 at
  n=80, high signal). Intervals too narrow → residual-variance posterior too tight
  in this regime. Follow-up on the R-structure prior / n dependence.

## Track A — Numerical-correctness hardening  (DONE 2026-07-11)

Cheap, high-value, and exactly the "accuracy + testing" emphasis. Hardened the
pipeline before spending compute on it (Track B). Full suite after: 1274 pass /
0 fail / 0 warn / 0 skip.

- [x] **A1. Stale `K=3` comments in `.predict_bace`** — fixed in all FIVE sample
  branches (gaussian, poisson, threshold/ordinal, categorical, and the OVR path in
  `.fit_predict_ovr`). Comments now say "single draw (K=1)"; the K-loop/median/
  `.cat_mode` scaffolding is left in place (correct, RNG-equivalent) and annotated as
  a no-op hook — deliberately NOT refactored, to avoid regression risk on validated
  numerical code.
- [x] **A2. Explicit `sample = FALSE` at `bace_imp.R:352`** with a comment explaining
  why it must never be TRUE (two-phase invariant).
- [x] **A3. Determinism regression tests** (`tests/testthat/test-track-a-hardening.R`):
  `.predict_bace(sample=FALSE)` identical across calls (no RNG); `sample=TRUE` varies.
- [x] **A4. Aligned oracle-comparator prior** in `dev/10` to BACE's own gaussian prior
  (R nu=2 + Gelman par-expanded G) so `ci_width_ratio` is apples-to-apples.
- [x] **A5. Stacking-variance test** — `pool_posteriors()` coefficient variance =
  within + between (mixture identity), validating it propagates imputation uncertainty.
  (Full nominal-coverage validation is Track B's multi-rep job.)

## Track B — Run the evidence base  (compute-heavy; the manuscript's Performance section)

- [ ] **B1. Local smoke run** of `dev/10` (1-2 reps) to confirm the Track-A-fixed
  pipeline end-to-end before committing hours of MCMC.
- [ ] **B2. Production reference-benchmark run** (sim_ideal/typical/heterogeneous/hard,
  180 reps) via the GitHub Actions `sim-reference-benchmark` workflow or overnight
  locally. Produces the accuracy / coverage / calibration / beta-recovery tables.
- [ ] **B3. Competitor benchmark arm.** Abstract promises comparison "to leading
  alternatives"; the eval currently only scores vs oracle + dialled truth. Add at
  minimum `missForest` + a phylo-aware baseline (`Rphylopars`); `mice` if feasible.
  Note: `pigauto` is a companion package, not a competitor — frame accordingly.
- [ ] **B4. Full 2000-species AVONET run** at production MCMC (April report only did
  the 300-species subset; flagged as "worth doing before release").

## Track C — Manuscript + release

- [ ] **C1. Restore missing bib/template assets** — `ms/ms.qmd` references
  `bib/refs.bib`, a `.csl`, and `template.docx` that don't exist (only `bib/BACE.bib`
  is present), so it won't render. Blocking for any write-up.
- [ ] **C2. Draft Methods + Results** around Track B outputs (currently a bullet skeleton).
- [ ] **C3. Pooling section** — write up stacking vs Rubin with the O(1/M) table, and
  the pigauto companion-package positioning (division of labour, not rivalry).
- [ ] **C4. Release hygiene** — `NEWS.md`; prune exported internals (dozens of
  `.`-prefixed functions are in NAMESPACE); tighten toward a CRAN-ready surface.

## Track D — Rubin's-rules pooling pathway  (D1/D2/D4 DONE 2026-07-11)

Model-agnostic MI layer + free choice of combiner. Plan: `.agents/plans/pool_mi_rubin.md`.

- [x] **D1. `with_imputations(object, .f, tree=)`** (`R/with_imputations.R`) — fits any
  model (frequentist or MCMCglmm) per imputed dataset; accepts bace_complete /
  bace_final / list; captures per-fit errors. Class `bace_mi_fits`.
- [x] **D2. `pool_mi()`** (`R/pool_mi.R`) — Rubin's rules; accepts MCMCglmm (posterior
  mean/var of fixed effects, with a NOTE steering to pool_posteriors + away from
  variance components), Barnard-Rubin df, reports `fmi`/`riv`. Mirrors pigauto columns.
- [ ] **D3. Generalize `pool_posteriors()`** to accept `with_imputations()` output —
  DEFERRED (keep this PR focused; existing pool_posteriors tests untouched).
- [x] **D4. Tests** (`tests/testthat/test-pool_mi.R`, 26 assertions) — **golden test vs
  `mice::pool` matches to machine precision** (estimate/SE/df/fmi, ~1e-15); MCMCglmm
  path extracts only fixed effects + warns; edge cases (M<2, term mismatch, dropped
  failures). Full suite: 1300 pass / 0 fail.

## Recovery simulation — parameter + value recovery under MCAR and MAR

`dev/12_recovery_simulation.R` (env-var scaled; smoke defaults tiny). Uses the Track D
pathway: bace_imp -> bace_final_imp -> with_imputations(gls) -> pool_mi(). Compares
oracle / complete-case / BACE on slope recovery (bias, 95% coverage), marginal-mean
recovery (the MAR signature), and hidden-cell correlation.

DONE 2026-07-11 — two production studies, full report in
`dev/simulation_results/SIMULATION_REPORT.qmd`:
- **Study A** (dev/12, 40 reps, gaussian, known b1): BACE slope 95% coverage 0.925 (MAR)
  / 1.00 (MCAR) = nominal; slope bias equals the oracle's (unbiased vs full data); under
  MAR complete-case mean bias -0.26 vs BACE +0.005 (MI corrects MAR selection).
- **Study B** (dev/13, 300 datasets, 5 response types x MCAR/MAR x 30 reps, 0 failures,
  18.2 min): recovery stable across MCAR/MAR. Ranking ordinal (acc 0.82, ord_mae 0.17) >
  binary (0.63) > gaussian/poisson (cor ~0.4) > **categorical (bal-acc 0.35-0.39, near
  chance — FLAGGED)**. Timing 15-18 s/run except categorical ~42 s. Cell-level PI
  coverage 0.73-0.89 (mild under-coverage at n=80 weak-signal).
- [ ] **Follow-up: categorical sub-study** — dial multinomial signal strength + n, test
  ovr on/off, to see if weak multinomial recovery is a DGP artifact or a real limitation
  before any paper claim.
- [ ] Optional: larger-n runs for tighter cell-PI coverage.

---

## Recommended order

1. **Track A** (all of it) — days, low risk, closes the correctness gaps.
2. **Track D** — build on the hardened pipeline; independent of compute.
3. **Track B** (B1 smoke → B2 production → B3 competitors → B4 AVONET) — the long pole.
4. **Track C** — needs B's numbers; C1 (bib fix) can happen anytime.

## Already done (context, not TODO)

- Coverage-calibration investigation + fix (`dev/benchmark_report_2026-04-19.md`).
- Reference-dataset generation pipeline (`dev/09`) + eval/aggregate scripts (`dev/10`,
  `dev/11`) + sharded GHA workflow. Datasets on disk (180 reps).
- `phylo_signal_summary()` (branch `feature/phylo-signal-summary`; tested).
- `sim_bace()` simulation engine; OVR categorical; `nitt_cat_mult`; posterior-
  predictive final imputation; parallel final runs; accessors.

## Open decisions still to confirm

- Track B: run production benchmark on GitHub Actions vs locally overnight?
- Track B3: which competitors are in scope (missForest + Rphylopars minimum; add mice?).
- Track C4: how aggressive on trimming the exported API for v1?
