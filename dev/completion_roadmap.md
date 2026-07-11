# BACE completion roadmap

Consolidated view of everything remaining to take BACE from "usable research code"
to "submittable manuscript + release-ready package". Created 2026-07-11.
Priorities weight the user's stated emphasis: **numerical accuracy + careful testing**.

Current health: `devtools::test()` = 1268 pass / 0 fail / 0 skip / 0 warn (R 4.4.2,
MCMCglmm 2.36). Coverage-calibration fix (R nu=2, par_expand, K=1) merged to main.
Reference datasets generated (180 reps); production evaluation NOT yet run.
Manuscript `ms/ms.qmd` is a 98-line skeleton whose bib/template assets are missing.

---

## Track A — Numerical-correctness hardening  (low compute; DO FIRST)

Cheap, high-value, and exactly the "accuracy + testing" emphasis. Hardens the
pipeline before we spend compute on it (Track B).

- [ ] **A1. Stale `K=3` comments in `.predict_bace`** (`R/model_functions.R` ~585-771).
  Every branch sets `K <- 1L` but comments still describe "median of K=3 draws".
  Fix comments; remove/annotate the now-vestigial median/`.cat_mode` machinery.
- [ ] **A2. Explicit `sample = FALSE` at `bace_imp.R:352`.** The only prediction call
  that relies on the default; it's the one site where `sample=TRUE` would *silently*
  break the convergence diagnostic. Make it explicit (its OVR sibling already is).
- [ ] **A3. Determinism regression test.** Call the convergence-chain prediction twice
  on one fit; assert identical imputed values. Locks the `sample=FALSE` contract.
- [ ] **A4. Align oracle-comparator prior** in `dev/10_evaluate_reference_datasets.R`
  (currently `nu=0.002` while BACE uses `nu=2`) — the `ci_width_ratio` is not
  apples-to-apples. Match BACE's prior, or add a second matched-prior oracle arm.
- [ ] **A5. Pooled-coverage test.** On simulated data with known beta, assert
  `pool_posteriors()` 95% intervals cover the truth near-nominally at `n_final=50`.
  Empirically validates the stacking combiner (currently only exercised structurally).

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

## Track D — Rubin's-rules pooling pathway  (planned: `dev/pool_mi_rubin_plan.md`)

Model-agnostic MI layer + free choice of combiner. Land after Track A.

- [ ] **D1. `with_imputations(object, .f)`** — fit any model (frequentist or MCMCglmm)
  per imputed dataset.
- [ ] **D2. `pool_mi()`** — Rubin's rules; accepts MCMCglmm (posterior mean/var),
  Barnard-Rubin df, reports `fmi`/`riv`. Mirrors pigauto's column set.
- [ ] **D3. Generalize `pool_posteriors()`** to accept `with_imputations()` output, so
  the same fits go to either combiner.
- [ ] **D4. Tests** — golden test vs `mice::pool`; Bayesian-fit tests; variance-
  component guard (never Rubin-pool variance components).

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
