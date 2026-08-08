# Project state (updated 2026-08-08)

Keep this current — it's the shared "where are we" so nobody re-derives it.

## Health
- `devtools::test()` = **1337 pass / 0 fail / 0 skip** (R 4.4.2, MCMCglmm 2.36).
  (One benign, stochastic warning from the intentional non-convergence fixture
  in test-bace_wrapper.R.)
- **2026-08: six root causes of the Study B performance failures found, fixed,
  and regression-tested.** Full write-up: [investigation-2026-08.md](investigation-2026-08.md).
  Headline: `bace_final_imp` never reinstated response NAs, so every final
  "posterior-predictive" model conditioned on the convergence chain's point
  imputations as observed data (improper MI → the gaussian 0.76 / poisson 0.69
  under-coverage); plus the poisson mean-of-exponentials point estimate &
  fixed 1e6 rate clip (the catastrophic counts), gaussian-copied poisson
  priors, R nu=2 (the 2026-04 "fix" — actually inflates residual variance and
  biases λ̂ down; reverted to nu=0.002), a degenerate multinomial sim DGP
  (Bayes ceiling ~0.50 — categorical was near-optimal all along), and a
  marginal-AIC poisson/gaussian type-detection race that silently imputed
  strong-signal count variables as gaussian.

## Done this cycle (2026-08)
- **Package fixes** (all with regression tests in `test-performance-fixes.R`):
  NA reinstatement in `.one_final_imp` runs; poisson `post_median` point
  estimate + data-adaptive rate guard (3·max_obs+5); family-specific priors
  (gaussian R V=1 nu=0.002; poisson R V=0.1 nu=1, PX alpha.V=1);
  `ovr_categorical` default → FALSE (true multinomial: better calibrated,
  ~2.6× faster); OVR models kept in `all_ovr_models`; deterministic
  non-negative-integer → poisson type rule.
- **Simulation machinery**: `sim_bace` non-degenerate multinomial DGP
  (`beta_resp_liab` per-liability slopes, independent per-liability random
  draws, returns `response_probs` for Bayes-ceiling reporting); dev/13 v3
  metrics separating model calibration from finite-draw artifacts
  (`cover95_t`, prob-based discrete sets, Spearman/log recovery,
  per-rep `bayes_bal_acc`).
- **Spot-checked on the exact worst pre-fix replicates**: poisson rep 587
  max imputation 1,039,177 → 17; poisson MAR rep 1 coverage 0.464 → 0.964;
  gaussian MAR rep 1 coverage → 0.895(q)/1.00(t).

## Re-validation (2026-08-09, full re-runs complete)
Study B = 9,997/10,000 replicates ok, 5.3 h wallclock. Same seeds/config as the
2026-07 run, so differences are pure fix effect. Report re-rendered
(`simulation-report.html`); pre-fix archive in
`dev/simulation_results/archive_2026-07_prebugfix/`.

| Type | Recovery pre→post | Cover(q) pre→post (ceiling 0.845) | Cover t/prob (target .95) |
|---|---|---|---|
| gaussian | 0.83/0.88 → same | 0.756/0.784 → **0.844/0.845 (= at ceiling)** | **t: 0.949/0.947** |
| poisson | 0.58/0.62 → **0.69/0.74** (log-scale 0.75/0.78) | 0.691/0.880 → 0.867/0.902 | t: 0.916/0.932 |
| binary | 0.86 → same | 0.944 → 0.98 (draw-set, inflated) | prob: 0.941/0.938 |
| ordinal | 0.88 → same | 0.968 → 0.99 | prob: 0.959/0.957 |
| categorical | acc 0.50/0.52 → **0.81/0.78**; bal-acc 0.63/0.65 vs Bayes ceiling 0.71/0.72 | 0.751/0.769 → 0.984 | prob: 0.984 |

Poisson catastrophes: P(|std bias|>1) 14.6%→6.5% (MAR), 8.8%→2.8% (MCAR); max
|std bias| **217,838 → 5.3**. Categorical runtime ~30-42s → ~22s (true
multinomial default). Study A: MAR slope coverage 0.925→**0.950**, MCAR 0.975;
complete-case MAR mean-bias still corrected (−0.259→+0.011); CI width 1.59→1.83
(restored between-imputation variance); cell correlation unchanged.

**Remaining honest limitations** (for the ms): poisson MAR t-interval coverage
0.916 (log-link extrapolation on MAR-hidden high-leverage cells is
intrinsically hard; mild positive mean bias ~0.15-0.18 SD, median ~0.04-0.08);
categorical prediction sets slightly conservative (0.984); 3/10,000 replicate
failures (2 length-zero-argument errors, 1 singular mixed-model equations) —
low-priority follow-up.

## Open / priority (details in roadmap.md)
- **Track B** — production reference benchmark + competitor arm
  (missForest / Rphylopars). **Track C** — manuscript (restore missing bib
  assets; draft Methods/Results). **Release** — NEWS.md, trim exported internals.
- Optional follow-ups from the investigation: multiple stochastic sweeps per
  final run for multi-variable-missing data; strengthen the convergence
  diagnostic (Geweke silently disabled); fit-once-draw-many optimization when
  only one variable is missing.
