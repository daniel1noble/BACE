# 2026-08 performance investigation: root causes and fixes

Why BACE under-performed in the 2026-07 Study B (gaussian coverage 0.76,
poisson 0.69 + catastrophic counts, categorical near chance), what was actually
wrong, and what was changed. Investigated by four independent agent
investigations with adversarial verification; every load-bearing claim below
was reproduced on the actual failing replicates (exact seeds) or verified by
direct code reading. Archived pre-fix outputs:
`dev/simulation_results/archive_2026-07_prebugfix/`.

## Root causes (ranked by impact)

### 1. Improper final-phase MI — THE bug (fixed)
`bace_imp` reinstates response NAs before each fit so MCMCglmm does data
augmentation (R/bace_imp.R). **`bace_final_imp` never did**: every final
"posterior-predictive" model was fit on the converged dataset with the
deterministic point imputations treated as *observed data*. Consequences
(all reproduced numerically): residual variance deflated (imputed cells sit on
the fitted surface: residual SD 0.16 vs 1.45 at real cells); BLUPs anchored to
the fills (poisson: cor 0.996 between point fills and final-draw means;
per-cell draw SD 5.7 vs 17.8 proper); between-imputation variance collapsed.
This is textbook improper MI (Rubin 1987 Ch. 4; van Buuren FIMD §3.2) and hurt
**every** variable type. Fix: reinstate the response NAs per formula in
`.one_final_run` (mirrors `bace_imp`); regression test asserts the fit sees
NAs at `miss_dat` rows. Single-fix effect on failing reps: poisson MAR rep-1
coverage 0.464→~1.0; gaussian per-cell predictive SD ratio 0.68→~1.0.

### 2. Finite-draw metric artifact (evaluation, not model)
The sims used n_final=15; the empirical 2.5–97.5% quantile interval of m draws
covers at most ~0.845 at m=15 **even for a perfectly calibrated imputer**
(order stats: the [min,max] of m draws covers (m−1)/(m+1)). Poisson MCAR
(0.880) was already *at* its ceiling. The discrete "95% set" from 15 draws
degenerates to the full draw support (0.95 unreachable at 1/15 resolution), so
binary/ordinal's "good" coverage was partly granularity-inflated. Fix (dev/13
v3): added `cover95_t` (mean ± t·sd·√(1+1/m), not ceiling-capped), prob-based
discrete sets from the models' averaged probability matrices, and
Spearman/log1p recovery for poisson (Pearson-on-raw-counts is tail-dominated).
The v2 quantile metrics are retained for before/after comparability.

### 3. Poisson: mean-of-exponentials point estimate + fixed 1e6 rate clip (fixed)
`.pred_count` point-imputed `colMeans(exp(Liab))` — one chain excursion
dominates the lognormal mean (rep 587: one exp(20.4) draw ⇒ point imputation
1,039,177 for truth 5), and bug #1 froze it into every dataset via
`pmin(rate, 1e6)` ⇒ deterministic `rpois(1e6)`. Fix: median-based point
estimate (`post_median` = exp of median liability) + data-adaptive rate guard
(3·max observed count + 5) in both branches. Post-fix rep 587: max imputation
17 (truth 5), coverage 0.935. PMM was evaluated head-to-head and **rejected**
(−0.3 to −0.5 SD bias under extrapolative MAR); NB/OLRE redundant (the units
variance already is an OLRE — its prior was the problem, see #4).

### 4. Priors: R nu=2 everywhere, gaussian-copied poisson priors (fixed)
Roadmap P2's hypothesis ("residual-variance posterior too tight") was
**backwards**: IW(V=1, nu=2) on z-scored data (true units var ~0.08 at λ=0.9)
*inflates* the residual variance ~3× and deflates λ̂ (0.755 vs 0.879 at
nu=0.002) — partially masking bug #1 by widening intervals. The 2026-04 nu=2
change was tuned against a benchmark whose real defects were bug #1 and the
K=3 median guard. Fix: gaussian R back to Hadfield's V=1, nu=0.002 (the
par-expanded G priors, which nu=0.002 predated, prevent the old near-zero
collapse); poisson gets log-link-scaled priors R = (V=0.1, nu=1) and PX
working prior alpha.V=1 (SD~1 on the log scale, not SD~100).

### 5. Categorical: mostly a benchmark artifact + a bad default (both fixed)
The sim DGP was **degenerate**: all K−1 multinomial liabilities shared one
linear predictor and one set of random draws, differing only by intercepts
(±0.25), so class-probability ratios were constant across observations and one
class could never be the Bayes argmax; dev/13 additionally gave categorical
random weak betas (`runif(−0.5, 0.5)`) instead of the strong known slopes every
other type got. **Bayes balanced-accuracy ceiling ≈ 0.50; observed 0.43–0.46
was near-optimal.** Fixes: `sim_bace` gained `beta_resp_liab` (per-liability
slopes) + independent per-liability phylo/species/residual draws + returns
`response_probs` (true class probabilities) so studies can report the ceiling;
dev/13 now uses strong per-liability betas and records per-rep
`bayes_bal_acc`. Separately, the OVR default was genuinely miscalibrated in
the 1-obs-per-species regime (quasi-separation, exploding phylo variance,
incoherent normalized probabilities) — **`ovr_categorical` default flipped to
FALSE** (true multinomial: better calibrated AND ~2.6× faster); OVR models are
now kept in `all_ovr_models` instead of being dropped.

### 6. Poisson-vs-gaussian type detection was a marginal-AIC race (fixed)
Found during post-fix validation: `.get_type` chose poisson for integer
columns only if plain-Poisson beat gaussian by ≥2 AIC on the *marginal*
intercept-only glm. A count trait with strong phylo signal is marginally
overdispersed, so poisson lost exactly when signal was strong → the variable
was silently imputed as gaussian (negative, non-integer "counts"), and the
family flip-flopped between replicates of the same design. New rule:
non-negative integers → poisson (the conditional model's units variance
handles overdispersion as an OLRE); integers with negatives → gaussian.

## Validation status (final, 2026-08-09)
- Full test suite green after fixes (1337 pass / 0 fail); new regression
  tests in `tests/testthat/test-performance-fixes.R`.
- Spot checks on the exact worst pre-fix replicates (Study B seeds):
  gaussian MAR rep 1 cover 0.895(q)/1.00(t); poisson MAR rep 1 cover 0.964;
  poisson rep 587 max imputation 1,039,177 → 17.
- **Full re-runs complete** (same seeds; 9,997/10,000 ok): gaussian
  t-interval coverage 0.949/0.947 (nominal; quantile metric at its 0.845
  ceiling); poisson recovery 0.58/0.62 → 0.69/0.74 with max |std bias|
  217,838 → 5.3 and catastrophe rate 14.6% → 6.5% (MAR); categorical accuracy
  0.50 → 0.81 (bal-acc at ~89% of the now-reported Bayes ceiling), set
  coverage 0.75 → 0.98, ~1.7× faster; binary/ordinal prob-set coverage
  0.94–0.96. Study A: MAR slope coverage 0.925 → 0.950, MAR mean-bias
  correction intact, CI width 1.59 → 1.83 (restored between-imputation
  variance). Remaining known gaps: poisson MAR t-coverage 0.916 + mild
  positive mean bias (intrinsic log-link extrapolation under MAR);
  categorical sets slightly conservative (0.984).

## Investigation provenance
4 parallel investigator agents (poisson pathway / interval calibration & MI
propriety / categorical & OVR / workflow & evaluation audit), findings
adversarially verified by independent agents instructed to refute them.
Literature anchors: Rubin (1987) ch.4; Little & Rubin (2002) §4.3; van Buuren
(2018) FIMD §3.2 §4.5; Gelman (2006) Bayesian Analysis 1:515; Hadfield (2010)
J. Stat. Softw. 33(2) + MCMCglmm Course Notes.
