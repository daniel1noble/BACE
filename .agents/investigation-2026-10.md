# Investigation 2026-10: Ayumi's review of BACE 0.1.0 (pre-release)

Status: fixes implemented and committed 2026-10-26; validation campaign in
progress. Plan with full detail: [plans/fixes-2026-10-ayumi.md](plans/fixes-2026-10-ayumi.md).

## Origin

While building the bird-migration tutorial (draft PR #15, 10,563 AVONET
species, ordinal + mixed-variable imputation with 5-fold CV), Ayumi Mizuno
reported four implementation concerns, each with empirical evidence. We
verified every claim against the code (HEAD 01e4f6d at the time): **all four
confirmed**, every line number she cited correct. The verification also found
a fifth, worse bug she had not raised. Her results commit (b321ec3) and HEAD
had identical R/ code, so all findings apply to her tutorial numbers.

Key context for "how did our sims miss this": both existing studies (dev/12,
dev/13) only ever made **one** variable (y) missing. With a single missing
variable the pre-fix final phase is actually exact, so Study A's 0.949
coverage was real — it just never tested the regime Ayumi hit (two traits
missing together).

## The five problems (verified numbers)

| # | Problem | Reproduced evidence | Fix (commit) |
|---|---------|---------------------|--------------|
| 5 (new, worst) | Default categorical path **ignores the phylogeny**: species-BLUP lookup built `trait.level.<tip>` but MCMCglmm names columns `trait.level.<REterm>.<tip>` → never matched → effect 0 | 40-spp fit: zeroing ALL species effects changed predictions by exactly **0**; correct lookup moves probabilities up to **0.65**, flips **32%** of predicted categories. Mean abs BLUP 1.83 | prefix+suffix matching via `.gather_blups(_mat)` + zero-match tripwire warning; also adds the missing `Species2` effects (species=TRUE) for categorical AND ordinal (d2ac904, tests b69f268) |
| 1 (Ayumi pt 1) | Final runs don't propagate predictor uncertainty: one pass per run, in formula order, from the same point-imputed start → first-listed variable conditions on identical point fills in all runs | Spy: x values seen by y's kept model had **SD = 0** across runs at every x-missing row (matches her 2,113-row AVONET observation); after fix SD 0.095–0.78 | `sweeps` argument (default 3, warm-up passes short-tailed ≈ 1.4x cost, auto-collapse to 1 when ≤ 1 variable missing) + per-run `seed` streams (serial ≡ parallel) + fallback redoes only failed runs (ace9ed5) |
| 2 (pt 2) | `prob_preds` plug-in (posterior-mean coefficients) → overconfident | Her data: mean max-prob 0.901 vs 0.848 honest, RPS 0.0647 vs 0.0622; our toy: 0.999 vs 0.917 | forward helpers refactored into 2 cores; probabilities averaged per draw; single-draw imputation path unchanged by construction (ad57001) |
| 3 (pt 3) | Diagnostics never check cutpoints, and never check final-phase models at all | Her Furnariidae pilot: cutpoint ESS 5–22, unreported | `.diagnose_one_model` core incl. CP; `.check_final_diagnostics` on every final model (incl. OVR sub-models), worst-case across runs, stored as `$diagnostics`; deliberately no cross-run R-hat (runs target different conditional posteriors) (8a34066) |
| 4 (pt 4) | Priors hard-coded; `prior=` silently swallowed; phases disagreed on gelman | `prior=` verified to vanish into `...`; `bace_final_imp` hard-coded gelman=2 vs option-read in `bace_imp` | per-variable `priors` argument on all three functions, blocks replaced wholesale, R-override refused for latent-scale types, OVR reuses the variable's prior, gelman unified, unknown-arg warning (840eaba) |

Defaults were deliberately NOT changed (Ayumi's suggested ordering): the knob
ships first; new defaults only with benchmarks.

## Statistical rationale for the sweeps fix (problem 1)

BACE has no joint model (one model per trait — that's what makes mixed types
feasible), so "sample the joint posterior predictive" must be manufactured by
Gibbs cycling of the per-trait conditionals (fully conditional specification;
van Buuren 2007, Stat Methods Med Res 16:219–242). One pass from a fixed
deterministic start is a Gibbs sampler stopped after half a step: the
first-visited trait never receives feedback. Expected damage where two+
variables are missing together: between-imputation variance of the
first-listed variable collapses to MC noise → Rubin's total variance
understated → ~0.90–0.93 actual coverage where 0.95 is claimed (fmi 0.25–0.35
regime), plus formula-order dependence. Warm-up passes need only one
approximately-stationary draw, so they keep full burn-in but a 50-draw tail
(~21% of a full fit at defaults).

## Test coverage added (suite: 1337 → 1418 pass / 0 fail, ~1.4 min)

- Hand-computed forward probabilities INCLUDING species effects (categorical +
  ordinal, averaged + single-draw), dotted-tip and node-column decoys,
  tripwire warning, degenerate-posterior plug-in equivalence, and a REAL-fit
  regression test (mocks cannot catch MCMCglmm naming drift).
- Sweeps: fit counts (runs x sweeps x variables), warm-tail nitt per sweep,
  the SD>0-across-runs regression, auto-collapse, seed reproducibility
  (serial == parallel), fallback-equals-reference-on-streams.
- Priors: validation, wholesale block replacement, dimension guards, spies
  proving the override reaches `.model_fit` in both phases and only for the
  targeted variable, gelman-option unification (failed pre-fix).
- Diagnostics: cutpoint rows + labels, worst-run aggregation, OVR coverage,
  end-to-end `$diagnostics` on a real pipeline.

## Validation campaign

- `dev/15_sweeps_calibration.R` (new): BOTH x and y 30% missing, MCAR + MAR,
  sweeps k in {1,2,3,5,8} paired within replicate (shared data + chained
  phase), 100 reps, n_final 15. Metrics: Rubin coverage/CI width (pool_mi,
  deliberately — understated B shows directly in T = W + (1+1/M)B),
  between-imputation variance B_m, between-run SD at hidden cells, joint
  y-x draw correlation at doubly-missing rows (≈ 0 at k=1 by construction),
  accuracy, wall time. k=1 column = published quantification of the old bug.
  **Running 2026-10-26** → `dev/simulation_results/sweeps_calibration.rds`.
- Study A re-run (dev/12, y-only; sweeps auto-collapse) as a no-regression
  guard, and dev/13 categorical cells to quantify problem 5's cost — queued
  after calibration.
- The 2026-08 categorical Study B numbers are **stale by construction** (they
  were produced with the phylogeny silently off) and must be regenerated for
  the manuscript.

## Knock-on effects

- CRAN: still on hold; tarball must be rebuilt and win-builder re-run after
  the validation campaign. The fixes land BEFORE first submission.
- Ayumi's tutorial (PR #15): her ordinal headline results are NOT affected by
  problem 5 (threshold path was correct); her mixed-variable section
  (problem 1) and plug-in probability figure (problem 2) need a re-run on the
  fixed package; her prior pilots become actionable via `priors`. Reply
  drafted for Dan to post.
- Tutorial/ms text that said "adjust priors" is now true; ms pseudo-code
  listing `priors` as an input is now true.
