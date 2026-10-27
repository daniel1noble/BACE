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
  y-x draw correlation at doubly-missing rows (predicted ≈ 0 at k=1 — wrong,
  see results: flat ~0.28 at every k, because the second-drawn trait always
  conditions on the first's fresh draw; the real k=1 damage is pairing bias),
  accuracy, wall time. k=1 column = published quantification of the old bug.
  **Done 2026-10-27** → results below; raw output in
  `dev/simulation_results/sweeps_calibration.{rds,log}` (currently gitignored).
- Study A re-run DONE 2026-10-27: numerically identical to the 2026-08 run
  (max |Δ slope| ≈ 1e-8; 240/240 ok) — no regression; the sweeps
  auto-collapse on y-only data confirmed empirically.
- dev/13 categorical cells DONE 2026-10-27 (1000 reps/cell, same seeds as
  2026-08; the stale cells are archived in
  `dev/simulation_results/archive_2026-08_pre-oct-fixes/`). Fixed vs stale:
  accuracy 0.783→0.799 (MCAR), 0.807→0.815 (MAR); balanced accuracy
  0.635→0.651, 0.628→0.642 — ~90% of the 0.722/0.709 Bayes ceilings (from
  ~88%); prob-set 95% coverage 0.984→0.992/0.993 (the posterior-averaged
  probabilities are slightly more conservative); class-proportion bias
  0.111→0.105 (MAR); runtime HALVED (22→12.5 s/rep — the vectorised
  prediction refactor). 1 failure/2000 (singular MME), same class as the
  three known 2026-08 failures. The modest accuracy gain is expected: this
  DGP's signal is mostly in the fixed effects; the fix matters most when
  phylogeny is the main source of information.

## Sweeps calibration results (2026-10-27) — default sweeps = 3 CONFIRMED

Setup: 100 replicates, 60 species, x and y each 30% missing, MCAR and MAR;
each replicate's data and chained phase shared across k ∈ {1,2,3,5,8}, so the
k contrast is paired (only the final phase differs). nitt 3000, n_final 15,
7 cores. 200/200 replicate rows ok; 20.5 min total. True slope b1 = 0.8.

| mech | k | slope bias | cover95 | CI width | B_m  | SD at hidden y |
|------|---|-----------:|--------:|---------:|-----:|---------------:|
| MAR  | 1 | +0.137     | 0.96    | 2.10     | 0.142| 0.634          |
| MAR  | 2 | +0.069     | 0.98    | 2.24     | 0.164| 0.671          |
| MAR  | 3 | +0.030     | 0.99    | 2.26     | 0.163| 0.678          |
| MAR  | 5 | +0.038     | 0.99    | 2.28     | 0.170| 0.677          |
| MAR  | 8 | +0.038     | 0.98    | 2.24     | 0.162| 0.678          |
| MCAR | 1 | +0.109     | 0.96    | 2.08     | 0.126| 0.599          |
| MCAR | 2 | +0.055     | 0.97    | 2.19     | 0.143| 0.637          |
| MCAR | 3 | +0.053     | 0.97    | 2.24     | 0.157| 0.645          |
| MCAR | 5 | +0.060     | 0.97    | 2.21     | 0.148| 0.641          |
| MCAR | 8 | +0.047     | 0.97    | 2.17     | 0.142| 0.640          |

Final-phase wall time per run by k: 4.6 / 5.8 / 7.1 / 9.8 / 13.7 s
(k=3 ≈ 1.54× k=1, ≈ half of k=8). Point accuracy (cor of imputed vs true y)
flat across k: ~0.71 MAR / ~0.77 MCAR. Joint y–x draw correlation at
doubly-missing rows flat ~0.27–0.29 across k.

**Verdict: k = 3 sits on the k = 8 plateau on every metric** — slope bias
(0.030 vs 0.038 MAR), B_m (0.163 vs 0.162 MAR; 0.157 vs 0.142 MCAR), SD at
hidden cells (0.678 vs 0.678; 0.645 vs 0.640), coverage within 0.01 — at
about half the k = 8 cost. Default stays `sweeps = 3`; no code change needed.

**The k = 1 arm (the old behaviour) shows the bug was worse than an
uncertainty problem.** The kept y draws were generated against x fills that
were then replaced by fresh x draws, so the shipped (y, x) pairs do not
belong together. That pairing mismatch biased the downstream slope upward by
+0.11 to +0.14 (14–17% of the true 0.8), on top of understating
between-imputation variance by ~13–20% and the spread at hidden cells by ~7%.
Coverage itself held at 0.96 in this DGP (Rubin's small-M t is conservative
here, masking the understatement), so the headline harm to report is the
bias, not coverage. One prediction was wrong and is corrected above: the
joint y–x correlation was NOT ≈0 at k=1 — the damage mechanism is incoherent
pairs, not missing correlation.

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
