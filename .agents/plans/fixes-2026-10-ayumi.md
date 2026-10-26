# Plan: Fix the problems Ayumi found in BACE (and one she didn't)

## Status (2026-10-26)

- [x] Problem 5 (categorical phylogeny lookup) — fixed d2ac904, tests b69f268
- [x] Problem 4 (priors argument + gelman unification) — 840eaba
- [x] Problem 1 (sweeps + seed + fallback) — ace9ed5
- [x] Problem 2 (posterior-averaged probabilities) — ad57001
- [x] Problem 3 (cutpoint + final-phase diagnostics) — 8a34066
- [x] NEWS.md, investigation note (.agents/investigation-2026-10.md),
      PR #15 reply draft (.agents/pr15-reply-draft.md — Dan to review/post)
- [ ] dev/15 sweeps calibration — RUNNING (started 2026-10-26, ~2-3 h;
      output dev/simulation_results/sweeps_calibration.{rds,log});
      confirm/revise sweeps default from its plateau
- [ ] Study A re-run (dev/12, no-regression guard) + dev/13 categorical cells
      (quantify problem 5) — queued after calibration
- [ ] Tutorial text updates (priors passages now true; mention sweeps) +
      full re-render; SIMULATION_REPORT update
- [ ] Rebuild tarball, fresh win-builder round; CRAN stays on hold
- Suite: 1418 pass / 0 fail (~1.4 min). Verified reproduction numbers are in
  the investigation note.

## Background

While building the bird-migration tutorial (PR #15), Ayumi reported four problems in
BACE 0.1.0. I verified every claim against the code: **all four are real**, and the
verification found a fifth, worse bug she didn't spot. None of them crash — they
silently give wrong or overconfident numbers, which is worse.

BACE is not on CRAN yet. All fixes land before submission (CRAN stays on hold).

Why our own simulations never caught problem 1: both existing studies (dev/12,
dev/13) only ever made **one** trait (y) missing — and with a single missing trait
the current code is actually exact. The problem only appears when two or more
traits are missing together, which is Ayumi's real-data situation.

---

## Problem 1 — Final imputations understate uncertainty (Ayumi's point 1)

**The problem.** BACE builds 50 final datasets. The traits depend on each other, so
they are imputed one at a time, in formula order. Whichever trait goes first is
imputed while the other traits' missing cells still hold the *fixed best guesses*
from the convergence phase — the same values in all 50 datasets (Ayumi measured
SD = 0). `sample = TRUE` is already on everywhere; it randomises the value being
drawn, not the values sitting in the other columns. So the first trait (usually the
analysis response y) never feels the other traits' uncertainty, and formula order
changes results.

**Why it matters.** Understated between-imputation variance → confidence intervals
too narrow (estimated ~0.90–0.93 actual coverage where 0.95 is claimed, in the
both-traits-missing regime). Capturing this uncertainty is the whole point of
multiple imputation. Dan: must fix.

**The fix.** Repeat the fill-in cycle a few times within each run (new `sweeps`
argument on `bace_final_imp()` and `bace()`, default 3). By pass 2 the first trait
is re-imputed against the *random draws* of the others from pass 1. Only the last
pass's values and models are kept — one value per missing cell per dataset. This is
a Gibbs sampler: cycling the per-trait draws is how a chained system produces one
draw from the joint distribution, cross-trait covariance included (van Buuren 2007,
Stat Methods Med Res 16:219–242). The convergence phase is untouched (`sample =
FALSE` there, by design).

Three refinements from the design review:

- **Auto-skip:** when only one trait has missing data, extra passes change nothing
  (mathematically exact), so BACE collapses to 1 pass — no cost for the most common
  use case, and all previously validated single-missing-trait behaviour is
  unchanged.
- **Cheap warm-up passes:** only the final pass feeds the pooled inference, so the
  earlier passes need just one valid draw, not a long chain — they run with a short
  post-burn-in tail (burn-in never shortened). Total cost ≈ 1.4× current, not 3×.
- **Reproducibility:** optional `seed` argument with one RNG stream per run, so
  serial and parallel give identical results and reruns reproduce. The parallel
  failure-fallback also gets fixed to redo only the failed runs instead of all 50.

Main code: wrap the variable loop at `R/bace_final_imp.R:131` in a sweep loop.
Default `sweeps = 3` is confirmed (or revised) by the calibration study below, not
guesswork.

---

## Problem 2 — Reported probabilities are overconfident (Ayumi's point 2)

**The problem.** For categorical/ordinal traits BACE reports each category's
probability (`prob_preds`). These are computed by plugging the posterior *average*
into the probability formula, which exaggerates certainty (Ayumi: "90% sure" where
honest is ~85%; our toy check: 0.999 vs 0.917). The imputed values themselves are
proper draws and are fine — only the reported probabilities are off. They are also
returned undocumented and untested.

**The fix.** Compute the probability for each posterior draw and average the
probabilities (correct order of operations). Implementation: merge the four
prediction helpers into two core functions where "one draw" and "average over
draws" share the same math — the single-draw path used for imputation stays
bit-identical in behaviour. Document `prob_preds` (per-run matrices; averaging
across runs gives the pooled probabilities). Must land after Problem 5's fix,
otherwise we'd be averaging species effects that are all zero.

---

## Problem 3 — The MCMC health check skips cutpoints and final models (Ayumi's point 3)

**The problem.** The health check looks at fixed effects and variances but never at
the *cutpoints* (the thresholds defining ordinal categories) — Ayumi had cutpoint
effective sample sizes of 5–22 and BACE said nothing. Worse: the check only ever
runs on the convergence-phase models; the 50 final models — the ones actually used
for inference — are never checked at all.

**The fix.** Add cutpoints to the check (`R/model_functions.R:1544`), run it on
every final-phase model (including the one-vs-rest sub-models, currently
invisible), and store a compact summary in the result: worst effective sample size
and worst Geweke score per parameter across runs, with only *flagged* variables
printed. One deliberate exclusion: no "R-hat across the 50 runs" — the runs are
*supposed* to differ (that spread IS imputation uncertainty), so a cross-run
convergence statistic would flag healthy behaviour as failure.

---

## Problem 4 — Users can't change the priors (Ayumi's point 4)

**The problem.** Every imputation model uses one hard-coded prior; `prior = ...` is
silently swallowed. Ayumi showed the default can let the phylogenetic variance
drift on small/sparse data, and a better prior fixed it. Our tutorial even tells
users to "adjust priors" — currently impossible. Bonus bug: the two phases disagree
on the categorical fixed-effect prior (`R/bace_final_imp.R:175` hard-codes what
`bace_imp` reads from an option).

**The fix.** New `priors` argument: a named list, one prior per trait, defaults
used for any trait not named. Guard rails from the design review:

- Users may replace the G (random-effect), R (residual) or B (fixed-effect) blocks
  whole — no element-wise merging (MCMCglmm's blocks are internally coupled and
  partial merges produce undiagnosable errors).
- **Overriding R is refused for ordinal/categorical traits** — the prediction math
  hard-codes residual variance 1 (probit scale); changing it would silently corrupt
  every probability. The error message says exactly that.
- A one-vs-rest categorical trait reuses its single prior for all its binary
  sub-models.
- Unify the phase inconsistency; warn on unknown arguments in `bace()` so typos
  stop vanishing.
- **Defaults are NOT changed now** — expose the knob first, revisit defaults later
  with benchmarks (Ayumi's suggested order; correct).

---

## Problem 5 (new, worst) — Categorical imputation silently ignores the phylogeny

**The problem.** For unordered categorical traits (the default path), the code
looks up each species' phylogenetic effect under the wrong column name — MCMCglmm
names them `trait.level.Species.tip`, the code builds `trait.level.tip` — so the
lookup never matches and every species effect becomes 0. Verified two ways (code
read + toy fit): predictions are identical to a no-phylogeny model to 2e-16. The
model *estimates* the effects fine; prediction just drops them. Sibling bug: with
`species = TRUE`, the second species effect is likewise never added
(ordinal AND categorical).

**Consequence that outlives the fix:** the 2026-08 categorical simulation numbers
were produced with the phylogeny silently off — stale by construction; the
categorical cells must be re-run for the manuscript.

**The fix.** Match the actual column names by prefix (`trait.level.Species.`) and
look species up against the suffix — exact for tip labels containing dots, immune
to MCMCglmm's node-effect columns, handles repeated species rows, and extends
naturally to the second species effect. Plus a tripwire: if effects exist but zero
species match, warn loudly — that's the condition that stayed silent for a year.
Locked in by tests comparing against hand-computed predictions AND one real-fit
test (mocks alone could never have caught a naming drift).

---

## Order of work

0. Copy this write-up into the repo as `.agents/plans/fixes-2026-10-ayumi.md` and
   maintain it there (Dan: plans always live in the repo). Save that preference to
   auto-memory.
1. **Problem 5** (phylogeny lookup) — worst error, self-contained, unblocks 2.
2. **Problem 4** (priors) — independent; unblocks Ayumi's prior pilots.
3. **Problem 1** (sweeps + seed + fallback fix).
4. **Problem 2** (posterior-averaged probabilities, on top of 5's corrected code).
5. **Problem 3** (diagnostics incl. cutpoints + final models).
6. Docs/NEWS/tutorial corrections; `.agents/` investigation note; draft PR #15
   reply for Dan to review and post (confirmed 2026-10-26).

## How we prove it

- **Reproduce first:** small script per problem showing the broken number before
  the fix (same discipline as the 2026-08 investigation).
- **New unit tests** (all tiny-MCMC or mock-based, ~15 s added): predictor values
  must differ across final runs with two missing traits (the exact SD = 0
  signature); fit count = runs × sweeps × traits, collapsing to 1 sweep when only
  one trait is missing; same-seed runs identical, serial ≡ parallel; categorical
  predictions must use the species effects (hand-computed + real-fit tests);
  cutpoint rows present in diagnostics; a supplied prior must reach the model;
  forbidden R override must error; final-phase prior option no longer hard-coded.
  Existing suites that test plumbing with several missing traits get `sweeps = 1`
  explicitly so the suite stays ~2 min.
- **Calibration study** (new `dev/15_sweeps_calibration.R`): both traits 30%
  missing, MCAR and MAR, sweeps ∈ {1, 2, 3, 5, 8}; accept the default k if
  coverage ≥ 0.93 and within 0.01 of the k = 8 reference, between-imputation
  variance ≥ 85% of reference. The k = 1 column doubles as the published
  quantification of the old bug. Also: warm-up-pass shortening must match
  full-length passes within Monte Carlo error, else lengthen it.
- **Targeted re-runs (Dan: confirmed 2026-10-26):** Study A full (cheap) as a
  no-regression guard; Study B categorical cells regenerated (stale by
  construction). Full Study B grid stays on the manuscript roadmap.
- **Gates:** `devtools::test()` all green, `R CMD check --as-cran` clean, rebuild
  tarball, fresh win-builder round before CRAN submission.

## Files touched

- `R/model_functions.R` (problems 2, 3, 5; prior guards)
- `R/bace_final_imp.R` (problems 1, 4; fallback; gelman unification)
- `R/bace_imp.R`, `R/bace.R` (thread `priors`, `sweeps`, `seed`; unknown-arg warning)
- `tests/testthat/` (new tests above; `sweeps = 1` in plumbing suites)
- `dev/15_sweeps_calibration.R` (new), `dev/simulation_results/SIMULATION_REPORT.qmd`
- `NEWS.md`, `man/` (re-document), `vignettes/bace.qmd` (prior/uncertainty text),
  `.agents/` (plan copy, investigation note, roadmap/project-state updates)

## What this means for Ayumi's tutorial

Her bird example is ordinal (threshold path), so the phylogeny bug does NOT affect
her headline results; but her mixed-variable section (problem 1) and plug-in
probability figure (problem 2) will need a re-run on the fixed package, and her
prior pilots become actionable once `priors` exists. All of this goes in the
drafted PR #15 reply.
