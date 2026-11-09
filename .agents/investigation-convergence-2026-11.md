# The convergence check was broken — investigation and repair (2026-11-09)

**Status: fixed.** The chained-equations convergence check in
`R/convergence_diagnostics.R` was replaced wholesale; `bace()`'s retry ladder
now fires only when retrying can help. This note is the full record: what we
saw, how we proved it, what the literature says, what the new check is, and
the calibration numbers behind its thresholds. Written in plain language; the
numbers are all reproducible.

## What we saw (the symptom)

Dogfooding the 2026-11 simulation, BACE said "not converged" on essentially
every dataset:

- 20-replicate test: **2 of 20** converged.
- Chain-length pilot: **0 of 17**, and tripling the chained runs (5 → 15) and
  tripling the final chain (50k → 150k) moved that *not at all*.
- Yet interval coverage was **nominal (0.95)** throughout, and across the 20
  replicates the verdict had **no relationship** with bias, coverage,
  interval width or imputation accuracy (|bias| 0.125 not-converged vs 0.178
  converged; width 0.637 vs 0.639; no significant correlation with the
  number of "converged" variables, signs inconsistent).

So the flag was both almost-always-on and uninformative. That smelled like a
broken test, not a broken sampler.

## The proof (mechanism, not assertion)

Feed the check a series that is **converged by construction** — pure noise
around a fixed value, no drift, no trend, at the Monte-Carlo noise floor of
a mean over 30 imputed cells. Pass rates at campaign settings (runs = 5,
4000 simulations per cell):

```
trait's mean:   0.00   0.25   0.50   1.00   2.00   5.00
%-change pass:  0.000  0.018  0.115  0.460  0.922  1.000
ACF pass:       0.369  0.378  0.393  0.381  0.381  0.391
trend pass:     0.954  0.947  0.954  0.943  0.951  0.949
whole dataset:  3.4%  (needs >= 5 of 6 variables)
```

Three root causes, each verified in code and numerically:

1. **The %-change criterion divided by the variable's own level**
   (`abs(diff(s) / (s[-n] + 1e-10))`). A mean-centred trait has a
   denominator near zero, so the "percentage" explodes and the criterion can
   never pass; a trait at level 20 drifting 0.5 units per run passes it 100%
   of the time (0.5/20 = 2.5% < 5%). The criterion measured *where a
   variable sits*, not *whether it settled* — wrong in both directions.
2. **The ACF criterion** correlated 4 points against a fixed |r| < 0.3
   cutoff: it rejected ~60% of perfectly settled chains (the sampling SD of
   a correlation over 4 points is ~0.5).
3. **Geweke was silently disabled** (hard-coded `NA`, and its own window
   logic needs ~20+ points anyway), leaving 3 live criteria of which one was
   structurally unpassable. The per-variable rule (>50% of live criteria)
   plus the overall rule (>70% of variables) then gave the ~3% dataset pass
   rate — matching the observed 2/20 and 0/17.

Worse findings from the full audit (4-agent review, 2026-11-09):

- **A zero-variance chain was special-cased as "converged".** A collapsing
  imputation spread is the *textbook* non-convergence (the locked-up sampler
  of van Buuren & Groothuis-Oudshoorn 2011, p. 37). Exactly backwards.
- **An unrecognised variable type fell through silently**, leaving an
  all-zero series that auto-passed every criterion.
- **The alternative methods were worse.** Energy and Wasserstein computed
  their distances correctly but thresholded the *coefficient of variation of
  the step size*: a steady drift has constant steps, so both methods passed
  drifting chains ~98–100% of the time, while coin-flipping (energy 54–62%,
  Wasserstein 20–24%) on settled ones. Energy also paired imputed cells
  across unrelated taxa and used an unseeded `sample()` (non-reproducible).
- The overall 70% rule was non-monotone in the number of variables
  (1–3 variables needed 100% agreement, 7 needed 71%).
- `bace()`'s retry ladder (5 → 8 → 11 runs ≈ 4.8× compute) therefore fired
  on nearly every real dataset, and no amount of chain could ever satisfy
  the check — the retries were pure waste, ending in a false warning.

## What real chains look like (why the fix is shaped the way it is)

We harvested real chained series (3 datasets × 15 runs at campaign MCMC
settings, all six trait types). Three facts drove the design:

1. **Warm-up is real.** Early runs can sit far from the settled region
   (a count trait's mean was 12–17 in runs 2–4, then ~2.9 ever after). The
   old check had no warm-up notion; trend/ACF saw the warm-up too.
2. **Mild alternation is normal.** Nearly every real series — including
   visibly healthy ones — has negative lag-1 autocorrelation: chained
   equations condition each run on the previous values of the other
   variables, and quantised proportions (1/30 jumps) can hit acf = −1 while
   moving only at the noise floor. On datasets with exactly this structure,
   measured coverage was still 0.95. So alternation *per se* must not fail
   a dataset; only **loud** flip-flops (movement well above the noise floor)
   are a pathology.
3. **The right yardstick is the imputation spread.** Under settledness, the
   mean over n missing cells moves by about spread × sqrt(2/n) between
   runs — a floor you can draw on a plot and test against.

## The new check (`method = "stability"`)

Per variable, watch the **mean and SD of the imputed values** across chained
runs (mice's chain means/variances; van Buuren 2018 ch. 6.5.2). Discard the
first half as warm-up. Measure everything in units of the variable's own
imputation spread. Flag three named pathologies, each with its remedy:

| verdict | meaning | remedy |
|---|---|---|
| drift | late chain still travelling (net movement > 0.75 spread units with Kendall trend support) | more chained runs — the one curable pathology; `bace()` retries only for this |
| oscillation | loud flip-flop: lag-1 acf < −0.6 AND mean step > 3 noise floors | more runs won't help; look for variables feeding back on each other |
| stuck | imputation spread collapsed below 25% of its early value | check that variable's model (priors, separation, near-constant values) |

Everything else is "ok", with advisory *notes* (mild alternation; identical
imputations every run) for the eye. **A pass means "no pathology detected",
not a certificate** — the output and docs say so explicitly, because no
clear-cut convergence test exists for chained equations (van Buuren 2018
ch. 6.5.2, quoting Cowles & Carlin 1996 that automated convergence
monitoring "is unsafe and should be avoided"), and the one empirical study
(Oberman, van Buuren & Vink 2021, arXiv:2110.11951) found valid inference
after 5–10 iterations while formal diagnostics kept "improving" to 20–30 —
exactly our coverage-vs-flag pattern. Multi-chain R-hat (as in the `mi`
package, Su, Gelman, Hill & Yajima 2011, JSS 45(2)) is the reviewer-friendly
upgrade path, but it needs several chained sequences from dispersed starts —
an architectural change, parked as a possible future feature.

## Calibration (thresholds are measured, not hand-picked)

Constructed chains, 1500–4000 per cell, fixed seeds; full harness in the
session scratchpad, operating characteristics pinned in
`tests/testthat/test-convergence_diagnostics.R`:

```
runs = 10 (bace() default)      pass/catch rate
  settled (mean 0 and 20)        0.995 - 0.999
  geometric approach (settles)   0.993 - 0.999
  constant                       1.000
  steady drift 0.5 spread/run    0.994 - 0.995   -> drift
  loud flip-flop (±0.75 spread)  0.996 - 0.999   -> oscillation
  collapsing spread              1.000            -> stuck
  scale invariance               verdicts identical at mean 0 vs 20

runs = 5 (campaign)              settled pass 0.99; drift power ~0.75;
                                 oscillation not assessable (too few kept
                                 points) - the output says so
runs = 15                        everything >= 0.999
```

Real chains (18 variable-series from the 3 harvested fits, whose downstream
coverage was nominal): **all pass**, with notes on 2 (quantised alternation;
a constant binary) — exactly the desired behaviour.

## The adversarial review (2026-11-09, same day) — 12 must-fix defects

The first version of the replacement was handed to three independent reviewers
briefed to **break** it (hostile inputs / hostile statistics / seam hunting).
They found twelve must-fix defects, all with reproductions. Every one is fixed
and pinned by a regression test. The important ones, because they say what a
convergence check is easy to get wrong:

1. **`sd()` without `na.rm`** → a single NA imputed cell in the last two runs
   made the spread read 0 → false "stuck", which `bace()` treats as uncurable.
2. **Categorical variables were watched on the first level only** → a chain
   flipping every cell between levels 2 and 3 left level 1 flat and scored a
   clean pass, with the note "identical imputations in every run". Now every
   level is monitored and the variable takes the worst level's verdict.
3. **`sqrt(p(1-p))` is not a spread measurement**, it is an arithmetic
   function of the mean. Reading it as one made (a) a settled binary over few
   cells false-"stuck" up to 16% of the time (one run's cells all landing in
   one level gives exactly 0), and (b) a genuinely *drifting* proportion
   report as "stuck" a quarter of the time — the wrong remedy, and it
   suppressed the retry that would have helped. Proportions are now measured
   in proportion units with no spread check, and their gates scale with the
   sampling noise of a share over that many cells.
4. **The drift gate's trend test could veto any amount of movement**: a chain
   jumping 500,000 spreads in its last run passed, because ties made Kendall's
   p = 0.14. There is now an unconditional escape at 3x the threshold, plus a
   "wander" catch for sustained large steps that go nowhere (a 6-spread wave
   returning home also passed).
5. **The noise floor shrank with the number of imputed cells but real
   between-run movement does not.** At >= ~350 cells the amplitude gate stopped
   screening and a 5-pair correlation decided alone: false oscillation 10-27%
   at n_cells 500-1000 — and the campaign's own n = 500 scenarios sit at 150.
   Fixed with an absolute amplitude floor; measured false rate now <= 3.8% at
   the most pessimistic noise level we have observed, <= 0.7% at the realistic
   one, with true oscillation still caught >= 98% at n_cells >= 150.
6. **Spread was only ever watched downward**: an exploding imputation variance
   passed silently. New verdict **diverging**.
7. **The constant-chain branch used an absolute 1e-12 tolerance**, so a clean
   0 -> 1e-13 drift was reported as "identical imputations". Now relative.
8. **"unassessable" counted as a pathology** → `converged = FALSE`, print said
   "Pathologies detected: YES" and `bace()` warned about an uncurable problem
   when nothing had been assessed. Now `NA`, as the docs always promised.
9. **Short chains passed silently** where the tests structurally cannot run
   (drift needs >= 5 runs, flip-flop >= 7). They now say so in a note.
10. **Four plot crashes** on states the check deliberately tolerates: a
    constant series (acf of zero variance → infinite ylim), an all-NaN series
    (the default trace plot, reachable from `bace(plot = TRUE)`), one NA run
    (lowess), and a non-matching `variables` argument.
11. **The bird_social case study read `stationarity_tests`**, removed with the
    old check — it would have failed at render. Rewritten against the new
    fields.
12. **Roxygen `\%` escaping truncated two help pages mid-sentence** under
    `markdown = TRUE` (the re-escaped `%` starts an Rd comment).

Also fixed: "pre-1.1" wording everywhere (this ships as **1.0.1**), the stale
README feature list, the main vignette's `skip_conv` description and its two
`method = "summary"` examples, a NEWS claim that the output says something it
did not, and `min_iterations` validation.

**What this says about the first attempt.** The design was right — scale-free,
warm-up aware, pathology-named — but five of the twelve defects were the same
mistake in a new place: assuming a quantity is a measurement when it is an
artefact of another one (sqrt(p(1-p)) for proportions, the first level as a
proxy for a categorical, a shrinking noise floor as a proxy for real movement).
That is exactly the family of error that produced the original bug, which is
worth remembering next time.

## API changes (all loud, pre-CRAN so cheap)

- `assess_convergence()`: default method `"stability"`; `"summary"` is a
  deprecated alias; the per-criterion and energy/wasserstein/all methods are
  **defunct with explanatory errors**; `pct_change_threshold`, `lag`,
  `alpha` removed; `drift_threshold` (0.75) and `stuck_ratio` (0.25) added;
  returns `verdicts`, `notes`, `sd_stats` alongside the old fields.
- `bace()`: retries **only on drift**; oscillation/stuck are reported with
  remedies and the run proceeds (MI propagates the remaining disagreement
  into wider intervals — that is what our 0.95 coverage under the old
  "failures" was showing).
- `plot()`: type `"steps"` replaces `"pct_change"` (steps vs noise floor,
  warm-up shaded); energy/wasserstein plot helpers defunct.
- Unknown variable types now error loudly in the assessment.
- Verdicts are `ok` / `drift` / `oscillation` / `stuck` / `diverging` /
  `unassessable`; `converged` is `NA` when nothing could be assessed.

## Consequences elsewhere

- **Simulation campaign**: the convergence-rate columns in older results
  (28%/51%/82% by n in the 2026-10 campaign; 2/20 and 0/17 this week)
  measured the bug, not the sampler. The chain-length pilot built to raise
  those rates is retired; its real question (is runs = 5, nitt = 50k
  enough?) was answered by coverage: yes.
- **Vignettes**: bace.qmd's callout advising users to route around false
  flags with `skip_conv` is rewritten. The bird_social vignette has a
  paragraph explaining away a false flag and cached fits made under the old
  check — needs a re-run, tracked separately.
- The old criteria produced no wrong *published* numbers: coverage was
  nominal throughout (the flag never fed inference; only compute and
  messaging).

## Test coverage added

`tests/testthat/test-convergence_diagnostics.R` now pins, with fixed seeds:
verdicts on every constructed truth state at means 0 and 20; measured
operating characteristics (settled ≥ 97% pass, drift/osc ≥ 95% catch,
stuck = 100%, verdict distribution identical under a +20 shift — the exact
property the old check violated); hand-computed statistics on a
deterministic series; and the break battery (all-NA series, Inf, single
imputed cell, one-bad-variable-names-it, all-one-level binary flips,
near-constant series just above the old 1e-10 variance cutoff). The old
suite's convergence tests were smoke tests — the one named for correctness
asserted only the S3 class, which is how this shipped.
