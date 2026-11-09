# Numerical accuracy audit — November 2026

**Date:** 2026-11-10
**Why:** after the convergence repair landed, the question was no longer "does it
run?" but "are the numbers right?". This is the record of a full numerical audit
of the package and what it found.

Short version: **eight silent bugs**, all of the same kind — code that returns an
ordinary-looking number that is wrong, with no error and no warning. Two of them
would have changed results in the manuscript. All are fixed, each with a test that
pins the actual number rather than just checking the code runs.

Companion notes: [investigation-convergence-2026-11.md](investigation-convergence-2026-11.md)
(the convergence rewrite), [investigation-2026-08.md](investigation-2026-08.md)
(the improper-MI defects).

---

## How the audit was run

Six independent reviewers, each given one lens and told to find numbers that are
wrong rather than code that crashes:

1. scales, links and back-transformation
2. the pooling mathematics
3. whether the multiple imputation is genuinely proper
4. whether the phylogeny is actually used, for every trait type
5. type detection and the prediction core
6. fresh eyes on the newly rewritten convergence code

Rules: read-only, and **no finding counts unless you ran code and printed the wrong
number**. Every claim below was then re-checked by hand before anything was changed.

Separately, the whole test suite was run with R's own silent-mismatch detectors
armed (`warnPartialMatchDollar`, `warnPartialMatchArgs`, `warnPartialMatchAttr`).

---

## What was found

### 1. A categorical level with a dot in its name could never be imputed

The worst one, because it is completely invisible.

MCMCglmm names a trait `trait<response>.<level>`. With levels
`forest | forest.edge | open` the trait is `traithabitat.forest.edge`. The code
labelled the probability column by taking **everything after the last dot**, so
that column came out labelled `edge`.

The probabilities were right. Only the label was wrong. But `.predict_bace()`
builds its sampling vocabulary with
`intersect(levels_var, colnames(probs))`, so the mislabelled column was **dropped**
and the remaining probabilities renormalised.

Measured: across 210 posterior-predictive draws, `forest.edge` was imputed **0
times**, although it was 35% of the observed data.

After the fix, dotted and plain level names give **identical** counts (68 / 28 / 44
either way), which is the right outcome — a level's *name* should never affect the
imputation.

Why nothing caught it: the `sample = FALSE` path escapes, because `.impute_levels()`
falls back to positional alignment when the names do not match. So the convergence
phase looked perfectly healthy.

Fix: strip the `trait<response>.` **prefix** (non-greedy), the rule `.pred_cat()`
already used. Plus a tripwire — a probability column whose name is not a declared
level now warns instead of silently deleting a category.

### 2. Multinomial heritability was inflated, and worse with more categories

`.get_multinom_diag_cols()` picked the residual diagonal with a pattern that
required a **numeric** trait index. MCMCglmm labels traits by factor level
(`traity.B`), never by number, so the match always failed on real fits and a silent
fallback took the first J columns. In MCMCglmm's column-major order those are
column 1 of R — one variance plus J−1 **covariances**.

R is fixed at `(I + 11')/K`, so that sum is exactly 1 for every K:

| K | true `tr(R) = 2(K−1)/K` | what BACE used |
|---|---|---|
| 3 | 1.333 | 1.000 |
| 4 | 1.500 | 1.000 |
| 5 | 1.600 | 1.000 |
| 6 | 1.667 | 1.000 |

So the residual trace was understated by 25% at K=3 and 40% at K=6, and reported
H2 grew too large with the number of categories. For a true H2 of 0.50, BACE
reported 0.571 (K=3), 0.600 (K=4), 0.615 (K=5).

**The package's own test passed because it fed in made-up numbered labels that
MCMCglmm never produces.** A test pinning a fiction is worse than no test.

Fix: match the diagonal by splitting the name on `:` and comparing the two trait
identifiers — which also copes with dotted level names. The no-columns fallback
changed from `J²/K` to `2J/K` (they are equal only at K=3, by coincidence).

### 3. The reported posterior SD was not a standard deviation

The gaussian back-transform multiplied the whole prediction frame by `sd` and added
`mean`, so the location shift reached the SD column too. But `sd(a + bX) = |b|·sd(X)`
— there is no `a` term.

On a variable with mean 1000 and sd 0.01 the reported SD was **1000.003** instead of
0.0031, and it was **negative** for any variable with a negative mean.

Imputed values were never affected (only the mean column is consumed downstream),
but `post_sd` is handed back to the user and nothing tested it.

### 4. `pool_mi()` reported missing information where there was none

`B = 0` (every imputation agreeing exactly — which happens for real whenever the
analysis model uses only fully observed variables) and `W = 0` shared a fallback of
`M − 1`. Their limits differ:

- as `r → ∞` (no within-imputation variance), df → `M − 1` ✔
- as `r → 0` (no between-imputation variance), df → **infinity** ✘

So `fmi = (r + 2/(ν+3))/(r+1)` reported `2/(M+2)` of the information missing when
none was. Now returns `df = Inf`, `fmi = 0`, `riv = 0` — the limit `mice` reaches by
clamping lambda (df ≈ 9e8, fmi ≈ 2e-9).

### 5. A character vector of formulas lost every response but the first

`fixformula = c("y ~ x1", "z ~ x1")` — using `c()` instead of `list()` — is not a
list, so it took the single-formula branch. There `stats::as.formula()` keeps only
the first element (R warns only that `formula()` is deprecated for length > 1) and
auto-expands from it. `z` disappeared entirely.

Worse: the guard that checks every variable with missing data has a formula is
inside `if (is.list(fixformula))`, so it was **skipped on the same path**. Nothing
complained.

A character vector of several formulas is now treated as the list the user meant, in
`bace()`, `bace_imp()` and `bace_final_imp()`.

### 6–8. Found earlier the same day

- **`with_imputations()` paired fits with the wrong tree.** When more trees are
  supplied than imputations — the normal case, e.g. 100 posterior trees and 50
  imputations — `bace_final_imp()` spreads across trees 1, 3, 5, … 99 while the
  downstream fits used 1, 2, 3, … 50. 49 of 50 models conditioned on a tree that
  had not generated their data. The guard compared only tree *counts*, so passing
  the same 100 trees matched on count and said nothing. Now read from the recorded
  `$tree_index` and checked by content hash.
- **`pool_mi(df_fun=)` silently skipped the small-sample correction** for `gls`,
  `phylolm` and `MCMCglmm`, whose `df.residual()` is `NULL` (issue #17).
- **`$` partial matching** in the new `.complete_data_df()`: `list(nitt = 13000)$n`
  returns `13000`, so a chain length could be read as a sample size.

---

## What was checked and found correct

Worth recording, so nobody re-audits it:

- Gaussian values round-trip **bit-identically** through standardisation
  (max error 0.000e+00 at scales from sd 0.01 to sd 250).
- Imputations match the closed-form conditional mean of the true model
  (r = 0.976–0.987, regression slope 1.00, RMSE at the theoretical optimum).
- The threshold, one-vs-rest and multinomial probability formulas each reproduce an
  independent re-derivation to **0.000e+00**; probabilities sum to 1 within 2e-15.
- Phylogenetic effects genuinely enter predictions — zeroing the species BLUPs moves
  probabilities by up to 0.67. (The 2026-10 silent-zero bug has **not** returned.)
- Across ~370 known-truth replicates, **no trait type shows detectable bias**:
  gaussian −0.018 SD units (95% CI −0.051 to +0.016), poisson +2.30% (95% CI −2.74%
  to +7.33%). Per-cell coverage 0.9545 against a nominal 0.95.
- Reported class probabilities are calibrated: binary slope 0.994 (95% CI
  0.908–1.079), categorical 0.969 (0.888–1.049) — both indistinguishable from 1.
- Poisson draws are genuine posterior-predictive draws, not `round(exp(eta))`.
- Factor round-trip preserves class, level set, level **order** and the ordered flag.
- **No `$` partial-match sites anywhere in BACE's own code.** All 4,239 warnings
  raised by R's detector come from inside MCMCglmm (`GRprior[[i]]$n` matching `nu`,
  and `print.summary` matching `Gcovariances`). `coda` and `ape`: none.
- `.complete_data_df()` agrees with `mice:::get.dfcom()` on `lm`, `glm`, `gls`,
  Brownian `gls` and `phylolm` — and both equal n − p.

---

## The recurring pattern

Five of the eight are the same mistake in different clothes: **a derived quantity
treated as if it were a measurement.** A label re-derived by string surgery instead
of carried forward. A tree pairing re-derived instead of read from `$tree_index`. A
degrees-of-freedom re-derived from whatever field happened to match. The data was
already there and correct; the code went and worked it out again, differently.

The second pattern is **tests that pin a fiction.** The multinomial test fabricated
label formats MCMCglmm never produces. A `with_imputations` test asserted the exact
silent-drop behaviour that was the bug. A green suite is not evidence that the
inputs in the test resemble the inputs in the wild.

---

## Second pass: six more lenses, two skeptics per finding

A second audit (37 agents; every claim re-checked by two independent reviewers
who had to reproduce a wrong number before it counted) returned **12 confirmed
findings**. One mattered, and it was a regression in the convergence rewrite
itself.

### The stability check called settled variables "drift" ~85% of the time

To judge whether a variable has stopped moving you need a ruler.
`.conv_assess_series()` used the imputation spread — fine — but when the
imputed cells had **no** within-run spread it fell back to `sd()` of the very
series it was judging. That is circular, and it has a closed form:

```
step_mean = mean|diff(m)| / sd(m)  ->  2/sqrt(pi) = 1.128   for ANY iid wobble
wander gate                         =  1.00
```

So pure noise cleared the gate. Measured on a settled series over 2000
replicates: drift 0.85 at `runs = 5`, 0.81 at 10, 0.88 at **50** — *more runs
made it worse*, because the statistic carries no scale.

Two ordinary situations land there: a variable with exactly **one** missing
cell (SD of a length-1 vector is `NA`, recorded as 0), and any **poisson**
variable whose rounded point imputations land on the same integer. The live
campaign has a count trait, so this was not hypothetical.

Fixed by taking the yardstick from **outside** the series — the spread of the
variable's observed values — and returning `"unassessable"` when no such scale
exists rather than inventing one. Settled series now return `"ok"` 100% of the
time at every `runs`; a genuine 3-data-SD drift is still caught 100% of the
time; series with a real within-run spread are unchanged.

**This is the third appearance of the same defect class in one cycle** — a
derived quantity used as an independent measurement — and this time it was
introduced *by the fix for the previous instance*. Worth remembering: fixing a
defect class does not immunise the fix against it.

### The other eleven, triaged

Scored by whether the live campaign can reach them:

| | finding | reaches the campaign? |
|---|---|---|
| 1 | `phylo_signal_summary(prior=)` un-fixes the probit residual while H² still divides by 1 (H² 0.85 reported vs 0.00 true) | no — campaign never calls it |
| 2 | random slopes in `ran_phylo_form` dropped when `species = FALSE` | no — campaign uses `~ 1 \| Species` |
| 3 | `.pred_threshold()`/`.pred_cat()` add a second residual to an already-realised latent | no — forward path is used, and its error vs reference is exactly 0 |
| 4 | `pool_posteriors()` stacks effects each in that imputation's own units | no — campaign pools with `pool_mi()` on data-scale datasets |
| 5 | standardisation constants not stored on `bace_final`, so no exact back-transform | no — completed datasets are already on the data scale |
| 6 | threshold/ordinal forward prediction drops the phylogenetic effect when the RE carries slopes | no — no slopes in the campaign |
| 7 | a declared-but-unobserved factor level shifts threshold imputations onto the wrong label (accuracy 0.73 -> 0.07) | **possible** — `threshold4` under clade missingness could lose a level |
| 8 | `phylo_signal_summary()` re-derives types differently from `bace()`; a whole-number double is fitted as Poisson | no — campaign never calls it |
| 10 | `plot_steps_convergence()` re-derives the yardstick; overstates steps 5-18x | cosmetic (plot only) |
| 11 | the `wander` gate is unreachable for proportion series (step <= 1 by construction, gate > 1.00) | missed alarm, not a wrong number |
| 12 | `thresholds` reported `osc_floor = 0.50` while 0.60 was applied | reporting only — **fixed** |

Fixed in this pass: **9** (the drift regression) and **12**. **11** is documented
as a known gap rather than guessed at — setting a proportion-scale gate needs
calibration against the harvested chains, and an uncalibrated gate would
reintroduce exactly the false-alarm failure just removed.

**Finding 7 is the one to look at next** for the campaign: if a `threshold4`
trait loses a level entirely under clade missingness, imputations may be
assigned to the wrong category label. It carried one of two skeptic votes, so
it needs reproducing against the campaign's own generator before acting.

## State after the audit

- `devtools::test()` — **1895 pass / 0 fail / 0 error / 0 warnings**.
- `R CMD check` — clean (0 errors, 0 warnings, 0 notes).
- Every fix has a regression test asserting the **number**, not just that it runs.

## Still open — for Dan and the co-authors

These are in the **simulation**, not the package, and they are judgement calls:

1. **The pgls fallback is backwards.** When the optimiser returns a *negative*
   lambda (less structure than zero), `pgls_slope()` refits with lambda **fixed at
   1** — the maximum-structure extreme. It fires in roughly two of every three BACE
   replicates at ps = 0.3.
2. **The arms are not compared like with like.** BACE's intervals use uncapped
   Rubin df; the oracle and complete-case arms use an ordinary n−2 t interval. The
   headline coverage comparison is therefore apples to oranges. Related: `mice`
   applies the Barnard–Rubin correction **by default** and BACE does not — worth
   about 0.5 coverage points at n=100, but **1–2 points on a 30-tip tree** (0.93
   actual against a claimed 0.95).
3. **`collect.R`'s `n` column is always 1**, so no summary row reveals how many
   replicates actually survived.
4. **Per-cell coverage is computed 54,000 times and never read** — stored in each
   `.rds`, used by no summary.
5. **The complete-case comparator deletes on all six traits** although the analysis
   model is only `x2 ~ x1`. Documented design choice, but the cross-mechanism
   pattern in Claim 1 is substantially an artefact of that rule.
6. **The calibration guard-rail is too loose** — `collect.R` warns only outside
   [0.90, 0.99], and 0.90 sits 7.3 sigma below nominal at 1000 replicates.
