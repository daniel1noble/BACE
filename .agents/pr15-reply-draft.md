# Draft reply for PR #15 (for Dan to review, edit, and post as a comment)

Audience: Ayumi, on the PR thread. Paste everything below the line.

---

Ayumi — thank you for this. The bird-migration example is a serious piece of
work, and all four of your points were specific and reproducible. We checked
each one against the code before deciding anything. The short version: **all
four are confirmed, and all four are now fixed on `main`.** Your report also
led us to a fifth bug none of us had caught.

**Point 1 — final runs reuse the same point imputations.** Confirmed. We
reproduced your SD = 0 observation directly: the first-listed variable's model
saw identical predictor values in every final run. The fix: each final run now
cycles through the variables several times (new `sweeps` argument, default 3),
drawing random values on every pass. By the second pass, every variable is
imputed against fresh draws of the others, so the uncertainty propagates and
formula order no longer matters. The early passes use shorter chains, so the
default costs about 1.4x the old single pass. When only one variable is
missing, BACE skips the extra passes (one pass is exact there — which is also
why our own single-missing-variable simulations never caught this). There is
also a new `seed` argument that makes the final runs exactly reproducible, in
serial or parallel.

We also ran a calibration study to choose the number of passes rather than
guess it: 100 simulated datasets, two variables each 30% missing, the same
data imputed with 1, 2, 3, 5 and 8 passes. Three passes matches an
eight-pass reference on every metric at about half its cost, so 3 stays the
default. The one-pass arm measured the old behaviour, and it turned out to
be worse than an uncertainty problem: the kept draws of the first variable
were generated against other-variable values that were then replaced, so the
imputed pairs did not belong together. In this setting that biased the
downstream regression slope upward by about 15% of its true value, on top of
understating the between-imputation spread. Three passes removes both
problems.

**Point 2 — plug-in probabilities.** Confirmed, and it affected binary and
multinomial too, not just ordinal. Your 0.901-vs-0.848 overconfidence is
Jensen's inequality. `prob_preds` are now computed the honest way: probability
per posterior draw, then averaged. The imputed values themselves were always
proper draws and are unchanged. `prob_preds` is now documented, including how
to pool it across runs.

**Point 3 — cutpoints missing from diagnostics.** Confirmed — and it was
worse: the health check never ran on the final-phase models at all, and those
are the models the results come from. Your Furnariidae cutpoint ESS of 5–22
is exactly what it should have flagged. `bace_final_imp()` output now carries
`$diagnostics`: every final model is checked (cutpoints and one-vs-rest
sub-models included), summarised as the worst case across runs, and only
flagged parameters are printed. We deliberately did not add a cross-run
R-hat: the final runs are supposed to differ — that spread IS the imputation
uncertainty — so a cross-run convergence statistic would flag healthy
behaviour.

**Point 4 — priors.** Confirmed: `prior = ...` was silently discarded, and
the two phases even used different categorical fixed-effect priors. There is
now a `priors` argument on `bace()`, `bace_imp()` and `bace_final_imp()`: a
named list, one element per variable, giving any of the MCMCglmm `G`/`R`/`B`
blocks. It is validated before any model runs. Your half-Cauchy pilot becomes:

```r
priors = list(habitat = list(G = list(list(V = 1, nu = 1,
                                           alpha.mu = 0, alpha.V = 25))))
```

One guard rail: `R` cannot be changed for ordinal/categorical variables,
because the probability calculations rely on the fixed latent scale. We
followed your suggested order — expose the knob first, revisit the default
only with benchmarks. Your prior-sensitivity pilots would be a great basis
for that follow-up.

**The fifth bug.** While verifying your points we found that the default
multinomial path looked species effects up under the wrong column name, never
found them, and silently used zero — categorical imputation was using no
phylogenetic information at all. Your example is ordinal, so **your headline
results are not affected** (the ordinal lookup was correct). Fixed, with a
warning that fires if the naming ever drifts again, and a real-fit regression
test. With `species = TRUE`, the second species effect is now also included
in ordinal and categorical predictions (it was silently dropped too).

**What this means for the tutorial.** Two sections will want a re-run on the
fixed package: the mixed-variable migration + HWI section (point 1 — expect
the 0.061 vs 0.054 RPS gap to shrink, and your SD = 0 observation to
disappear) and the plug-in probability figure (point 2 — `prob_preds` are now
the honest numbers, so your manual averaging workaround can go). The neutral
descriptions of both behaviours can then become "fixed in 0.1.0" notes. Happy
to coordinate on the re-run.

Full verification write-up: `.agents/investigation-2026-10.md`. Test suite:
1418 passing. CRAN submission stays on hold until the validation runs are
green.

Thank you again — this is exactly the kind of review that makes the package
defensible.
