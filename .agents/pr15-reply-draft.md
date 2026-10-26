# Draft reply for PR #15 (for Dan to review, edit, and post as a comment)

Audience: Ayumi, on the PR thread. Paste everything below the line.

---

Ayumi — thank you for this. The bird-migration example is a serious piece of
work in its own right, and the four implementation points you raised were all
specific, evidenced, and reproducible. We went through each one against the
code before deciding anything, and the short version is: **all four are
confirmed, and all four are now fixed on `main`**, together with a fifth bug
your report led us to that none of us had caught.

**Your point 1 (final runs condition on fixed point imputations).** Confirmed
exactly as you described — we reproduced your SD = 0 observation with a spy on
the model fits: the first-listed variable's model saw identical predictor
values in every final run. The cause: one pass per run, in formula order, from
the chained phase's point-imputed dataset. The fix is the one you hinted at:
each final run now performs `sweeps` stochastic passes over the variables
(default 3, `sample = TRUE` throughout), so every variable — including the
first-listed one — conditions on fresh posterior-predictive draws of the
others. Warm-up passes use a shortened post-burn-in tail, so the default costs
~1.4x the old single pass rather than 3x, and when only one variable is
missing it collapses to a single pass (exact there, and why our y-only
simulations never caught this). Formula-order dependence disappears with the
passes. There's also a new `seed` argument giving per-run RNG streams, so
results are reproducible and identical in serial and parallel. A calibration
study (`dev/15_sweeps_calibration.R`, both variables 30% missing, k in
{1,2,3,5,8}) is running to confirm the default; the k = 1 arm doubles as the
quantification of the old behaviour.

**Your point 2 (plug-in `prob_preds`).** Confirmed, and broader than ordinal —
binary, multinomial and OVR all plugged posterior-mean coefficients into the
probability formulas. Your 0.901-vs-0.848 overconfidence is exactly Jensen's
inequality at work. The forward prediction paths now average the per-draw
probabilities over the stored posterior (the imputed category draws themselves
were already proper single-draw samples and are unchanged). `prob_preds` is
now documented, including how to pool it across runs.

**Your point 3 (cutpoints missing from diagnostics).** Confirmed — and worse
than you reported: the health check also never ran on the final-phase models
at all (the ones actually pooled for inference; your Furnariidae cutpoint ESS
of 5–22 is precisely what it would have flagged). `bace_final_imp()` results
now carry a `$diagnostics` summary covering fixed effects, variances,
cutpoints, and one-vs-rest sub-models, aggregated as worst-case across the
runs, with only flagged parameters printed. We deliberately did NOT add a
cross-run R-hat: the runs target different conditional posteriors, so
between-run spread is imputation uncertainty doing its job, not
non-convergence.

**Your point 4 (priors not settable).** Confirmed (`prior = ...` silently
vanished into `...`; the two phases even disagreed on the categorical
fixed-effect prior). There is now a `priors` argument on `bace()`,
`bace_imp()` and `bace_final_imp()`: a named list, one element per variable,
supplying any of the MCMCglmm `G`/`R`/`B` blocks, which replace the defaults
wholesale, validated before any MCMC runs. Your half-Cauchy pilot becomes:
`priors = list(habitat = list(G = list(list(V = 1, nu = 1, alpha.mu = 0,
alpha.V = 25))))`. One guard rail: `R` cannot be overridden for
ordinal/categorical variables, because the prediction math hard-codes the
fixed latent scale. We followed your suggested ordering — expose the knob
first, revisit the default only with benchmarks; your prior-sensitivity
pilots would be a great basis for that follow-up, and `bace()` now warns on
unknown arguments instead of swallowing them.

**The fifth bug (found while verifying yours).** The default multinomial path
was looking species effects up under the wrong MCMCglmm column name
(`trait.level.<tip>` instead of `trait.level.Species.<tip>`), never matching,
and silently using zero — i.e., non-phylogenetic imputation for categorical
traits. Your bird example is ordinal, so **your headline results are not
affected** (the threshold lookup was correct). It's fixed with exact prefix
matching, a tripwire warning if naming ever drifts again, and a real-fit
regression test; with `species = TRUE` the second species effect is now also
included in ordinal/categorical predictions (it was silently dropped too).

**What this means for the tutorial.** Two sections will need a re-run on the
fixed package once the calibration lands: the mixed-variable migration + HWI
section (point 1 — expect your 0.061 vs 0.054 RPS gap to shrink and the HWI
SD = 0 observation to disappear) and the plug-in probability figure (point 2 —
`prob_preds` are now the honest probabilities, so the plug-in/averaged
contrast becomes historical). The neutral descriptions of both behaviours can
then become "fixed in 0.1.0" notes, which is a nicer story to publish. Happy
to coordinate on the re-run — the scripts in the PR look like they'll run
unchanged apart from dropping your manual posterior-averaging workaround.

The full write-up (verification evidence, line references, test list) is in
`.agents/investigation-2026-10.md`, and the running plan in
`.agents/plans/fixes-2026-10-ayumi.md`. Suite is at 1418 passing tests; CRAN
submission stays on hold until the validation campaign is green.

Thank you again — this is exactly the kind of review that makes the package
defensible.
