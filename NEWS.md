# BACE 1.0.1

* **The chained-equations convergence check was defective and has been
  replaced.** The old check compared each variable's movement between runs to
  the variable's own *level* (`%` change) and its own mean (the energy and
  Wasserstein variants), so its verdict tracked where a variable sat rather
  than whether it had settled: a mean-centred trait could *never* pass
  (division by ~0), while a trait far from zero passed even while steadily
  drifting. Fed chains that are converged by construction, the old check
  passed ~3% of them at default settings; fed steady drifts at level 20, it
  passed essentially all of them. In practice this meant `bace()` silently
  retried (5 → 8 → 11 runs, ~4.8x the compute) on nearly every real dataset
  and then warned that convergence had failed when it had not. Separately,
  a chain whose imputation variance collapsed to zero — the textbook
  non-convergence of van Buuren & Groothuis-Oudshoorn (2011, p. 37) — was
  special-cased as *converged*.

  The new check (`assess_convergence()`, method `"stability"`) watches the
  mean *and* the SD of each variable's imputed values across chained runs,
  treats the first half as warm-up, measures all movement in units of the
  variable's imputation spread, and flags three named pathologies: **drift**
  (still travelling; the one thing more runs can cure — `bace()` now retries
  only for this), **oscillation** (variables flip-flopping in step; retrying
  is pointless and `bace()` now says so and proceeds), and **stuck**
  (collapsing imputation variance) and **diverging** (the spread growing
  without bound). Categorical and ordinal variables are monitored on every
  level and take the worst level's verdict — watching only the first level
  is blind to a chain that flips its cells between the others.

  Calibration on constructed chains (pinned in the test suite): settled
  chains pass >= 99% at `runs = 10` wherever their mean sits and however
  many cells are imputed; drift, loud oscillation, collapse and divergence
  are each caught >= 99%; false oscillation stays under 5% even at 1000
  imputed cells and the most pessimistic between-run noise we have measured.
  Short chains say what they could not do: below `runs = 5` the drift test
  cannot run and below `runs = 7` the flip-flop test cannot, and the output
  carries a note saying so. A pass means *no pathology detected*, not a
  certificate — no clear-cut convergence test exists for chained equations
  (van Buuren 2018, ch. 6.5.2; Oberman, van Buuren & Vink 2021,
  arXiv:2110.11951) — so the trace plots remain the primary diagnostic.

  Breaking changes, all loud: the per-criterion methods
  (`"summary.acf"`, `"summary.percentage"`, `"summary.trend"`,
  `"summary.geweke"`) and the distribution methods (`"energy"`,
  `"wasserstein"`, `"all"`) are defunct with explanatory errors
  (`"summary"` itself is a deprecated alias for `"stability"`); the
  `pct_change_threshold`, `lag` and `alpha` arguments are gone (replaced by
  `drift_threshold` and `stuck_ratio`); `plot()` type `"pct_change"` is
  replaced by `"steps"` (movement against the Monte-Carlo noise floor) and
  the energy/Wasserstein plot helpers are defunct. Fits made with earlier
  versions re-assess cleanly with the new check.

* Fixed a crash in `sim_bace()` when a multinomial predictor's liabilities
  were large (softmax overflow; GitHub issue #19).

* Fixed a crash in the final-model diagnostics when a chain was too short or
  degenerate for `effectiveSize()`/`geweke.diag()`/`acf()`; each is now
  reported as NA instead of aborting the fit (GitHub issue #20).

# BACE 1.0.0

First CRAN release. Earlier versions were distributed from GitHub only. The
1.0.0 number marks a stable user-facing API, not a change in method. Every
change below is additive, with one small exception noted under the per-phase
MCMC settings (`fixformula` can no longer be abbreviated as far as `f`/`fi`).

* **Long runs can now be resumed instead of restarted.** Give `bace()` (or
  `bace_imp()` / `bace_final_imp()`) a `checkpoint_dir` and the convergence
  chain is saved when it finishes, and each final imputation run as soon as it
  completes. Calling again with the same arguments and the same directory
  loads that work and computes only what is missing — so a session that dies
  during the final runs picks up where it left off rather than starting over.
  (The convergence chain is saved once it finishes, so a crash during that
  phase still restarts it.) A resumed run gives
  exactly the same answer as one that was never interrupted, which is why
  `seed` is required. If any setting that would change the results differs
  from the saved work, the resume is refused and the differing settings are
  named, so results from two different analyses are never mixed. Default
  `NULL`, i.e. off, and nothing else changes.

* **The final imputation runs can now have their own MCMC settings.**
  `bace()` gains `final_nitt`, `final_thin` and `final_burnin` (default
  `NULL` = reuse the chained-phase settings, exactly as before). The two
  phases do different jobs — the chained phase only needs point imputations
  to settle, while the final runs supply the draws that are pooled for
  inference — so they often want very different chains, for example
  `nitt = 50000, thin = 40` with `final_nitt = 5e5, final_thin = 450` (scale
  the thinning with the chain: every retained draw is stored for all
  `n_final` runs). Previously this required bypassing `bace()` and calling
  `bace_imp()` and `bace_final_imp()` separately, which is what the bird case
  study had to do. (The amphibian study splits the phases for a different
  reason — memory on an 18,000-row dataset — which this does not change.)
  The new arguments sit at the end of `bace()`'s argument list, so existing
  calls are unaffected; the only change is that `fixformula` can no longer be
  abbreviated as far as `f` or `fi` (`fix` and longer still work).

* **New `bace_validate()`: check how well BACE imputes YOUR data.** The
  natural call is on a finished run — `bace_validate(fit)` reuses the run's
  data, formulas, tree(s) and settings (now stored on the result as
  `$settings` and `$trees_used`) and re-runs the pipeline on freshly masked
  observed cells, scoring the imputations against the hidden truth; the
  direct data-and-formulas form works without a prior run. zRMSE (with its
  mean-filling baseline), 95% interval coverage (with a simulated fair
  expectation for counts) and accuracy (with fractional credit for tied
  votes) match the manuscript's definitions exactly; balanced accuracy and
  the Brier score are reported in addition. Repetitions (`reps`) give an
  honest between-repetition uncertainty. Inspired by Amelia's `overimpute()`
  (Honaker, King & Blackwell 2011, *J Stat Softw* 45(7)), but deliberately
  harsher (a share of cells is hidden at once). The docs state plainly what
  the check cannot do: it validates prediction where data were observed and
  cannot detect missing-not-at-random bias.

* **New `get_imputation_intervals()`: calibrated per-cell intervals for
  imputed values.** The standard multiple-imputation t-interval (Rubin 1987),
  which holds its 95% level for any `n_final` of about 20 or more — unlike
  the raw 2.5/97.5% quantiles of the imputed values (≈92% at `n_final = 50`),
  which the documentation previously suggested. The `n_final` documentation
  now points here.

* **Phylogenetic (tree) uncertainty can now be propagated through the
  imputations.** `bace()`, `bace_imp()` and `bace_final_imp()` accept a
  `multiPhylo` (or list of trees) as `phylo`: the convergence chain conditions
  on the first tree (its stabilisation diagnostic needs a fixed target), and
  each of the `n_final` final imputation runs is fit under its own tree —
  tree `m` for run `m`, selected deterministically from the run index (an
  evenly spaced subsequence when more trees than runs are supplied, balanced
  recycling otherwise), so `seed` reproducibility and serial/parallel equality
  are untouched. Between-imputation variance then includes tree uncertainty:
  proper multiple imputation with the phylogeny as a nuisance parameter drawn
  from its posterior (Rubin 1987; Nakagawa & de Villemereuil 2019, *Systematic
  Biology* 68:632–641). `bace_final_imp()` records `$tree_index` and
  `$tree_hashes`; `pool_posteriors()` detects differing trees, aligns
  species-level BLUPs by name across runs and drops internal-node BLUPs (node
  identities differ between trees), with `BACE_pooling$trees_vary` set;
  `with_imputations(tree = <multiPhylo>)` passes tree `m` to the downstream
  fit on imputation `m`. Runtime cost is negligible (`inverseA` is ~0.05 s at
  10,000 tips; the MCMC fits are indifferent to which tree they receive).

* **Trees are validated up front, with the offending tree named.** Tip
  coverage of the data, negative or zero-length branch lengths (with the
  standard `1e-8` remedy in the message), and every condition
  `MCMCglmm::inverseA()` enforces — most importantly ultrametricity, which
  previously surfaced as an opaque error deep inside a fit (and for a
  `multiPhylo`, only after hours) — are now checked at entry. Node labels
  (e.g. support values) are stripped internally: `inverseA()` rejects
  duplicated ones, and they fragment internal bookkeeping. A warning flags
  trees whose root depths suggest mixed branch-length scalings.

* **`inverseA()` is computed once per distinct tree per phase** instead of
  once per model fit (it was recomputed for every variable × sweep × run),
  via an internal per-tree cache cleared at each phase entry. Results are
  bitwise unchanged; the saving is seconds on long runs, but the cache is
  what makes per-run trees cheap.

# BACE 0.1.0

First public release.

### Fixes before first release (2026-10)

Found while verifying a co-author's review of the pre-release builds (four
reported issues, all confirmed, plus one found during verification); earlier
development builds are affected. None of these crash - they silently degrade
results, which is why they are listed prominently.

* **Warm-up sweeps are robust to transient fit failures.** With several
  variables missing, a warm-up sweep can resample a rare factor level out of
  existence, making the next conditional model's equations singular (most
  likely at strong phylogenetic signal). Such a failure now skips that
  variable's update for that sweep (the previous sweep's values carry
  forward - a skipped Gibbs update) with a warning, instead of aborting the
  run; a failure in the final, kept sweep still errors, so returned results
  are never silently degraded.

* **The default categorical prior no longer warns on every model build.**
  With the default `BACE.gelman = 2`, every categorical prior construction
  emitted "Not enough complete cases to compute Gelman prior or gelman=2"
  even though the pseudo-Gelman prior *is* the documented default - 60+
  spurious warnings per `bace()` call on a dataset with one categorical
  variable, enough to bury real warnings behind R's 50-warning cap. The
  warning now fires only on a genuine unrequested fallback (`gelman = 1`
  with fewer complete cases than J). `bace_final_imp()` also no longer
  builds (and then discards) a categorical prior when
  `ovr_categorical = TRUE` delegates prior construction to the one-vs-rest
  binary sub-models - matching `bace_imp()`, which never built one there.
  Prior values, imputations and RNG streams are unchanged.

* **Categorical imputation silently ignored the phylogeny.** The default
  (multinomial) prediction path looked species effects up under the wrong
  MCMCglmm column name, never matched, and used 0 - predictions equalled a
  non-phylogenetic model's. Fixed (exact prefix matching, with a tripwire
  warning should naming ever drift again); with `species = TRUE`, the second
  species effect is now also included for categorical and ordinal predictions.
* **Final imputations understated uncertainty when two or more variables were
  missing.** Each final run filled the variables once, in formula order, from
  the same point-imputed start, so the first-listed variable conditioned on
  identical point estimates of the others in every run (between-imputation
  variance collapsed; formula order mattered). The final phase now performs
  `sweeps` stochastic passes per run (default 3; warm-up passes use a short
  post-burn-in tail, so the default costs about 1.4x a single pass; collapses
  to 1 pass when at most one variable is missing, where a single pass is
  exact). New optional `seed` gives per-run RNG streams: results reproduce
  and are identical whether runs execute serially or in parallel; the
  parallel fallback now redoes only failed runs.
* **Reported category probabilities were overconfident.** `prob_preds` (and
  the chained phase's argmax fills) evaluated probabilities at posterior-mean
  coefficients; they now average the per-draw probabilities over the
  posterior (imputed values themselves were always proper draws and are
  unchanged in distribution). `prob_preds` is now documented.
* **MCMC diagnostics now cover ordinal cutpoints and the final-phase
  models.** Previously only chained-phase fixed effects and variances were
  checked; the models actually pooled for inference were never diagnosed.
  `bace_final_imp()` results carry a `diagnostics` summary (worst-case ESS /
  Geweke across runs, including one-vs-rest sub-models); printing shows
  flagged parameters only.
* **Priors are now user-settable.** New `priors` argument (named list, one
  element per variable) on `bace()`, `bace_imp()` and `bace_final_imp()`;
  blocks `G`/`R`/`B` replace the defaults wholesale, with validation before
  any MCMC runs. Overriding `R` is refused for threshold/ordinal/categorical
  variables (fixed latent scale). Defaults are unchanged. `bace()` now warns
  about unknown arguments instead of silently ignoring them, and both phases
  read the same `BACE.gelman` option (the final phase used to hard-code a
  different categorical fixed-effect prior than the chained phase).

* Core imputation pipeline: `bace()` orchestrates chained-equations imputation
  (`bace_imp()`), convergence assessment of the chained loop
  (`assess_convergence()`, with a family of `plot_*_convergence()`
  visualisations), multiple posterior-predictive final imputations
  (`bace_final_imp()`, parallelisable via `n_cores`), and posterior pooling
  (`pool_posteriors()`), with automatic retry of the chained phase on
  non-convergence.
* Every per-variable imputation model is a phylogenetic mixed model fitted
  with MCMCglmm; supported response types are gaussian, poisson (count),
  binary and ordinal (threshold scale), and unordered categorical
  (multinomial probit by default; one-vs-rest binary threshold models
  available via `ovr_categorical = TRUE`).
* Deterministic variable-type detection from column storage: doubles are
  gaussian, non-negative integers are poisson counts, factors map to
  threshold/multinomial models by level count and orderedness.
* Rubin's-rules pathway for arbitrary downstream analysis models:
  `with_imputations()` fits any model to each imputed dataset (with optional
  phylogeny pass-through) and `pool_mi()` combines estimates with Rubin's
  rules (Barnard-Rubin small-sample degrees of freedom, `fmi`/`riv`
  diagnostics; matches `mice::pool` to machine precision on identical inputs).
* Accessors for fitted objects: `get_pooled_model()`, `get_imputed_data()`.
* Pre-flight phylogenetic signal screening with `phylo_signal_summary()`
  (latent-scale H2 for all variable types, plus lambda, K and D where
  applicable, with MCMC reliability flags).
* Simulation engine `sim_bace()` for phylogenetically structured mixed-type
  datasets with known ground truth, missingness mechanisms, and optional
  within-species replication (`sim_tree()` simulates the phylogeny alone).
* The public API is deliberately small: imputation pipeline, pooling,
  accessors, diagnostics/plots, signal screening, and the two simulation
  entry points. All other helpers are internal and unexported.
