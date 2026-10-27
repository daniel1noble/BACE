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
