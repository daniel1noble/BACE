# BACE 0.1.0

First public release.

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
