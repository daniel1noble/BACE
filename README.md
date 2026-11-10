<img src="vignettes/figs/bace_hex_sticker.png" align="right" width="150" alt="BACE hex sticker" />

# **B**ayesian **A**ugmentation by **C**hained **E**quations (**BACE**) for Comparative Biology

<!-- badges: start -->
[![R-CMD-check](https://github.com/daniel1noble/BACE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daniel1noble/BACE/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/daniel1noble/BACE/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/daniel1noble/BACE/actions/workflows/test-coverage.yaml)
[![Coverage Status](https://img.shields.io/coveralls/github/daniel1noble/BACE/main)](https://coveralls.io/github/daniel1noble/BACE?branch=main)
[![Ask Us Anything\ !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://github.com/daniel1noble/BACE/issues/new)
![Open Source Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)
<!-- badges: end -->

**B**ayesian **A**ugmentation by **C**hained **E**quations (**BACE**) is an R package for imputing missing data in comparative (species-level) datasets using Bayesian chained equations with a phylogenetic random effect. Each per-variable imputation model is a phylogenetic mixed model fit with [`MCMCglmm`](https://cran.r-project.org/package=MCMCglmm), so closely related taxa inform each other's imputations — something general-purpose MI tools (e.g. `mice`, `Amelia`) do not do. 

BACE handles mixed variable types in a single workflow:

- **Continuous** (Gaussian)
- **Count** (Poisson)
- **Binary** (modelled on the threshold / probit scale)
- **Multinomial categorical**
- **Ordinal** (threshold models with ≥ 3 ordered levels)

Variable types are detected automatically from column classes; users can also pass per-variable formulas when they want explicit control over which predictors enter each imputation model.

> ⚠️ *BACE* is under active development and the API may still change. It is usable for research workflows today, but please read the vignette carefully and sanity-check both the per-imputation MCMC convergence and the chained-equations convergence before relying on results for publication.

## Installation
You can install the development version of BACE from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("daniel1noble/BACE")
```

## Quick start

The example below is fully self-contained: `sim_bace()` simulates a comparative
dataset (two gaussian predictors by default) with missing values and a matching
phylogeny, so you can paste and run it as is — then swap in your own data and tree.

```R
library(BACE)

# Simulate a dataset with 20% missingness in each variable, moderate
# phylogenetic signal, and a phylogeny (returned alongside the data)
set.seed(123)
sim <- sim_bace(
  phylo_signal = c(0.6, 0.6, 0.6),
  missingness  = c(0.2, 0.2, 0.2),
  n_cases      = 100,
  n_species    = 100
)

# Full workflow: initial chained imputation -> convergence check -> final
# independent imputations -> pooled posterior.
result <- bace(
  fixformula     = "y ~ x1 + x2",
  ran_phylo_form = "~ 1 | Species",
  phylo          = sim$tree,
  data           = sim$data,
  runs           = 15,
  n_final        = 10,
  nitt           = 100000,
  burnin         = 10000,
  thin           = 10
)

# Pooled posterior for the analysis model
summary(get_pooled_model(result, "y"))

# Imputed datasets (one per final run) for downstream analyses
imputed_list <- get_imputed_data(result)

# Each imputed cell with its mean and a calibrated 95% interval
# (the error bars to report alongside imputed values)
intervals <- get_imputation_intervals(result)

# How well does the imputation work on THIS dataset? Hide some known
# cells, re-impute them, and compare against the hidden truth —
# reusing the finished run's data, tree and settings:
validation <- bace_validate(result, prop = 0.1, seed = 1)
validation
```

See `?bace`, `?bace_imp`, `?bace_final_imp`, `?pool_posteriors`, `?pool_mi`, `?assess_convergence`, `?get_pooled_model`, `?get_imputed_data`, `?get_imputation_intervals`, `?get_scaling`, `?phylo_signal_summary` and `?bace_validate` for details.

## Tutorial

A full tutorial — including when BACE is appropriate relative to other MI tools, the four-step workflow, interpreting pooled output, and the distinction between within-chain MCMC convergence and chained-equations convergence — is at <https://daniel1noble.github.io/BACE/>.

## Currently implemented

- [x] Chained-equations imputation for continuous, count, binary, multinomial and ordinal variables
- [x] Phylogenetic random effect in every per-variable model via `MCMCglmm`
- [x] Optional decomposition of phylogenetic and non-phylogenetic species effects (`species = TRUE`) for datasets with replicated observations per species
- [x] Per-variable MCMC settings — `nitt`, `burnin`, `thin` can be supplied as lists to use different values for different response models
- [x] **One-vs-rest (OVR) categorical imputation** — as an alternative to the default multinomial probit, unordered categorical variables can be modelled as J independent binary threshold models (one per level) via `ovr_categorical = TRUE`
- [x] **`nitt_cat_mult` parameter** — integer multiplier applied to `nitt` and `burnin` for categorical and ordinal variables only, for cases where harder-to-mix models need longer chains
- [x] Gelman / pseudo-Gelman prior options for categorical models to improve mixing, and fully user-settable priors per variable via the `priors` argument
- [x] Stability check on the chained-equations loop (`assess_convergence()`): watches the mean and SD of each variable's imputed values across runs, discards the warm-up half, and names what it finds — drift, oscillation, a collapsing or diverging spread. Movement is measured against the variable's imputation spread, or, when the imputed cells have none (a single missing cell, or counts that all round to the same integer), against the spread of its **observed** values. Where no scale independent of the series exists the verdict is `unassessable` rather than a guess. A pass means no pathology detected, not a certificate (van Buuren 2018, ch. 6.5.2)
- [x] Auto-restart of the chained imputation when — and only when — the stability check finds drift, the one pathology more runs can cure (`bace()` with `max_attempts`)
- [x] Final independent imputations with posterior-predictive sampling (`bace_final_imp()`), so each of the `n_final` imputed datasets is a proper draw rather than a repeated point estimate — continuous variables draw from the posterior predictive via the fixed-effect design matrix, categorical/ordinal variables sample from the posterior probability distribution
- [x] Posterior pooling across imputations by stacking per-imputation MCMC chains — a Monte Carlo approximation to the marginal posterior integrating over the imputation distribution; see `?pool_posteriors` for the references and assumptions
- [x] **Rubin's-rules pathway** — `with_imputations()` fits any downstream analysis model to each imputed dataset and `pool_mi()` combines the estimates with Rubin's rules (Barnard–Rubin small-sample degrees of freedom, `fmi`/`riv` diagnostics; matches `mice::pool` to machine precision on identical inputs). Where `df.residual()` is undefined — `gls`, `phylolm`, `MCMCglmm` — the complete-data degrees of freedom are derived automatically rather than the correction being skipped
- [x] Parallel execution of final imputation runs (`bace_final_imp(n_cores = …)`) via `parallel::mclapply`
- [x] Accessor helpers `get_pooled_model()` and `get_imputed_data()` for extracting pooled models and imputed datasets from `bace_complete` / `bace_pooled` / `bace_final` objects
- [x] `get_imputation_intervals()`: calibrated per-cell 95% intervals for imputed values (Rubin t), for reporting and plotting
- [x] Checkpoint/resume for long runs (`checkpoint_dir`): saved progress is reloaded, so an interrupted analysis continues instead of restarting
- [x] Separate MCMC settings for the final imputation runs (`final_nitt`, `final_thin`, `final_burnin` in `bace()`), so the chained phase and the pooled-inference phase can have different chain lengths
- [x] **`bace_validate()`: cross-validation on your own data** — hide a proportion of the cells you *do* observe, re-impute them with the finished run's own tree and settings, and compare against the held-out truth. Reports zRMSE against a baseline, correlation, and per-cell interval coverage, per variable. For counts it simulates the attainable coverage ceiling rather than assuming 0.95, because whole numbers make exactly 0.95 unreachable
- [x] **Phylogenetic (tree) uncertainty** — pass a `multiPhylo` (or a list of `phylo`) and each of the `n_final` imputation runs is fit under its own tree, so between-imputation variance includes tree uncertainty as well as missing-data uncertainty. The assignment is deterministic and recorded: an evenly spaced subsequence of the supplied trees when there are at least as many trees as runs, balanced recycling otherwise, with `$tree_index` and a content hash per tree kept on the result. The chained convergence phase conditions on the first tree only, by necessity — its stabilisation diagnostic reads imputed values settling toward a fixed target, and cycling trees would inject tree-to-tree jitter into exactly those trajectories. `with_imputations()` **reads that recorded pairing**, so a downstream phylogenetic model (e.g. a PGLS) is fitted against the same tree that generated each imputation, verified by hash rather than assumed
- [x] **Phylogenetic signal pre-flight (`phylo_signal_summary()`)** — fits a univariate phylogenetic mixed model per variable and reports posterior-mean heritability H² with 95% HPD, and optionally Pagel's λ, Blomberg's K and Fritz–Purvis D. A cheap way to see how much the tree can actually contribute before committing to a full run
- [x] **`get_scaling()`** — BACE z-scores continuous variables before fitting, so fitted models carry coefficients in SD units. The constants used are recorded per imputation run *and per model*, so any coefficient can be put back on the data scale exactly. A variable's divisor is not a single number: as the response of its own model it is the observed-only SD, but as a predictor elsewhere it is recomputed from that run's imputed values
- [x] Convergence plots: `plot_trace_convergence()`, `plot_steps_convergence()`, `plot_density_convergence()`, `plot_acf_convergence()` and `plot_convergence_summary()`, plus `plot_bace_imputation_comparison()` for observed-versus-imputed distributions
- [x] Session-wide defaults via `bace_options()` / `bace_option_defaults()`
- [x] Simulation engine `sim_bace()` for generating phylogenetically-structured comparative datasets with arbitrary response types, interaction terms and random slopes (`rr`, `rr_form`), plus `sim_tree()` for the tree and species assignment on its own — useful for method validation

## Not currently supported

- **Random slopes in the imputation models.** `ran_phylo_form` takes a random-intercept
  formula only. A slope is **refused with an error** rather than silently fitted as an
  intercept: BACE cannot yet build a correct prior for one (the inverse-Wishart `nu` is not
  widened with `V`, so a 2×2 structure gets an improper prior). Tracked as
  [#22](https://github.com/daniel1noble/BACE/issues/22). Note this is separate from
  `sim_bace()`, which *can* simulate data with random slopes.

## Future development

- [ ] Random slopes and more flexible random-effect structures per response variable (accept a *list* of random-effect formulas so each imputation model can have its own random-effect specification) — [#22](https://github.com/daniel1noble/BACE/issues/22)
- [ ] Parallelisation within a single chained-equations iteration (currently only the `n_final` independent runs are parallelised)

Bug reports, feature requests and worked examples are very welcome at <https://github.com/daniel1noble/BACE/issues>.
