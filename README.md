# **B**ayesian **A**ugmentation by **C**hained **E**quations (**BACE**) or BACE with phylogeny (Phylo-BACE)

<!-- badges: start -->
[![R-CMD-check](https://github.com/daniel1noble/BACE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daniel1noble/BACE/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/daniel1noble/BACE/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/daniel1noble/BACE/actions/workflows/test-coverage.yaml)
[![Coverage Status](https://img.shields.io/coveralls/github/daniel1noble/BACE/main)](https://coveralls.io/github/daniel1noble/BACE?branch=main)
[![Ask Us Anything\ !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://github.com/daniel1noble/BACE/issues/new)
![Open Source Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)
<!-- badges: end -->

**B**ayesian **A**ugmentation by **C**hained **E**quations (**BACE**) is an R package for imputing missing data in comparative (species-level) datasets using Bayesian chained equations with a phylogenetic random effect. Each per-variable imputation model is a phylogenetic mixed model fit with [`MCMCglmm`](https://cran.r-project.org/package=MCMCglmm), so closely related taxa inform each other's imputations — something general-purpose MI tools (e.g. `mice`, `Amelia`) do not do. The package is sometimes also referred to as **Phylo-BACE**.

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

```R
library(BACE)

# Full workflow: initial chained imputation -> convergence check -> final
# independent imputations -> pooled posterior.
result <- bace(
  fixformula     = "y ~ x1 + x2",
  ran_phylo_form = "~ 1 | species",
  phylo          = tree,
  data           = my_data,
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
```

See `?bace`, `?bace_imp`, `?bace_final_imp`, `?pool_posteriors`, `?assess_convergence`, `?get_pooled_model` and `?get_imputed_data` for details.

## Tutorial

A full tutorial — including when Phylo-BACE is appropriate relative to other MI tools, the four-step workflow, interpreting pooled output, and the distinction between within-chain MCMC convergence and chained-equations convergence — is at <https://daniel1noble.github.io/BACE/>.

## Currently implemented

- [x] Chained-equations imputation for continuous, count, binary, multinomial and ordinal variables
- [x] Phylogenetic random effect in every per-variable model via `MCMCglmm`
- [x] Optional decomposition of phylogenetic and non-phylogenetic species effects (`species = TRUE`) for datasets with replicated observations per species
- [x] Per-variable MCMC settings — `nitt`, `burnin`, `thin` can be supplied as lists to use different values for different response models
- [x] **One-vs-rest (OVR) categorical imputation** — unordered categorical variables are modelled as J independent binary threshold models (one per level) rather than a single multinomial probit. OVR chains mix more reliably and is the default (`ovr_categorical = TRUE`); the multinomial probit remains available via `ovr_categorical = FALSE`
- [x] **`nitt_cat_mult` parameter** — integer multiplier applied to `nitt` and `burnin` for categorical and ordinal variables only, for cases where harder-to-mix models need longer chains
- [x] Gelman / pseudo-Gelman prior options for categorical models to improve mixing
- [x] Convergence diagnostics across the chained-equations loop (`assess_convergence()`), with summary, Wasserstein and energy-distance methods
- [x] Auto-restart of initial imputation when chained-equations convergence is not reached (`bace()` with `max_attempts`)
- [x] Final independent imputations with posterior-predictive sampling (`bace_final_imp()`), so each of the `n_final` imputed datasets is a proper draw rather than a repeated point estimate — continuous variables draw from the posterior predictive via the fixed-effect design matrix, categorical/ordinal variables sample from the posterior probability distribution
- [x] Posterior pooling across imputations by stacking per-imputation MCMC chains — a Monte Carlo approximation to the marginal posterior integrating over the imputation distribution; see `?pool_posteriors` for the references and assumptions
- [x] Parallel execution of final imputation runs (`bace_final_imp(n_cores = …)`) via `parallel::mclapply`
- [x] Accessor helpers `get_pooled_model()` and `get_imputed_data()` for extracting pooled models and imputed datasets from `bace_complete` / `bace_pooled` / `bace_final` objects
- [x] Simulation engine `sim_bace()` for generating phylogenetically-structured comparative datasets with arbitrary response types, interaction terms and random slopes — useful for method validation

## Future development

- [ ] Random slopes and more flexible random-effect structures per response variable (accept a *list* of random-effect formulas so each imputation model can have its own random-effect specification)
- [ ] Parallelisation within a single chained-equations iteration (currently only the `n_final` independent runs are parallelised)
- [ ] A Rubin's-rules combiner as an alternative scalar-summary path alongside the current stacked-chain combiner
- [ ] Tests for `pool_posteriors()`, `assess_convergence()` and the parallel branch of `bace_final_imp()`
- [ ] Tighter CRAN-ready API surface (fewer exported internals) and a `NEWS.md`

Bug reports, feature requests and worked examples are very welcome at <https://github.com/daniel1noble/BACE/issues>.
