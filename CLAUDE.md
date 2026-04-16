# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and test commands

```r
# In R (after opening BACE.Rproj or setwd to repo root):
devtools::load_all()          # load package for interactive development
devtools::document()          # regenerate NAMESPACE and man/ from roxygen2
devtools::test()              # run all tests
devtools::test(filter = "bace_imp")  # run a single test file by pattern
devtools::check()             # full R CMD check (warnings = errors)
devtools::install()           # install locally
```

From the shell:
```bash
cd /path/to/BACE
Rscript -e "devtools::check()"
Rscript -e "devtools::test()"
```

Benchmark script lives outside the package at `/Users/z3437171/Dropbox/Github Local/bace_test/avonet_mixed_benchmark.R`. Run it via a wrapper like `/tmp/long_run_v3.R` (see that file for the pattern), redirecting stdout to a log file.

## Architecture

BACE imputes missing data in phylogenetically-structured comparative datasets using Bayesian chained equations. Every per-variable imputation model is a phylogenetic mixed model fit with `MCMCglmm`.

### Two-phase design — the central invariant

```
bace_imp()          → convergence chain (runs = N iterations)
                       uses argmax / posterior mean for imputed values
                       DO NOT pass sample=TRUE here — it breaks Gelman-Rubin

assess_convergence() → chained-equations convergence check

bace_final_imp()    → n_final independent posterior-predictive draws
                       uses sample=TRUE — each run draws from the posterior
                       predictive so datasets differ between runs

pool_posteriors()   → stack all n_final MCMC chains into one posterior
```

`bace()` is the user-facing wrapper that calls all four in sequence, with auto-retry on non-convergence.

### Variable type detection

Types are inferred automatically from column class. The mapping is in `R/prep_functions.R` (`.detect_type()`):

| Column class | BACE type | MCMCglmm family |
|---|---|---|
| `numeric` | `gaussian` | `gaussian` |
| `integer` | `poisson` | `poisson` |
| `factor`, 2 levels | `binary` | `threshold` |
| `factor`, ≥3 levels | `categorical` | `categorical` (multinomial probit) |
| `ordered` | `threshold` | `threshold` (ordinal probit) |

Override with `ovr_categorical = TRUE` to use J one-vs-rest binary threshold models instead of multinomial probit for categorical variables (better MCMC mixing but slower).

### Key internal functions (all `.`-prefixed, not exported)

| Function | File | Purpose |
|---|---|---|
| `.model_fit()` | `model_functions.R` | Fit one MCMCglmm model; handles single/dual random effects and `slice` sampling |
| `.predict_bace()` | `model_functions.R` | Generate predictions; `sample=TRUE` draws a random MCMC iteration (gaussian: from `Sol %*% t(W)`; poisson: from `model$Liab`; categorical/threshold: via `.impute_levels()`) |
| `.fit_predict_ovr()` | `model_functions.R` | OVR approach: J binary threshold models, aggregate normalised probabilities |
| `.impute_levels()` | `model_functions.R` | Sample a category from probability matrix rows |
| `.data_prep()` | `prep_functions.R` | Scale continuous vars (stores mean/sd for back-transform), detect types, handle factors |
| `.build_formula_string()` | `build_functions.R` | From one formula `y ~ x1 + x2`, generate all variable-swapped formulas (one per variable) |
| `.make_prior()` | `model_functions.R` | Build MCMCglmm prior list; Gelman/pseudo-Gelman options for categorical |
| `.standardize_mcmc_params()` | `bace_imp.R` | Normalise `nitt`/`thin`/`burnin` to per-formula lists |

### Posterior pooling

`pool_posteriors()` stacks `Sol`, `VCV`, `CP`, `Liab` and `Deviance` chains across all `n_final` imputed fits into a single `mcmc` object. This is a Monte Carlo approximation to the marginal posterior integrating over the imputation distribution — **not** Rubin's rules. See the `@details` section of `pool_posteriors()` for the theoretical justification and key references.

### Probability pooling for categorical/threshold variables

`bace_final_imp` stores per-run probability matrices (`all_prob_preds`) with **animal names as rownames**. Downstream pooling must join by name, not by row position — the rows of probability matrices may be in different order across runs. See `pool_cat()` in the benchmark script for the pattern.

### MCMC parameters for categorical/threshold models

- `nitt_cat_mult` (integer, default 1): multiplier applied to `nitt` and `burnin` for categorical and threshold response variables only. Set to 2–3 for hard-to-mix models.
- `slice = TRUE` is set automatically inside `.model_fit()` for threshold/ordinal/categorical types (Neal slice sampling for fixed effects).

### Output object classes

| Class | Created by | Key slots |
|---|---|---|
| `bace` | `bace_imp()` | `$data` (list of imputed datasets per run), `$types`, `$miss_dat`, `$phylo_ran` |
| `bace_final` | `bace_final_imp()` | `$all_models`, `$all_datasets`, `$all_prob_preds` |
| `bace_pooled` | `pool_posteriors()` | `$models` (named list of MCMCglmm objects), `$variables`, `$n_imputations` |
| `bace_complete` | `bace()` | `$pooled_models`, `$final_results`, `$imputed_datasets`, `$prob_preds`, `$convergence`, `$converged` |

## Active development branch

`feature/categorical-improvements` — contains:
- Posterior-predictive sampling in `bace_final_imp` (`sample=TRUE`)
- OVR categorical models (`ovr_categorical` parameter)
- `nitt_cat_mult` parameter
- Probability matrix pooling with animal-name rownames
- Gaussian/Poisson sampling paths in `.predict_bace()`
