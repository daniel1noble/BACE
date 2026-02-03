# **B**ayesian **A**ugmentation by **C**hained **E**quations (**BACE**) or BACE with phylogeny (Phylo-BACE)

<!-- badges: start -->
[![R-CMD-check](https://github.com/daniel1noble/BACE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daniel1noble/BACE/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/daniel1noble/BACE/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/daniel1noble/BACE/actions/workflows/test-coverage.yaml)
[![Ask Us Anything\ !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://github.com/daniel1noble/BACE/issues/new)
![Open Source Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)
<!-- badges: end -->

**B**ayesian **A**ugmentation by **C**hained **E**quations (**BACE**) is a package that allows users to used chained equations to impute missing data for incomplete comparative datasets using phylogenetic information (or more correctly, **Phylo-BACE**). 

It is currently able to impute missing data for continuous, categorical, and binary variables with variable imputation model structures. Phylogenetic information is included as a random effect in all the imputation models. It uses Bayesian models to estimate missing values across the data. 

*Phylo-BACE* is currently under active development and should not be used until further testing has been done.

## Installation
You can install the development version of BACE from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("daniel1noble/BACE")
```

## Tutorial

A tutorial for using BACE can be found at: https://daniel1noble.github.io/BACE/

## Future Development
- [ ] multiple random effect levels (phylogeny and species) to be used to imputing missing data, including random slopes. 
- [ ] more flexibility in the random effect structure for different variables. Feed list in the `ran_phylo_form` and evaluate it differently. Probably easier to re-construct function as a list with length based on types. Then, if only one formula provided, all elements get the same, if a list, they are all unique. Name by `response_var` then add within functions and loops a call to `ran_phylo`. Then, turn off `ginverse` or evaluate so it is an identity matrix so that argument can remain consistent in `MCMCglmm`. We'll need to adjust the formula construction functions too because MCMCGlmm sets up random slopes differently (e.g., ~us(1 + Time):Chick - from coursenotes.) To make things more general in the random effect formula build we could always use `us()` functions in MCMCglmm to set up unstructured variance-covariance matrices, even for intercept only models. 
- [ ] `nitt`, `burnin`, and `thin` arguments to be variable for different variable types. Probably provide a list or vector. Need checks to ensure these are the same length as the number of variable types.
- [ ] end model simulations to account for imputation uncertainty - sampling from the posterior
- [ ] options for parallelization of imputation of variables during each run.
- [X] Prior B options for categorical variables - Gelman prior or flat prior. Improve mixing.
- [ ] Convergence checks across runs - need to store full predictions for each run, much like data.
- [ ] auto-restart if convergence not reached or, user can specify number of runs.

list("~1|phylo + ~1 + x1|species", "~1|species")
