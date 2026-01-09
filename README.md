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
- multiple random effect levels (phylogeny and species) to be used to imputing missing data, including random slopes. 