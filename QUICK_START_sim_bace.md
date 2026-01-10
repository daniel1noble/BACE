# sim_bace Quick Start Guide

## Basic Examples

### 1. Simplest Simulation (All Defaults)

```r
# Gaussian response, 2 gaussian predictors, auto-generated everything
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian"),
  n_cases = 200,
  n_species = 50
)


# View the simulation
print(sim)
head(sim$data)
```

### 2. Specify Random Effect Structure

```r
# Define variance components
my_raneff <- make_raneff(
  var_names = c("y", "x1", "x2"),
  phylo_frac = c(0.4, 0.3, 0.2),      # 40%, 30%, 20% phylogenetic
  species_frac = c(0.2, 0.3, 0.3)     # 20%, 30%, 30% species
)

sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian"),
  raneff = my_raneff,
  n_cases = 200,
  n_species = 50
)
```

### 3. Specify Fixed Effects (Dependencies)

```r
# Define dependencies and effect sizes
my_fixeff <- make_fixeff(
  var_names = c("y", "x1", "x2", "x3"),
  predictor_types = c("gaussian", "gaussian", "gaussian"),
  formulas = list(
    y ~ x1 + x2 + x3,    # y depends on all predictors
    x2 ~ x1,             # x2 depends on x1
    x3 ~ x1 + x2         # x3 depends on x1 and x2
  ),
  betas = list(
    y = c(x1 = 0.5, x2 = 0.3, x3 = 0.2),
    x2 = c(x1 = 0.4),
    x3 = c(x1 = 0.3, x2 = 0.25)
  )
)

sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian", "gaussian"),
  fixeff = my_fixeff,
  n_cases = 200,
  n_species = 50
)
```

### 4. Poisson Response (Count Data)

```r
sim <- sim_bace(
  response_type = "poisson",
  predictor_types = c("gaussian", "gaussian"),
  n_cases = 200,
  n_species = 50
)

# Intercepts should be on log scale for poisson
my_fixeff <- make_fixeff(
  var_names = c("y", "x1", "x2"),
  predictor_types = c("gaussian", "gaussian"),
  formulas = list(y ~ x1 + x2),
  betas = list(y = c(0.3, 0.2)),
  intercepts = c(y = log(5), x1 = 0, x2 = 0)  # log(5) ≈ 1.6 for mean count of 5
)
```

### 5. Binary Response

```r
sim <- sim_bace(
  response_type = "binary",
  predictor_types = c("gaussian", "gaussian"),
  n_cases = 200,
  n_species = 50
)

# View distribution
table(sim$data$y)
```

### 6. Categorical Predictors

```r
# Ordered categories (threshold)
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "threshold3", "binary"),
  n_cases = 200,
  n_species = 50
)

# Unordered categories (multinomial)
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "multinomial4"),
  n_cases = 200,
  n_species = 50
)

# Categorical response
sim <- sim_bace(
  response_type = "threshold4",
  predictor_types = c("gaussian", "gaussian"),
  n_cases = 200,
  n_species = 50
)
```

### 7. Add Interactions

```r
my_fixeff <- make_fixeff(
  var_names = c("y", "x1", "x2", "x3"),
  predictor_types = c("gaussian", "gaussian", "gaussian"),
  formulas = list(y ~ x1 + x2 + x3),
  betas = list(y = c(0.5, 0.3, 0.2)),
  interactions = list(
    formulas = list(y ~ x1:x2 + x1:x3),
    strengths = list(y = c("x1:x2" = 0.5, "x1:x3" = 0.8))
  )
)

sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian", "gaussian"),
  fixeff = my_fixeff,
  n_cases = 200,
  n_species = 50
)
```

### 8. Random Slopes

```r
# Both phylogenetic and species random effects, with random slopes
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian"),
  random_formula = "~phylo+species",
  rr = TRUE,
  rr_form = list(
    phylo = c("x1"),           # x1 has phylo random slope
    species = c("x1", "x2")    # x1 and x2 have species random slopes
  ),
  n_cases = 200,
  n_species = 50
)
```

### 9. Missing Data

```r
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian"),
  missingness = c(0.1, 0, 0.2),  # 10% missing in y, 0% in x1, 20% in x2
  n_cases = 200,
  n_species = 50
)

# Check missing data
sapply(sim$data, function(x) sum(is.na(x)))
```

### 10. Complex Mixed Design

```r
# Define everything
my_raneff <- make_raneff(
  var_names = c("y", "x1", "x2", "x3"),
  phylo_frac = c(0.3, 0.2, 0.1, 0.15),
  species_frac = c(0.3, 0.4, 0.4, 0.35)
)

my_fixeff <- make_fixeff(
  var_names = c("y", "x1", "x2", "x3"),
  predictor_types = c("gaussian", "binary", "threshold3"),
  formulas = list(
    y ~ x1 + x2 + x3,
    x2 ~ x1,
    x3 ~ x1
  ),
  betas = list(
    y = c(x1 = 0.5, x2 = 0.3, x3 = 0.2),
    x2 = c(x1 = 0.4),
    x3 = c(x1 = 0.3)
  ),
  intercepts = c(y = 0, x1 = 0, x2 = 0, x3 = 0)
)

sim <- sim_bace(
  response_type = "poisson",
  predictor_types = c("gaussian", "binary", "threshold3"),
  raneff = my_raneff,
  fixeff = my_fixeff,
  random_formula = "~phylo+species",
  rr = TRUE,
  rr_form = list(phylo = c("x1"), species = c("x1")),
  missingness = c(0.05, 0, 0.1, 0.15),
  n_cases = 300,
  n_species = 100,
  birth = 0.8,
  death = 0.4
)

print(sim)
```

## Convenience Wrappers

For quick simulations:

```r
# Gaussian
sim <- sim_bace_gaussian(n_predictors = 3, n_cases = 200, n_species = 50)

# Poisson
sim <- sim_bace_poisson(n_predictors = 2, n_cases = 200, n_species = 50)

# Binary
sim <- sim_bace_binary(n_predictors = 2, n_cases = 200, n_species = 50)
```

## Accessing Output

```r
# Simulated data
data <- sim$data
head(data)

# Phylogenetic tree
tree <- sim$tree
plot(tree)

# Parameters used
params <- sim$params
names(params)

# Random effects
random_effects <- sim$random_effects
names(random_effects)

# Underlying liabilities (new!)
liabilities <- sim$random_effects$liabilities
head(liabilities$y)  # Continuous liability for response
```

## Key stuff


1. **Check Intercepts**: For non-gaussian responses, ensure intercepts are on appropriate scales
   - Poisson: log scale (e.g., log(5) ≈ 1.6)
   - Binary/Categorical: logit scale (e.g., logit(0.7) ≈ 0.85)
2. **Categorical Betas**: For multinomialK predictors, you can provide:
   - Single value (auto-spreads across K-1 levels)
   - K-1 values (one per level contrast)
3. **Random Effect Structure**: Choose based on your needs
   - "~phylo": Phylogenetic signal only
   - "~species": Between-species variance only (no phylogeny)
   - "~phylo+species": Both sources of variance
4. **Simulation Order**: Variables are simulated in dependency order automatically
5. **Random Slopes**: Only for continuous predictors (gaussian, poisson)
6. **Normalized Liabilities**: All variables start at mean=0, sd=1 on liability scale
7. **Missing Data**: Applied after simulation, so doesn't affect dependencies

## Common Issues

**Issue**: "Variables X not yet simulated"
**Solution**: Check your formulas - variables can only depend on those simulated earlier

**Issue**: "Beta length mismatch"
**Solution**: Factor predictors create multiple dummy variables. Let auto-expansion handle it or provide correct number of betas

**Issue**: "Random slopes ignored for non-continuous covariates"
**Solution**: Random slopes only work for gaussian and poisson predictors, not categorical

## Getting Help

```r
# Print simulation summary
print(sim)
print_sim_bace_summary(sim)
```