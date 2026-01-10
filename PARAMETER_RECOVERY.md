# Parameter Recovery in sim_bace

## The Question

**After normalizing variables on the liability scale, are the parameters used to simulate them still recoverable from analyses?**

## Short Answer

**Yes, with caveats** - the parameters are recoverable, but the recovery depends on the variable type:

- **Categorical variables** (binary, threshold, multinomial): ✅ **YES** - parameters fully recoverable
- **Gaussian variables**: ✅ **YES** - betas and variance components recoverable (fixed in latest version)
- **Poisson variables**: ✅ **YES** - betas recoverable (fixed in latest version)

## Detailed Explanation

### How Normalization Works (Updated)

The simulation now uses **different normalization strategies** for different variable types:

#### 1. Categorical Variables (binary, threshold, multinomial)
```r
# Normalized to mean=0, sd=1
linear_pred <- scale(linear_pred, center = TRUE, scale = TRUE)[, 1]
```

**Why sd=1?** In latent variable models (probit, logit, ordinal), the liability is **unobserved**. Only the categories are observed. The variance of the liability is not identified from data alone - we need to fix it for identification. Standard practice is to set it to 1.

**Parameter recovery:** When you fit a probit/logit/ordinal model to the observed categories, you recover the **original betas** (up to the residual variance scale). The model estimates:
```
β_estimated ≈ β_true / σ_residual
```

Since we set σ_residual = 1 in the latent liability, β_estimated ≈ β_true.

#### 2. Gaussian Variables
```r
# Only centered to mean=0, variance preserved
linear_pred <- scale(linear_pred, center = TRUE, scale = FALSE)[, 1]
```

**Why preserve variance?** For gaussian variables, the observed data IS the liability. The variance contains real information about the random effects and residuals.

**Parameter recovery:** When you fit a linear mixed model to the observed data, you recover:
- **Betas**: Original fixed effect coefficients
- **Variance components**: Original phylogenetic, species, and residual variances
- **Random slopes**: Original random slope variances

#### 3. Poisson Variables
```r
# Only centered to mean=0, scale preserved
linear_pred <- scale(linear_pred, center = TRUE, scale = FALSE)[, 1]
```

**Why preserve scale?** For poisson variables, the liability is on the log scale (log link function). The scale affects the implied count distribution.

**Parameter recovery:** When you fit a poisson mixed model with log link, you recover:
- **Betas**: Original regression coefficients (on log scale)
- **Variance components**: The variance structure is preserved

## Example: Parameter Recovery Test

### Simulate Data
```r
library(ape)
library(MASS)
source("R/simulate_auxiliary.r")
source("R/simulate_simBACE.R")

# Define known parameters
my_fixeff <- make_fixeff(
  var_names = c("y", "x1", "x2"),
  predictor_types = c("gaussian", "gaussian"),
  formulas = list(y ~ x1 + x2, x2 ~ x1),
  betas = list(
    y = c(x1 = 0.5, x2 = 0.3),
    x2 = c(x1 = 0.4)
  )
)

my_raneff <- make_raneff(
  var_names = c("y", "x1", "x2"),
  phylo_frac = c(0.3, 0.2, 0.1),
  species_frac = c(0.2, 0.3, 0.4)
)

# Simulate
set.seed(123)
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian"),
  fixeff = my_fixeff,
  raneff = my_raneff,
  n_cases = 500,
  n_species = 50
)

# True parameters
true_beta_x1 <- 0.5
true_beta_x2 <- 0.3
```

### Recover Parameters (Example with simple linear model)
```r
# Fit linear model (ignoring phylogeny for simplicity)
fit <- lm(y ~ x1 + x2, data = sim$data)
summary(fit)

# Compare estimated vs true betas
coef(fit)
# Expected: x1 ≈ 0.5, x2 ≈ 0.3

# For proper analysis including phylogeny, use MCMCglmm, brms, or similar
```

### Recover Parameters (Proper mixed model)
```r
# Using MCMCglmm (example)
library(MCMCglmm)

# Prepare inverse phylogenetic covariance matrix
inv_phylo <- inverseA(sim$tree, nodes = "TIPS")

# Fit model with phylogenetic random effect
prior <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)

fit_phylo <- MCMCglmm(
  y ~ x1 + x2,
  random = ~ species,
  ginverse = list(species = inv_phylo$Ainv),
  data = sim$data,
  prior = prior,
  nitt = 13000,
  burnin = 3000,
  thin = 10,
  verbose = FALSE
)

# Check beta recovery
summary(fit_phylo$Sol)
# Expected: x1 ≈ 0.5, x2 ≈ 0.3

# Check variance component recovery
summary(fit_phylo$VCV)
# Compare to true values:
# Phylo variance: my_raneff$phylo_var["y"]
# Residual variance: my_raneff$residual_var["y"]
```

## Why This Matters

The difference in normalization strategy is crucial:

### Before Fix (ALL variables normalized to sd=1)
- **Categorical**: ✅ Correct (needed for identification)
- **Gaussian**: ❌ Wrong (converts to standardized coefficients)
- **Poisson**: ❌ Wrong (changes count distribution)

### After Fix (differential normalization)
- **Categorical**: ✅ Correct (sd=1 for identification)
- **Gaussian**: ✅ Correct (preserves variance for recovery)
- **Poisson**: ✅ Correct (preserves scale for counts)

## Variance Decomposition

For gaussian variables, you can now recover the full variance decomposition:

```r
# True variance components (from simulation)
total_var <- my_raneff$total_var["y"]
phylo_var <- my_raneff$phylo_var["y"]
species_var <- my_raneff$species_var["y"]
residual_var <- my_raneff$residual_var["y"]

# These should sum to the observed variance in the data
var(sim$data$y, na.rm = TRUE)

# Phylogenetic signal (recoverable)
lambda_true <- phylo_var / total_var
# Can estimate with phylolm, MCMCglmm, brms, etc.
```

## Technical Details

### Liability Scale Mathematics

For **categorical variables**, the probit model assumes:
```
liability_i = X_i * β + u_phylo + u_species + ε
y_i = I(liability_i > threshold)
```

Where ε ~ N(0, σ²_residual). For identification, we set σ²_residual = 1.

The estimated β from a probit model = true_β / σ_residual.
Since σ_residual = 1, estimated_β = true_β ✓

For **gaussian variables**, no latent variable:
```
y_i = X_i * β + u_phylo + u_species + ε
```

All variances are identified. If we standardized y to sd=1, we'd lose this information.
By only centering (mean=0), we preserve:
- True β values
- True variance components

For **poisson variables**, log link:
```
log(E[y_i]) = X_i * β + u_phylo + u_species
```

The scale of the linear predictor matters for the count distribution.
By only centering, we preserve the original β values and variance structure.

## Summary Table

| Variable Type | Normalization | β Recovery | Variance Recovery | Notes |
|---------------|---------------|------------|-------------------|-------|
| Binary | mean=0, **sd=1** | ✅ Yes | ✅ Yes (relative)* | *Residual var fixed at 1 |
| ThresholdK | mean=0, **sd=1** | ✅ Yes | ✅ Yes (relative)* | *Residual var fixed at 1 |
| MultinomialK | mean=0, **sd=1** | ✅ Yes | ✅ Yes (relative)* | *Residual var fixed at 1 |
| Gaussian | **mean=0 only** | ✅ Yes | ✅ Yes (absolute) | Full variance recovered |
| Poisson | **mean=0 only** | ✅ Yes | ✅ Yes (absolute) | On log scale |

## Recommendation

For simulation studies testing imputation or analysis methods:

1. **Use the updated version** - it preserves variance appropriately
2. **For gaussian**: You'll recover exact β and variance components
3. **For categorical**: You'll recover β relative to residual variance (standard practice)
4. **For poisson**: You'll recover β on the log scale
5. **Always check**: Run parameter recovery tests to validate your analysis pipeline

## Code Change Summary

**Before:**
```r
# All variables
linear_pred <- scale(linear_pred, center = TRUE, scale = TRUE)[, 1]  # ❌
```

**After:**
```r
# Categorical variables
if (is_categorical) {
  linear_pred <- scale(linear_pred, center = TRUE, scale = TRUE)[, 1]  # ✅
} else {
  # Gaussian/Poisson
  linear_pred <- scale(linear_pred, center = TRUE, scale = FALSE)[, 1]  # ✅
}
```

This ensures parameters are recoverable while maintaining appropriate identification for each variable type.
