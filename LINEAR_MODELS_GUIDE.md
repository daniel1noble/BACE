# Linear Models Component in sim_bace Output

## Overview

The `sim_bace()` function now returns a `linear_models` element that provides detailed information about the actual linear model used to simulate each variable. This allows you to:

1. Understand the exact data generating process
2. Access design matrices (X, Z) and parameter vectors (beta, u, e)
3. Verify parameter recovery in your analyses
4. Decompose variance contributions
5. Understand link functions for each variable type

## Structure

For each simulated variable (response and predictors), `linear_models` contains:

```r
sim$linear_models$variable_name <- list(
  link = "link_function",    # Character: "identity", "log", or "probit"
  X = design_matrix,         # Matrix: Fixed effects design matrix
  beta = coefficient_vector, # Numeric vector: Fixed effect coefficients
  Z = random_matrix,         # Matrix: Random effects design matrix
  u = random_effects,        # Numeric vector: Combined random effects (phylo + species)
  e = residuals,             # Numeric vector: Residual errors
  note = NULL                # Character or NULL: Special notes (e.g., for multinomial)
)
```

## Link Functions by Variable Type

| Variable Type | Link Function | Notes |
|---------------|---------------|-------|
| `gaussian` | `identity` | y = X*β + Z*u + e |
| `poisson` | `log` | log(E[y]) = X*β + Z*u |
| `binary` | `probit` | Pr(y=1) = Φ(X*β + Z*u + e) |
| `thresholdK` | `probit` | Latent liability with K-1 thresholds |
| `multinomialK` | `probit` | K-1 latent liabilities (matrix structure)* |

*Note: Multinomial variables require K-1 latent variables. The current implementation stores a note about this but doesn't fully capture the matrix structure.

## Design Matrices

### Fixed Effects Design Matrix (X)
- **Rows**: Number of observations (n_cases)
- **Columns**: Number of predictors + intercept
- **Structure**: Standard model matrix with dummy coding for categorical predictors
- **Example**: For `y ~ x1 + x2` with continuous predictors:
  ```
  X = [1  x1_1  x2_1]
      [1  x1_2  x2_2]
      [...         ]
      [1  x1_n  x2_n]
  ```

### Random Effects Design Matrix (Z)
- **Rows**: Number of observations (n_cases)
- **Columns**: Number of species (n_species)
- **Structure**: Indicator matrix mapping observations to species
- **Example**: Each row has a single 1 indicating which species that observation belongs to

## Coefficients and Effects

### Beta (β) - Fixed Effect Coefficients
- **Length**: Number of columns in X
- **Order**: [intercept, predictor1, predictor2, ...]
- **Interpretation**: Effect sizes specified in `fixeff$betas`

### u - Combined Random Effects
- **Length**: Number of species
- **Composition**: `u = u_phylo + u_species`
- **Interpretation**: Total random effect for each species (phylogenetic + species-specific)
  - Access separate components via `sim$random_effects$u_phylo` and `sim$random_effects$u_species`

### e - Residual Errors
- **Length**: Number of observations
- **Distribution**: Independent N(0, σ²_residual) errors
- **Interpretation**: Observation-level noise

## Important Notes

### Gaussian Variables: Centering

For **gaussian** variables, `sim_bace()` centers the linear predictor to mean=0 before using it as the observed data. This means:

```r
# The stored components represent the TRUE generative model:
linear_pred_raw = X %*% beta + Z %*% u + e

# The observed data is:
observed = scale(linear_pred_raw, center=TRUE, scale=FALSE)
```

**Why?** Centering preserves variance structure while standardizing the location, ensuring parameters can be recovered from analyses.

**For categorical variables** (binary, threshold, multinomial), the linear predictor is normalized to mean=0, sd=1 for identification (standard practice in latent variable models).

### Categorical Predictors

When a predictor is categorical (e.g., `multinomial3`), the design matrix X includes dummy-coded columns:
- `multinomial3` becomes 2 columns in X (K-1 parameterization)
- The beta vector has length matching the number of columns in X

### Independent Predictors

Predictors with no dependencies have:
- X with single column (intercept only)
- beta with single value (intercept)
- Example: `sim$linear_models$x1$X = matrix(1, nrow=n_cases, ncol=1)`

## Example Usage

### Basic Access

```r
sim <- sim_bace(
  response_type = 'gaussian',
  predictor_types = c('gaussian', 'binary'),
  n_cases = 100,
  n_species = 20
)

# View all available linear models
names(sim$linear_models)
# [1] "y"  "x1" "x2"

# Access response linear model
sim$linear_models$y$link        # "identity"
sim$linear_models$y$X           # 100 x 3 matrix
sim$linear_models$y$beta        # c(intercept, x1_coef, x2_coef)
```

### Reconstruct Linear Predictor

```r
# Reconstruct the raw (before centering) linear predictor
linear_pred <- as.numeric(sim$linear_models$y$X %*% sim$linear_models$y$beta) +
               as.numeric(sim$linear_models$y$Z %*% sim$linear_models$y$u) +
               sim$linear_models$y$e

# For gaussian response, the observed data is:
# observed_y = scale(linear_pred, center=TRUE, scale=FALSE)
```

### Variance Decomposition

```r
# Separate phylogenetic and species effects
u_phylo <- sim$random_effects$u_phylo$y
u_species <- sim$random_effects$u_species$y

# Calculate variance contributions
phylo_var <- var(as.numeric(sim$linear_models$y$Z %*% u_phylo))
species_var <- var(as.numeric(sim$linear_models$y$Z %*% u_species))
residual_var <- var(sim$linear_models$y$e)

# Total variance
total_var <- phylo_var + species_var + residual_var

# Phylogenetic signal
lambda <- phylo_var / total_var
```

### Parameter Recovery Test

```r
# Use stored true parameters to test if your analysis recovers them
true_beta_x1 <- sim$linear_models$y$beta[2]  # Coefficient for x1

# Fit your model
library(MCMCglmm)
# ... fit model ...

# Compare estimated vs true
estimated_beta_x1 <- posterior.mode(fit$Sol[, "x1"])
recovery_error <- abs(estimated_beta_x1 - true_beta_x1)
```

### Examine Design Matrix

```r
# View the actual design matrix used
head(sim$linear_models$y$X)
#   (Intercept)         x1         x2
# 1           1 -0.5060455  1.1034216
# 2           1  0.5743192  0.7077346
# ...

# Check for categorical expansions
colnames(sim$linear_models$y$X)
# May include dummy variables like "x2cat2", "x2cat3" if x2 is multinomial
```

## Comparison with random_effects

The `linear_models` component differs from `random_effects`:

| Component | Contains | Purpose |
|-----------|----------|---------|
| `linear_models` | X, beta, Z, u (combined), e | Full model structure for each variable |
| `random_effects` | u_phylo, u_species (separate), residuals | Separate random effect components |

Both are provided for different use cases:
- Use `linear_models` when you need the complete model specification
- Use `random_effects` when you need separate variance components

## Advanced: Working with Multiple Variables

```r
# Iterate through all variables
for (var_name in names(sim$linear_models)) {
  lm_info <- sim$linear_models[[var_name]]
  
  cat("\nVariable:", var_name, "\n")
  cat("  Link:", lm_info$link, "\n")
  cat("  Number of predictors:", ncol(lm_info$X) - 1, "\n")  # -1 for intercept
  cat("  Beta coefficients:", lm_info$beta, "\n")
  
  if (!is.null(lm_info$note)) {
    cat("  Note:", lm_info$note, "\n")
  }
}
```

## Limitations

1. **Multinomial variables**: Full K-1 liability matrix structure not yet implemented. A note is provided indicating this limitation.

2. **Random slopes**: When random slopes are included, the stored `u` represents the base random effects. Random slope contributions are not separately stored in `linear_models` (but are included in the simulated data).

3. **Interaction terms**: Interaction effects are included in the simulation but are not explicitly broken out in the `linear_models` structure. They affect the final linear predictor before it's stored.

4. **Normalization**: Remember that for gaussian/poisson variables, the observed data is centered (mean=0) compared to the raw linear predictor. For categorical variables, it's normalized (mean=0, sd=1).

## See Also

- `?sim_bace` - Main simulation function
- `?make_fixeff` - Specify fixed effects structure
- `?make_raneff` - Specify random effects variance
- `example_linear_models.R` - Comprehensive example script
- `PARAMETER_RECOVERY.md` - Guide to parameter recovery and normalization
