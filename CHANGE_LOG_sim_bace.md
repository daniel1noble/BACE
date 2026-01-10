# sim_bace Function Rewrite - Summary

## Overview

The `sim_bace` function and its auxiliary functions have been completely rewritten to provide more deliberate control over data simulation while maintaining compatibility with the BACE imputation framework.

## Key Changes

### 1. New Helper Functions

#### `make_raneff()`
- Defines random effect variance components
- Parameters:
  - `var_names`: variable names (response + predictors)
  - `phylo_frac`: fraction of variance explained by phylogeny (0-1)
  - `species_frac`: fraction of variance explained by non-phylo species effects (0-1)
  - `total_var`: total variance for each variable
  - `rr_var`: random slope variances (optional)
- Auto-generation: defaults to 0% phylo, 30% species, residual = remainder
- Returns: object of class "raneff" with variance components

#### `make_fixeff()`
- Defines fixed effect structure
- Parameters:
  - `var_names`: variable names
  - `predictor_types`: types of predictors (gaussian, poisson, binary, thresholdK, multinomialK)
  - `formulas`: dependency formulas (list or single formula)
  - `betas`: regression coefficients (auto-expands for categorical variables)
  - `interactions`: interaction specifications (formulas + strengths)
  - `intercepts`: intercept values (with scale warnings)
  - `sparsity`: for auto-generation (proportion of zero coefficients)
- Auto-generation: creates weak to moderate effects with sparse dependencies
- Returns: object of class "fixeff" with all fixed effect specifications

### 2. Updated Main Function Arguments

**New `sim_bace()` signature:**
```r
sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian"),
  n_cases = 200,
  n_species = 75,
  birth = 0.8,
  death = 0.4,
  missingness = NULL,
  var_names = NULL,
  random_formula = "~phylo",
  rr = FALSE,
  rr_form = NULL,
  raneff = NULL,
  fixeff = NULL,
  sparsity = 0.7
)
```

**Key argument changes:**
- `response_type`: unchanged - supports gaussian, poisson, binary, thresholdK, multinomialK
- `predictor_types`: unchanged - same types as response
- `random_formula`: NEW - formula specifying random effects ("~phylo", "~species", or "~phylo+species")
- `raneff`: NEW - output from `make_raneff()` (auto-generated if NULL)
- `fixeff`: NEW - output from `make_fixeff()` (auto-generated if NULL)
- `rr`: unchanged - logical for random slopes
- `rr_form`: updated format - list(raneffect = c("predictor1", "predictor2"))

**Removed arguments (replaced by `raneff` and `fixeff`):**
- `beta_matrix` → now part of `fixeff$formulas` and `fixeff$betas`
- `beta_resp` → now part of `fixeff$betas`
- `beta_sparsity` → replaced by `sparsity` parameter
- `ix_matrix` → replaced by `fixeff$interactions`
- `beta_ix` → replaced by `fixeff$interactions$strengths`
- `intercepts` → now part of `fixeff$intercepts`
- `phylo_signal` → replaced by `raneff$phylo_frac`
- `sigmas` → replaced by `raneff` variance components

### 3. Treatment of Variables

**All variables simulated on liability scale:**
- All responses and predictors start as normalized liabilities (mean=0, sd=1)
- Liabilities are transformed based on variable type:
  - **Gaussian**: remains as continuous liability
  - **Poisson**: not transformed on liability scale (uses exp link on linear predictor)
  - **Binary**: liability > 0 → category 1, else category 0
  - **ThresholdK**: liability binned into K ordered categories using evenly-spaced thresholds
  - **MultinomialK**: K-1 liabilities transformed via softmax to K categories

**Categorical variables properly treated:**
- After liability simulation, categorical variables are converted to factors
- These factors are then used as predictors in subsequent variable simulations
- Model matrices automatically handle factor expansion for multinomial predictors

### 4. New Auxiliary Functions

**Added to simulate_auxiliary.r:**
- `.build_dep_matrix_from_formulas()`: converts formulas to dependency matrix
- `.determine_sim_order()`: topological sort for simulation order
- `.liability_to_categories()`: converts liabilities to appropriate variable types
- `.parse_random_formula()`: parses random effect specifications
- `.auto_generate_formulas()`: creates default dependency structure
- `.auto_generate_betas()`: creates default regression coefficients
- `.expand_categorical_betas()`: expands single beta to K-1 betas for categorical predictors
- `.validate_intercepts()`: checks intercepts and issues warnings

### 5. Simulation Workflow

1. **Setup**: Parse arguments, auto-generate `raneff`/`fixeff` if needed
2. **Tree**: Simulate phylogenetic tree using birth-death process
3. **Random Effects**: Sample phylo, species, and residual random effects
4. **Determine Order**: Use topological sort to determine simulation order based on dependencies
5. **Simulate Predictors**: In dependency order, simulate each predictor
   - Build linear predictor from dependencies
   - Add random effects (phylo + species + residual)
   - Add random slopes if specified
   - Add interactions if specified
   - Normalize to mean=0, sd=1
   - Convert liability to appropriate type
6. **Simulate Response**: Similar to predictors but for response variable
7. **Apply Missingness**: Introduce missing data as specified
8. **Output**: Return data, tree, parameters, and random effects

### 6. Output Structure

```r
list(
  data = data.frame,          # Simulated data
  tree = phylo object,        # Phylogenetic tree
  params = list(              # All simulation parameters
    response_type, predictor_types, var_names,
    n_cases, n_species, n_species_actual,
    birth, death, missingness,
    random_formula, random_structure,
    rr, rr_form, raneff, fixeff, sparsity,
    sim_order
  ),
  random_effects = list(      # All random components
    u_phylo,                  # Phylogenetic random effects
    u_species,                # Species random effects
    u_slopes,                 # Random slopes
    residuals,                # Residual errors
    liabilities               # Underlying liabilities (NEW)
  )
)
```

## Examples

### Basic Usage (Auto-generated)
```r
# Simplest call - all parameters auto-generated
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian"),
  n_cases = 200,
  n_species = 50
)
```

### Custom Fixed and Random Effects
```r
# Define random effect structure
my_raneff <- make_raneff(
  var_names = c("y", "x1", "x2"),
  phylo_frac = c(0.5, 0.3, 0.1),
  species_frac = c(0.2, 0.3, 0.4)
)

# Define fixed effect structure
my_fixeff <- make_fixeff(
  var_names = c("y", "x1", "x2"),
  predictor_types = c("gaussian", "gaussian"),
  formulas = list(
    y ~ x1 + x2,
    x2 ~ x1
  ),
  betas = list(
    y = c(x1 = 0.5, x2 = 0.3),
    x2 = c(x1 = 0.4)
  )
)

# Simulate
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian"),
  raneff = my_raneff,
  fixeff = my_fixeff,
  n_cases = 200,
  n_species = 50
)
```

### With Interactions
```r
my_fixeff <- make_fixeff(
  var_names = c("y", "x1", "x2", "x3"),
  predictor_types = c("gaussian", "gaussian", "gaussian"),
  formulas = list(y ~ x1 + x2 + x3),
  betas = list(y = c(0.5, 0.3, 0.2)),
  interactions = list(
    formulas = list(y ~ x1:x2 + x1:x3),
    strengths = list(y = c("x1:x2" = 0.5, "x1:x3" = 1.0))
  )
)

sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian", "gaussian"),
  fixeff = my_fixeff,
  n_cases = 200
)
```

### With Random Slopes
```r
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian"),
  random_formula = "~phylo+species",
  rr = TRUE,
  rr_form = list(
    phylo = c("x1"),
    species = c("x1", "x2")
  ),
  n_cases = 200,
  n_species = 50
)
```

### Mixed Variable Types
```r
sim <- sim_bace(
  response_type = "poisson",
  predictor_types = c("gaussian", "binary", "threshold3", "multinomial4"),
  random_formula = "~phylo+species",
  n_cases = 300,
  n_species = 100
)
```

## Convenience Wrappers

Updated to use new interface:
```r
sim_bace_gaussian(n_predictors = 3, n_cases = 200, n_species = 75, sparsity = 0.7)
sim_bace_poisson(n_predictors = 3, n_cases = 200, n_species = 75, sparsity = 0.7)
sim_bace_binary(n_predictors = 3, n_cases = 200, n_species = 75, sparsity = 0.7)
```

## Benefits of Rewrite

1. **More Deliberate Control**: Users can precisely specify all aspects of simulation
2. **Better Categorical Handling**: Proper K-1 parameterization with automatic expansion
3. **Cleaner Interface**: Separate `make_raneff()` and `make_fixeff()` functions
4. **Formula-Based Dependencies**: More intuitive than beta matrices
5. **Normalized Liabilities**: All variables start from same scale
6. **Flexible Random Effects**: Easy to specify phylo-only, species-only, or both
7. **Better Warnings**: Intercept scale warnings, random slope validation
8. **Auto-Generation**: Sensible defaults when parameters not specified
9. **Stored Liabilities**: Access to underlying continuous scale for all variables

## Backward Compatibility Notes

The old interface (with `beta_matrix`, `beta_resp`, `phylo_signal`, etc.) is **not supported** in the rewritten version. Users must migrate to the new interface using `make_raneff()` and `make_fixeff()`.

Migration is straightforward:
- Old `phylo_signal` → New `raneff` with `phylo_frac`
- Old `beta_matrix` → New `fixeff` with `formulas` and `betas`
- Old `beta_resp` → New `fixeff$betas$y`
- Old `ix_matrix` → New `fixeff$interactions`
- Old `sigmas` → New `raneff` variance components

## Testing

All functionality tested and validated:
- ✓ Basic gaussian simulation with auto-generation
- ✓ Custom raneff and fixeff specifications
- ✓ Mixed variable types (gaussian, binary, threshold, multinomial, poisson)
- ✓ Interaction terms
- ✓ Random slopes
- ✓ Convenience wrappers

See `test_new_sim_bace.R` for complete test suite.
