# =============================================================================
# AUXILIARY FUNCTIONS
# =============================================================================

#' @title Generate default variable names
#' @description Generates default variable names if not provided
#' @param n_predictors Number of predictor variables
#' @return Character vector of variable names (y, x1, x2, ...)
var_name_gen <- function(n_predictors) {
  c("y", paste0("x", seq_len(n_predictors)))
}

#' @title Create random effect specification
#' @description Defines random effect variance components for simulation.
#'   Specifies fractions of variance explained by phylogenetic effect,
#'   non-phylogenetic species effect, and residual. Also handles random slopes.
#' @param var_names Character vector of variable names (response + predictors).
#'   If NULL, auto-generation requires n_vars to be specified.
#' @param n_vars Total number of variables (response + predictors). Only used
#'   if var_names is NULL.
#' @param phylo_frac Numeric vector of fractions of variance explained by
#'   phylogenetic effect for each variable. Values between 0 and 1.
#'   Default is 0 (no phylogenetic signal).
#' @param species_frac Numeric vector of fractions of variance explained by
#'   non-phylogenetic species random effect. Values between 0 and 1.
#'   Default is 0.3 (30% between-species variance).
#' @param total_var Numeric vector of total variance for each variable.
#'   Residual variance = total_var * (1 - phylo_frac - species_frac).
#'   Default is 1 for all variables.
#' @param rr_var Named list specifying random slope variances. Each element
#'   is named by the random effect type (e.g., "phylo", "species") and contains
#'   a named vector of variances for each covariate with random slopes.
#'   Default is 0.25 * corresponding random intercept variance.
#'   Example: list(phylo = c(x1 = 0.15, x2 = 0.20), species = c(x4 = 0.10))
#' @return List with class "raneff" containing:
#'   - var_names: variable names
#'   - phylo_frac: phylogenetic variance fractions
#'   - species_frac: species variance fractions
#'   - total_var: total variances
#'   - phylo_var: phylogenetic variances (phylo_frac * total_var)
#'   - species_var: species variances (species_frac * total_var)
#'   - residual_var: residual variances
#'   - rr_var: random slope variances (if provided)
#' @examples
#' # Auto-generate with defaults
#' raneff <- make_raneff(n_vars = 4)
#'
#' # Custom specification
#' raneff <- make_raneff(
#'   var_names = c("y", "x1", "x2", "x3"),
#'   phylo_frac = c(0.5, 0.3, 0.1, 0.0),
#'   species_frac = c(0.2, 0.3, 0.4, 0.3),
#'   total_var = c(1, 1, 1, 1)
#' )
#'
#' # With random slopes
#' raneff <- make_raneff(
#'   var_names = c("y", "x1", "x2"),
#'   phylo_frac = c(0.5, 0.3, 0.0),
#'   species_frac = c(0.2, 0.3, 0.3),
#'   rr_var = list(phylo = c(x1 = 0.15), species = c(x1 = 0.10))
#' )
#' @export
make_raneff <- function(var_names = NULL, n_vars = NULL,
                        phylo_frac = NULL, species_frac = NULL,
                        total_var = NULL, rr_var = NULL) {
  
  # Determine number of variables
  if (!is.null(var_names)) {
    n_vars <- length(var_names)
  } else if (!is.null(n_vars)) {
    var_names <- var_name_gen(n_vars - 1)
  } else {
    stop("Either var_names or n_vars must be provided")
  }
  
  # Set defaults
  if (is.null(phylo_frac)) {
    phylo_frac <- rep(0, n_vars)
    message("Using default phylo_frac: 0 (no phylogenetic signal)")
  }
  if (is.null(species_frac)) {
    species_frac <- rep(0.3, n_vars)
    message("Using default species_frac: 0.3 (30% between-species variance)")
  }
  if (is.null(total_var)) {
    total_var <- rep(1, n_vars)
  }
  
  # Validate inputs
  if (length(phylo_frac) != n_vars) {
    stop("phylo_frac must have length equal to n_vars")
  }
  if (length(species_frac) != n_vars) {
    stop("species_frac must have length equal to n_vars")
  }
  if (length(total_var) != n_vars) {
    stop("total_var must have length equal to n_vars")
  }
  
  # Check bounds
  if (any(phylo_frac < 0) || any(phylo_frac >= 1)) {
    stop("phylo_frac values must be in [0, 1)")
  }
  if (any(species_frac < 0) || any(species_frac >= 1)) {
    stop("species_frac values must be in [0, 1)")
  }
  if (any(phylo_frac + species_frac >= 1)) {
    stop("phylo_frac + species_frac must be < 1 to leave room for residual variance")
  }
  if (any(total_var <= 0)) {
    stop("total_var values must be positive")
  }
  
  # Calculate variance components
  phylo_var <- phylo_frac * total_var
  species_var <- species_frac * total_var
  residual_var <- total_var * (1 - phylo_frac - species_frac)
  
  # Name vectors
  names(phylo_frac) <- var_names
  names(species_frac) <- var_names
  names(total_var) <- var_names
  names(phylo_var) <- var_names
  names(species_var) <- var_names
  names(residual_var) <- var_names
  
  # Create output structure
  out <- list(
    var_names = var_names,
    phylo_frac = phylo_frac,
    species_frac = species_frac,
    total_var = total_var,
    phylo_var = phylo_var,
    species_var = species_var,
    residual_var = residual_var,
    rr_var = rr_var
  )
  
  class(out) <- c("raneff", "list")
  return(out)
}

#' @title Create fixed effect specification
#' @description Defines fixed effect structure for simulation including
#'   dependency formulas, beta coefficients, interaction terms, and intercepts.
#' @param var_names Character vector of variable names (response + predictors).
#'   Required unless using auto-generation.
#' @param predictor_types Character vector of predictor types.
#'   Required for proper handling of categorical variables.
#' @param formulas List of formulas or single formula defining dependencies.
#'   If single formula (e.g., y ~ x1 + x2), assumes all variables depend on
#'   all other variables. If list, each element defines dependencies for one
#'   variable (e.g., list(y ~ x1 + x2, x2 ~ x1)).
#'   If NULL, auto-generated based on sparsity parameter.
#' @param betas Named list of beta coefficients. Each element corresponds to
#'   a variable (e.g., y, x1, x2) and contains either:
#'   - A numeric vector of coefficients (one per predictor in formula)
#'   - For categorical predictors: K-1 coefficients for multinomialK, or
#'     single coefficient that will be spread across K-1 levels
#'   If NULL, auto-generated with weak to moderate effects.
#' @param interactions List with two elements:
#'   - formulas: List of formulas specifying interactions (e.g., y ~ x1:x2 + x1:x3)
#'   - strengths: Named list of interaction strengths as relative to main effects
#'     (e.g., list(y = c("x1:x2" = 0.5, "x1:x3" = 1.0)))
#'     0 = no interaction, 1 = interaction as strong as main effect
#'   If NULL, no interactions are added.
#' @param intercepts Named numeric vector of intercepts for each variable.
#'   Warnings issued if intercepts are large relative to variance or on
#'   inappropriate scales (should be: identity for gaussian, log for poisson,
#'   logit for categorical). Default is 0 for all variables.
#' @param sparsity Proportion of potential dependencies to set to zero when
#'   auto-generating formulas. Default is 0.7 (70% of edges removed).
#' @param n_vars Total number of variables. Only used for auto-generation when
#'   var_names is NULL.
#' @return List with class "fixeff" containing:
#'   - var_names: variable names
#'   - predictor_types: predictor types
#'   - formulas: dependency formulas
#'   - betas: beta coefficients (expanded for categorical variables)
#'   - interactions: interaction specification
#'   - intercepts: intercept values
#' @examples
#' # Manual specification
#' fixeff <- make_fixeff(
#'   var_names = c("y", "x1", "x2", "x3"),
#'   predictor_types = c("gaussian", "gaussian", "multinomial3"),
#'   formulas = list(
#'     y ~ x1 + x2 + x3,
#'     x2 ~ x1,
#'     x3 ~ x1 + x2
#'   ),
#'   betas = list(
#'     y = c(0.5, 0.3, 0.4),      # One per predictor
#'     x2 = 0.4,
#'     x3 = c(0.3, 0.2)            # Two for multinomial3: K-1 = 2
#'   ),
#'   intercepts = c(y = 0, x1 = 0, x2 = 0, x3 = 0)
#' )
#'
#' # Auto-generation
#' fixeff <- make_fixeff(
#'   var_names = c("y", "x1", "x2"),
#'   predictor_types = c("gaussian", "gaussian"),
#'   sparsity = 0.5
#' )
#'
#' # With interactions
#' fixeff <- make_fixeff(
#'   var_names = c("y", "x1", "x2", "x3"),
#'   predictor_types = c("gaussian", "gaussian", "gaussian"),
#'   formulas = list(y ~ x1 + x2 + x3),
#'   betas = list(y = c(0.5, 0.3, 0.2)),
#'   interactions = list(
#'     formulas = list(y ~ x1:x2),
#'     strengths = list(y = c("x1:x2" = 0.5))
#'   )
#' )
#' @export
make_fixeff <- function(var_names = NULL, predictor_types = NULL,
                        formulas = NULL, betas = NULL,
                        interactions = NULL, intercepts = NULL,
                        sparsity = 0.7, n_vars = NULL) {
  
  # Determine number of variables
  if (!is.null(var_names)) {
    n_vars <- length(var_names)
    n_predictors <- n_vars - 1
  } else if (!is.null(n_vars)) {
    n_predictors <- n_vars - 1
    var_names <- var_name_gen(n_predictors)
  } else {
    stop("Either var_names or n_vars must be provided")
  }
  
  # Validate predictor_types
  if (!is.null(predictor_types)) {
    if (length(predictor_types) != n_predictors) {
      stop("predictor_types must have length equal to n_predictors")
    }
  } else {
    # Default to all gaussian
    predictor_types <- rep("gaussian", n_predictors)
    message("Using default predictor_types: all gaussian")
  }
  
  # Auto-generate formulas if not provided
  if (is.null(formulas)) {
    formulas <- .auto_generate_formulas(var_names, sparsity)
    message("Auto-generated formulas with sparsity = ", sparsity)
  }
  
  # Normalize formulas to list format
  if (!is.list(formulas)) {
    formulas <- list(formulas)
  }
  
  # Auto-generate betas if not provided
  if (is.null(betas)) {
    betas <- .auto_generate_betas(formulas, var_names, predictor_types)
    message("Auto-generated betas with weak to moderate effects")
  }
  
  # Expand betas for categorical variables
  betas <- .expand_categorical_betas(betas, formulas, var_names, predictor_types)
  
  # Set default intercepts
  if (is.null(intercepts)) {
    intercepts <- rep(0, n_vars)
    names(intercepts) <- var_names
  } else {
    if (is.null(names(intercepts))) {
      names(intercepts) <- var_names
    }
    # Validate intercepts
    .validate_intercepts(intercepts, var_names, predictor_types)
  }
  
  # Create output structure
  out <- list(
    var_names = var_names,
    predictor_types = predictor_types,
    formulas = formulas,
    betas = betas,
    interactions = interactions,
    intercepts = intercepts
  )
  
  class(out) <- c("fixeff", "list")
  return(out)
}

# Internal helper functions for make_fixeff

.auto_generate_formulas <- function(var_names, sparsity) {
  n_vars <- length(var_names)
  formulas <- list()
  
  # First variable (response 'y') should depend on predictors
  # Predictors depend on each other based on sparsity
  
  for (i in seq_len(n_vars)) {
    focal_var <- var_names[i]
    
    if (i == 1) {
      # Response variable - depends on all predictors
      if (n_vars > 1) {
        all_preds <- var_names[2:n_vars]
        formula_string <- paste(focal_var, "~", paste(all_preds, collapse = " + "))
        formulas[[length(formulas) + 1]] <- as.formula(formula_string)
      }
    } else if (i == 2) {
      # First predictor - no dependencies (independent)
      next
    } else {
      # Other predictors - can depend on earlier predictors (not response)
      potential_preds <- var_names[2:(i - 1)]  # Exclude response (var 1)
      
      # Randomly select predictors based on sparsity
      keep <- runif(length(potential_preds)) > sparsity
      if (!any(keep)) {
        # Ensure at least one predictor
        keep[sample(length(potential_preds), 1)] <- TRUE
      }
      
      selected_preds <- potential_preds[keep]
      
      # Create formula
      formula_string <- paste(focal_var, "~", paste(selected_preds, collapse = " + "))
      formulas[[length(formulas) + 1]] <- as.formula(formula_string)
    }
  }
  
  return(formulas)
}

.auto_generate_betas <- function(formulas, var_names, predictor_types) {
  betas <- list()
  
  for (form in formulas) {
    lhs <- all.vars(form)[1]  # Response variable
    rhs <- all.vars(form)[-1] # Predictor variables
    
    # Generate weak to moderate effects: uniform(-0.5, 0.5)
    beta_vals <- runif(length(rhs), -0.5, 0.5)
    names(beta_vals) <- rhs
    
    betas[[lhs]] <- beta_vals
  }
  
  return(betas)
}

.expand_categorical_betas <- function(betas, formulas, var_names, predictor_types) {
  predictor_names <- var_names[-1]
  
  for (form in formulas) {
    lhs <- all.vars(form)[1]
    rhs <- all.vars(form)[-1]
    
    if (is.null(betas[[lhs]])) next
    
    beta_vals <- betas[[lhs]]
    expanded_betas <- list()
    
    for (pred_name in rhs) {
      # Find predictor index
      pred_idx <- which(predictor_names == pred_name)
      if (length(pred_idx) == 0) next
      
      pred_type <- predictor_types[pred_idx]
      beta_val <- beta_vals[pred_name]
      
      if (is.na(beta_val)) next
      
      # Check if categorical
      if (grepl("^multinomial", pred_type)) {
        n_cats <- as.numeric(gsub("multinomial", "", pred_type))
        n_params <- n_cats - 1
        
        # If single beta provided, spread across levels
        if (length(beta_val) == 1) {
          expanded <- seq(-abs(beta_val) / 2, abs(beta_val) / 2, 
                          length.out = n_params) * sign(beta_val)
          expanded_betas[[pred_name]] <- expanded
        } else if (length(beta_val) == n_params) {
          expanded_betas[[pred_name]] <- beta_val
        } else {
          stop("Beta for ", pred_name, " (multinomial", n_cats, 
               ") should be length 1 or ", n_params)
        }
      } else {
        # Gaussian, poisson, binary, thresholdK all use single beta
        expanded_betas[[pred_name]] <- beta_val
      }
    }
    
    betas[[lhs]] <- expanded_betas
  }
  
  return(betas)
}

.validate_intercepts <- function(intercepts, var_names, predictor_types) {
  # Check for large intercepts (warning if |intercept| > 5)
  large_intercepts <- abs(intercepts) > 5
  if (any(large_intercepts)) {
    warning("Large intercepts detected (|value| > 5) for: ",
            paste(names(intercepts)[large_intercepts], collapse = ", "),
            "\nEnsure intercepts are on appropriate scales:",
            "\n  - Identity for gaussian",
            "\n  - Log for poisson (e.g., log(5) ≈ 1.6 for mean count of 5)",
            "\n  - Logit for categorical (e.g., logit(0.7) ≈ 0.85 for 70% probability)")
  }
  
  return(invisible(NULL))
}

#' @title Simulate phylogenetic tree and assign cases to species
#' @description Simulates a phylogenetic tree using birth-death process
#'   and assigns observations to species
#' @param n_species Number of species on the tree
#' @param birth Birth rate for tree simulation
#' @param death Death rate for tree simulation
#' @param n_cases Number of observations to generate
#' @param str_len Length of random species name strings
#' @return List containing tree and case-to-species assignment
#' @export
sim_tree <- function(n_species, birth, death, n_cases, str_len = 5) {
  # Generate random species names
  mock_names <- sort(stringi::stri_rand_strings(n_species, str_len))

  # Assign cases to species (with replacement if n_cases > n_species)
  case_species <- sample(mock_names, n_cases, replace = TRUE)
  species_in_data <- sort(unique(case_species))


  # Simulate tree
  tree <- ape::rphylo(
    n = length(species_in_data),
    birth = birth, death = death, T0 = 100,
    fossils = FALSE
  )
  tree$tip.label <- species_in_data

  return(list(tree = tree, cases = case_species, species = species_in_data))
}

#' @title Calculate design size and simulation order
#' @description Determines the number of parameters needed per predictor
#'   and the order in which predictors should be simulated based on dependencies.
#'   To avoid numerical issues - the first predictor should be continuous gaussian.
#' @param predictor_types Character vector of predictor types
#' @param beta_matrix Matrix indicating dependencies between predictors
#' @return List with ns (number of parameters per predictor) and sim_order
design_size <- function(predictor_types, beta_matrix) {
  n_pred <- length(predictor_types)

  if (!(nrow(beta_matrix) == ncol(beta_matrix) &&
          nrow(beta_matrix) == n_pred)) {
    stop("Beta matrix dimensions must match number of predictors")
  }

  # Calculate number of parameters per predictor
  ns <- sapply(predictor_types, function(type) {
    if (type %in% c("gaussian", "binary", "poisson")) {
      1
    } else if (grepl("^multinomial", type)) {
      n_cats <- as.numeric(gsub("multinomial", "", type))
      n_cats - 1 # K-1 liabilities for K categories
    } else if (grepl("^threshold", type)) {
      # Threshold type uses single latent variable with thresholds
      1
    } else {
      stop("Unknown predictor type: ", type)
    }
  })

  # Create dependency matrix (1 where beta != 0)
  dep_matrix <- (beta_matrix != 0) * 1

  # Find independent predictors (no dependencies)
  if (any(rowSums(dep_matrix) == 0)) {
    sim_order <- which(rowSums(dep_matrix) == 0)
  } else {
    stop("At least one predictor must be independent (no dependencies)")
  }

  # Determine simulation order using topological sort
  while (!all(seq_len(n_pred) %in% sim_order)) {
    next_ix <- apply(dep_matrix, 1, function(x, ix) {
      # Check if all dependencies are already in sim_order
      deps <- which(x == 1)
      if (length(deps) > 0 && all(deps %in% ix)) {
        return(TRUE)
      } else if (length(deps) == 0) {
        return(TRUE)
      }
      return(FALSE)
    }, ix = sim_order)

    new_vars <- which(next_ix & !(seq_len(n_pred) %in% sim_order))
    if (length(new_vars) == 0) {
      stop("Circular dependency detected in beta_matrix")
    }
    sim_order <- c(sim_order, new_vars)
  }

  # Check first independent predictor is gaussian
  if (predictor_types[sim_order[1]] != "gaussian") {
    warning("First independent predictor should ideally be gaussian for stability")
  }

  return(list(ns = ns, sim_order = sim_order))
}

#' @title Generate beta coefficients from beta matrix row
#' @description Extracts and expands beta coefficients for a variable. Uses the
#'   beta matrix and response betas to generate the full set of coefficients
#'   needed for predictors with multiple parameters. For cases of multiple categories,
#'   the current versions is just spreading the effect evenly across levels.
#' @param beta_row Row of beta matrix for the focal variable
#' @param ns Vector of parameter counts per predictor
#' @param default_beta Default beta value if using placeholders
#' @return Vector of beta coefficients
beta_generator <- function(beta_row, ns, default_beta = 0.5) {
  vars <- which(beta_row != 0)
  betas <- numeric()

  for (i in vars) {
    if (ns[i] == 1) {
      betas <- c(betas, beta_row[i])
    } else {
      # For multinomial predictors, spread betas across levels
      betas <- c(betas, seq(
        from = -abs(beta_row[i]) / 2,
        to = abs(beta_row[i]) / 2,
        length.out = ns[i]
      ) * sign(beta_row[i]))
    }
  }
  return(betas)
}

#' @title Generate default beta matrix
#' @description Creates a default beta matrix with dependency structure set
#'   through the sparsity parameter. it produces a single beta coefficient per
#'   predictor - in cases multiple values are needed the other function
#'   beta_generator() will expand them accordingly.
#' @param n_predictors Number of predictors
#' @param sparsity Proportion of zero entries (default 0.7)
#' @param beta_range Range for non-zero beta values
#' @return Square matrix of beta coefficients
#' @export
generate_default_beta_matrix <- function(n_predictors, sparsity = 0.7,
                                         beta_range = c(-0.5, 0.5)) {
  beta_matrix <- matrix(0, nrow = n_predictors, ncol = n_predictors)

  # Create lower triangular structure (predictors depend on earlier ones)
  for (i in 2:n_predictors) {
    for (j in 1:(i - 1)) {
      if (runif(1) > sparsity) {
        beta_matrix[i, j] <- runif(1, beta_range[1], beta_range[2])
      }
    }
  }

  return(beta_matrix)
}

#' @title Multinomial liability to category
#' @description Converts multivariate liability scale to categorical outcomes
#' @param liability Matrix of liabilities (n_cases x K-1 for K categories)
#' @param categories Vector of category labels
#' @return Character vector of sampled categories
#' @export
mnom_liab2cat <- function(liability, categories) {
  n_cats <- length(categories)
  n_cases <- nrow(liability)

  # Softmax transformation with reference category
  # Add column of zeros for reference category
  liab_full <- cbind(0, liability)

  # Convert to probabilities via softmax
  probs <- exp(liab_full) / rowSums(exp(liab_full))

  # Sample categories
  result <- apply(probs, 1, function(p) {
    sample(categories, 1, prob = p)
  })

  return(result)
}

#' @title Ordinal liability to category
#' @description Converts latent liability to ordered categorical outcomes (1, 2, 3, ...)
#'   using cumulative probit/logit model with evenly spaced thresholds
#' @param liability Vector of latent liabilities (n_cases)
#' @param n_cats Number of ordered categories
#' @param threshold_spread Controls spread of thresholds (default 1.5)
#' @return Integer vector of ordered categories (1, 2, ..., n_cats)
#' @export
ordinal_liab2cat <- function(liability, n_cats, threshold_spread = 1.5) {
  n_cases <- length(liability)

  # Create evenly spaced thresholds
  # For K categories, we need K-1 thresholds
  thresholds <- seq(
    from = -threshold_spread * (n_cats - 1) / 2,
    to = threshold_spread * (n_cats - 1) / 2,
    length.out = n_cats - 1
  )

  # Assign categories based on thresholds
  # Category k if threshold[k-1] < liability <= threshold[k]
  result <- rep(1L, n_cases)
  for (k in seq_along(thresholds)) {
    result[liability > thresholds[k]] <- as.integer(k + 1)
  }

  return(result)
}

#' @title Sample random effects
#' @description Generates random effect coefficients for a given variance
#' @param sigma2 Variance of the random effect
#' @param n Number of levels
#' @param cor_matrix Correlation matrix (for phylogenetic effects)
#' @return Vector of random effect values
sample_random_effects <- function(sigma2, n, cor_matrix = NULL) {
  if (is.null(cor_matrix)) {
    # Independent random effects
    return(rnorm(n, mean = 0, sd = sqrt(sigma2)))
  } else {
    # Correlated random effects (e.g., phylogenetic)
    return(MASS::mvrnorm(n = 1, mu = rep(0, n), Sigma = sigma2 * cor_matrix))
  }

  # TODO: add later option to generate a RE matrix with given structure
  # e.g., block-diagonal for multiple REs for subsequent categories
  # in categorical (multinomial) variables (for now we assume
  # homogeneous RE variance across categories)
  
}

#' @title Apply missingness to data
#' @description Introduces missing values (NA) at random
#' @param x Vector or column to apply missingness to
#' @param prop Proportion of values to set as missing
#' @return Vector with NA values
apply_missingness <- function(x, prop) {
  if (prop <= 0) {
    return(x)
  }
  n <- length(x)
  missing_idx <- sample(seq_len(n), size = floor(n * prop), replace = FALSE)
  x[missing_idx] <- NA
  return(x)
}

#' @title Parse interaction matrix to extract interaction pairs
#' @description Decodes the interaction matrix to identify which variables interact.
#'   The matrix uses integer codes where variables sharing the same non-zero integer
#'   in a row form an interaction. Returns a list of interaction pairs for each row.
#' @param ix_matrix Square matrix with interaction codes (lower triangular)
#' @param var_names Character vector of variable names (predictors + response)
#' @return Named list where each element contains interaction pairs for that variable
#' @examples
#' # ix_matrix:
#' # 0 0 0 0
#' # 0 0 0 0
#' # 1 1 0 0
#' # 1 1 2 0
#' # For x3: x1:x2 (both have code 1)
#' # For y: x1:x2 (code 1), x1:x3 (code 2 - wait, 1 and 2 share code? No.)
#' # Correction: For row 4 (y): columns with same integer interact
#' # 1,2 share '1' -> x1:x2; columns 1,3 have 1,2 -> no shared code
#' # Actually: row 4 is [1, 1, 2, 0] -> x1 and x2 share '1' -> x1:x2
#' # x1 has 1, x3 has 2 -> different codes, no interaction
#' # Hmm, user example says "1 interacts with 3 (share a 2)" but col1=1, col3=2
#' # Re-reading: "12" in position [4,1] means digit 1 AND digit 2
#' # So [12, 1, 2, 0] means x1 has codes {1,2}, x2 has {1}, x3 has {2}
#' # x1 shares 1 with x2 -> x1:x2; x1 shares 2 with x3 -> x1:x3
parse_ix_matrix <- function(ix_matrix, var_names) {
  n_vars <- nrow(ix_matrix)
  interactions <- list()
  
  for (row_idx in seq_len(n_vars)) {
    row_name <- var_names[row_idx]
    row_interactions <- list()
    
    # Get the row values (only look at columns before current row - lower triangular)
    row_vals <- ix_matrix[row_idx, seq_len(row_idx - 1), drop = TRUE]
    
    if (length(row_vals) == 0 || all(row_vals == 0)) {
      interactions[[row_name]] <- list()
      next
    }
    
    # Parse each cell to extract digit codes
    # E.g., 12 -> c(1, 2); 1 -> c(1); 123 -> c(1, 2, 3)
    col_codes <- lapply(row_vals, function(val) {
      if (val == 0) return(integer(0))
      digits <- as.integer(strsplit(as.character(as.integer(val)), "")[[1]])
      return(digits)
    })
    names(col_codes) <- var_names[seq_len(row_idx - 1)]
    
    # Find all unique codes in this row
    all_codes <- unique(unlist(col_codes))
    all_codes <- all_codes[all_codes != 0]  # Remove zeros
    
    # For each code, find which columns share it -> they interact
    for (code in all_codes) {
      vars_with_code <- names(col_codes)[sapply(col_codes, function(x) code %in% x)]
      if (length(vars_with_code) >= 2) {
        # Create all pairwise interactions
        pairs <- combn(vars_with_code, 2, simplify = FALSE)
        for (pair in pairs) {
          # Sort pair for consistency
          pair <- sort(pair)
          pair_name <- paste(pair, collapse = ":")
          if (!pair_name %in% names(row_interactions)) {
            row_interactions[[pair_name]] <- pair
          }
        }
      }
    }
    
    interactions[[row_name]] <- row_interactions
  }
  
  return(interactions)
}

#' @title Generate interaction coefficients
#' @description Creates coefficients for interaction terms. Can be specified
#'   directly or generated randomly.
#' @param interactions List of interaction pairs (from parse_ix_matrix)
#' @param beta_ix User-specified interaction coefficients (named list or NULL)
#' @param beta_range Range for random coefficient generation
#' @return Named list of interaction coefficients
generate_ix_betas <- function(interactions, beta_ix = NULL, beta_range = c(-0.3, 0.3)) {
  ix_betas <- list()
  
  for (var_name in names(interactions)) {
    var_interactions <- interactions[[var_name]]
    if (length(var_interactions) == 0) {
      ix_betas[[var_name]] <- numeric(0)
      next
    }
    
    var_ix_betas <- numeric(length(var_interactions))
    names(var_ix_betas) <- names(var_interactions)
    
    for (ix_name in names(var_interactions)) {
      # Check if user provided this coefficient
      if (!is.null(beta_ix) && !is.null(beta_ix[[var_name]]) && 
          !is.null(beta_ix[[var_name]][[ix_name]])) {
        var_ix_betas[ix_name] <- beta_ix[[var_name]][[ix_name]]
      } else {
        # Generate random coefficient
        var_ix_betas[ix_name] <- runif(1, beta_range[1], beta_range[2])
      }
    }
    
    ix_betas[[var_name]] <- var_ix_betas
  }
  
  return(ix_betas)
}

#' @title Calculate interaction term values
#' @description Computes the product of interacting variables for use in 
#'   linear predictors. Handles numeric and categorical variables appropriately.
#' @param covars Data frame of covariate values
#' @param interaction_pair Character vector of length 2 with variable names
#' @return Numeric vector of interaction values (or matrix for categorical)
calculate_ix_term <- function(covars, interaction_pair) {
  var1 <- interaction_pair[1]
  var2 <- interaction_pair[2]
  
  # Get values, converting factors to numeric if needed
  val1 <- covars[[var1]]
  val2 <- covars[[var2]]
  
  if (is.factor(val1)) val1 <- as.numeric(val1) - 1  # 0-based for factors
  if (is.factor(val2)) val2 <- as.numeric(val2) - 1
  if (is.character(val1)) val1 <- as.numeric(as.factor(val1)) - 1
  if (is.character(val2)) val2 <- as.numeric(as.factor(val2)) - 1
  
  return(val1 * val2)
}

#' @title Expand beta_resp to full coefficient vector
#' @description Converts beta_resp (vector or list) to a full coefficient vector
#'   that matches the model matrix columns. Supports both simple vectors (where
#'   categorical effects are spread automatically) and lists (where users can
#'   specify exact coefficients for each category level).
#' @param beta_resp Either a numeric vector of length n_predictors, or a list
#'   where each element is either a single value or a vector of K-1 values for
#'   multinomial predictors with K categories
#' @param predictor_types Character vector of predictor types
#' @param var_names Character vector of variable names (including response)
#' @param intercept Response intercept value
#' @return Named list with:
#'   - beta_full: full coefficient vector including intercept
#'   - beta_resp_stored: standardized list format for storage
expand_beta_resp <- function(beta_resp, predictor_types, var_names, intercept = 0) {
  n_predictors <- length(predictor_types)
  predictor_names <- var_names[-1]
  

  # Convert vector to list if needed

if (is.numeric(beta_resp) && !is.list(beta_resp)) {
    if (length(beta_resp) != n_predictors) {
      stop("beta_resp vector must have length equal to number of predictors")
    }
    beta_resp_list <- as.list(beta_resp)
    names(beta_resp_list) <- predictor_names
  } else if (is.list(beta_resp)) {
    # Validate list has correct names or length
    if (is.null(names(beta_resp))) {
      if (length(beta_resp) != n_predictors) {
        stop("Unnamed beta_resp list must have length equal to number of predictors")
      }
      names(beta_resp) <- predictor_names
    }
    beta_resp_list <- beta_resp
  } else {
    stop("beta_resp must be a numeric vector or a list")
  }
  
  # Build full beta vector
  beta_full <- intercept
  beta_resp_stored <- list()
  
  for (i in seq_len(n_predictors)) {
    pred_name <- predictor_names[i]
    pred_type <- predictor_types[i]
    
    # Get beta value(s) for this predictor
    beta_val <- beta_resp_list[[pred_name]]
    if (is.null(beta_val)) beta_val <- beta_resp_list[[i]]
    
    if (grepl("^multinomial", pred_type)) {
      n_cats <- as.numeric(gsub("multinomial", "", pred_type))
      n_levels <- n_cats - 1  # K-1 dummy variables
      
      if (length(beta_val) == 1) {
        # Spread single value across levels
        expanded <- seq(-abs(beta_val) / 2, abs(beta_val) / 2, length.out = n_levels) * sign(beta_val)
        beta_resp_stored[[pred_name]] <- expanded
      } else if (length(beta_val) == n_levels) {
        # Use provided values directly
        expanded <- beta_val
        beta_resp_stored[[pred_name]] <- beta_val
      } else {
        stop("For multinomial", n_cats, " predictor '", pred_name, 
             "', beta_resp must be length 1 or ", n_levels, ", got ", length(beta_val))
      }
      beta_full <- c(beta_full, expanded)
    } else {
      # gaussian, binary, poisson, ordinal - single coefficient
      if (length(beta_val) != 1) {
        stop("For ", pred_type, " predictor '", pred_name, "', beta_resp must be length 1")
      }
      beta_full <- c(beta_full, beta_val)
      beta_resp_stored[[pred_name]] <- beta_val
    }
  }
  
  return(list(
    beta_full = beta_full,
    beta_resp_stored = beta_resp_stored
  ))
}

#' @title Build dependency matrix from formulas
#' @description Converts formula-based dependency specification to a beta matrix
#'   structure for simulation order determination.
#' @param formulas List of formulas defining dependencies
#' @param var_names Character vector of variable names
#' @param betas Named list of beta coefficients
#' @return Square matrix (n_predictors x n_predictors) with dependency structure
.build_dep_matrix_from_formulas <- function(formulas, var_names, betas) {
  predictor_names <- var_names[-1]  # Exclude response
  n_predictors <- length(predictor_names)
  
  # Initialize zero matrix
  dep_matrix <- matrix(0, nrow = n_predictors, ncol = n_predictors)
  rownames(dep_matrix) <- predictor_names
  colnames(dep_matrix) <- predictor_names
  
  for (form in formulas) {
    lhs <- all.vars(form)[1]  # Response variable
    rhs <- all.vars(form)[-1] # Predictor variables
    
    # Find indices in predictor matrix (skip if lhs is response 'y')
    lhs_idx <- which(predictor_names == lhs)
    if (length(lhs_idx) == 0) next  # Response variable, not a predictor
    
    for (pred in rhs) {
      pred_idx <- which(predictor_names == pred)
      if (length(pred_idx) > 0) {
        # Get beta value (use first element if it's a vector for categorical)
        beta_val <- betas[[lhs]][[pred]]
        if (is.null(beta_val)) beta_val <- betas[[lhs]][pred]
        if (!is.null(beta_val) && length(beta_val) > 0) {
          # Use first beta if it's a vector
          dep_matrix[lhs_idx, pred_idx] <- beta_val[1]
        }
      }
    }
  }
  
  return(dep_matrix)
}

#' @title Determine simulation order from formulas
#' @description Determines the order in which variables should be simulated
#'   based on their dependencies specified in formulas. Uses topological sort.
#' @param formulas List of formulas
#' @param var_names Character vector of variable names
#' @param predictor_types Character vector of predictor types
#' @return Integer vector of simulation order (indices into predictor_types)
.determine_sim_order <- function(formulas, var_names, predictor_types) {
  predictor_names <- var_names[-1]
  n_predictors <- length(predictor_names)
  
  # Build dependency graph
  depends_on <- vector("list", n_predictors)
  names(depends_on) <- predictor_names
  
  for (i in seq_len(n_predictors)) {
    depends_on[[i]] <- character(0)
  }
  
  for (form in formulas) {
    lhs <- all.vars(form)[1]
    rhs <- all.vars(form)[-1]
    
    # Find lhs in predictor_names
    lhs_idx <- which(predictor_names == lhs)
    if (length(lhs_idx) == 0) next  # Response variable
    
    # Add dependencies
    for (pred in rhs) {
      if (pred %in% predictor_names) {
        depends_on[[lhs_idx]] <- unique(c(depends_on[[lhs_idx]], pred))
      }
    }
  }
  
  # Topological sort
  sim_order <- integer(0)
  remaining <- seq_len(n_predictors)
  
  while (length(remaining) > 0) {
    # Find variables with all dependencies satisfied
    ready <- sapply(remaining, function(idx) {
      deps <- depends_on[[idx]]
      all(deps %in% predictor_names[sim_order])
    })
    
    if (!any(ready)) {
      stop("Circular dependency detected in formulas")
    }
    
    # Add ready variables to sim_order
    new_vars <- remaining[ready]
    sim_order <- c(sim_order, new_vars)
    remaining <- remaining[!ready]
  }
  
  # Check first predictor is gaussian
  if (predictor_types[sim_order[1]] != "gaussian") {
    warning("First independent predictor should ideally be gaussian for numerical stability")
  }
  
  return(sim_order)
}

#' @title Convert categorical liability to actual categories
#' @description After simulating liability for categorical variables, converts
#'   them to actual factor levels so they can be used as categorical predictors
#'   in subsequent simulations.
#' @param liability Vector or matrix of liabilities
#' @param var_type Character string indicating variable type
#' @param var_name Variable name (for error messages)
#' @return Factor or integer vector of categories
.liability_to_categories <- function(liability, var_type, var_name) {
  if (var_type == "binary") {
    # Binary: threshold at 0
    categories <- as.integer(liability > 0)
    return(factor(categories, levels = c(0, 1)))
  } else if (grepl("^threshold", var_type)) {
    # Ordered categorical
    n_cats <- as.numeric(gsub("threshold", "", var_type))
    categories <- ordinal_liab2cat(liability, n_cats)
    return(factor(categories, levels = seq_len(n_cats), ordered = TRUE))
  } else if (grepl("^multinomial", var_type)) {
    # Unordered categorical (multinomial)
    n_cats <- as.numeric(gsub("multinomial", "", var_type))
    # liability should be a matrix with K-1 columns
    if (!is.matrix(liability)) {
      stop("For multinomial variables, liability must be a matrix")
    }
    cat_labels <- paste0("cat", seq_len(n_cats))
    categories <- mnom_liab2cat(liability, cat_labels)
    return(factor(categories, levels = cat_labels))
  } else {
    # For gaussian/poisson, return as-is
    return(liability)
  }
}

#' @title Parse random effect formula
#' @description Parses random effect specification (e.g., "~phylo", "~species",
#'   "~phylo+species") to determine which random effects to include.
#' @param random_formula Formula or character string specifying random effects
#' @return List with logical flags: has_phylo, has_species
.parse_random_formula <- function(random_formula) {
  if (is.null(random_formula)) {
    random_formula <- "~phylo"
  }
  
  if (is.character(random_formula)) {
    random_formula <- as.formula(random_formula)
  }
  
  if (!inherits(random_formula, "formula")) {
    stop("random_formula must be a formula or character string")
  }
  
  # Extract terms
  terms_str <- as.character(random_formula)[2]  # RHS of formula
  terms_vec <- strsplit(terms_str, "\\+")[[1]]
  terms_vec <- trimws(terms_vec)
  
  has_phylo <- "phylo" %in% terms_vec
  has_species <- "species" %in% terms_vec
  
  if (!has_phylo && !has_species) {
    warning("random_formula contains neither 'phylo' nor 'species'. Using ~phylo as default.")
    has_phylo <- TRUE
  }
  
  return(list(has_phylo = has_phylo, has_species = has_species))
}