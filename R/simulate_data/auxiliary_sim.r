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

#' @title Simulate phylogenetic tree and assign cases to species
#' @description Simulates a phylogenetic tree using birth-death process
#'   and assigns observations to species
#' @param n_species Number of species on the tree
#' @param birth Birth rate for tree simulation
#' @param death Death rate for tree simulation
#' @param n_cases Number of observations to generate
#' @param str_len Length of random species name strings
#' @return List containing tree and case-to-species assignment
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
    } else if (grepl("^ordinal", type)) {
      # Ordinal uses single latent variable with thresholds
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