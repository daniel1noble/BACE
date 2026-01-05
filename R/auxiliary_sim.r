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
#'   and the order in which predictors should be simulated based on dependencies
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
      return(1)
    } else if (grepl("^multinomial", type)) {
      n_cats <- as.numeric(gsub("multinomial", "", type))
      return(n_cats - 1)  # K-1 liabilities for K categories
    } else if (grepl("^ordinal", type)) {
      # Ordinal uses single latent variable with thresholds
      return(1)
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
#' @description Extracts and expands beta coefficients for a variable
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
#' @description Creates a default beta matrix with sparse dependencies
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
  if (prop <= 0) return(x)
  n <- length(x)
  missing_idx <- sample(seq_len(n), size = floor(n * prop), replace = FALSE)
  x[missing_idx] <- NA
  return(x)
}