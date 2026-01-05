#' @title simBACE: Comprehensive simulator of mock data for BACE
#' @description simBACE is a function to simulate mock data for BACE.
#'   The function simulates data with defined phylogenetic structure and
#'   with multiple dependencies on additional covariates. The data can be
#'   gaussian, binary, poissonian, multinomial (unordered categories), 
#'   or ordinal (ordered categories coded as integers 1,2,3,...).
#' @author BACE Development Team
#' @importFrom ape rphylo vcv Ntip
#' @importFrom MASS mvrnorm ginv
#' @importFrom stats rnorm rbinom rpois plogis model.matrix as.formula
#' @importFrom stringi stri_rand_strings

# Source auxiliary functions
library(here)
source(here::here("R", "auxiliary_sim.r"))

# =============================================================================
# MAIN SIMULATION FUNCTION
# =============================================================================

#' @param response_type Type of response variable: "gaussian", "poisson",
#'   "binary", "multinomialK" (unordered categories), or "ordinalK" 
#'   (ordered categories coded as integers 1,2,...,K) where K is number of categories
#' @param predictor_types Character vector specifying type of each predictor
#' @param var_names Character vector of variable names (response + predictors).
#'   If NULL, default names (y, x1, x2, ...) are generated
#' @param beta_matrix Square matrix of regression coefficients. Entry [i,j]
#'   indicates effect of predictor j on predictor i. If NULL, a default
#'   sparse matrix is generated
#' @param beta_resp Vector of regression coefficients for response on predictors
#' @param intercepts Named list with 'predictors' (vector) and 'response' (scalar)
#' @param birth Birth rate for phylogenetic tree simulation (default 0.8)
#' @param death Death rate for phylogenetic tree simulation (default 0.4)
#' @param phylo_signal Numeric vector of phylogenetic signals for response and
#'   each predictor. Values between 0 and 1 (default 0 = no phylogenetic signal)
#' @param n_cases Total number of observations (default 200)
#' @param n_species Number of species on the tree (default 75)
#' @param missingness Numeric vector indicating proportion of missing values
#'   for response and each predictor (default all 0s = no missing data)
#' @param sigmas Named list of variance components:
#'   - sigma_species: variances for non-phylogenetic species effects
#'   - sigma_residual: residual variances
#'   Each element is a vector for response + predictors
#' @param rr Logical, if TRUE allow random slopes (default FALSE
#' @param rr_form Named list specifying random slopes. Each element is named
#'   by random effect and contains character vector of covariate names
#'
#' @return Named list containing:
#'   - data: data.frame with simulated data
#'   - tree: phylogenetic tree
#'   - params: list of simulation parameters used
#'   - random_effects: list of simulated random effect coefficients
#'
#' @examples
#' # Basic gaussian simulation
#' sim <- simBACE(
#'   response_type = "gaussian",
#'   predictor_types = c("gaussian", "gaussian", "binary"),
#'   n_cases = 100,
#'   n_species = 30
#' )
#'
#' # With phylogenetic signal and random slopes
#' sim <- simBACE(
#'   response_type = "poisson",
#'   predictor_types = c("gaussian", "binary"),
#'   phylo_signal = c(0.5, 0.3, 0.1),
#'   rr = TRUE,
#'   rr_form = list(species = c("x1"))
#' )
#'
#' @export

simBACE <- function(
    response_type = "gaussian",
    predictor_types = c("gaussian", "gaussian"),
    var_names = NULL,
    beta_matrix = NULL,
    beta_resp = NULL,
    beta_sparsity = 0.7,
    intercepts = NULL,
    birth = 0.8,
    death = 0.4,
    phylo_signal = NULL,
    n_cases = 200,
    n_species = 75,
    missingness = NULL,
    sigmas = NULL,
    rr = FALSE,
    rr_form = NULL) {
  # -------------------------------------------------------------------------
  # SETUP AND VALIDATION
  # -------------------------------------------------------------------------

  n_predictors <- length(predictor_types)
  n_vars <- n_predictors + 1 # response + predictors

  # Generate default variable names if not provided
  if (is.null(var_names)) {
    var_names <- var_name_gen(n_predictors)
  }
  if (length(var_names) != n_vars) {
    stop("var_names must have length equal to number of predictors + 1 (response)")
  }

  # Generate default beta matrix if not provided
  if (is.null(beta_matrix)) {
    beta_matrix <- generate_default_beta_matrix(n_predictors, sparsity = beta_sparsity)
    message("Generated default beta_matrix with sparse dependencies")
  }
  if (!is.matrix(beta_matrix) || nrow(beta_matrix) != n_predictors ||
    ncol(beta_matrix) != n_predictors) {
    stop("beta_matrix must be a square matrix with dimensions equal to n_predictors")
  }

  # Name beta_matrix rows and columns according to predictor names
  predictor_names <- var_names[-1]  # Exclude response name
  rownames(beta_matrix) <- predictor_names
  colnames(beta_matrix) <- predictor_names

  # Generate default beta_resp if not provided
  if (is.null(beta_resp)) {
    beta_resp <- runif(n_predictors, -0.5, 0.5)
    message("Generated default beta_resp: ", paste(round(beta_resp, 3), collapse = ", "))
  }
  if (length(beta_resp) != n_predictors) {
    stop("beta_resp must have length equal to number of predictors")
  }

  # Setup intercepts
  if (is.null(intercepts)) {
    intercepts <- list(
      predictors = rep(0, n_predictors),
      response = 0
    )
  }

  # Setup phylogenetic signal (default: no signal)
  if (is.null(phylo_signal)) {
    phylo_signal <- rep(0, n_vars)
  }
  if (length(phylo_signal) != n_vars) {
    stop("phylo_signal must have length equal to n_vars (response + predictors)")
  }
  if (any(phylo_signal < 0) || any(phylo_signal >= 1)) {
    stop("phylo_signal values must be in [0, 1)")
  }

  # Setup missingness (default: no missing data)
  if (is.null(missingness)) {
    missingness <- rep(0, n_vars)
  }
  if (length(missingness) != n_vars) {
    stop("missingness must have length equal to n_vars (response + predictors)")
  }

  # Setup variance components
  if (is.null(sigmas)) {
    sigmas <- list(
      sigma_species = rep(1, n_vars),
      sigma_residual = rep(1, n_vars)
    )
  }

  # Validate random slopes specification
  if (rr && is.null(rr_form)) {
    warning("rr=TRUE but rr_form not specified. No random slopes will be generated.")
    rr <- FALSE
  }

  # -------------------------------------------------------------------------
  # GENERATE PHYLOGENETIC TREE
  # -------------------------------------------------------------------------

  taxa <- sim_tree(n_species, birth, death, n_cases)
  tree <- taxa$tree
  case_species <- taxa$cases
  species_list <- taxa$species
  n_species_actual <- length(species_list)

  # Phylogenetic correlation matrix
  cor_phylo <- ape::vcv(tree, corr = TRUE)

  # -------------------------------------------------------------------------
  # CALCULATE VARIANCE COMPONENTS
  # -------------------------------------------------------------------------

  sigma2_species <- sigmas$sigma_species^2
  sigma2_residual <- sigmas$sigma_residual^2

  # Calculate phylogenetic variance from heritability
  # phylo_signal = sigma2_phylo / (sigma2_phylo + sigma2_species + sigma2_residual)
  sigma2_phylo <- (phylo_signal * (sigma2_species + sigma2_residual)) /
    pmax(1 - phylo_signal, 1e-6)

  # -------------------------------------------------------------------------
  # SAMPLE RANDOM EFFECTS
  # -------------------------------------------------------------------------

  # Species random effects (non-phylogenetic)
  u_species <- lapply(sigma2_species, function(s2) {
    sample_random_effects(s2, n_species_actual, cor_matrix = NULL)
  })
  names(u_species) <- var_names

  # Phylogenetic random effects
  u_phylo <- lapply(sigma2_phylo, function(s2) {
    if (s2 > 0) {
      sample_random_effects(s2, n_species_actual, cor_matrix = cor_phylo)
    } else {
      rep(0, n_species_actual)
    }
  })
  names(u_phylo) <- var_names

  # Residuals
  residuals <- lapply(sigma2_residual, function(s2) {
    rnorm(n_cases, mean = 0, sd = sqrt(s2))
  })
  names(residuals) <- var_names

  # Random slopes (if specified)
  u_slopes <- list()
  if (rr && !is.null(rr_form)) {
    for (re_name in names(rr_form)) {
      covars_with_slopes <- rr_form[[re_name]]
      u_slopes[[re_name]] <- list()

      for (cov in covars_with_slopes) {
        if (re_name == "species") {
          # Random slopes at species level
          u_slopes[[re_name]][[cov]] <- rnorm(
            n_species_actual,
            mean = 0,
            sd = sqrt(sigma2_species[1] * 0.5) # Half the species variance
          )
        }
      }
    }
  }

  # -------------------------------------------------------------------------
  # CREATE DESIGN MATRICES
  # -------------------------------------------------------------------------

  # Z-matrix for species (maps observations to species)
  Z <- model.matrix(
    ~ 0 + factor(species),
    data.frame(species = case_species)
  )
  colnames(Z) <- species_list

  # Get design information
  mydesign <- design_size(predictor_types, beta_matrix)

  # -------------------------------------------------------------------------
  # SIMULATE PREDICTORS
  # -------------------------------------------------------------------------

  # Initialize predictor data frame
  covars <- as.data.frame(matrix(
    0,
    nrow = n_cases,
    ncol = n_predictors,
    dimnames = list(NULL, var_names[-1])
  ))

  # Simulate predictors in dependency order
  for (vari in mydesign$sim_order) {
    pred_type <- predictor_types[vari]
    pred_name <- var_names[vari + 1] # +1 because var_names includes response
    var_idx <- vari + 1 # Index in full variable list (response = 1)

    # Get dependencies
    which_deps <- which(beta_matrix[vari, ] != 0)

    # Build linear predictor
    if (length(which_deps) == 0) {
      # Independent predictor
      linear_pred <- intercepts$predictors[vari]
    } else {
      # Dependent predictor
      dep_names <- var_names[which_deps + 1]

      X_formula <- paste("~ 1 +", paste(dep_names, collapse = " + "))
      X <- model.matrix(as.formula(X_formula), covars[, dep_names, drop = FALSE])

      betas <- c(
        intercepts$predictors[vari],
        beta_generator(beta_matrix[vari, ], mydesign$ns)
      )

      linear_pred <- as.numeric(X %*% betas)
    }

    # Add random effects
    linear_pred <- linear_pred +
      as.numeric(Z %*% u_species[[var_idx]]) +
      as.numeric(Z %*% u_phylo[[var_idx]]) +
      residuals[[var_idx]]

    # Add random slopes if applicable
    if (rr && "species" %in% names(u_slopes)) {
      for (cov in names(u_slopes$species)) {
        if (cov %in% names(covars) && !all(covars[[cov]] == 0)) {
          slope_contribution <- as.numeric(Z %*% u_slopes$species[[cov]]) * covars[[cov]]
          linear_pred <- linear_pred + slope_contribution
        }
      }
    }

    # Generate response based on predictor type
    if (pred_type == "gaussian") {
      covars[[pred_name]] <- linear_pred
    } else if (pred_type == "binary") {
      probs <- plogis(linear_pred)
      covars[[pred_name]] <- rbinom(n_cases, 1, probs)
    } else if (pred_type == "poisson") {
      rates <- exp(linear_pred)
      covars[[pred_name]] <- rpois(n_cases, rates)
    } else if (grepl("^multinomial", pred_type)) {
      n_cats <- as.numeric(gsub("multinomial", "", pred_type))
      categories <- LETTERS[1:n_cats]

      # Generate K-1 liabilities
      n_liab <- n_cats - 1
      liabilities <- matrix(0, nrow = n_cases, ncol = n_liab)

      for (k in 1:n_liab) {
        # Each liability has slightly different intercept
        liab_intercept <- intercepts$predictors[vari] + (k - n_liab / 2) * 0.5

        if (length(which_deps) == 0) {
          liab_linear <- liab_intercept
        } else {
          liab_linear <- liab_intercept + linear_pred - intercepts$predictors[vari]
        }

        # Add liability-specific noise
        liabilities[, k] <- liab_linear + rnorm(n_cases, 0, 0.5)
      }

      covars[[pred_name]] <- mnom_liab2cat(liabilities, categories)
    } else if (grepl("^ordinal", pred_type)) {
      n_cats <- as.numeric(gsub("ordinal", "", pred_type))
      
      # Ordinal uses single latent variable with thresholds
      covars[[pred_name]] <- ordinal_liab2cat(linear_pred, n_cats)
    }
  }

  # -------------------------------------------------------------------------
  # SIMULATE RESPONSE
  # -------------------------------------------------------------------------

  resp_name <- var_names[1]

  # Build design matrix for response
  X_resp_formula <- paste("~ 1 +", paste(var_names[-1], collapse = " + "))

  # Handle factor predictors for model matrix
  covars_for_X <- covars
  for (i in seq_len(n_predictors)) {
    if (grepl("^multinomial", predictor_types[i])) {
      covars_for_X[[var_names[i + 1]]] <- factor(covars_for_X[[var_names[i + 1]]])
    }
    # Note: ordinal predictors are kept as integers (treated as numeric in regression)
  }

  X_resp <- model.matrix(as.formula(X_resp_formula), covars_for_X)

  # Adjust beta_resp for factor expansion
  beta_resp_full <- intercepts$response
  for (i in seq_len(n_predictors)) {
    if (grepl("^multinomial", predictor_types[i])) {
      n_cats <- as.numeric(gsub("multinomial", "", predictor_types[i]))
      # Spread the effect across categories
      beta_resp_full <- c(
        beta_resp_full,
        seq(-abs(beta_resp[i]) / 2, abs(beta_resp[i]) / 2,
          length.out = n_cats - 1
        ) * sign(beta_resp[i])
      )
    } else {
      beta_resp_full <- c(beta_resp_full, beta_resp[i])
    }
  }

  # Ensure correct length
  if (length(beta_resp_full) != ncol(X_resp)) {
    beta_resp_full <- c(
      intercepts$response,
      rep(beta_resp, length.out = ncol(X_resp) - 1)
    )
  }

  # Linear predictor for response
  linear_pred_resp <- as.numeric(X_resp %*% beta_resp_full) +
    as.numeric(Z %*% u_species[[1]]) +
    as.numeric(Z %*% u_phylo[[1]]) +
    residuals[[1]]

  # Add random slopes for response
  if (rr && "species" %in% names(u_slopes)) {
    for (cov in names(u_slopes$species)) {
      if (cov %in% names(covars)) {
        cov_vals <- if (is.numeric(covars[[cov]])) covars[[cov]] else as.numeric(factor(covars[[cov]]))
        slope_contribution <- as.numeric(Z %*% u_slopes$species[[cov]]) * cov_vals
        linear_pred_resp <- linear_pred_resp + slope_contribution
      }
    }
  }

  # Generate response based on type
  if (response_type == "gaussian") {
    response <- linear_pred_resp
  } else if (response_type == "binary") {
    probs <- plogis(linear_pred_resp)
    response <- rbinom(n_cases, 1, probs)
  } else if (response_type == "poisson") {
    rates <- exp(linear_pred_resp)
    response <- rpois(n_cases, rates)
  } else if (grepl("^multinomial", response_type)) {
    n_cats <- as.numeric(gsub("multinomial", "", response_type))
    categories <- LETTERS[1:n_cats]

    n_liab <- n_cats - 1
    liabilities <- matrix(0, nrow = n_cases, ncol = n_liab)

    for (k in 1:n_liab) {
      liab_intercept <- intercepts$response + (k - n_liab / 2) * 0.5
      liabilities[, k] <- liab_intercept + linear_pred_resp - intercepts$response +
        rnorm(n_cases, 0, 0.5)
    }

    response <- mnom_liab2cat(liabilities, categories)
  } else if (grepl("^ordinal", response_type)) {
    n_cats <- as.numeric(gsub("ordinal", "", response_type))
    
    # Ordinal uses single latent variable with thresholds
    response <- ordinal_liab2cat(linear_pred_resp, n_cats)
  } else {
    stop("Unknown response_type: ", response_type)
  }

  # -------------------------------------------------------------------------
  # APPLY MISSINGNESS
  # -------------------------------------------------------------------------

  response <- apply_missingness(response, missingness[1])

  for (i in seq_len(n_predictors)) {
    covars[[var_names[i + 1]]] <- apply_missingness(
      covars[[var_names[i + 1]]],
      missingness[i + 1]
    )
  }

  # -------------------------------------------------------------------------
  # ASSEMBLE OUTPUT
  # -------------------------------------------------------------------------

  # Create output data frame
  out_data <- data.frame(
    species = case_species,
    response = response,
    covars,
    stringsAsFactors = FALSE
  )
  names(out_data)[2] <- resp_name

  # Store parameters used
  params <- list(
    response_type = response_type,
    predictor_types = predictor_types,
    var_names = var_names,
    beta_matrix = beta_matrix,
    beta_resp = beta_resp,
    intercepts = intercepts,
    phylo_signal = phylo_signal,
    sigmas = sigmas,
    n_cases = n_cases,
    n_species = n_species,
    n_species_actual = n_species_actual,
    birth = birth,
    death = death,
    missingness = missingness,
    rr = rr,
    rr_form = rr_form,
    design = mydesign
  )

  # Store random effects
  random_effects <- list(
    u_species = u_species,
    u_phylo = u_phylo,
    u_slopes = u_slopes,
    residuals = residuals
  )

  # Warn about integer representation of categorical variables
  int_vars <- c()
  if (response_type == "binary" || grepl("^ordinal", response_type)) {
    int_vars <- c(int_vars, var_names[1])
  }
  for (i in seq_len(n_predictors)) {
    if (predictor_types[i] == "binary" || grepl("^ordinal", predictor_types[i])) {
      int_vars <- c(int_vars, var_names[i + 1])
    }
  }
  if (length(int_vars) > 0) {
    message("Note: Variables '", paste(int_vars, collapse = "', '"), 
            "' are categorical but represented as integers. ",
            "Consider converting to factor if needed for analysis.")
  }

  return(list(
    data = out_data,
    tree = tree,
    params = params,
    random_effects = random_effects
  ))
}


# =============================================================================
# CONVENIENCE WRAPPER FUNCTIONS
# =============================================================================

#' @title Quick Gaussian simulation
#' @description Simplified wrapper for gaussian response simulation
#' @param n_predictors Number of gaussian predictors
#' @param n_cases Number of observations
#' @param n_species Number of species
#' @param phylo_signal Single phylogenetic signal value applied to all variables
#' @return simBACE output list
#' @export
simBACE_gaussian <- function(n_predictors = 3, n_cases = 200, n_species = 75,
                             phylo_signal = 0, beta_sparsity = 0.7) {
  simBACE(
    response_type = "gaussian",
    predictor_types = rep("gaussian", n_predictors),
    phylo_signal = rep(phylo_signal, n_predictors + 1),
    n_cases = n_cases,
    n_species = n_species,
    beta_sparsity = beta_sparsity
  )
}

#' @title Quick Poisson simulation
#' @description Simplified wrapper for count response simulation
#' @param n_predictors Number of gaussian predictors
#' @param n_cases Number of observations
#' @param n_species Number of species
#' @param phylo_signal Single phylogenetic signal value applied to all variables
#' @return simBACE output list
#' @export
simBACE_poisson <- function(n_predictors = 3, n_cases = 200, n_species = 75,
                            phylo_signal = 0) {
  simBACE(
    response_type = "poisson",
    predictor_types = rep("gaussian", n_predictors),
    phylo_signal = rep(phylo_signal, n_predictors + 1),
    n_cases = n_cases,
    n_species = n_species
  )
}

#' @title Quick Binary simulation
#' @description Simplified wrapper for binary response simulation
#' @param n_predictors Number of gaussian predictors
#' @param n_cases Number of observations
#' @param n_species Number of species
#' @param phylo_signal Single phylogenetic signal value applied to all variables
#' @return simBACE output list
#' @export
simBACE_binary <- function(n_predictors = 3, n_cases = 200, n_species = 75,
                           phylo_signal = 0) {
  simBACE(
    response_type = "binary",
    predictor_types = rep("gaussian", n_predictors),
    phylo_signal = rep(phylo_signal, n_predictors + 1),
    n_cases = n_cases,
    n_species = n_species
  )
}

#' @title Print summary of simBACE output
#' @description Prints a summary of the simulated data
#' @param sim_output Output from simBACE function
#' @export
print_simBACE_summary <- function(sim_output) {
  cat("=== simBACE Simulation Summary ===\n\n")

  cat("Response type:", sim_output$params$response_type, "\n")
  cat("Predictor types:", paste(sim_output$params$predictor_types, collapse = ", "), "\n")
  cat("Variable names:", paste(sim_output$params$var_names, collapse = ", "), "\n\n")

  cat("Sample size:", sim_output$params$n_cases, "observations\n")
  cat(
    "Species:", sim_output$params$n_species_actual,
    "(requested:", sim_output$params$n_species, ")\n"
  )
  cat("Tree tips:", ape::Ntip(sim_output$tree), "\n\n")

  cat("Phylogenetic signal:", paste(round(sim_output$params$phylo_signal, 3),
    collapse = ", "
  ), "\n")
  cat("Missingness:", paste(round(sim_output$params$missingness * 100, 1), "%",
    collapse = ", "
  ), "\n\n")

  cat("Beta matrix:\n")
  print(sim_output$params$beta_matrix)
  cat("\n")

  cat("\nBeta response coefficients:\n")
  print(sim_output$params$beta_resp)
  cat("\n")

  if (sim_output$params$rr) {
    cat("Random slopes: YES\n")
    for (re in names(sim_output$params$rr_form)) {
      cat("  -", re, ":", paste(sim_output$params$rr_form[[re]], collapse = ", "), "\n")
    }
  } else {
    cat("Random slopes: NO\n")
  }

  cat("\n--- Data preview ---\n")
  print(head(sim_output$data))

  cat("\n--- Data summary ---\n")
  print(summary(sim_output$data))

  # Warn about integer representation of categorical variables
  int_vars <- c()
  if (sim_output$params$response_type == "binary" || 
      grepl("^ordinal", sim_output$params$response_type)) {
    int_vars <- c(int_vars, sim_output$params$var_names[1])
  }
  for (i in seq_along(sim_output$params$predictor_types)) {
    ptype <- sim_output$params$predictor_types[i]
    if (ptype == "binary" || grepl("^ordinal", ptype)) {
      int_vars <- c(int_vars, sim_output$params$var_names[i + 1])
    }
  }
  if (length(int_vars) > 0) {
    cat("\n*** NOTE ***\n")
    cat("Variables '", paste(int_vars, collapse = "', '"), 
        "' are categorical but represented as integers.\n", sep = "")
    cat("Consider converting to factor if needed for analysis.\n")
  }

  invisible(sim_output)
}
