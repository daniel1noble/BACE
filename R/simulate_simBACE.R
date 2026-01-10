#' @title sim_bace: Comprehensive simulator of mock data for BACE
#' @description sim_bace is a function to simulate mock data for BACE with
#'   deliberate control over all aspects of data generation. The function
#'   simulates data with phylogenetic structure, between-species variance,
#'   and complex dependencies between variables. All response types are
#'   simulated on the liability scale (normalized, mean=0) and then
#'   transformed as appropriate.
#' @author Szymek Drobniak
#' @author Daniel Noble
#' @author Shinichi Nakagawa
#' @importFrom ape rphylo vcv Ntip
#' @importFrom MASS mvrnorm ginv
#' @importFrom stats rnorm rbinom rpois plogis model.matrix as.formula
#' @importFrom stringi stri_rand_strings

# =============================================================================
# MAIN SIMULATION FUNCTION (REWRITTEN)
# =============================================================================

#' @param response_type Type of response variable: "gaussian", "poisson" (integer count),
#'   "binary", "thresholdK" (ordered factor with K categories), or "multinomialK"
#'   (unordered factor with K categories), where K is number of categories (e.g., "threshold3")
#' @param predictor_types Character vector specifying type of each predictor.
#'   Options: "gaussian", "poisson", "binary", "thresholdK", "multinomialK"
#' @param n_cases Total number of observations (default 200)
#' @param n_species Number of species for tree simulation (default 75)
#' @param birth Birth rate for phylogenetic tree simulation (default 0.8)
#' @param death Death rate for phylogenetic tree simulation (default 0.4)
#' @param missingness Numeric vector of missingness proportions for each variable
#'   (response + predictors). 0 = no missing data. If NULL, defaults to 0 for all.
#'   Example: c(0, 0.1, 0, 0.05) for 0% missing in y, 10% in x1, 0% in x2, 5% in x3
#' @param var_names Character vector of variable names (response + predictors).
#'   If NULL, defaults to c("y", "x1", "x2", ...)
#' @param random_formula Formula specifying random effects structure.
#'   Options: "~phylo" (default - phylogenetic random effect only),
#'   "~species" (non-phylogenetic between-species variance only),
#'   "~phylo+species" (both phylogenetic and non-phylogenetic species effects).
#'   Can also be a formula object.
#' @param rr Logical, whether to include random slopes for numerical/integer
#'   (gaussian or poisson) predictors (default FALSE)
#' @param rr_form List specifying which predictors have random slopes for each
#'   random effect. Format: list(raneffect = c("predictor1", "predictor2")).
#'   Example: list(phylo = c("x1", "x2"), species = c("x4")).
#'   Only used if rr = TRUE.
#' @param raneff Output from make_raneff() function defining random effect
#'   variance components. If NULL, auto-generated with default parameters
#'   (0% phylo, 30% species, remainder residual)
#' @param fixeff Output from make_fixeff() function defining fixed effect
#'   structure (formulas, betas, interactions, intercepts). If NULL,
#'   auto-generated based on sparsity parameter with weak to moderate effects.
#' @param sparsity For auto-generation only: proportion of potential dependencies
#'   to set to zero (default 0.7 = sparse dependencies)
#'
#' @return Named list containing:
#'   - data: data.frame with simulated data (species, response, predictors)
#'   - tree: phylogenetic tree (ape phylo object)
#'   - params: list of all simulation parameters used
#'   - random_effects: list of simulated random effect values
#'
#' @examples
#' # Basic gaussian simulation with defaults
#' sim <- sim_bace(
#'   response_type = "gaussian",
#'   predictor_types = c("gaussian", "gaussian"),
#'   n_cases = 100,
#'   n_species = 30
#' )
#'
#' # Custom fixed and random effects
#' my_raneff <- make_raneff(
#'   var_names = c("y", "x1", "x2"),
#'   phylo_frac = c(0.5, 0.3, 0.1),
#'   species_frac = c(0.2, 0.3, 0.4)
#' )
#'
#' my_fixeff <- make_fixeff(
#'   var_names = c("y", "x1", "x2"),
#'   predictor_types = c("gaussian", "gaussian"),
#'   formulas = list(y ~ x1 + x2, x2 ~ x1),
#'   betas = list(y = c(0.5, 0.3), x2 = 0.4)
#' )
#'
#' sim <- sim_bace(
#'   response_type = "gaussian",
#'   predictor_types = c("gaussian", "gaussian"),
#'   raneff = my_raneff,
#'   fixeff = my_fixeff,
#'   n_cases = 200,
#'   n_species = 50
#' )
#'
#' # With interactions
#' my_fixeff_ix <- make_fixeff(
#'   var_names = c("y", "x1", "x2", "x3"),
#'   predictor_types = c("gaussian", "gaussian", "gaussian"),
#'   formulas = list(y ~ x1 + x2 + x3),
#'   betas = list(y = c(0.5, 0.3, 0.2)),
#'   interactions = list(
#'     formulas = list(y ~ x1:x2 + x1:x3),
#'     strengths = list(y = c("x1:x2" = 0.5, "x1:x3" = 1.0))
#'   )
#' )
#'
#' sim <- sim_bace(
#'   response_type = "gaussian",
#'   predictor_types = c("gaussian", "gaussian", "gaussian"),
#'   fixeff = my_fixeff_ix,
#'   n_cases = 200
#' )
#'
#' # Mixed variable types with random slopes
#' sim <- sim_bace(
#'   response_type = "poisson",
#'   predictor_types = c("gaussian", "binary", "threshold3"),
#'   random_formula = "~phylo+species",
#'   rr = TRUE,
#'   rr_form = list(phylo = c("x1"), species = c("x1")),
#'   n_cases = 300,
#'   n_species = 100
#' )
#'
#' @export

sim_bace <- function(
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
    sparsity = 0.7) {

  # -------------------------------------------------------------------------
  # SETUP AND VALIDATION
  # -------------------------------------------------------------------------

  n_predictors <- length(predictor_types)
  n_vars <- n_predictors + 1  # response + predictors

  # Generate default variable names if not provided
  if (is.null(var_names)) {
    var_names <- var_name_gen(n_predictors)
    message("Generated variable names: ", paste(var_names, collapse = ", "))
  }
  if (length(var_names) != n_vars) {
    stop("var_names must have length equal to number of predictors + 1 (response)")
  }

  # Parse random effect structure
  random_structure <- .parse_random_formula(random_formula)

  # Auto-generate raneff if not provided
  if (is.null(raneff)) {
    # First create basic raneff without random slopes
    phylo_frac_vals <- if (random_structure$has_phylo) rep(0, n_vars) else rep(0, n_vars)
    species_frac_vals <- if (random_structure$has_species) rep(0.3, n_vars) else rep(0, n_vars)
    
    raneff <- make_raneff(
      var_names = var_names,
      phylo_frac = phylo_frac_vals,
      species_frac = species_frac_vals
    )
    
    # Now add random slope variances if needed
    if (rr && !is.null(rr_form)) {
      rr_var_list <- list()
      for (re_name in names(rr_form)) {
        if (re_name %in% c("phylo", "species") && 
            ((re_name == "phylo" && random_structure$has_phylo) || 
             (re_name == "species" && random_structure$has_species))) {
          # Set random slope variance to 0.25 * mean random intercept variance
          re_vars <- if (re_name == "phylo") raneff$phylo_var else raneff$species_var
          mean_var <- mean(re_vars[re_vars > 0])
          if (is.na(mean_var) || mean_var == 0) mean_var <- 0.3 * 0.25  # Default
          
          rr_var_list[[re_name]] <- setNames(
            rep(0.25 * mean_var, length(rr_form[[re_name]])),
            rr_form[[re_name]]
          )
        }
      }
      if (length(rr_var_list) > 0) {
        raneff$rr_var <- rr_var_list
      }
    }
  }

  # Validate raneff structure
  if (!inherits(raneff, "raneff")) {
    stop("raneff must be output from make_raneff() function")
  }
  if (!identical(raneff$var_names, var_names)) {
    stop("raneff$var_names does not match var_names")
  }

  # Auto-generate fixeff if not provided
  if (is.null(fixeff)) {
    fixeff <- make_fixeff(
      var_names = var_names,
      predictor_types = predictor_types,
      sparsity = sparsity
    )
  }

  # Validate fixeff structure
  if (!inherits(fixeff, "fixeff")) {
    stop("fixeff must be output from make_fixeff() function")
  }
  if (!identical(fixeff$var_names, var_names)) {
    stop("fixeff$var_names does not match var_names")
  }

  # Setup missingness (default: no missing data)
  if (is.null(missingness)) {
    missingness <- rep(0, n_vars)
    names(missingness) <- var_names
  } else if (length(missingness) != n_vars) {
    stop("missingness must have length equal to n_vars (response + predictors)")
  } else {
    names(missingness) <- var_names
  }

  # Validate random slopes specification
  if (rr && is.null(rr_form)) {
    warning("rr=TRUE but rr_form not specified. No random slopes will be generated.")
    rr <- FALSE
  }
  
  # Validate rr_form against random_structure
  if (rr && !is.null(rr_form)) {
    for (re_name in names(rr_form)) {
      if (re_name == "phylo" && !random_structure$has_phylo) {
        warning("rr_form includes 'phylo' but random_formula does not. Ignoring phylo random slopes.")
        rr_form[[re_name]] <- NULL
      }
      if (re_name == "species" && !random_structure$has_species) {
        warning("rr_form includes 'species' but random_formula does not. Ignoring species random slopes.")
        rr_form[[re_name]] <- NULL
      }
    }
    # Remove empty elements
    rr_form <- rr_form[sapply(rr_form, length) > 0]
    if (length(rr_form) == 0) {
      rr <- FALSE
    }
  }

  # -------------------------------------------------------------------------
  # GENERATE PHYLOGENETIC TREE
  # -------------------------------------------------------------------------

  taxa <- sim_tree(n_species, birth, death, n_cases)
  tree <- taxa$tree
  case_species <- taxa$cases
  species_list <- taxa$species
  n_species_actual <- length(species_list)

  # Phylogenetic correlation matrix (only if needed)
  cor_phylo <- if (random_structure$has_phylo) {
    ape::vcv(tree, corr = TRUE)
  } else {
    NULL
  }

  # -------------------------------------------------------------------------
  # SAMPLE RANDOM EFFECTS
  # -------------------------------------------------------------------------

  # Z-matrix for species (maps observations to species)
  Z <- model.matrix(
    ~ 0 + factor(species),
    data.frame(species = case_species)
  )
  colnames(Z) <- species_list

  # Sample phylogenetic random effects
  u_phylo <- if (random_structure$has_phylo) {
    lapply(raneff$phylo_var, function(var) {
      if (var > 0) {
        sample_random_effects(var, n_species_actual, cor_matrix = cor_phylo)
      } else {
        rep(0, n_species_actual)
      }
    })
  } else {
    lapply(seq_len(n_vars), function(i) rep(0, n_species_actual))
  }
  names(u_phylo) <- var_names

  # Sample non-phylogenetic species random effects
  u_species <- if (random_structure$has_species) {
    lapply(raneff$species_var, function(var) {
      if (var > 0) {
        sample_random_effects(var, n_species_actual, cor_matrix = NULL)
      } else {
        rep(0, n_species_actual)
      }
    })
  } else {
    lapply(seq_len(n_vars), function(i) rep(0, n_species_actual))
  }
  names(u_species) <- var_names

  # Sample residuals
  residuals <- lapply(raneff$residual_var, function(var) {
    rnorm(n_cases, mean = 0, sd = sqrt(var))
  })
  names(residuals) <- var_names

  # Sample random slopes (if specified)
  u_slopes <- list()
  ignored_rr_covars <- c()

  if (rr && !is.null(rr_form)) {
    # Identify continuous predictors (gaussian, poisson) - valid for random slopes
    predictor_names <- var_names[-1]
    continuous_predictors <- predictor_names[predictor_types %in% c("gaussian", "poisson")]

    for (re_name in names(rr_form)) {
      covars_with_slopes <- rr_form[[re_name]]
      u_slopes[[re_name]] <- list()

      for (cov in covars_with_slopes) {
        # Check if covariate is continuous
        if (cov %in% continuous_predictors) {
          # Get variance from raneff$rr_var if available
          if (!is.null(raneff$rr_var) && !is.null(raneff$rr_var[[re_name]]) &&
              !is.null(raneff$rr_var[[re_name]][[cov]])) {
            slope_var <- raneff$rr_var[[re_name]][[cov]]
          } else {
            # Default: 0.25 * corresponding random intercept variance
            if (re_name == "phylo" && random_structure$has_phylo) {
              slope_var <- 0.25 * mean(raneff$phylo_var[raneff$phylo_var > 0])
            } else if (re_name == "species" && random_structure$has_species) {
              slope_var <- 0.25 * mean(raneff$species_var[raneff$species_var > 0])
            } else {
              slope_var <- 0
            }
          }

          if (slope_var > 0) {
            if (re_name == "phylo" && random_structure$has_phylo) {
              u_slopes[[re_name]][[cov]] <- sample_random_effects(
                slope_var, n_species_actual, cor_matrix = cor_phylo
              )
            } else if (re_name == "species" && random_structure$has_species) {
              u_slopes[[re_name]][[cov]] <- sample_random_effects(
                slope_var, n_species_actual, cor_matrix = NULL
              )
            }
          }
        } else {
          # Track ignored non-continuous covariates
          ignored_rr_covars <- c(ignored_rr_covars, cov)
        }
      }
    }

    # Issue warning for ignored covariates
    if (length(ignored_rr_covars) > 0) {
      warning(
        "Random slopes ignored for non-continuous covariates: '",
        paste(unique(ignored_rr_covars), collapse = "', '"),
        "'. Random slopes are only applied to gaussian and poisson predictors."
      )
    }
  }

  # -------------------------------------------------------------------------
  # DETERMINE SIMULATION ORDER
  # -------------------------------------------------------------------------

  sim_order <- .determine_sim_order(fixeff$formulas, var_names, predictor_types)

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

  # Store liabilities for categorical variables (needed for subsequent dependencies)
  liabilities <- list()

  # Store linear model components for each variable
  linear_models <- list()

  # Simulate predictors in dependency order
  for (pred_idx in sim_order) {
    pred_type <- predictor_types[pred_idx]
    pred_name <- var_names[pred_idx + 1]  # +1 because var_names includes response
    var_idx <- pred_idx + 1  # Index in full variable list (response = 1)

    # Find formula for this predictor
    pred_formula <- NULL
    for (form in fixeff$formulas) {
      if (all.vars(form)[1] == pred_name) {
        pred_formula <- form
        break
      }
    }

    # Build linear predictor
    if (is.null(pred_formula)) {
      # Independent predictor (no dependencies)
      # Create intercept-only design matrix
      X <- matrix(1, nrow = n_cases, ncol = 1)
      colnames(X) <- "(Intercept)"
      beta_vec <- c(fixeff$intercepts[pred_name])
      linear_pred <- rep(fixeff$intercepts[pred_name], n_cases)
    } else {
      # Dependent predictor
      rhs_vars <- all.vars(pred_formula)[-1]  # Predictor variables

      # Check if all rhs_vars exist in covars
      missing_vars <- setdiff(rhs_vars, names(covars))
      if (length(missing_vars) > 0) {
        stop("Variables ", paste(missing_vars, collapse = ", "), 
             " not yet simulated but required for ", pred_name,
             ". Check simulation order or formula dependencies.")
      }

      # Build design matrix
      X_formula_str <- paste("~ 1 +", paste(rhs_vars, collapse = " + "))
      
      # Prepare data for model matrix (convert factors properly)
      covars_for_X <- covars[, rhs_vars, drop = FALSE]
      for (rv in rhs_vars) {
        if (is.factor(covars_for_X[[rv]])) {
          # Already a factor, keep as-is
        } else {
          # Check if should be factor based on predictor type
          rv_idx <- which(var_names[-1] == rv)
          if (length(rv_idx) > 0) {
            rv_type <- predictor_types[rv_idx]
            if (grepl("^multinomial", rv_type)) {
              covars_for_X[[rv]] <- factor(covars_for_X[[rv]])
            }
          }
        }
      }

      X <- model.matrix(as.formula(X_formula_str), covars_for_X)

      # Get betas for this predictor
      pred_betas <- fixeff$betas[[pred_name]]
      
      # Flatten betas list to vector matching X columns
      beta_vec <- c(fixeff$intercepts[pred_name])  # Intercept
      for (rv in rhs_vars) {
        if (is.list(pred_betas[[rv]])) {
          beta_vec <- c(beta_vec, unlist(pred_betas[[rv]]))
        } else if (!is.null(pred_betas[[rv]])) {
          beta_vec <- c(beta_vec, pred_betas[[rv]])
        }
      }

      # Verify dimensions match
      if (length(beta_vec) != ncol(X)) {
        warning("Beta length mismatch for ", pred_name, ": expected ", ncol(X),
                ", got ", length(beta_vec), ". Padding/truncating.")
        if (length(beta_vec) < ncol(X)) {
          beta_vec <- c(beta_vec, rep(0, ncol(X) - length(beta_vec)))
        } else {
          beta_vec <- beta_vec[seq_len(ncol(X))]
        }
      }

      linear_pred <- as.numeric(X %*% beta_vec)
    }

    # Add random effects
    linear_pred <- linear_pred +
      as.numeric(Z %*% u_phylo[[var_idx]]) +
      as.numeric(Z %*% u_species[[var_idx]]) +
      residuals[[var_idx]]

    # Add random slopes if applicable
    if (rr && !is.null(u_slopes)) {
      for (re_name in names(u_slopes)) {
        if (!is.null(u_slopes[[re_name]])) {
          for (cov in names(u_slopes[[re_name]])) {
            if (cov %in% names(covars) && !all(covars[[cov]] == 0)) {
              cov_vals <- if (is.numeric(covars[[cov]])) covars[[cov]] else as.numeric(covars[[cov]])
              slope_contribution <- as.numeric(Z %*% u_slopes[[re_name]][[cov]]) * cov_vals
              linear_pred <- linear_pred + slope_contribution
            }
          }
        }
      }
    }

    # Add interaction terms if applicable
    if (!is.null(fixeff$interactions) && !is.null(fixeff$interactions$formulas)) {
      for (ix_form in fixeff$interactions$formulas) {
        if (all.vars(ix_form)[1] == pred_name) {
          # Extract interaction terms from formula
          ix_terms <- attr(terms(ix_form), "term.labels")
          ix_terms <- ix_terms[grepl(":", ix_terms)]  # Only interaction terms

          for (ix_term in ix_terms) {
            ix_vars <- strsplit(ix_term, ":")[[1]]
            if (all(ix_vars %in% names(covars)) && !all(covars[[ix_vars[1]]] == 0) &&
                !all(covars[[ix_vars[2]]] == 0)) {
              ix_val <- calculate_ix_term(covars, ix_vars)
              ix_strength <- if (!is.null(fixeff$interactions$strengths) &&
                                  !is.null(fixeff$interactions$strengths[[pred_name]]) &&
                                  !is.null(fixeff$interactions$strengths[[pred_name]][[ix_term]])) {
                fixeff$interactions$strengths[[pred_name]][[ix_term]]
              } else {
                0.5  # Default interaction strength
              }

              # Interaction effect as fraction of main effect
              # Get numeric beta values
              beta_vals <- fixeff$betas[[pred_name]]
              if (is.list(beta_vals)) {
                beta_nums <- unlist(lapply(beta_vals, function(b) {
                  if (is.numeric(b)) b else NULL
                }))
              } else {
                beta_nums <- beta_vals
              }
              main_effect <- if (length(beta_nums) > 0) mean(abs(beta_nums)) else 0.3
              ix_coef <- ix_strength * main_effect

              linear_pred <- linear_pred + ix_coef * ix_val
            }
          }
        }
      }
    }

    # Normalize liability
    # For categorical variables: mean=0, sd=1 (needed for identification)
    # For gaussian/poisson: only center to mean=0 (preserve variance for parameter recovery)
    if (pred_type %in% c("binary") || grepl("^threshold", pred_type) || grepl("^multinomial", pred_type)) {
      # Categorical: normalize to sd=1 (standard for latent variable models)
      linear_pred <- scale(linear_pred, center = TRUE, scale = TRUE)[, 1]
    } else {
      # Gaussian/Poisson: only center, preserve variance
      linear_pred <- scale(linear_pred, center = TRUE, scale = FALSE)[, 1]
    }

    # Store liability
    liabilities[[pred_name]] <- linear_pred

    # Determine link function for this variable type
    link_function <- if (pred_type == "gaussian") {
      "identity"
    } else if (pred_type == "poisson") {
      "log"
    } else if (pred_type == "binary" || grepl("^threshold", pred_type) || grepl("^multinomial", pred_type)) {
      "probit"  # Using probit link for liability threshold models
    } else {
      "identity"
    }

    # Store linear model components
    # Note: For multinomial variables (K>2), the full model involves K-1 latent
    # liabilities stored as a matrix. Currently we only store scalar model components.
    # For multinomial predictors, only the structural information is stored.
    linear_models[[pred_name]] <- list(
      link = link_function,
      X = X,
      beta = beta_vec,
      Z = Z,
      u = u_phylo[[pred_idx]] + u_species[[pred_idx]],  # Combined random effects
      e = residuals[[pred_idx]],  # Residual error
      note = if (grepl("^multinomial", pred_type)) {
        paste0("Multinomial variable with K=", gsub("multinomial", "", pred_type),
               " requires K-1 latent liabilities (matrix structure)")
      } else NULL
    )

    # Convert liability to appropriate scale/type
    covars[[pred_name]] <- .liability_to_categories(linear_pred, pred_type, pred_name)
  }

  # -------------------------------------------------------------------------
  # SIMULATE RESPONSE
  # -------------------------------------------------------------------------

  resp_name <- var_names[1]

  # Find formula for response
  resp_formula <- NULL
  for (form in fixeff$formulas) {
    if (all.vars(form)[1] == resp_name) {
      resp_formula <- form
      break
    }
  }

  if (is.null(resp_formula)) {
    stop("No formula found for response variable '", resp_name, 
         "'. fixeff$formulas must include a formula for the response.")
  }

  # Build design matrix for response
  rhs_vars <- all.vars(resp_formula)[-1]
  X_formula_str <- paste("~ 1 +", paste(rhs_vars, collapse = " + "))

  # Prepare data for model matrix
  covars_for_X_resp <- covars[, rhs_vars, drop = FALSE]
  for (rv in rhs_vars) {
    if (is.factor(covars_for_X_resp[[rv]])) {
      # Already a factor
    } else {
      rv_idx <- which(var_names[-1] == rv)
      if (length(rv_idx) > 0) {
        rv_type <- predictor_types[rv_idx]
        if (grepl("^multinomial", rv_type)) {
          covars_for_X_resp[[rv]] <- factor(covars_for_X_resp[[rv]])
        }
      }
    }
  }

  X_resp <- model.matrix(as.formula(X_formula_str), covars_for_X_resp)

  # Get betas for response
  resp_betas <- fixeff$betas[[resp_name]]

  # Flatten betas to vector
  beta_vec_resp <- c(fixeff$intercepts[resp_name])
  for (rv in rhs_vars) {
    if (is.list(resp_betas[[rv]])) {
      beta_vec_resp <- c(beta_vec_resp, unlist(resp_betas[[rv]]))
    } else {
      beta_vec_resp <- c(beta_vec_resp, resp_betas[[rv]])
    }
  }

  # Verify dimensions
  if (length(beta_vec_resp) != ncol(X_resp)) {
    warning("Beta length mismatch for response: expected ", ncol(X_resp),
            ", got ", length(beta_vec_resp), ". Padding/truncating.")
    if (length(beta_vec_resp) < ncol(X_resp)) {
      beta_vec_resp <- c(beta_vec_resp, rep(0, ncol(X_resp) - length(beta_vec_resp)))
    } else {
      beta_vec_resp <- beta_vec_resp[seq_len(ncol(X_resp))]
    }
  }

  linear_pred_resp <- as.numeric(X_resp %*% beta_vec_resp) +
    as.numeric(Z %*% u_phylo[[1]]) +
    as.numeric(Z %*% u_species[[1]]) +
    residuals[[1]]

  # Add random slopes for response
  if (rr && !is.null(u_slopes)) {
    for (re_name in names(u_slopes)) {
      if (!is.null(u_slopes[[re_name]])) {
        for (cov in names(u_slopes[[re_name]])) {
          if (cov %in% rhs_vars) {
            cov_vals <- if (is.numeric(covars[[cov]])) covars[[cov]] else as.numeric(covars[[cov]])
            slope_contribution <- as.numeric(Z %*% u_slopes[[re_name]][[cov]]) * cov_vals
            linear_pred_resp <- linear_pred_resp + slope_contribution
          }
        }
      }
    }
  }

  # Add interaction terms for response
  if (!is.null(fixeff$interactions) && !is.null(fixeff$interactions$formulas)) {
    for (ix_form in fixeff$interactions$formulas) {
      if (all.vars(ix_form)[1] == resp_name) {
        ix_terms <- attr(terms(ix_form), "term.labels")
        ix_terms <- ix_terms[grepl(":", ix_terms)]

        for (ix_term in ix_terms) {
          ix_vars <- strsplit(ix_term, ":")[[1]]
          if (all(ix_vars %in% names(covars))) {
            ix_val <- calculate_ix_term(covars, ix_vars)
            ix_strength <- if (!is.null(fixeff$interactions$strengths) &&
                                !is.null(fixeff$interactions$strengths[[resp_name]]) &&
                                !is.null(fixeff$interactions$strengths[[resp_name]][[ix_term]])) {
              fixeff$interactions$strengths[[resp_name]][[ix_term]]
            } else {
              0.5
            }

            # Get numeric beta values
            beta_vals <- fixeff$betas[[resp_name]]
            if (is.list(beta_vals)) {
              beta_nums <- unlist(lapply(beta_vals, function(b) {
                if (is.numeric(b)) b else NULL
              }))
            } else {
              beta_nums <- beta_vals
            }
            main_effect <- if (length(beta_nums) > 0) mean(abs(beta_nums)) else 0.3
            ix_coef <- ix_strength * main_effect

            linear_pred_resp <- linear_pred_resp + ix_coef * ix_val
          }
        }
      }
    }
  }

  # Normalize response liability
  # For categorical: mean=0, sd=1 (needed for identification)
  # For gaussian/poisson: only center to mean=0 (preserve variance for parameter recovery)
  if (response_type %in% c("binary") || grepl("^threshold", response_type) || grepl("^multinomial", response_type)) {
    # Categorical: normalize to sd=1 (standard for latent variable models)
    linear_pred_resp <- scale(linear_pred_resp, center = TRUE, scale = TRUE)[, 1]
  } else {
    # Gaussian/Poisson: only center, preserve variance
    linear_pred_resp <- scale(linear_pred_resp, center = TRUE, scale = FALSE)[, 1]
  }

  # Store response liability
  liabilities[[resp_name]] <- linear_pred_resp

  # Determine link function for response
  link_function_resp <- if (response_type == "gaussian") {
    "identity"
  } else if (response_type == "poisson") {
    "log"
  } else if (response_type == "binary" || grepl("^threshold", response_type) || grepl("^multinomial", response_type)) {
    "probit"  # Using probit link for liability threshold models
  } else {
    "identity"
  }

  # Store linear model components for response
  # Note: For multinomial responses (K>2), the full model involves K-1 latent
  # liabilities stored as a matrix. Currently we only store scalar model components.
  linear_models[[resp_name]] <- list(
    link = link_function_resp,
    X = X_resp,
    beta = beta_vec_resp,
    Z = Z,
    u = u_phylo[[1]] + u_species[[1]],  # Combined random effects
    e = residuals[[1]],  # Residual error
    note = if (grepl("^multinomial", response_type)) {
      paste0("Multinomial variable with K=", gsub("multinomial", "", response_type),
             " requires K-1 latent liabilities (matrix structure)")
    } else NULL
  )

  # Convert liability to appropriate scale/type
  response <- .liability_to_categories(linear_pred_resp, response_type, resp_name)

  # -------------------------------------------------------------------------
  # APPLY MISSINGNESS
  # -------------------------------------------------------------------------

  response <- apply_missingness(response, missingness[resp_name])

  for (pred_name in var_names[-1]) {
    covars[[pred_name]] <- apply_missingness(
      covars[[pred_name]],
      missingness[pred_name]
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
    n_cases = n_cases,
    n_species = n_species,
    n_species_actual = n_species_actual,
    birth = birth,
    death = death,
    missingness = missingness,
    random_formula = random_formula,
    random_structure = random_structure,
    rr = rr,
    rr_form = rr_form,
    raneff = raneff,
    fixeff = fixeff,
    sparsity = sparsity,
    sim_order = sim_order
  )

  # Store random effects
  random_effects <- list(
    u_phylo = u_phylo,
    u_species = u_species,
    u_slopes = u_slopes,
    residuals = residuals,
    liabilities = liabilities
  )

  output <- list(
    data = out_data,
    tree = tree,
    params = params,
    random_effects = random_effects,
    linear_models = linear_models
  )

  class(output) <- c("sim_bace", "list")
  return(output)
}


# =============================================================================
# CONVENIENCE WRAPPER FUNCTIONS
# =============================================================================

#' @title Quick Gaussian simulation
#' @description Simplified wrapper for gaussian response simulation
#' @param n_predictors Number of gaussian predictors
#' @param n_cases Number of observations
#' @param n_species Number of species
#' @param sparsity Proportion of dependencies to set to zero (default 0.7)
#' @return sim_bace output list
#' @export
sim_bace_gaussian <- function(n_predictors = 3, n_cases = 200, n_species = 75,
                               sparsity = 0.7) {
  sim_bace(
    response_type = "gaussian",
    predictor_types = rep("gaussian", n_predictors),
    n_cases = n_cases,
    n_species = n_species,
    sparsity = sparsity
  )
}

#' @title Quick Poisson simulation
#' @description Simplified wrapper for count response simulation
#' @param n_predictors Number of gaussian predictors
#' @param n_cases Number of observations
#' @param n_species Number of species
#' @param sparsity Proportion of dependencies to set to zero (default 0.7)
#' @return sim_bace output list
#' @export
sim_bace_poisson <- function(n_predictors = 3, n_cases = 200, n_species = 75,
                              sparsity = 0.7) {
  sim_bace(
    response_type = "poisson",
    predictor_types = rep("gaussian", n_predictors),
    n_cases = n_cases,
    n_species = n_species,
    sparsity = sparsity
  )
}

#' @title Quick Binary simulation
#' @description Simplified wrapper for binary response simulation
#' @param n_predictors Number of gaussian predictors
#' @param n_cases Number of observations
#' @param n_species Number of species
#' @param sparsity Proportion of dependencies to set to zero (default 0.7)
#' @return sim_bace output list
#' @export
sim_bace_binary <- function(n_predictors = 3, n_cases = 200, n_species = 75,
                             sparsity = 0.7) {
  sim_bace(
    response_type = "binary",
    predictor_types = rep("gaussian", n_predictors),
    n_cases = n_cases,
    n_species = n_species,
    sparsity = sparsity
  )
}

#' @title Print summary of sim_bace output
#' @description Prints a summary of the simulated data
#' @param x Output from sim_bace function
#' @param ... Additional arguments (unused)
#' @export
print.sim_bace <- function(x, ...) {
  cat("=== sim_bace Simulation Output ===\n\n")

  cat("Response type:", x$params$response_type, "\n")
  cat("Predictor types:", paste(x$params$predictor_types, collapse = ", "), "\n")
  cat("Variable names:", paste(x$params$var_names, collapse = ", "), "\n\n")

  cat("Sample size:", x$params$n_cases, "observations\n")
  cat("Species:", x$params$n_species_actual, 
      "(requested:", x$params$n_species, ")\n")
  cat("Tree tips:", ape::Ntip(x$tree), "\n\n")

  cat("Random effects structure:", as.character(x$params$random_formula)[2], "\n")
  if (x$params$random_structure$has_phylo) {
    cat("  - Phylogenetic variance fractions:",
        paste(round(x$params$raneff$phylo_frac, 3), collapse = ", "), "\n")
  }
  if (x$params$random_structure$has_species) {
    cat("  - Species variance fractions:",
        paste(round(x$params$raneff$species_frac, 3), collapse = ", "), "\n")
  }
  cat("  - Residual variance fractions:",
      paste(round(x$params$raneff$residual_var / x$params$raneff$total_var, 3),
            collapse = ", "), "\n\n")

  cat("Missingness:", paste(round(x$params$missingness * 100, 1), "%",
                             collapse = ", "), "\n\n")

  if (x$params$rr && !is.null(x$params$rr_form) && length(x$params$rr_form) > 0) {
    cat("Random slopes: YES\n")
    for (re in names(x$params$rr_form)) {
      cat("  -", re, ":", paste(x$params$rr_form[[re]], collapse = ", "), "\n")
    }
  } else {
    cat("Random slopes: NO\n")
  }

  # Display formulas
  cat("\nFixed effect formulas:\n")
  for (form in x$params$fixeff$formulas) {
    cat("  ", deparse(form), "\n")
  }

  # Display interaction terms
  if (!is.null(x$params$fixeff$interactions) && 
      !is.null(x$params$fixeff$interactions$formulas) &&
      length(x$params$fixeff$interactions$formulas) > 0) {
    cat("\nInteraction formulas:\n")
    for (ix_form in x$params$fixeff$interactions$formulas) {
      cat("  ", deparse(ix_form), "\n")
    }
  } else {
    cat("\nInteraction terms: NONE\n")
  }

  cat("\n--- Data preview ---\n")
  print(head(x$data))

  invisible(x)
}

#' @title Print summary of sim_bace output (alternative name)
#' @description Prints a summary of the simulated data
#' @param sim_output Output from sim_bace function
#' @export
print_sim_bace_summary <- function(sim_output) {
  print.sim_bace(sim_output)
  invisible(sim_output)
}
