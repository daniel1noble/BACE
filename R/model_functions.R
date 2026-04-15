#' @title .model_fit
#' @description Function takes incomplete data and a formula string and fits a model using MCMCglmm
#' @param data A dataframe containing missing data.
#' @param tree A phylogenetic tree of class 'phylo' from the ape package.
#' @param fixformula A string that specifies the fixed effect structure in the model.
#' @param randformula A formula or list of formulas that specifies the random effect structure in the model. If a list, should have elements named 'phylo' and 'species'.
#' @param type A string that specifies the type of model to fit. Options are "normal", "binary", "count", "categorical", "ordinal" and the appropriate family will be used in MCMCglmm.
#' @param nitt An integer specifying the number of iterations to run the MCMC algorithm for this specific model. Default is 6000. Note: when called from bace_imp, this can be model-specific if a list was provided.
#' @param thin An integer specifying the thinning rate for the MCMC algorithm for this specific model. Default is 5. Note: when called from bace_imp, this can be model-specific if a list was provided.
#' @param burnin An integer specifying the number of iterations to discard as burnin for this specific model. Default is 1000. Note: when called from bace_imp, this can be model-specific if a list was provided.
#' @param prior A list specifying the prior distributions for the MCMCglmm model.
#' @return A list of draws from the posterior distribution of the model parameters.
#' @importFrom methods as
#' @export

.model_fit <- function(data, tree, fixformula, randformula, type, prior, nitt = 6000, thin = 5, burnin = 1000) {

	# Create sparse matrix of phylogeny. Make sure to include all nodes because the algorithm used for nodes = "ALL" is more stable and orders of magnitude faster. Very important for large trees.
		   A  <- MCMCglmm::inverseA(tree, nodes = "ALL")$Ainv
	
	# Handle random formula structure - can be a single formula or a list with phylo and species
	if (is.list(randformula) && !inherits(randformula, "formula")) {
		# Dual random effects: phylo + species
		if (!all(c("phylo", "species") %in% names(randformula))) {
			stop("When randformula is a list, it must have elements named 'phylo' and 'species'")
		}
		
		# Extract grouping variable from formulas
		# For formulas like ~Species or ~us(1 + x):Species, need to get the grouping variable
		phylo_vars <- all.vars(randformula$phylo)
		species_vars <- all.vars(randformula$species)
		
		# The grouping variable is the last variable in the formula
		# For ~Species, it's Species
		# For ~us(1 + x):Species, it's Species (last in all.vars)
		phylo_name <- phylo_vars[length(phylo_vars)]
		species_name <- species_vars[length(species_vars)]
		
		# The second copy of the species column should already exist in data (added by bace_imp)
		species2_name <- paste0(species_name, "2")
		if (!species2_name %in% colnames(data)) {
			stop("When using dual random effects, data must contain a column named '", species2_name, "'")
		}
		
		# Create identity matrix for non-phylogenetic species effect
		# Use tree tip labels to ensure matrix dimensions match phylogeny
		species_levels <- tree$tip.label
		n_species <- length(species_levels)
		I_species <- diag(n_species)
		rownames(I_species) <- colnames(I_species) <- species_levels
		
		# Convert to sparse matrix format (dgCMatrix) as required by MCMCglmm
		I_species_sparse <- as(I_species, "dgCMatrix")
		
		# Combine random formulas
		# Need to handle different cases: intercept only vs random slopes
		phylo_formula_char <- deparse(randformula$phylo)
		species_formula_char <- deparse(randformula$species)
		
		# Check if phylo has random slopes (contains : or us() or idh())
		has_slopes <- grepl(":", phylo_formula_char) || grepl("us\\(", phylo_formula_char) || grepl("idh\\(", phylo_formula_char)
		
		if (has_slopes) {
			# Phylo has random slopes, species only intercept
			# Use original name for phylo, species2_name for non-phylo
			combined_formula <- stats::as.formula(paste(phylo_formula_char, "+", species2_name))
		} else {
			# Both are random intercepts only
			# Use original name for phylo, species2_name for non-phylo
			combined_formula <- stats::as.formula(paste("~", phylo_name, "+", species2_name))
		}
		
		# Set up ginverse with both matrices
		ginv_list <- list()
		ginv_list[[phylo_name]] <- A
		ginv_list[[species2_name]] <- I_species_sparse
		
		# Check if phylo has random slopes that need wider prior matrices
		# For us(1 + x):Species, we need a 2x2 matrix, for us(1 + x + y):Species, 3x3, etc.
		phylo_formula_char <- deparse(randformula$phylo)
		if (has_slopes) {
			# Extract the number of random effects from the formula
			# For us(1 + x):Species, we have intercept (1) + slope (x) = 2 dimensions
			# Extract variables inside us()
			us_match <- regexpr("us\\(([^)]+)\\)", phylo_formula_char, perl = TRUE)
			if (us_match > 0) {
				us_content <- regmatches(phylo_formula_char, us_match)
				# Extract content between us( and )
				us_vars <- gsub("us\\(|\\)", "", us_content)
				# Count terms: "1 + x" has 2 terms, "1 + x + y" has 3, etc.
				# Split by + and count
				terms <- strsplit(us_vars, "\\s*\\+\\s*")[[1]]
				n_phylo_dims <- length(terms)
			} else {
				# Fallback: assume 2 dimensions (intercept + 1 slope)
				n_phylo_dims <- 2
			}
			
			# Adjust the prior for the phylo random effect (first G element)
			# Need to expand V to n_phylo_dims x n_phylo_dims
			if (length(prior$G) >= 1 && nrow(prior$G[[1]]$V) == 1) {
				prior$G[[1]]$V <- diag(n_phylo_dims)
				if (!is.null(prior$G[[1]]$alpha.mu)) {
					# Parameter expansion case
					prior$G[[1]]$alpha.mu <- rep(0, n_phylo_dims)
					prior$G[[1]]$alpha.V <- 1e4 * diag(n_phylo_dims)
				}
			}
		}
		
		
	} else{
		# Single random effect (phylo only)
		combined_formula <- randformula
		phylo_name <- all.vars(randformula)[1]
		ginv_list <- setNames(list(A), phylo_name)
	}
	
	# Fit the model using MCMCglmm
  # slice = TRUE improves MCMC mixing for threshold/ordinal/categorical types
  use_slice <- type %in% c("threshold", "ordinal", "categorical")

  if(type != "categorical"){
  	model <- MCMCglmm::MCMCglmm(fixed = fixformula,
                               random = combined_formula,
                                 data = data,
                               family = type,
							               ginverse = ginv_list,
                              verbose = FALSE, pr = TRUE, pl = TRUE,
                                saveX = TRUE, saveZ = TRUE,
                                 nitt = nitt,
                                 thin = thin,
                               burnin = burnin,
							    prior = prior, singular.ok=TRUE,
							    slice = use_slice)
        } else {

      # Categorical model needs special treatment. Append -1 to the right side of ~ formula to remove intercept
      fixformula_cat <- as.formula(paste0(as.character(fixformula)[2], "~ trait - 1 + trait:(", as.character(fixformula)[3], ")"))
      
      # Handle random formula for categorical models
      # For categorical models, each random effect needs idh(trait): prefix
      if (is.list(randformula) && !inherits(randformula, "formula")) {
      	# Dual random effects: both need idh(trait) expansion
      	
      	# Check if phylo has random slopes - not supported for categorical with dual random effects
      	phylo_formula_char <- deparse(randformula$phylo)
      	has_slopes_cat <- grepl("us\\(", phylo_formula_char) || grepl("idh\\(", phylo_formula_char) || grepl("\\+", phylo_formula_char)
      	if (has_slopes_cat) {
      		stop("Categorical models with random slopes and dual random effects (species=TRUE) are not currently supported. Please use either:\n",
      		     "  1. Categorical model with phylo + species random intercepts (no random slopes), OR\n",
      		     "  2. Categorical model with phylo random slopes only (species=FALSE)")
      	}
      	
      	# Extract the random effect terms from combined_formula
      	combined_rhs <- as.character(combined_formula)[2]
      	
      	# Smart split: only split on + that are NOT inside parentheses
      	# Use a simple approach: identify the main grouping variables
      	# For "us(1 + x):Species + Species2", we want ["us(1 + x):Species", "Species2"]
      	# For "Species + Species2", we want ["Species", "Species2"]
      	
      	# Count parentheses to track depth
      	re_terms <- character()
      	current_term <- ""
      	paren_depth <- 0
      	
      	for (i in seq_len(nchar(combined_rhs))) {
      		char <- substr(combined_rhs, i, i)
      		if (char == "(") {
      			paren_depth <- paren_depth + 1
      			current_term <- paste0(current_term, char)
      		} else if (char == ")") {
      			paren_depth <- paren_depth - 1
      			current_term <- paste0(current_term, char)
      		} else if (char == "+" && paren_depth == 0) {
      			# This is a top-level +, so split here
      			re_terms <- c(re_terms, trimws(current_term))
      			current_term <- ""
      		} else {
      			current_term <- paste0(current_term, char)
      		}
      	}
      	# Add the last term
      	if (nchar(trimws(current_term)) > 0) {
      		re_terms <- c(re_terms, trimws(current_term))
      	}
      	
      	# Add idh(trait): to each term
      	re_terms_expanded <- paste0("idh(trait):", re_terms)
      	# Combine back
      	ranformula_cat <- as.formula(paste0("~", paste(re_terms_expanded, collapse = " + ")))
      } else {
      	# Single random effect: original behavior
      	ranformula_cat <- as.formula(paste0("~", "idh(trait):",as.character(combined_formula)[2]))
      }

      model <- MCMCglmm::MCMCglmm(fixed = fixformula_cat,
                                 random = ranformula_cat,
                                 data = data,
                               family = type,
                               rcov = ~us(trait):units,
							               ginverse = ginv_list,
                              verbose = FALSE, pr = TRUE, pl = TRUE,
                                saveX = TRUE, saveZ = TRUE,
                                 nitt = nitt,
                                 thin = thin,
                               burnin = burnin,
							    prior = prior, singular.ok=TRUE,
							    slice = use_slice)
                  }								
  	return(model)
}

#' @title .count_categorical_fixef
#' @description Helper function to count the number of fixed effect coefficients 
#'              in a categorical MCMCglmm model with trait expansion
#' @param fixform Formula for fixed effects (without trait expansion)
#' @param data Data frame containing predictor variables (but not 'trait')
#' @param n_levels Number of levels in the categorical response variable
#' @return Integer giving the total number of beta coefficients after trait expansion
#' @details For categorical models, MCMCglmm internally expands the formula to 
#'          J = n_levels - 1 traits. This function creates a temporary 'trait' variable,
#'          builds the expanded formula, and counts resulting coefficients.
#' @export
.count_categorical_fixef <- function(fixform, data, n_levels) {

  J <- n_levels - 1  # Number of non-baseline traits

  # Compute analytically: avoids the NA / missing-factor-level problem that
  # arises when model.matrix is called on a tiny subset of data that does not
  # include all factor levels.
  #
  # MCMCglmm categorical expansion:
  #   original:  response ~ x1 + x2 + ... (with intercept)
  #   expanded:  response ~ trait:(x1 + x2 + ...) - 1   (J traits, no overall intercept)
  #
  # Each trait gets:  1 (intercept) + sum(nlevels(factor) - 1 for factors, 1 for continuous)
  # Total = J * n_cols_per_trait

  # Predictor names only (drop the response, left-hand side)
  pred_vars <- all.vars(update(fixform, NULL ~ .))

  n_cols_per_trait <- 1L  # intercept term contributed by each trait level

  for (nm in pred_vars) {
    col <- data[[nm]]
    if (is.factor(col) || is.ordered(col)) {
      # Treatment / poly contrasts: nlevels - 1 columns
      n_cols_per_trait <- n_cols_per_trait + nlevels(col) - 1L
    } else {
      # Continuous or integer: 1 column
      n_cols_per_trait <- n_cols_per_trait + 1L
    }
  }

  return(J * n_cols_per_trait)
}

#' @title .list_of_G
#' @description Function creates a list of G priors for the MCMCglmm model
#' @param n_rand An integer specifying the number of random effects in the model.
#' @param nu A numeric specifying the nu parameter for the prior.
#' @param par_expand A logical indicating whether to use parameter expansion.
#' @param diag An integer specifying the dimension of the variance-covariance matrix.
#' @return A list of G priors for the MCMCglmm model.

.list_of_G <- function(n_rand, nu = NULL, par_expand = FALSE, diag = 1) {
  
  if(is.null(nu)) {nu <- 0.002}
  
  prior_G <- list()

  if (par_expand) {
    for (i in 1:n_rand) {
      prior_G[[i]] <- list(
        V = diag(diag), nu = nu,
        alpha.mu = rep(0, diag), alpha.V = 1e4 * diag(diag)
      )
    }
  } else {
    for (i in 1:n_rand) {
      prior_G[[i]] <- list(V = diag(diag), nu = nu)
    }
  }
  
  return(prior_G)
}

#' @title .make_prior
#' @description Function creates the prior for the MCMCglmm model
#' @param n_rand An integer specifying the number of random effects in the model (1 for phylo only, 2 for phylo + species).
#' @param type A string that specifies the type of model to fit.
#' @param nu A numeric specifying the nu parameter for the prior.
#' @param n_levels An integer specifying the number of levels for categorical or ordinal data.
#' @param par_expand A logical indicating whether to use parameter expansion.
#' @param fixform A formula object specifying the fixed effects structure (required for categorical models with Gelman prior).
#'    This is formula without trait expansion.
#' @param data A data frame containing the data (required for categorical models with Gelman prior).
#' @param gelman A logical indicating whether to use Gelman prior for fixed effects in categorical models. Default is TRUE.
#' @return A list of priors for the MCMCglmm model.
#' @importFrom stats complete.cases
#' @export

# TODO: Fix the gelamn prior definition for categorical models
.make_prior <- function(
    n_rand,
    type,
    nu = NULL,
    n_levels = NULL,
    par_expand = FALSE,
    fixform = NULL, # formula without trait expansion!
    data = NULL,
    gelman = 0) {
  
  # Validate n_rand
  if (!n_rand %in% c(1, 2)) {
    stop("n_rand must be 1 (phylo only) or 2 (phylo + species)")
  }
  
  if (type == "gaussian") {
    if (is.null(nu)) {
      nu <- 0.002
    }

    prior_G <- .list_of_G(n_rand, nu, par_expand)
    prior <- list(
      R = list(V = 1, nu = nu),
      G = prior_G
    )
  }

  if (type == "poisson") {
    if (is.null(nu)) {
      nu <- 0.02
    }
    prior_G <- .list_of_G(n_rand, nu, par_expand)
    prior <- list(
      R = list(V = 1, nu = nu),
      G = prior_G
    )
  }

  if (type == "categorical") {
    stopifnot(!is.null(n_levels))

    

    J <- n_levels - 1

    J_matrix <- array(1, dim = c(J, J)) # matrix of ones
    I_matrix <- diag(J) # identity matrix
    IJ <- (I_matrix + J_matrix) / n_levels # see Hadfield's Course notes p. 97
    prior_G <- .list_of_G(n_rand, nu = nu, par_expand, diag = J)

    prior <- list(
      R = list(V = IJ, fix = 1),
      G = prior_G
    )

    if(gelman) {
      # Gelman prior for fixed effects
      # For categorical models, MCMCglmm expands the formula internally to J traits
      # So we need J times the number of fixed effects (including intercept terms per trait)
      stopifnot(!is.null(fixform) | !is.null(data))
      # Remove rows with NA in variables used in the formula
      formula_vars <- all.vars(fixform)

      complete_data <- data[complete.cases(data[, formula_vars, drop = FALSE]), ]
      
      # Check if we have enough complete cases
      if(nrow(complete_data) < J | gelman == 2) {
        warning("Not enough complete cases to compute Gelman prior or gelman=2. Using pseudo-Gelman prior.")
        n_fixef_total <- .count_categorical_fixef(fixform, data, n_levels)
        prior_B <- list(
          mu = rep(0, n_fixef_total),
          V = diag(n_fixef_total) * (1 + pi^2 / 3)
        )
        
      } else {
        # Compute Gelman prior with complete data
        # Note: gelman.prior returns only V (variance-covariance matrix), not mu
        # For categorical with J traits, we need to expand this
        complete_data$trait <- factor(rep(1:J, length.out = nrow(complete_data)))
        fixform_cat <- as.formula(paste0(
          as.character(fixform)[2], 
          "~ trait:(", 
          as.character(fixform)[3], 
          ") - 1"
        ))
        gelman_V <- MCMCglmm::gelman.prior(fixform_cat, data = complete_data, scale = 1 + pi^2 / 3)
        
        # Expand to J traits: block diagonal structure
        # Each trait gets the same prior structure
        gelman_size <- nrow(gelman_V)  # Number of coefficients per trait
        prior_B <- list(
          mu = rep(0, gelman_size),  # mu is always 0 (centered prior)
          V = gelman_V  # Block diagonal for J traits
        )
      }
      prior$B <- prior_B
    }
  }

  if (type == "threshold" || type == "ordinal") {
    # stopifnot(!is.null(n_levels)) # this check not needed for ordinal prior

    if (is.null(nu)) {
      nu <- 0.002
    }

    # Fix residual variance R at 1
    # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm
    # priors from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q2/018194.html
    # The residual variance has to be fixed, because the latent variable is without scale

    # J <- n_levels # number of categories not needed for ordinal prior

    prior_G <- .list_of_G(n_rand, nu = nu, par_expand)
    prior <- list(
      R = list(V = 1, fix = 1),
      G = prior_G
    )
  }

  return(prior)
}

#' @title .fit_predict_ovr
#' @description One-vs-rest binary fitting for categorical variables.
#'   Fits J binary threshold MCMCglmm models (one per category level) and
#'   aggregates predicted probabilities into a normalised probability matrix.
#'   Binary threshold models mix far better than multinomial probit and are
#'   the recommended path when ovr_categorical = TRUE.
#' @param data_i Prepared data frame (output of .data_prep)
#' @param response_var Character name of the categorical response variable
#' @param fixformula Fixed-effects formula
#' @param phylo Phylogenetic tree (phylo object)
#' @param ran_phylo_form Random-effects formula
#' @param nitt MCMC iterations
#' @param thin MCMC thinning interval
#' @param burnin MCMC burn-in
#' @param n_rand_eff Number of random effects (1 or 2)
#' @param cluster_col Name of the phylogenetic grouping column
#' @param dat_prep Output of .data_prep (used to extract factor levels)
#' @param sample Logical; if TRUE sample from probability distribution, else argmax
#' @return Named list with full_prediction (n_obs x J probability matrix) and pred_values
.fit_predict_ovr <- function(data_i, response_var, fixformula, phylo,
                              ran_phylo_form, nitt, thin, burnin,
                              n_rand_eff, cluster_col, dat_prep, sample = FALSE) {

  lv          <- dat_prep[[1]][[response_var]]
  level_names <- if (is.factor(lv)) levels(lv) else sort(unique(na.omit(as.character(lv))))
  n_obs       <- nrow(data_i)

  # Initialise probability matrix with equal priors (overwritten per level)
  prob_mat <- matrix(0.5, nrow = n_obs, ncol = length(level_names),
                     dimnames = list(NULL, level_names))

  for (j in seq_along(level_names)) {
    lv_j     <- level_names[j]
    data_bin <- data_i
    # Binarize: level_j -> "yes", all others -> "no", NA stays NA
    data_bin[[response_var]] <- factor(
      ifelse(is.na(data_i[[response_var]]), NA_character_,
             ifelse(as.character(data_i[[response_var]]) == lv_j, "yes", "no")),
      levels = c("no", "yes"))

    prior_j <- .make_prior(n_rand = n_rand_eff, n_levels = 2L, type = "threshold",
                            fixform = fixformula, data = data_bin, gelman = 0)

    model_j <- .model_fit(data = data_bin, tree = phylo,
                          fixformula = fixformula, randformula = ran_phylo_form,
                          type = "threshold", prior = prior_j,
                          nitt = nitt, thin = thin, burnin = burnin)

    pred_j <- .pred_threshold_forward(model_j, fixformula, data_bin,
                                       cluster_col = cluster_col,
                                       level_names = c("no", "yes"))
    prob_mat[, j] <- pred_j[, "yes"]
  }

  # Normalize rows to sum to 1
  rs            <- rowSums(prob_mat)
  rs[rs == 0]   <- 1L
  prob_mat      <- prob_mat / rs

  list(full_prediction = prob_mat,
       pred_values     = .impute_levels(prob_mat, level_names, sample = sample))
}


#' @title .predict_bace
#' @description Function creates a predcition from MCMCglmm model
#' @param model A MCMCglmm model object
#' @param dat_prep A list containing the prepared data frame and attributes for continuous variables.
#' @param response_var A string specifying the name of the response variable.
#' @param type A string that specifies the type of model to fit.
#' @param sample A logical indicating whether to sample from the distribution for categorical/ordinal variables.
#' @param formula Optional: the fixed-effects formula used to fit the model. When supplied together
#'   with data_full, enables forward prediction for categorical/threshold models that covers ALL
#'   rows (including those with NA response) rather than only the complete cases stored in model$Liab.
#' @param data_full Optional: the full data frame passed to MCMCglmm (all rows, including NA response).
#'   Required together with formula to activate forward prediction for categorical/threshold types.
#' @param cluster_col Name of the random-effect grouping column (default "animal").
#' @param ... Additional arguments (not used).
#' @return A vector of predicted values from the MCMCglmm model.
#' @export
.predict_bace <- function(model, dat_prep, response_var, type = NULL, sample = FALSE,
                           formula = NULL, data_full = NULL, cluster_col = "animal", ...) {

  # Initialise so return() is always valid even if a type branch is skipped
  pred_prob   <- NULL
  pred_values <- NULL

				if(type == "gaussian"){
					# z-transformed so need to back-transform
					mean_val <- dat_prep[[2]][[response_var]]$mean
					sd_val   <- dat_prep[[2]][[response_var]]$sd

					if (sample) {
					  # Draw one random MCMC iteration from the posterior
					  X   <- as.matrix(model$X)
					  Sol <- as.matrix(model$Sol)
					  W   <- if (!is.null(model$Z)) cbind(X, as.matrix(model$Z)) else X
					  common <- intersect(colnames(W), colnames(Sol))
					  eta <- Sol[, common, drop = FALSE] %*% t(W[, common, drop = FALSE])
					  i_samp <- sample.int(nrow(eta), 1L)
					  pred_values <- as.numeric(eta[i_samp, ]) * sd_val + mean_val
					} else {
					  # Predict from model and back-transform
					  pred_prob   <- .pred_cont(model) * sd_val + mean_val # Full prediction
					  pred_values <- pred_prob[, 1]                        # Extract posterior mean
					}
			     }

				 if(type == "poisson"){
					if (sample) {
					  # Draw one random MCMC iteration from the liability scale
					  liab <- as.matrix(model$Liab)
					  i_samp <- sample.int(nrow(liab), 1L)
					  pred_values <- round(exp(as.numeric(liab[i_samp, ])), digits = 0)
					} else {
					  # Predict from model and round to nearest integer to retain count data
					  pred_prob   <- .pred_count(model)                    # Full prediction
					  pred_values <- round(pred_prob[, 1], digits = 0)    # Extract posterior mean and round
					}
				 }

				 if(type == "threshold" || type == "ordinal"){
					# Identify number of categories and their levels from the data
			      	     lv  <- dat_prep[[1]][[response_var]]
				  levels_var <- if (is.factor(lv)) levels(lv) else sort(unique(na.omit(as.character(lv))))

				   # Use forward prediction (all rows) when formula + full data supplied;
				   # otherwise fall back to Liab-based prediction (complete cases only).
				   if (!is.null(formula) && !is.null(data_full)) {
				     pred_prob <- .pred_threshold_forward(model, formula, data_full,
				                                          cluster_col = cluster_col,
				                                          level_names = levels_var)
				   } else {
				     pred_prob <- .pred_threshold(model, level_names = levels_var)
				   }

				   # For each observation, sample from the categorical distribution based on the predicted probabilities. TO DO: Note we could also just take the max probability for baseline level
             pred_values <- .impute_levels(pred_prob, levels_var, sample = sample)
          }


				 if(type == "categorical"){
					# Identify number of categories and their levels from the data
			      	     lv  <- dat_prep[[1]][[response_var]]
				  levels_var <- if (is.factor(lv)) levels(lv) else sort(unique(na.omit(as.character(lv))))

					# Use forward prediction (all rows) when formula + full data supplied;
					# otherwise fall back to Liab-based prediction (complete cases only).
					if (!is.null(formula) && !is.null(data_full)) {
					  pred_prob <- .pred_cat_forward(model, formula, data_full,
					                                 cluster_col = cluster_col,
					                                 baseline_name = levels_var[1])
					} else {
					  pred_prob <- .pred_cat(model, baseline_name = levels_var[1])
					}

					# For each observation, sample from the categorical distribution based on the predicted probabilities
				   pred_values <- .impute_levels(pred_prob, levels_var, sample = sample)
				 }

	return(list(full_prediction = pred_prob,
                  pred_values = pred_values))
}

#' @title .impute_levels
#' @description Function samples predicted levels for categorical or ordinal variables based on predicted probabilities
#' @param pred_prob A data frame of predicted probabilities for each category
#' @param levels_var A character vector specifying the names of the levels/categories
#' @param sample A logical indicating whether to sample from the distribution or take the maximum probability
#' @return A vector of predicted levels for each observation
.impute_levels <- function(pred_prob, levels_var, sample = FALSE) {
  if (sample) {
    pred_values <- apply(pred_prob, 1, function(probs) {
      sample(levels_var, size = 1, prob = probs)
    })
  } else {
    pred_values <- levels_var[apply(pred_prob, 1, which.max)]
  }
  return(as.character(pred_values))
}

#' @title .pred_cat
#' @description Function calculates predicted probabilities for each category from a categorical MCMCglmm model
#' @param model A MCMCglmm model object
#' @param baseline_name A string specifying the name of the baseline category
#' @return A data frame of predicted probabilities for each category
#' @export
.pred_cat <- function(model, baseline_name = "Baseline") {

  # 1. Basic Dimensions
     n_obs <- model$Residual$nrl
  n_traits <- (ncol(model$Sol) - model$Fixed$nll) / (model$Fixed$nfl / (ncol(model$Sol) / (ncol(model$Liab)/n_obs)))

  # Simpler way to get n_traits for categorical:
  n_traits <- ncol(model$Liab) / n_obs

  # 2. Calculate the c2 Scaling Factor
  # c = 16 * sqrt(3) / (15 * pi)
        c2 <- (16 * sqrt(3) / (15 * pi))^2

  # Get Total Variance (G + R) for each sample in the chain
  # We sum the Phylo and Units variances.
  # Note: MCMCglmm categorical models usually have fixed residual variance (e.g. 1)

  # since residual variance is fixed - we can use directly the IJ matrix
  IJ <- IJ <- 1 / (n_traits+1) * (diag(n_traits) + matrix(1, nrow = n_traits, ncol = n_traits))
  scaling_factor <- sqrt(1 + c2 * diag(IJ))

  # 3. Extract and Scale Liabilities
  exp_liab_list <- list()

  # Initialize exp_sum with 1 (which is exp(0) for the baseline)
  exp_sum <- 1

  for (i in 1:n_traits) {
    cols <- ((i - 1) * n_obs + 1):(i * n_obs)

    # Apply c2 scaling to the liabilities before exponentiating
    # We divide each sample's liabilities by that sample's scaling factor
    scaled_liab <- model$Liab[, cols] / scaling_factor[i]

    exp_liab_list[[i]] <- exp(scaled_liab)
    exp_sum <- exp_sum + exp_liab_list[[i]]
  }

  # 4. Calculate Mean Probabilities (%)
  prop_results <- list()
  for (i in 1:n_traits) {
    prop_results[[i]] <- colMeans(exp_liab_list[[i]] / exp_sum)
  }

  # Calculate Baseline level %
  prop_results[[n_traits + 1]] <- colMeans(1 / exp_sum)

  # 5. Extract Names Generically
  # Look at Fixed effects to get trait/variable names
  raw_names <- colnames(model$Sol)[1:n_traits]

  # Removes "trait", the variable name, and the following dot
  clean_names <- gsub("^trait.*?\\.", "", raw_names)

  # 6. Final Data Frame Assembly
  df_final <- as.data.frame(do.call(cbind, prop_results))
  colnames(df_final) <- c(clean_names, baseline_name)

  # Rearrange columns alphabetically (Aquatic, Insessorial, Terrestrial)
  df_ordered <- df_final[, order(colnames(df_final))]
  rownames(df_ordered) <- paste0("Obs_", 1:n_obs)

  return(df_ordered)
}


#' @title .pred_cat_forward
#' @description Forward prediction of categorical probabilities for ALL rows (including
#'   those with NA response) using posterior mean fixed + random effects from model$Sol.
#'   Called instead of .pred_cat() when the response variable has missing values, because
#'   MCMCglmm drops NA-response rows from model$Liab/model$X, making .pred_cat() return
#'   fewer rows than the full dataset.
#' @param model MCMCglmm model object (fit with pr=TRUE, saveX=TRUE)
#' @param formula The fixed-effects formula used to fit the model
#' @param data_i Full data frame including rows with NA response (all predictors complete)
#' @param cluster_col Name of the random-effect grouping column (e.g. "animal")
#' @param baseline_name Name for the baseline (reference) category
#' @return Data frame with nrow(data_i) rows of predicted probabilities per category
#' @export
.pred_cat_forward <- function(model, formula, data_i, cluster_col = "animal",
                               baseline_name = "Baseline") {

  beta_mean <- colMeans(as.matrix(model$Sol))
  sol_names <- names(beta_mean)

  # MCMCglmm naming for categorical:
  #   intercept for trait k  : "trait{Response}.{Level}"           (no suffix)
  #   slope for predictor j  : "trait{Response}.{Level}:{j}"       (colon separator)
  #   species BLUP            : "trait{Response}.{Level}.{species}" (dot separator)
  #
  # Extract trait identifiers from slope columns (those containing ":").
  # This excludes BLUP entries (no colon) which also start with "trait",
  # e.g. "traitTrophic.Level.Herbivore.Accipiter_gentilis".
  trait_slope_cols <- sol_names[grep("^trait.*:", sol_names)]
  trait_names      <- unique(sub(":.*$", "", trait_slope_cols))
  n_traits         <- length(trait_names)
  n_obs          <- nrow(data_i)

  # c2 scaling (same as .pred_cat)
  c2             <- (16 * sqrt(3) / (15 * pi))^2
  IJ             <- 1 / (n_traits + 1) * (diag(n_traits) + matrix(1, n_traits, n_traits))
  scaling_factor <- sqrt(1 + c2 * diag(IJ))

  # Build fixed-effects design matrix for ALL rows
  rhs   <- formula[-2]   # drop LHS response
  X_all <- tryCatch(model.matrix(rhs, data_i), error = function(e) NULL)

  if (is.null(X_all)) {
    # Fallback: return uniform probabilities
    k  <- n_traits + 1
    df <- as.data.frame(matrix(1 / k, nrow = n_obs, ncol = k))
    # Extract level names from trait_names: last dot-separated component
    non_base <- sub(".*\\.", "", trait_names)
    colnames(df) <- c(baseline_name, non_base)
    return(df)
  }

  x_cols  <- colnames(X_all)
  exp_list <- list()
  exp_sum  <- 1

  for (ki in seq_len(n_traits)) {
    tk <- trait_names[ki]

    # MCMCglmm Sol column names:
    #   intercept (x_col == "(Intercept)") -> tk
    #   other predictors                   -> paste0(tk, ":", x_col)
    sol_fixed <- ifelse(x_cols == "(Intercept)", tk, paste0(tk, ":", x_cols))
    valid     <- sol_fixed %in% sol_names
    beta_k    <- beta_mean[sol_fixed[valid]]

    # Linear predictor from fixed effects
    eta_k <- as.numeric(X_all[, x_cols[valid], drop = FALSE] %*% beta_k)

    # Species BLUPs: "trait{Response}.{Level}.{species_id}" (dot separator)
    if (cluster_col %in% colnames(data_i)) {
      species_vals <- as.character(data_i[[cluster_col]])
      blup_cols    <- paste0(tk, ".", species_vals)
      blups        <- vapply(blup_cols, function(bc)
                               if (bc %in% sol_names) beta_mean[[bc]] else 0,
                             numeric(1))
      eta_k <- eta_k + blups
    }

    exp_list[[ki]] <- exp(eta_k / scaling_factor[ki])
    exp_sum        <- exp_sum + exp_list[[ki]]
  }

  # Softmax probabilities
  prop_results <- vector("list", n_traits + 1)
  for (ki in seq_len(n_traits)) {
    prop_results[[ki]] <- exp_list[[ki]] / exp_sum
  }
  prop_results[[n_traits + 1]] <- 1 / exp_sum

  # Extract category (level) names: last dot-separated component of trait identifier
  # e.g. "traitTrophic.Level.Herbivore" -> "Herbivore"
  non_baseline_names <- sub(".*\\.", "", trait_names)
  all_level_names    <- c(baseline_name, non_baseline_names)

  # Output columns: [baseline, trait1, ..., trait_{n_traits}] = levels_var order
  df_final  <- as.data.frame(do.call(cbind, prop_results[c(n_traits + 1L, seq_len(n_traits))]))
  colnames(df_final)  <- all_level_names
  rownames(df_final)  <- paste0("Obs_", seq_len(n_obs))
  return(df_final)
}


#' @title .pred_threshold
#' @description Function calculates predicted probabilities for each category from a threshold MCMCglmm model
#' @param model A MCMCglmm model object
#' @param level_names A character vector specifying the names of the levels/categories
#' @return A data frame of predicted probabilities for each category
#' @export
.pred_threshold <- function(model, level_names = NULL) {

  n_obs <- model$Residual$nrl

  if (!is.null(model$CP)) {
    cp_samples <- cbind(0, model$CP)
    n_levels   <- ncol(cp_samples) + 1
  } else {
    cp_samples <- matrix(0, nrow = nrow(model$Sol), ncol = 1)
    n_levels   <- 2
  }

  liab <- as.matrix(model$Liab)

  all_probs <- vector("list", n_levels)
  for (j in seq_len(n_levels)) {
    t_low  <- if (j == 1)        -Inf else cp_samples[, j - 1]
    t_high <- if (j == n_levels)  Inf else cp_samples[, j]
    all_probs[[j]] <- colMeans(pnorm(t_high, mean = liab, sd = 1) -
                                 pnorm(t_low,  mean = liab, sd = 1))
  }

  df_final <- as.data.frame(do.call(cbind, all_probs))
  colnames(df_final) <- if (is.null(level_names)) paste0("Level_", seq_len(n_levels)) else level_names
  return(df_final)
}


#' @title .pred_threshold_forward
#' @description Forward prediction of ordinal/threshold probabilities for ALL rows (including
#'   those with NA response) using posterior mean fixed + random effects from model$Sol.
#'   Called instead of .pred_threshold() when the response variable has missing values.
#' @param model MCMCglmm model object (fit with pr=TRUE)
#' @param formula The fixed-effects formula used to fit the model
#' @param data_i Full data frame including rows with NA response (all predictors complete)
#' @param cluster_col Name of the random-effect grouping column (e.g. "animal")
#' @param level_names Character vector of ordered level names
#' @return Data frame with nrow(data_i) rows of predicted probabilities per category
#' @export
.pred_threshold_forward <- function(model, formula, data_i,
                                     cluster_col = "animal",
                                     level_names = NULL) {

  beta_mean <- colMeans(as.matrix(model$Sol))
  sol_names <- names(beta_mean)
  n_obs     <- nrow(data_i)

  # Build fixed-effects design matrix for ALL rows
  rhs   <- formula[-2]
  X_all <- tryCatch(model.matrix(rhs, data_i), error = function(e) NULL)

  if (is.null(X_all)) {
    n_lv <- if (!is.null(model$CP)) ncol(as.matrix(model$CP)) + 2 else 2
    df <- as.data.frame(matrix(1 / n_lv, nrow = n_obs, ncol = n_lv))
    colnames(df) <- if (!is.null(level_names) && length(level_names) == n_lv)
                      level_names else paste0("Level_", seq_len(n_lv))
    return(df)
  }

  # Fixed effect columns: X_all column names that exist in Sol (no cluster prefix)
  x_cols <- colnames(X_all)
  valid  <- x_cols %in% sol_names
  beta_k <- beta_mean[x_cols[valid]]

  # Linear predictor from fixed effects
  eta <- as.numeric(X_all[, x_cols[valid], drop = FALSE] %*% beta_k)

  # Add species BLUPs (posterior mean)
  if (cluster_col %in% colnames(data_i)) {
    species_vals <- as.character(data_i[[cluster_col]])
    blup_cols    <- paste0(cluster_col, ".", species_vals)
    blups        <- vapply(blup_cols, function(bc)
                             if (bc %in% sol_names) beta_mean[[bc]] else 0,
                           numeric(1))
    eta <- eta + blups
  }

  # Cut-points: posterior mean
  # MCMCglmm fixes the first threshold at 0; additional ones are in model$CP.
  # Prepend 0 so cp_all = c(0, CP2_mean, CP3_mean, ...).
  if (!is.null(model$CP) && length(model$CP) > 0) {
    cp_all   <- c(0, colMeans(as.matrix(model$CP)))
    n_levels <- length(cp_all) + 1
  } else {
    cp_all   <- 0          # binary: single threshold at 0
    n_levels <- 2
  }

  # P(category j | eta_i) = Phi(CP_j - eta_i) - Phi(CP_{j-1} - eta_i)
  all_probs <- vector("list", n_levels)
  for (j in seq_len(n_levels)) {
    t_low  <- if (j == 1)        -Inf else cp_all[j - 1]
    t_high <- if (j == n_levels)  Inf else cp_all[j]
    all_probs[[j]] <- pnorm(t_high - eta) - pnorm(t_low - eta)
  }

  df_final <- as.data.frame(do.call(cbind, all_probs))
  if (!is.null(level_names) && length(level_names) == n_levels) {
    colnames(df_final) <- level_names
  } else {
    colnames(df_final) <- paste0("Level_", seq_len(n_levels))
  }
  return(df_final)
}


#' @title .pred_cont
#' @description Posterior mean, posterior SD, and 95% credible interval (2.5%, 97.5%)
#'              of fitted values for a Gaussian (identity-link) MCMCglmm model
#' @param model A MCMCglmm model object
#' @return A data frame with columns: post_mean, post_sd, ci_lower, ci_upper
#'         (one row per observation used in the fit)
#' @export
.pred_cont <- function(model) {

  # Need X and Sol (and Z if you want conditional fitted values)
  if (is.null(model$X))   stop("model$X is missing: fit with saveX=TRUE.")
  if (is.null(model$Sol)) stop("model$Sol is missing.")

  X   <- as.matrix(model$X)
  Sol <- as.matrix(model$Sol)

  # Use Z if present (pr=TRUE typically required to have RE columns in Sol)
  if (!is.null(model$Z)) {
    W <- cbind(X, as.matrix(model$Z))
  } else {
    W <- X
  }

  # Align columns (most robust)
  if (!is.null(colnames(W)) && !is.null(colnames(Sol))) {
    common <- intersect(colnames(W), colnames(Sol))
    if (length(common) == 0) {
      stop("No matching coefficient names between design matrix (X/Z) and Sol.")
    }
    W   <- W[, common, drop = FALSE]
    Sol <- Sol[, common, drop = FALSE]
  } else {
    # Fallback: assume ordering matches for first min columns
    p <- min(ncol(W), ncol(Sol))
    W   <- W[, seq_len(p), drop = FALSE]
    Sol <- Sol[, seq_len(p), drop = FALSE]
  }

  # eta draws: [n_iter x n_obs]
  eta <- Sol %*% t(W)

  ci <- t(apply(eta, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE))

  out <- data.frame(
    post_mean = as.numeric(colMeans(eta, na.rm = TRUE)),
    post_sd   = as.numeric(apply(eta, 2, stats::sd, na.rm = TRUE)),
    ci_lower  = as.numeric(ci[, 1]),
    ci_upper  = as.numeric(ci[, 2])
  )

  rownames(out) <- paste0("Obs_", seq_len(nrow(out)))
  return(out)
}


#' @title .pred_count
#' @description Posterior mean, posterior SD, and 95% credible interval (2.5%, 97.5%)
#'              of fitted values for a Poisson (log-link) MCMCglmm model
#' @param model A MCMCglmm model object
#' @return A data frame with columns: post_mean, post_sd, ci_lower, ci_upper
#'         (one row per observation)
#' @export
.pred_count <- function(model) {

  if (is.null(model$Liab)) stop("model$Liab is missing.")
  liab <- as.matrix(model$Liab)

  mu <- exp(liab)

  ci <- t(apply(mu, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE))

  out <- data.frame(
    post_mean = as.numeric(colMeans(mu, na.rm = TRUE)),
    post_sd   = as.numeric(apply(mu, 2, stats::sd, na.rm = TRUE)),
    ci_lower  = as.numeric(ci[, 1]),
    ci_upper  = as.numeric(ci[, 2])
  )

  rownames(out) <- paste0("Obs_", seq_len(nrow(out)))
  return(out)
}


#' @title .check_mcmc_diagnostics
#' @description Function checks MCMC convergence and mixing for all parameters in MCMCglmm models
#' @param bace_output Output from bace_imp function containing models_last_run
#' @return A list containing diagnostic statistics and a summary data frame
#' @details This function performs the following diagnostics:
#' \itemize{
#'   \item Effective Sample Size (ESS) for all parameters
#'   \item Geweke convergence diagnostic
#'   \item Autocorrelation at lag 1
#'   \item Summary statistics for fixed and random effects
#' }
#' @export
.check_mcmc_diagnostics <- function(bace_output) {
  
  if (!inherits(bace_output, "bace")) {
    stop("Input must be a 'bace' object from bace_imp()")
  }
  
  if (is.null(bace_output$models_last_run)) {
    stop("No models found in bace_output. Make sure models were saved during imputation.")
  }
  
  models <- bace_output$models_last_run
  var_types <- bace_output$types
  diagnostics <- list()
  
  for (var_name in names(models)) {
    model <- models[[var_name]]
    
    # Get variable type for this model
    var_type <- var_types[[var_name]]
    
    # Extract ONLY fixed effects and variance components (exclude BLUPs)
    # Fixed effects: First nfl columns of Sol (number of fixed effects levels)
    n_fixed <- model$Fixed$nfl
    fixed_effects <- model$Sol[, 1:n_fixed, drop = FALSE]
    
    # Random effect variance components
    random_variances <- model$VCV
    
    # Combine fixed effects and variance components only
    all_params <- cbind(fixed_effects, random_variances)
    param_names <- colnames(all_params)
    
    # Initialize results for this model
    var_diagnostics <- data.frame(
      parameter = param_names,
      mean = numeric(length(param_names)),
      sd = numeric(length(param_names)),
      ess = numeric(length(param_names)),
      geweke_z = numeric(length(param_names)),
      geweke_pval = numeric(length(param_names)),
      autocorr_lag1 = numeric(length(param_names)),
      convergence = character(length(param_names)),
      stringsAsFactors = FALSE
    )
    
    # Calculate diagnostics for each parameter
    for (i in seq_along(param_names)) {
      param <- all_params[, i]
      param_name <- param_names[i]
      
      # Check if this is a "units" parameter for threshold or categorical models
      # For these models, units variance is fixed and doesn't need convergence testing
      is_fixed_units <- (var_type %in% c("threshold", "categorical")) && 
                        grepl("units", param_name, ignore.case = TRUE)
      
      # Basic statistics
      var_diagnostics$mean[i] <- mean(param)
      var_diagnostics$sd[i] <- sd(param)
      
      if (is_fixed_units) {
        # For fixed units parameters, skip convergence tests
        var_diagnostics$ess[i] <- NA
        var_diagnostics$geweke_z[i] <- NA
        var_diagnostics$geweke_pval[i] <- NA
        var_diagnostics$autocorr_lag1[i] <- NA
        var_diagnostics$convergence[i] <- "Fixed"
      } else {
        # Effective Sample Size
        var_diagnostics$ess[i] <- coda::effectiveSize(param)
        
        # Geweke diagnostic (tests for equality of means in first 10% and last 50% of chain)
        geweke_result <- coda::geweke.diag(param)
        var_diagnostics$geweke_z[i] <- geweke_result$z
        var_diagnostics$geweke_pval[i] <- 2 * pnorm(-abs(geweke_result$z))
        
        # Autocorrelation at lag 1
        acf_result <- stats::acf(param, lag.max = 1, plot = FALSE)
        var_diagnostics$autocorr_lag1[i] <- acf_result$acf[2]
        
        # Convergence assessment
        # ESS > 200 is generally considered adequate
        # Geweke p-value > 0.05 suggests convergence
        # Low autocorrelation (< 0.4) suggests good mixing
        
        ess_ok <- var_diagnostics$ess[i] > 200
        geweke_ok <- var_diagnostics$geweke_pval[i] > 0.05
        autocorr_ok <- abs(var_diagnostics$autocorr_lag1[i]) < 0.4
        
        if (ess_ok && geweke_ok && autocorr_ok) {
          var_diagnostics$convergence[i] <- "Good"
        } else if (ess_ok && geweke_ok) {
          var_diagnostics$convergence[i] <- "Acceptable"
        } else {
          var_diagnostics$convergence[i] <- "Check"
        }
      }
    }
    
    # Identify fixed vs random variance components
    # Fixed effects are the first n_fixed parameters
    var_diagnostics$effect_type <- ifelse(
      seq_along(param_names) <= n_fixed,
      "Fixed",
      "Variance"
    )
    
    diagnostics[[var_name]] <- var_diagnostics
  }
  
  # Create overall summary
  all_diagnostics <- do.call(rbind, lapply(names(diagnostics), function(var) {
    df <- diagnostics[[var]]
    df$response_variable <- var
    df
  }))
  
  # Summary statistics (excluding Fixed units from convergence calculations)
  not_fixed <- all_diagnostics$convergence != "Fixed"
  summary_stats <- list(
    n_models = length(models),
    n_parameters_total = nrow(all_diagnostics),
    n_fixed_effects = sum(all_diagnostics$effect_type == "Fixed"),
    n_variance_components = sum(all_diagnostics$effect_type == "Variance"),
    n_fixed_units = sum(all_diagnostics$convergence == "Fixed"),
    convergence_summary = table(all_diagnostics$convergence),
    mean_ess_fixed = mean(all_diagnostics$ess[all_diagnostics$effect_type == "Fixed" & not_fixed], na.rm = TRUE),
    mean_ess_variance = mean(all_diagnostics$ess[all_diagnostics$effect_type == "Variance" & not_fixed], na.rm = TRUE),
    prop_good_convergence = mean(all_diagnostics$convergence[not_fixed] == "Good"),
    prop_poor_convergence = mean(all_diagnostics$convergence[not_fixed] == "Check")
  )
  
  result <- list(
    diagnostics_by_model = diagnostics,
    all_diagnostics = all_diagnostics,
    summary = summary_stats
  )
  
  class(result) <- c("bace_diagnostics", "list")
  return(result)
}


#' @title print.bace_diagnostics
#' @description Print method for MCMC diagnostics from bace_imp models
#' @param x A bace_diagnostics object from .check_mcmc_diagnostics()
#' @param ... Additional arguments (not used)
#' @export
print.bace_diagnostics <- function(x, ...) {
  
  cat("\n========================================\n")
  cat("BACE MCMC Diagnostics Summary\n")
  cat("========================================\n\n")
  
  # Overall summary
  cat("Number of models:", x$summary$n_models, "\n")
  cat("Total parameters:", x$summary$n_parameters_total, "\n")
  cat("  Fixed effects:", x$summary$n_fixed_effects, "\n")
  cat("  Variance components:", x$summary$n_variance_components, "\n")
  if (x$summary$n_fixed_units > 0) {
    cat("  Fixed units (threshold/categorical):", x$summary$n_fixed_units, "\n")
  }
  cat("\n")
  
  # Effective Sample Size
  cat("Mean Effective Sample Size:\n")
  cat(sprintf("  Fixed effects:      %.1f\n", x$summary$mean_ess_fixed))
  cat(sprintf("  Variance components: %.1f\n", x$summary$mean_ess_variance))
  cat("\n")
  
  # Convergence summary
  cat("Convergence Assessment:\n")
  conv_table <- x$summary$convergence_summary
  for (status in names(conv_table)) {
    prop <- conv_table[status] / sum(conv_table) * 100
    cat(sprintf("  %-12s: %3d parameters (%.1f%%)\n", 
                status, conv_table[status], prop))
  }
  cat("\n")
  
  # Overall assessment
  if (x$summary$prop_good_convergence >= 0.8) {
    cat("Overall: GOOD - Most parameters show good convergence\n")
  } else if (x$summary$prop_poor_convergence <= 0.2) {
    cat("Overall: ACCEPTABLE - Most parameters converged adequately\n")
  } else {
    cat("Overall: CHECK - Many parameters show convergence issues\n")
    cat("  Consider increasing nitt, adjusting priors, or checking model specification\n")
  }
  cat("\n")
  
  # Model-by-model summary
  cat("Convergence by Response Variable:\n")
  cat("----------------------------------------\n")
  
  for (var_name in names(x$diagnostics_by_model)) {
    diag <- x$diagnostics_by_model[[var_name]]
    
    # Exclude Fixed units from convergence counts
    not_fixed <- diag$convergence != "Fixed"
    n_good <- sum(diag$convergence == "Good")
    n_poor <- sum(diag$convergence == "Check")
    n_fixed_units <- sum(diag$convergence == "Fixed")
    n_total <- sum(not_fixed)  # Total excluding fixed units
    n_total_all <- nrow(diag)  # Including fixed units
    
    # Calculate status based on non-fixed parameters
    if (n_total > 0) {
      status_symbol <- if (n_good / n_total >= 0.8) "[OK]" else if (n_poor / n_total > 0.3) "[CH]" else "[~]"
    } else {
      status_symbol <- "[--]"  # No testable parameters
    }
    
    if (n_fixed_units > 0) {
      cat(sprintf("%s %-15s: %2d/%2d parameters good (%d fixed)\n", 
                  status_symbol, var_name, n_good, n_total, n_fixed_units))
    } else {
      cat(sprintf("%s %-15s: %2d/%2d parameters good\n", 
                  status_symbol, var_name, n_good, n_total))
    }
    
    # Highlight poorly converged parameters
    poor_params <- diag$parameter[diag$convergence == "Check"]
    if (length(poor_params) > 0) {
      cat(sprintf("    Check convergence: %s\n", 
                  paste(poor_params, collapse=", ")))
    }
  }
  
  cat("========================================\n\n")
  
  invisible(x)
}


#' @title summary.bace_diagnostics
#' @description Detailed summary method for MCMC diagnostics
#' @param object A bace_diagnostics object from .check_mcmc_diagnostics()
#' @param ... Additional arguments (not used)
#' @export
summary.bace_diagnostics <- function(object, ...) {
  
  cat("\n========================================\n")
  cat("Detailed MCMC Diagnostics\n")
  cat("========================================\n\n")
  
  for (var_name in names(object$diagnostics_by_model)) {
    cat("\nResponse Variable:", var_name, "\n")
    cat("----------------------------------------\n")
    
    diag <- object$diagnostics_by_model[[var_name]]
    
    # Print table
    print(diag[, c("parameter", "effect_type", "mean", "sd", "ess", 
                   "geweke_pval", "autocorr_lag1", "convergence")], 
          row.names = FALSE, digits = 3)
    
    cat("\n")
  }
  
  invisible(object)
}
