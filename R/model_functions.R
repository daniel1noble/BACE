#' @title .model_fit
#' @description Function takes incomplete data and a formula string and fits a model using MCMCglmm
#' @param data A dataframe containing missing data.
#' @param tree A phylogenetic tree of class 'phylo' from the ape package.
#' @param fixformula A string that specifies the fixed effect structure in the model.
#' @param randformula A string that specifies the random effect structure in the model.
#' @param type A string that specifies the type of model to fit. Options are "normal", "binary", "count", "categorical", "ordinal" and the appropriate family will be used in MCMCglmm.
#' @param nitt An integer specifying the number of iterations to run the MCMC algorithm. Default is 6000
#' @param thin An integer specifying the thinning rate for the MCMC algorithm. Default is 5.
#' @param burnin An integer specifying the number of iterations to discard as burnin. Default is 1000.
#' @param prior A list specifying the prior distributions for the MCMCglmm model.
#' @return A list of draws from the posterior distribution of the model parameters.
#' @export

.model_fit <- function(data, tree, fixformula, randformula, type, prior, nitt = 6000, thin = 5, burnin = 1000) {

	# Create sparse matrix of phylogeny
		   A  <- MCMCglmm::inverseA(tree, nodes = "TIPS")$Ainv
	
	# Name of the column in the data corresponding to the phylogeny
		name  <- all.vars(randformula)
		
	# Fit the model using MCMCglmm
  if(type != "categorical"){
  	model <- MCMCglmm::MCMCglmm(fixed = fixformula,
                               random = randformula,
                                 data = data,
                               family = type,
							               ginverse = setNames(list(A), name),
                              verbose = FALSE, pr = TRUE, pl = TRUE,
                                saveX = TRUE, saveZ = TRUE,
                                 nitt = nitt,
                                 thin = thin,
                               burnin = burnin, 
							    prior = prior, singular.ok=TRUE)	
        } else {

      # Categorical model needs special treatment. Append -1 to the right size of ~ formula to remove intercept
      fixformula_cat <- as.formula(paste0(as.character(fixformula)[2], "~ -1 + ", as.character(fixformula)[3]))
      ranformula_cat <- as.formula(paste0("~", "idh(trait):",as.character(randformula)[2]))

      model <- MCMCglmm::MCMCglmm(fixed = fixformula_cat,
                                 random = ranformula_cat,
                                 data = data,
                               family = type,
                               rcov = ~us(trait):units,
							               ginverse = setNames(list(A), name),
                              verbose = FALSE, pr = TRUE, pl = TRUE,
                                saveX = TRUE, saveZ = TRUE,
                                 nitt = nitt,
                                 thin = thin,
                               burnin = burnin, 
							    prior = prior, singular.ok=TRUE)	     
                  }								
  	return(model)
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
      prior_G[[paste0("G", i)]] <- list(
        V = diag(diag), nu = nu,
        alpha.mu = rep(0, diag), alpha.V = 1e4 * diag(diag)
      )
    }
  } else {
    for (i in 1:n_rand) {
      prior_G[[paste0("G", i)]] <- list(V = diag(diag), nu = nu)
    }
  }
  
  return(prior_G)
}

#' @title .make_prior
#' @description Function creates the prior for the MCMCglmm model
#' @param n_rand An integer specifying the number of random effects in the model.
#' @param type A string that specifies the type of model to fit.
#' @param nu A numeric specifying the nu parameter for the prior.
#' @param n_levels An integer specifying the number of levels for categorical or ordinal data.
#' @param par_expand A logical indicating whether to use parameter expansion.
#' @return A list of priors for the MCMCglmm model.
#' @export

.make_prior <- function(n_rand, type, nu = NULL, n_levels = NULL, par_expand = FALSE) {
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

  if (type  == "categorical") {

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

#' @title .predict_bace
#' @description Function creates a predcition from MCMCglmm model
#' @param model A MCMCglmm model object
#' @param dat_prep A list containing the prepared data frame and attributes for continuous variables.
#' @param response_var A string specifying the name of the response variable.
#' @param type A string that specifies the type of model to fit.
#' @param sample A logical indicating whether to sample from the distribution for categorical/ordinal variables.
#' @param ... Additional arguments (not used).
#' @return A vector of predicted values from the MCMCglmm model.
#' @export
.predict_bace <- function(model, dat_prep, response_var, type = NULL, sample = FALSE, ...) {				

				if(type == "gaussian"){
					# z-transformed so need to back-transform
					mean_val <- dat_prep[[2]][[response_var]]$mean
					sd_val   <- dat_prep[[2]][[response_var]]$sd
					
					# Predict from model and back-transform
				 		pred_values <- predict(model, marginal = NULL, type = "response", posterior = "all") * sd_val + mean_val
			     }

				 if(type == "poisson"){  
					# Predict from model and round to nearest integer to retain count data
						pred_values <- round(predict(model, marginal = NULL, type = "response", posterior = "all"), digits = 0)
				 }

				 if(type == "threshold" || type == "ordinal"){
					# Identify number of categories and their levels from the data
			      	     lv  <- dat_prep[[1]][[response_var]]
				  levels_var <- sort(unique(as.character(lv)))
				  
				   # Predicts probabilities for each category
				 pred_prob <- .pred_threshold(model, level_names = levels_var)

				   # For each observation, sample from the categorical distribution based on the predicted probabilities. TO DO: Note we could also just take the max probability for baseline level
             pred_values <- .impute_levels(pred_prob, levels_var, sample = sample)
          }
				 
				 
				 if(type == "categorical"){
					# Identify number of categories and their levels from the data
			      	     lv  <- dat_prep[[1]][[response_var]]
				  levels_var <- sort(unique(as.character(lv)))
					
					# Predict category probabilities
					pred_prob <- .pred_cat(model, baseline_name = levels_var[1])

					# For each observation, sample from the categorical distribution based on the predicted probabilities
				   pred_values <- .impute_levels(pred_prob, levels_var, sample = sample)
				 }

	return(pred_values)
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
  return(pred_values)
}

#' @title .pred_cat
#' @description Function calculates predicted probabilities for each category from a categorical MCMCglmm model
#' 	@param model A MCMCglmm model object
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
  prob_results <- list()
  for (i in 1:n_traits) {
    prob_results[[i]] <- colMeans(exp_liab_list[[i]] / exp_sum) * 100
  }
  
  # Calculate Baseline level %
  prob_results[[n_traits + 1]] <- colMeans(1 / exp_sum) * 100
  
  # 5. Extract Names Generically
  # Look at Fixed effects to get trait/variable names
  raw_names <- colnames(model$Sol)[1:n_traits]
  
  # Removes "trait", the variable name, and the following dot
  clean_names <- gsub("^trait.*?\\.", "", raw_names)
  
  # 6. Final Data Frame Assembly
  df_final <- as.data.frame(do.call(cbind, prob_results))
  colnames(df_final) <- c(clean_names, baseline_name)
  
  # Rearrange columns alphabetically (Aquatic, Insessorial, Terrestrial)
  df_ordered <- df_final[, order(colnames(df_final))]
  rownames(df_ordered) <- paste0("Obs_", 1:n_obs)
  
  return(round(df_ordered, 2))
}


#' @title .pred_threshold
#' @description Function calculates predicted probabilities for each category from a threshold MCMCglmm model
#' @param model A MCMCglmm model object
#' @param level_names A character vector specifying the names of the levels/categories	
#' @return A data frame of predicted probabilities for each category
#' @export
.pred_threshold <- function(model, level_names = NULL) {
  # 1. Dimensions and Data
           n_obs <- model$Residual$nrl
            liab <- model$Liab
         v_total <- rowSums(model$VCV)
  scaling_factor <- 1
  
  # 2. Thresholds (Cut-points)
  # MCMCglmm fixes the first threshold at 0.
  # Additional cut-points are in model$CP.
  if (!is.null(model$CP)) {
    # Ordinal case: J > 2 levels
    cp_samples <- cbind(0, model$CP) # First CP is 0
      n_levels <- ncol(cp_samples) + 1
  } else {
    # Binary case: 2 levels
    cp_samples <- matrix(0, nrow = nrow(liab), ncol = 1)
      n_levels <- 2
  }
  
  # 3. Calculate Probabilities
  # We iterate through iterations and calculate area under normal curve
  # Prob(category j) = Phi((CP_j - Liab) / scale) - Phi((CP_j-1 - Liab) / scale)
  all_probs <- list()
  for (j in 1:n_levels) {
    # Define upper and lower bounds for the current category
    # Lower bound (T_low)
    if (j == 1) {
      t_low <- -Inf
    } else {
      t_low <- cp_samples[, j - 1]
    }
    # Upper bound (T_high)
    if (j == n_levels) {
      t_high <- Inf
    } else {
      t_high <- cp_samples[, j]
    }
    
	# Calculate probability for this category across all MCMC samples
    # Probability = pnorm(Upper) - pnorm(Lower)
    # We use mapply or row-wise math to handle the scaling per iteration
    # Prob per sample: pnorm(t_high, mean=liab, sd=scaling_factor) - pnorm(t_low, ...)
             p_cat <- pnorm(t_high, mean = liab, sd = scaling_factor) -
                      pnorm(t_low, mean = liab, sd = scaling_factor)
    all_probs[[j]] <- colMeans(p_cat) 
  }
  
  # 4. Final Data Frame
  df_final <- as.data.frame(do.call(cbind, all_probs))
  
  # 5. Naming and Ordering
  if (is.null(level_names)) {
    colnames(df_final) <- paste0("Level_", 1:n_levels)
  } else {
    colnames(df_final) <- level_names
  }
  return(df_final)
}
