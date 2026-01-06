#' @title model_fit
#' @description Function takes incomplete data and a formula string and fits a model using MCMCglmm
#' @param data A dataframe containing missing data.
#' @param tree A phylogenetic tree of class 'phylo' from the ape package.
#' @param fixformula A string that specifies the fixed effect structure in the model.
#' @param randformula A string that specifies the random effect structure in the model.
#' @param type A string that specifies the type of model to fit. Options are "normal", "binary", "count", "categorical", "ordinal" and the appropriate family will be used in MCMCglmm.
#' @param nitt An integer specifying the number of iterations to run the MCMC algorithm. Default is 50000.
#' @param thin An integer specifying the thinning rate for the MCMC algorithm. Default is 10.
#' @param burnin An integer specifying the number of iterations to discard as burnin. Default is 1000.
#' @return A list of draws from the posterior distribution of the model parameters.
#' @export

model_fit <- function(data, tree, fixformula, randformula, type, prior, nitt = 50000, thin = 10, burnin = 1000) {

	# Create sparse matrix of phylogeny
		   A  <- MCMCglmm::inverseA(tree, nodes = "TIPS")$Ainv
	
	# Name of the column in the data corresponding to the phylogeny
		name  <- all.vars(randformula)
		
	# Fit the model using MCMCglmm
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
  	return(model)
}

#' @title make_prior
#' @description Function creates the prior for the MCMCglmm model
#' @param n_rand An integer specifying the number of random effects in the model.
#' @param type A string that specifies the type of model to fit.
#' @param nu A numeric specifying the nu parameter for the prior.
#' @return A list of priors for the MCMCglmm model.
#' @export
 
make_prior <- function(n_rand, type, nu = NULL) {

	if(type == "gaussian") {
		if(is.null(nu)) {nu <- 0.002}
		prior <- list(R = list(V = 1, nu = nu),
					  G = list(G1 = list(V = diag(n_rand), nu = nu)))
	}

	if(type == "poisson") {
		if(is.null(nu)) {nu <- 0.02}
		prior <- list(R = list(V = 1, nu = nu),
                      G = list(G1 = list(V = diag(n_rand), nu = nu)))
	}

	if(type == "categorical") {
		J <- length(unique(ph_sub)) #number of categories
  		number_fix_parameters <- ncol(X_model_matrix_1_sub) * (J-1)

		# Get the number of random effects variables
		number_random_parameters <- n_rand * (J - 1)

		J_matrix <- array(1, dim = c(J, J) - 1) # matrix of ones
		I_matrix <- diag(J - 1) #identity matrix

		IJ <- (I_matrix + J_matrix)/J # see Hadfield's Course notes p. 97
		prior <- list(R = list(V = IJ, fix = 1),
					  G = list(G1 = list(V = diag(number_random_parameters), nu = J+number_random_effects)))
	}

	if(type == "threshold") {
		if(is.null(nu)) {nu <- 0.02}
		prior <- list(R = list(V = 1, fix = TRUE),
                      G = list(G1 = list(V = diag(n_rand), nu = nu)))
	}

	if(type == "ordinal") {
		# Get the number of random effects variables
		number_random_effects <- length(znames_1)
		number_random_parameters <- number_random_effects
		#Fix residual variance R at 1
		# cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm

		J <- length(table(y_imp)) #number of categories
		#priors from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q2/018194.html

		#The residual variance has to be fixed, because the latent variable is without scale
		prior <- list(R = list(V = 1, fix = TRUE),
			G = list(G1 = list(V = diag(number_random_effects), nu = 0.002)))
	}

	return(prior)
}


#' @title predict_bace
#' @description Function creates a predcition from MCMCglmm model
#' @param model A MCMCglmm model object
#' @param dat_prep A list containing the prepared data frame and attributes for continuous variables.
#' @param type A string that specifies the type of model to fit.
#' @return A vector of predicted values from the MCMCglmm model.
#' @export
predict_bace <- function(model, dat_prep, type = NULL, ...) {				

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

				  # Check there are only two levels for threshold model
				  if(type == "threshold" && length(levels_var) != 2){
					  stop("Threshold models can only be used for binary data with two levels.")
				  }
				   
				   # Predicts probabilities for each category
				 pred_prob <- predict(model, marginal = NULL, type = "response", posterior = "all")

				   # For each observation, sample from the categorical distribution based on the predicted probabilities
				   pred_values <- apply(pred_prob, 1, function(probs) {
					   sample(levels_var, size = 1, prob = c(probs, 1-probs))
				   })
				 }
				 
				 if(type == "categorical"){
					# Identify number of categories and their levels from the data
			      	     lv  <- dat_prep[[1]][[response_var]]
				  levels_var <- sort(unique(as.character(lv)))
					
					# Predict category probabilities
					pred_prob <- pred_cat(model, baseline_name = levels_var[1])

					# For each observation, sample from the categorical distribution based on the predicted probabilities
				   pred_values <- apply(pred_prob, 1, function(probs) {
					   sample(levels_var, size = 1, prob = probs)
				   })
				 }

	return(pred_values)
}

#' @title pred_cat
#' @description Function calculates predicted probabilities for each category from a categorical MCMCglmm model
#' 	@param model A MCMCglmm model object
#' @param baseline_name A string specifying the name of the baseline category
#' @return A data frame of predicted probabilities for each category
#' @export
pred_cat <- function(model, baseline_name = "Baseline") {
  
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
  v_total <- rowSums(model$VCV)
  scaling_factor <- sqrt(1 + c2 * v_total)
  
  # 3. Extract and Scale Liabilities
  exp_liab_list <- list()
  # Initialize exp_sum with 1 (which is exp(0) for the baseline)
  exp_sum <- 1 
  
  for (i in 1:n_traits) {
    cols <- ((i - 1) * n_obs + 1):(i * n_obs)
    
    # Apply c2 scaling to the liabilities before exponentiating
    # We divide each sample's liabilities by that sample's scaling factor
    scaled_liab <- model$Liab[, cols] / scaling_factor
    
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

#' @title pred_threshold
#' @description Function calculates predicted probabilities for each category from a threshold MCMCglmm model
#' @param model A MCMCglmm model object
#' @param level_names A vector of strings specifying the names of the levels/categories
#' @return A data frame of predicted probabilities for each category
#' @export	
pred_threshold <- function(model, level_names = NULL) {
  # 1. Dimensions and Data
  n_obs <- model$Residual$nrl
  liab <- model$Liab
  v_total <- rowSums(model$VCV)
  scaling_factor <- sqrt(1 + v_total)
  
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
    
    all_probs[[j]] <- colMeans(p_cat) * 100
  }
  
  # 4. Final Data Frame
  df_final <- as.data.frame(do.call(cbind, all_probs))
  
  # 5. Naming and Ordering
  if (is.null(level_names)) {
    colnames(df_final) <- paste0("Level_", 1:n_levels)
  } else {
    colnames(df_final) <- level_names
  }
  
  # Sort alphabetically (e.g., Aquatic, Insessorial, Terrestrial)
  df_final <- df_final[, order(colnames(df_final))]
  rownames(df_final) <- paste0("Obs_", 1:n_obs)
  
  return(round(df_final, 2))
}

#' @title get_imputed
#' @description Function extracts the prior for the MCMCglmm model
#' @param model An integer specifying the number of random effects in the model.
#' @param type A string that specifies the type of model to fit.
#' @return imputated data based the MCMCglmm model.
#' @export
 
get_imputed <- function(model, type) {
    
    # the number of estimated random effects without Node



	# Extract the imputed data from the model

  posterior.mean <- summary(model$Sol)$statistics[,1]
  # get posterior mean of random effects
  zpost.mean.extra <- posterior.mean[(ncol(model$X) + 1):length(posterior.mean)]
  # remove the .Node columns
  zpost.mean <- zpost.mean.extra[!grepl("\\.Node", names(zpost.mean.extra))] 
  # get posteior mean of fixed effects
  xpost.mean <- posterior.mean[1:ncol(model$X)]
  # get residual variance  
  varaince.e <- as.numeric(posterior.mode(model$VCV)[dim(summary(model$VCV)$statistics)[1]]) 
  # the last column contains the variance (not standard deviation) of the residuals

  #number_of_draws <- nrow(pointdraws)
  #select.record <- sample(1:number_of_draws, size = 1)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  ###start imputation
  rand.eff.imp <- matrix(zpost.mean, ncol = n.par.rand)

  fix.eff.imp <- matrix(xpost.mean, nrow = ncol(model$X))

  sigma.y.imp <- sqrt(varaince.e)

	if(type == "normal") {
		imputed_y <- fix.eff.imp + X_model_matrix_1_sub %*% t(rand.eff.imp) + rnorm(nrow(X_model_matrix_1_sub), 0, sigma.y.imp)
	}

	if(type == "poisson") {
		prior <- list(R = list(V = 1e-07, nu = -2),
                      G = list(G1 = list(V = diag(ncol(Z)), nu = nu)))
	}

	if(type == "categorical") {
		J <- length(unique(ph_sub)) #number of categories
  				number_fix_parameters <- ncol(X_model_matrix_1_sub) * (J-1)

		# Get the number of random effects variables
		number_random_effects <- length(znames_1)
		number_random_parameters <- number_random_effects * (J - 1)
		#Fix residual variance R at 1
		# cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm

		J_matrix <- array(1, dim = c(J, J) - 1) # matrix of ones
		I_matrix <- diag(J - 1) #identiy matrix

		IJ <- (I_matrix + J_matrix)/J # see Hadfield's Course notes p. 97
		prior <- list(R = list(V = IJ, fix = 1),
					  G = list(G1 = list(V = diag(number_random_parameters), nu = J+number_random_effects)))
	}

	if(type == "threshold") {
		prior <- list(R = list(V = 1, fix = TRUE),
                      G = list(G1 = list(V = diag(n_rand), nu = nu)))
	}

	if(type == "ordinal") {
		# Get the number of random effects variables
		number_random_effects <- length(znames_1)
		number_random_parameters <- number_random_effects
		#Fix residual variance R at 1
		# cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm

		J <- length(table(y_imp)) #number of categories
		#priors from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q2/018194.html

		#The residual variance has to be fixed, because the latent variable is without scale
		prior <- list(R = list(V = 1, fix = TRUE),
			G = list(G1 = list(V = diag(number_random_effects), nu = 0.002)))
	}

	return(imputed_y)
}


#' @title predict.bace
#' @description Function creates a predcition from MCMCglmm model
#' @param model A MCMCglmm model object
#' @return A list of priors for the MCMCglmm model.
#' @export
#' 
#' 
predict.bace <- function(model, type = NULL, ...) {

	if(is.null(type)) {
		pred <- as.vector(stats::predict(model, marginal = NULL))
	}

	if(type == "categorical") {
		
	}

return(pred)
}

