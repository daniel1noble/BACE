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
                              verbose = FALSE, pr = TRUE,
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
						pred_values <- round(predict(model, marginal = NULL, type = "response", posterior = "all"), digits = 0)
				 }

				 if(type == "threshold" || type == "categorical" || type == "ordinal"){
					# Identify number of categories and their levels from the data
			      levels_var <- levels(dat_prep[[1]][[response_var]])
				   
				   # Predicts probabilities for each category
				 pred_prob <- predict(model, marginal = NULL, type = "response", posterior = "all")

				   # For each observation, sample from the categorical distribution based on the predicted probabilities
				   pred_values <- apply(pred_prob, 1, function(probs) {
					   sample(levels_var, size = 1, prob = probs)
				   })
				 }
				 
				return(pred_values)
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

