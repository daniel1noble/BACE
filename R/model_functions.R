#' @title model_fit
#' @description Function takes incomplete data and a formula string and fits a model using MCMCglmm
#' @param data A dataframe containing missing data.
#' @param fixformula A string that specifies the fixed effects in the model.
#' @param randformula A string that specifies the random effects in the model.
#' @param type A string that specifies the type of model to fit.
#' @param nitt An integer specifying the number of iterations to run the MCMC algorithm.
#' @param thin An integer specifying the thinning rate for the MCMC algorithm.
#' @param burnin An integer specifying the number of iterations to discard as burnin.
#' @param n.par.rand An integer specifying the number of random effects in the model.
#' @return A list of draws from the posterior distribution of the model parameters.
#' @export


model_fit <- function(data, fixformula, randformula, type, nitt, thin, burnin, n_rand) {

	# Create a prior list for residual and random effects
             prior <- make_prior(n_rand, type)

	# Fit the model using MCMCglmm
  	model <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                        random = randformula,
                                          data = data,
                                        family = type,
                                       verbose = FALSE, pr = TRUE, prior = prior,
                                         saveX = TRUE, saveZ = TRUE,
                                          nitt = nitt,
                                          thin = thin,
                                        burnin = burnin)
	class(model) <- "bace"									
  	return(model)
}

#' @title make_prior
#' @description Function creates the prior for the MCMCglmm model
#' @param n_rand An integer specifying the number of random effects in the model.
#' @param type A string that specifies the type of model to fit.
#' @param nu A numeric specifying the nu parameter for the prior.
#' @return A list of priors for the MCMCglmm model.
#' @export
 
make_prior <- function(n_rand, type, nu = 0.002) {

	if(type == "normal") {
		prior <- list(R = list(V = 1, nu = nu),
					  G = list(G1 = list(V = diag(n_rand), nu = nu)))
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

	if(type == "binary") {
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


#' @title predict.bace
#' @description Function creates a predcition from MCMCglmm model
#' @param model A MCMCglmm model object
#' @return A list of priors for the MCMCglmm model.
#' @export
#' 
#' 
predict.bace <- function(model, type = NULL, ...) {

	if(is.null(type)) {
		pred <- as.vector(predict(model, marginal = NULL))
	}

	if(type == "categorical") {
		pred <- 
	}

return(pred)
}
