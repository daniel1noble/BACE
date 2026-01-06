
#' @title bace_imp
#' @description Function to perform Bayesian imputation using BACE for missing data in a dataset
#' @param fixformula A character string specifying the fixed effects formula that is of major interest or a list of formulas that one wishes to estimate. This should be of the form: y ~ x. Providing a list gives users complex flexibility in the types of models fit for the different variables. Note that users can also include interaction terms using the standard R formula syntax (e.g., x1*x2). 
#' @param ran_phylo_form A character string specifying the random effects and phylogenetic structure formula used in the model.
#' @param species A logical indicating whether species is included as a random effect.
#' @param phylo A phylogenetic tree of class 'phylo' from the ape package.
#' @param data A data frame containing the dataset with missing values to be imputed.
#' @param runs An integer specifying the number of imputation iterations to perform. Default is 10.
#' @param ... Additional arguments to be passed to the underlying modeling functions.
#' @return A list containing imputed datasets and model summaries.
#' @examples \dontrun{
#' phylo <- force.ultrametric(ape::rtree(5)) # Example phylogenetic tree with 5 tips
#' data <- data.frame(y = c(1,2,3,NA,5), x1 = factor(c("A","B","A","B","A")), x2 = c(10,20,30,40,50), Species = phylo$tip.label)
#' bace_imp(fixformula = "y ~ x1 + x2", ran_phylo_form = "~ 1 |Species", phylo = phylo, data = data)
#' bace_imp(fixformula = list("y ~ x1 + x2", "x2 ~ x1"), ran_phylo_form = "~ 1 |Species", phylo = phylo, data = data, runs = 20)
#' }
#' @export

bace_imp <- function(fixformula, ran_phylo_form, phylo, data, runs = 10L, nitt = 50000, thin = 10, burnin = 1000, ...){

	# First, we need to get the variables from the formulas that are used to subset out of the data. Returns a list of variablesCheck if a list of formulas is provided and if not then create formulas
	
	if(!is.list(fixformula)){
	    fix <- as.vector(unlist(get_variables(fixformula)))
	} else {
		fix <- unique(as.vector(unlist(lapply(fixformula, get_variables))))
	}
	
	# Get the random effect and phylogenetic variables
		phylo_ran <- get_variables(ran_phylo_form, fix = FALSE)

	# Subset the data for fixed and random effects. Need to clean up column names for randdata
		 data_sub <- data[, c(fix, phylo_ran[["cluster"]])]

	# We need to obtain the class/type of the variables. list.
	       types <- get_type(data_sub)[fix]

	# Now, create the formulas for all the fixed effect variables that will change as we move through the chained equations. list.
	if(!is.list(fixformula)){
		formulas <- build_formula_string(fixformula) 
	} else {
		formulas <- lapply(fixformula, as.formula)
	}
	
	# Random effect formula
	   ran_phylo_form <- build_formula_string_random(ran_phylo_form)

	# How we need to use the data, type of variable class to fit the models
	complete_data_list  <- list()

	# For iteration 1 we want to impute missing data as a rough approximation 

		for(j in 2:runs){
			for(i in 1:length(fix$fix)){
			prior_i <- make_prior(n_rand = length(phylo_ran$ran), type = types[[i]])
			 
			  model <-  model_fit(data = data, 
			                      tree = phylo,
			                fixformula = formulas[[i]],  # Only simple structure for now
						   randformula = ran_phylo_form, # Only phylogeny, but need to add species too
							      type = types[[i]], 
							     prior = prior_i, 
								  nitt = nitt, 
								  thin = thin, 
								burnin = burnin)        # Prior Not working for all types yet


			complete_data_list[[j]]  <-  predict(model, data_list[[i]], marginal = NULL, posterior = "all")
		}
		}
}