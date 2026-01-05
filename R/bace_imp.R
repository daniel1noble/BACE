
#' @title bace_imp
#' @description Function to perform Bayesian imputation using BACE for missing data in a dataset
#' @param fixformula A character string specifying the fixed effects formula used in the model. This should be of the form: y ~ x.
#' @param ran_phylo_form A character string specifying the random effects and phylogenetic structure formula used in the model.
#' @param species A logical indicating whether species is included as a random effect.
#' @param phylo A phylogenetic tree of class 'phylo' from the ape package.
#' @param data A data frame containing the dataset with missing values to be imputed.
#' @param ... Additional arguments to be passed to the underlying modeling functions.
#' @return A list containing imputed datasets and model summaries.
#' @export

bace_imp <- function(fixformula, ran_phylo_form, species = TRUE, phylo, data, ...){

	# First, we need to get the variables from the formulas that are used to subset out of the data. Returns a list of variables
		  fix <- get_variables(fixformula)
		phylo <- get_variables(ran_phylo_form, fix = FALSE)

	# Subset the data for fixed and random effects. Need to clean up column names for randdata
		 fixdata <- data[, fix[["fix"]]]
		phylodata <- data.frame(data[, unlist(phylo)])
			if(ncol(phylodata) == 1){
				names(phylodata) <- phylo[["cluster"]]
			}

	# We need to obtain the class/type of the variables. list.
	       types <- get_types(fixdata)

	# Now, create the formulas for all the fixed effect variables that will change as we move through the chained equations. list.
		formulas <- build_formula_string(fixformula)

	# How we need to use the data, type of variable class to fit the models
		for(i in var){
			prior <- make_prior(i, type = types)
			model <- model_fit(i, type = types)
		}
}