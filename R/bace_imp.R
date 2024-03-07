



bace_imp <- function(fixformula,ran_phylo_form, species = TRUE, phylo, data, ...){

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