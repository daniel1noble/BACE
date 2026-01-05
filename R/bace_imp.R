
#' @title bace_imp
#' @description Function to perform Bayesian imputation using BACE for missing data in a dataset
#' @param fixformula A character string specifying the fixed effects formula used in the model. This should be of the form: y ~ x.
#' @param ran_phylo_form A character string specifying the random effects and phylogenetic structure formula used in the model.
#' @param species A logical indicating whether species is included as a random effect.
#' @param phylo A phylogenetic tree of class 'phylo' from the ape package.
#' @param data A data frame containing the dataset with missing values to be imputed.
#' @param ... Additional arguments to be passed to the underlying modeling functions.
#' @return A list containing imputed datasets and model summaries.
#' @examples \dontrun{
#' phylo <- ape::rtree(5) # Example phylogenetic tree with 5 tips
#' data <- data.frame(y = c(1,2,3,NA,5), x1 = factor(c("A","B","A","B","A")), x2 = c(10,20,30,40,50), Species = phylo$tip.label)
#' bace_imp(fixformula = "y ~ x1 + x2", ran_phylo_form = "~ 1 |Species", phylo = phylo, data = data)
#' }
#' @export

bace_imp <- function(fixformula, ran_phylo_form, phylo, data, ...){

	# First, we need to get the variables from the formulas that are used to subset out of the data. Returns a list of variables
		      fix <- get_variables(fixformula)
		phylo_ran <- get_variables(ran_phylo_form, fix = FALSE)

	# Subset the data for fixed and random effects. Need to clean up column names for randdata
		 fixdata <- data[, fix[["fix"]]]
		phylodata <- data.frame(data[, phylo_ran[["cluster"]]])
			if(ncol(phylodata) == 1){
				names(phylodata) <- phylo[["cluster"]]
			}

	# We need to obtain the class/type of the variables. list.
	       types <- get_type(fixdata)

	# Now, create the formulas for all the fixed effect variables that will change as we move through the chained equations. list.
		formulas <- build_formula_string(fixformula)

	# How we need to use the data, type of variable class to fit the models
		for(i in 1:length(fix$fix)){
			prior <- make_prior(n_rand = length(phylo_ran$ran), type = types[[i]])
			model <-  model_fit(data = data_list[[i]], 
			              fixformula = formulas[[i]],  # Only simple structure for now
						 randformula = ran_phylo_form, # Only phylogeny, but need to add species too
							    type = types[[i]], 
							   prior = prior)          # Prior Not working for all types yet
		}
}