
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
#' set.seed(123)
#' phylo <- force.ultrametric(ape::rtree(30)) # Example phylogenetic tree with 30 tips
#' phylo  <- compute.brlen(phylo, method = "Grafen")
#' data <- data.frame(y = rpois(30, lambda = 5), x1 = factor(rep(c("A", "B","A"), length.out = 30)), x2 = rnorm(30, 10, 2), Species = phylo$tip.label)
#' # Introduce some missing data
#' missing_indices <- sample(1:30, 10)
#' data$y[missing_indices] <- NA
#' data$x1[sample(1:30, 5)] <- NA
#' data$x2[sample(1:30, 5)] <- NA	
#' # Run BACE imputation
#' bace_imp(fixformula = "y ~ x1 + x2", ran_phylo_form = "~ 1 |Species", phylo = phylo, data = data)
#' bace_imp(fixformula = list("y ~ x1 + x2", "x2 ~ x1"), ran_phylo_form = "~ 1 |Species", phylo = phylo, data = data, runs = 20)
#' }
#' @export

bace_imp <- function(fixformula, ran_phylo_form, phylo, data, runs = 10L, nitt = 50000, thin = 10, burnin = 1000, runs = 10, ...){
	#---------------------------------------------#
	# Preparation steps & Checks
	#---------------------------------------------#
	# First, we need to get the variables from the formulas that are used to subset out of the data. Returns a list of variables. Check if a list of formulas is provided and if not then create formulas
	
	if(!is.list(fixformula)){
	    fix <- as.vector(unlist(get_variables(fixformula)))
	} else {
		fix <- unique(as.vector(unlist(lapply(fixformula, get_variables))))
	}
	
	# Get the random effect and phylogenetic variables
		phylo_ran <- get_variables(ran_phylo_form, fix = FALSE)

	# Check that all variable names are in the dataframe and if not stop
		all_vars <- c(fix, phylo_ran[["cluster"]])
			if(!all(all_vars %in% colnames(data))){
				missing_vars <- all_vars[!all_vars %in% colnames(data)]
				stop(paste("The following variables are not in the dataframe: ", paste(missing_vars, collapse = ", ")))
			}

	# Subset the data for fixed and random effects to keep data focused. Need to clean up column names for randdata
		 data_sub <- data[, c(fix, phylo_ran[["cluster"]])]

	#---------------------------------------------#
	# Build formulas if needed otherwise user specified
	#---------------------------------------------#
	# Now, create the formulas for all the fixed effect variables that will change as we move through the chained equations. list.
	if(!is.list(fixformula)){
		formulas <- build_formula_string(fixformula) 
	} else {
		formulas <- lapply(fixformula, as.formula)
	}
	
	# Random effect formula. TO DO: Need to be more sophisticated here to allow for multiple random effects
	   ran_phylo_form <- build_formula_string_random(ran_phylo_form)

	#---------------------------------------------#
	# Create indicators for where missing data are
	#---------------------------------------------#
	# Identify all cells indicators within a the data that are missing 
	   missing_matrix <- which(is.na(data_sub), arr.ind = TRUE)

	# Add a new column to indicate the variable (column) names for easier tracking
	           missing_matrix <- data.frame(missing_matrix)
	   missing_matrix$colname <- colnames(data_sub)[missing_matrix[,2]]

	#---------------------------------------------#
	# BACE Imputation
	#---------------------------------------------#
	# How we need to use the data, type of variable class to fit the models
	
		# List to hold runs
		pred_missing_run <- list()

		for(r in 1:runs){

			# List to hold variables predictions for this run
			pred_missing_vars <- matrix(0, nrow = nrow(data_sub), ncol = length(fix))
			
			for(i in 1:length(fix)){

			# We need to obtain the class/type of the variables. list.
	          types <- get_type(data_sub)
			
			# identify the response variable for the current formula			
			response_var <- all.vars(formulas[[i]][[2]])
			
			# Prepare the data 
			dat_prep <- data_prep(data = data_sub, formula = formulas[[i]], types = types)
			
			 if(r == 1){
				# For iteration 1 we want to impute missing data as a rough approximation z-transform for continuous data, and random sampling from the empirical distribution for categorical data.
				data_i <- dat_prep[[1]]

			 } else {
				# If not the first iteration then we use the predicted values for the missing cells of the predictors from the previous iteration to fill in missing data for the response
				data_i <- data_sub 

				# remove response variable from missing matrix
				col_response <- which(colnames(data_i) == response_var)
				
				# Filter out response variable column from missing matrix
				missing_matrix_pred <- data.frame(missing_matrix)  %>% dplyr::filter(col != col_response)

				# Fill in missing values with predicted values from previous run
				         idx <- cbind(missing_matrix_pred$row, missing_matrix_pred$col)
				 data_i[idx] <- pred_missing_run[[r-1]][idx]
			 }

			# Set up prior for the specific variable type
			     prior_i <- make_prior(n_rand = length(phylo_ran$ran), type = types[[response_var]])
			
			# Fit the model and predict missing data
			  model <-  model_fit(data = data_i, 
			                      tree = phylo,
			                fixformula = formulas[[i]],  # Only simple structure for now
						   randformula = ran_phylo_form, # Only phylogeny, but need to add species too
							      type = types[[response_var]], 
							     prior = prior_i, 
								  nitt = nitt, 
								  thin = thin, 
								burnin = burnin)         # Prior Not working for all types yet

			# Predict missing data and store in list to keep track across runs, if the variable was z-transformed then transform back using the attributes from data_i preparation which is done automatically for gaussian variables
				pred_values <- predict_bace(model, dat_prep, type = types[[response_var]])
			 
			# Store predicted values for only the missing data points
			                       id <- missing_matrix[with(missing_matrix, colname == response_var), "row"]
			 pred_missing_vars[id,i]  <- pred_values[id]
			 
		}
		
		     pred_missing_run[[r]] <- pred_missing_vars
		names(pred_missing_run)[r] <- paste0("Iter_", r)
}
return(list(iterations = pred_missing_run))
}	
	

