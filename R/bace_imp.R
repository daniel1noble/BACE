
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
#' data <- data.frame(y = rpois(30, lambda = 5), x1 = factor(rep(c("A", "B","A"), length.out = 30)), x2 = rnorm(30, 10, 2), x3 = factor(rep(c("A", "B", "C", "D", "E"), length.out = 30)), x4 = factor(rep(c("A", "B", "C", "D", "E"), length.out = 30), levels = c("B", "A", "C", "D", "E"), ordered = TRUE), Species = phylo$tip.label)
#' # Introduce some missing data
#' missing_indices <- sample(1:30, 10)
#' data$y[missing_indices] <- NA
#' data$x1[sample(1:30, 5)] <- NA
#' data$x2[sample(1:30, 5)] <- NA	
#' data$x3[sample(1:30, 5)] <- NA
#' data$x4[sample(1:30, 5)] <- NA	
#' # Run BACE imputation
#' bace_imp(fixformula = "y ~ x1 + x2", ran_phylo_form = "~ 1 |Species", phylo = phylo, data = data)
#' bace_imp(fixformula = list("y ~ x1 + x2", "x2 ~ x1", "x1 ~ x2", "x3 ~ x1 + x2", "x4 ~ x1 + x2"), ran_phylo_form = "~ 1 |Species", phylo = phylo, data = data, runs = 5)
#' }
#' @export
#' options(error = recover)
bace_imp <- function(fixformula, ran_phylo_form, phylo, data, nitt = 5000, thin = 10, burnin = 1000, runs = 10, ...){
	#---------------------------------------------#
	# Preparation steps & Checks
	#---------------------------------------------#
	# First, we need to get the variables from the formulas that are used to subset out of the data. Returns a list of variables. Check if a list of formulas is provided and if not then create formulas
	
	if(!is.list(fixformula)){
	    fix <- as.vector(unlist(.get_variables(fixformula)))
	} else {
		fix <- unique(as.vector(unlist(lapply(fixformula, .get_variables))))
	}
	
	# Get the random effect and phylogenetic variables
		phylo_ran <- .get_variables(ran_phylo_form, fix = FALSE)

	# Check that all variable names are in the dataframe and if not stop
		all_vars <- c(fix, phylo_ran[["cluster"]])
			if(!all(all_vars %in% colnames(data))){
				missing_vars <- all_vars[!all_vars %in% colnames(data)]
				stop(paste("The following variables are not in the dataframe: ", paste(missing_vars, collapse = ", ")))
			}

	# Subset the data for fixed and random effects to keep data focused. Need to clean up column names for randdata
		 data_sub <- data[, c(fix, phylo_ran[["cluster"]])]
		 data_sub2 <- data_sub  # Create a copy to hold imputed data across runs

	# Get a summary of the data types, classes etc
	   data_summary <- .summarise_var_types(data_sub)
	
	#---------------------------------------------#
	# Build formulas if needed otherwise user specified
	#---------------------------------------------#
	# Now, create the formulas for all the fixed effect variables that will change as we move through the chained equations. list.
	if(!is.list(fixformula)){
		formulas <- .build_formula_string(fixformula) 
	} else {
		formulas <- lapply(fixformula, as.formula)
	}
	
	# Random effect formula. TO DO: Need to be more sophisticated here to allow for multiple random effects
	   ran_phylo_form <- .build_formula_string_random(ran_phylo_form)

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
	
		# List to hold runs. 
		pred_missing_run <- list()

		# Store initial data as first run
		     pred_missing_run[[1]] <- data_sub
		names(pred_missing_run)[1] <- "Initial_Data"

		# We need to categorize the variable types for all fixed effect variables for modelling
	               types <- lapply(fix, function(var) .get_type(data_summary, var, data_sub))
			names(types) <- fix
	
	# Now, we can loop through the number of runs specified by the user
		for(r in 2:(runs+1)){
			for(i in 1:length(formulas)){

			# Identify the response variable for the current formula			
			response_var <- all.vars(formulas[[i]][[2]])
			
			# Prepare the data 
			    dat_prep <- .data_prep(data = pred_missing_run[[1]], formula = formulas[[i]], types = types)
			
			 if(r == 2){
				# For iteration 1 we want to impute missing data as a rough approximation z-transform for continuous data, and random sampling from the empirical distribution for categorical data.
				data_i <- dat_prep[[1]]

			 } else {
				data_i <- dat_prep[[1]]

					# predictors in *this* subset
					predictor_vars <- setdiff(names(data_i), response_var)

					# missing entries that correspond to predictors present in this subset
					mm <- missing_matrix[missing_matrix$colname %in% predictor_vars, , drop = FALSE]

					# predictions from previous run
					pred_full <- pred_missing_run[[r - 1]]

					# --- Safety checks ---
					stopifnot(is.data.frame(mm))
					stopifnot(all(c("row", "colname") %in% names(mm)))

					# Ensure rows are valid for data_i
					mm <- mm[mm$row >= 1 & mm$row <= nrow(data_i), , drop = FALSE]

					# Ensure predicted object has the needed columns (by name)
					if (!all(mm$colname %in% colnames(pred_full))) {
					missing_cols <- setdiff(unique(mm$colname), colnames(pred_full))
					stop("pred_full is missing columns: ", paste(missing_cols, collapse = ", "))
					}

					# Split rows by column name, then fill
					rows_by_col <- split(mm$row, mm$colname)

					for (cn in names(rows_by_col)) {
					rows <- rows_by_col[[cn]]

					# Use [[ ]] to avoid data.frame coercion issues
					data_i[[cn]][rows] <- pred_full[[cn]][rows]
					}
			 }

			# Set up prior for the specific variable type
			   levels_cat <- data_summary$n_levels[which(data_summary$variable == response_var)]
				  levels  <- if(anyNA(levels_cat)) NULL else levels_cat
			      prior_i <- .make_prior(n_rand = length(phylo_ran$ran), n_levels = levels, type = types[[response_var]])
			
			# Fit the model and predict missing data
			  model <-  .model_fit(data = data_i, 
			                       tree = phylo,
			                 fixformula = formulas[[i]],  # Only simple structure for now
						    randformula = ran_phylo_form, # Only phylogeny, but need to add species too
					 		       type = types[[response_var]], 
							      prior = prior_i, 
								   nitt = nitt, 
								   thin = thin, 
								 burnin = burnin)         # Prior Not working for all types yet

			# Predict missing data and store in list to keep track across runs, if the variable was z-transformed then transform back using the attributes from data_i preparation which is done automatically for gaussian variables
				pred_values <- .predict_bace(model, dat_prep, response_var = response_var, type = types[[response_var]])
			 
			# Store predicted values for only the missing data points
			                         id <- missing_matrix[with(missing_matrix, colname == response_var), "row"]
							    data_id <- which(colnames(data_sub2) == response_var)
			     data_sub2[id, data_id]  <- pred_values[id]
			 
		} # End of formula loop		
		     pred_missing_run[[r]] <- data_sub2
		names(pred_missing_run)[r] <- paste0("Iter_", r-1)
}

       out <- list(iterations = pred_missing_run)
class(out) <- "bace"

return(out)
}	
