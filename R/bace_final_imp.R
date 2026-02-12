#' @title bace_final_imp
#' @description Run final imputation iterations starting from converged data, saving all models
#' @param bace_object An object of class 'bace' from bace_imp
#' @param fixformula A character string or list specifying the fixed effects formulas
#' @param ran_phylo_form A character string specifying the random effects and phylogenetic structure formula
#' @param phylo A phylogenetic tree of class 'phylo' from the ape package
#' @param nitt An integer or list specifying the number of MCMC iterations. Default is 6000
#' @param thin An integer or list specifying the thinning rate. Default is 5
#' @param burnin An integer or list specifying the burn-in iterations. Default is 1000
#' @param n_final An integer specifying the number of final imputation runs. Default is 10
#' @param species A logical indicating whether to decompose phylogenetic and non-phylogenetic species effects. Default is FALSE
#' @param verbose A logical indicating whether to print progress messages. Default is TRUE
#' @param ... Additional arguments passed to modeling functions
#' @return A list of class 'bace_final' containing:
#'   - all_models: List of length n_final, each containing models for all variables
#'   - all_datasets: List of length n_final with imputed datasets
#'   - formulas: The formulas used
#'   - types: Variable types
#'   - phylo_ran: Phylogenetic random effects structure
#'   - call: The function call
#' @examples \dontrun{
#' # After running bace_imp and checking convergence
#' result <- bace_imp(fixformula = "y ~ x1 + x2", ...)
#' final <- bace_final_imp(result, fixformula = "y ~ x1 + x2", ...)
#' }
#' @export
bace_final_imp <- function(bace_object, fixformula, ran_phylo_form, phylo, 
                           nitt = 6000, thin = 5, burnin = 1000, 
                           n_final = 10, species = FALSE, verbose = TRUE, ...) {
  
  # Check inputs
  if (!inherits(bace_object, "bace")) {
    stop("bace_object must be an object of class 'bace' from bace_imp function")
  }
  
  # Extract components from bace_object
  last_data <- bace_object$data[[length(bace_object$data)]]
  miss_dat <- bace_object$miss_dat
  types <- bace_object$types
  phylo_ran <- bace_object$phylo_ran
  
  # Build formulas if needed
  if (!is.list(fixformula)) {
    formulas <- .build_formula_string(fixformula) 
  } else {
    formulas <- lapply(fixformula, as.formula)
  }
  
  # Random effect formula with species parameter
  ran_phylo_form <- .build_formula_string_random(ran_phylo_form, species = species)
  
  # Standardize MCMC parameters
  n_models <- length(formulas)
  nitt_list <- .standardize_mcmc_params(nitt, n_models, "nitt")
  thin_list <- .standardize_mcmc_params(thin, n_models, "thin")
  burnin_list <- .standardize_mcmc_params(burnin, n_models, "burnin")
  
  # Variables to impute
  fix <- names(types)
  
  # Storage for all runs
  all_models <- vector("list", n_final)
  all_datasets <- vector("list", n_final)
  
  # Starting dataset
  data_current <- last_data
  
  if (verbose) {
    cat("\n=== Running", n_final, "final imputation iterations ===\n\n")
  }
  
  # Run n_final imputation iterations
  for (run in 1:n_final) {
    
    # Storage for models in this run
    models_this_run <- list()
    
    # Loop through formulas to fit models and impute
    for (i in 1:length(formulas)) {
      
      # Identify response variable
      response_var <- all.vars(formulas[[i]][[2]])
      
      # Prepare data
      dat_prep <- .data_prep(
        data = data_current, 
        formula = formulas[[i]], 
        types = types, 
        ran_cluster = phylo_ran[["cluster"]]
      )
      
      data_i <- dat_prep[[1]]
      
      # If species decomposition is enabled, add a second copy of the species column
      # This is needed for MCMCglmm to distinguish between phylo and non-phylo effects
      # MUST happen BEFORE prior creation
      if (species) {
        species_col_name <- phylo_ran[["cluster"]]
        data_i[[paste0(species_col_name, "2")]] <- data_i[[species_col_name]]
      }
      
      # Check if variable needs imputation (has missing data)
      has_missing <- response_var %in% miss_dat$colname
      
      if (has_missing) {
        # Get prior
        fixform <- formulas[[i]]
        levels <- dat_prep$levels
        gelman <- NULL  # Use default prior for final runs
        
        n_rand_eff <- if (species) 2 else 1
        
        prior_i <- .make_prior(
          n_rand = n_rand_eff, 
          n_levels = levels,
          type = types[[response_var]], 
          fixform = fixform, 
          data = data_i, 
          gelman = gelman
        )
        
        # Fit model
        model <- .model_fit(
          data = data_i, 
          tree = phylo,
          fixformula = formulas[[i]],
          randformula = ran_phylo_form,
          type = types[[response_var]], 
          prior = prior_i, 
          nitt = nitt_list[[i]], 
          thin = thin_list[[i]], 
          burnin = burnin_list[[i]]
        )
        
        # Store model
        models_this_run[[response_var]] <- model
        
        # Predict missing data
        predictions <- .predict_bace(
          model, 
          dat_prep, 
          response_var = response_var, 
          type = types[[response_var]]
        )
        
        # Update data with new imputations
        id <- miss_dat[miss_dat$colname == response_var, "row"]
        data_id <- which(colnames(data_current) == response_var)
        data_current[id, data_id] <- predictions[["pred_values"]][id]
        
        if (verbose) {
          cat(paste0("Run ", run, "/", n_final, 
                     " - Imputed variable: ", response_var, "\n"))
        }
      } else {
        # Variable has no missing data, still fit model but don't impute
        if (verbose) {
          cat(paste0("Run ", run, "/", n_final, 
                     " - Variable ", response_var, 
                     " has no missing data (skipped imputation)\n"))
        }
      }
    }
    
    # Store results for this run
    all_models[[run]] <- models_this_run
    all_datasets[[run]] <- data_current
    
    if (verbose && run %% max(1, floor(n_final/4)) == 0) {
      cat(paste0("\nCompleted ", run, "/", n_final, " runs\n\n"))
    }
  }
  
  # Prepare output
  out <- list(
    all_models = all_models,
    all_datasets = all_datasets,
    formulas = formulas,
    types = types,
    phylo_ran = phylo_ran,
    miss_dat = miss_dat,
    n_final = n_final,
    call = match.call()
  )
  
  class(out) <- "bace_final"
  
  if (verbose) {
    cat("\n=== Final imputation complete ===\n")
    cat("Saved", n_final, "imputed datasets with corresponding models\n")
  }
  
  return(out)
}


#' @title print.bace_final
#' @description Print method for bace_final objects
#' @param x Object of class bace_final
#' @param ... Additional arguments
#' @export
print.bace_final <- function(x, ...) {
  cat("\n=== BACE Final Imputation Results ===\n\n")
  cat("Number of imputation runs:", x$n_final, "\n")
  cat("Number of variables with models:", 
      length(x$all_models[[1]]), "\n")
  cat("Variables:", paste(names(x$all_models[[1]]), collapse = ", "), "\n\n")
  cat("Use pool_posteriors() to combine model results across imputations\n")
}
