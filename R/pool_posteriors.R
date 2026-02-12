#' @title pool_posteriors
#' @description Pool posteriors from multiple imputation runs to account for imputation uncertainty
#' @param bace_final_object An object of class 'bace_final' from bace_final_imp
#' @param variable Character string specifying which variable's model to pool. 
#'   If NULL (default), pools all variables
#' @return A list of class 'bace_pooled' containing pooled models for each variable.
#'   Each pooled model is an MCMCglmm object with class c("MCMCglmm", "bace_pooled_MCMCglmm") 
#'   containing:
#'   - Sol: Pooled fixed effect posteriors (combined across imputations)
#'   - VCV: Pooled variance component posteriors (combined across imputations)
#'   - CP: Pooled cutpoint posteriors (for categorical/ordinal models, if present)
#'   - Liab: Pooled latent variables (for categorical models, if present)
#'   - DIC: Mean DIC across imputations
#'   - Deviance: Combined deviance samples
#'   - BACE_pooling: Metadata about pooling (n_imputations, n_samples_per_imputation, etc.)
#'   - All other MCMCglmm components (Fixed, Random, etc.)
#' @details This function implements posterior pooling by concatenating MCMC samples 
#'   from all imputation runs. This naturally accounts for both within-imputation and 
#'   between-imputation variance. The resulting MCMCglmm object can be used with 
#'   standard MCMCglmm methods like print() and summary(), which will display results
#'   that properly reflect uncertainty from both the model and the imputation process.
#'   
#'   The pooled models retain the MCMCglmm class so all standard MCMCglmm methods work,
#'   but also include class "bace_pooled_MCMCglmm" for custom behavior when needed.
#' @importFrom coda as.mcmc
#' @examples \dontrun{
#' # After running bace_final_imp
#' final <- bace_final_imp(bace_obj, ...)
#' pooled <- pool_posteriors(final)
#' 
#' # Extract pooled model for specific variable - works like MCMCglmm!
#' pooled_y <- pooled$y
#' summary(pooled_y)  # Standard MCMCglmm summary
#' print(pooled_y)    # Standard MCMCglmm print
#' plot(pooled_y)     # Standard MCMCglmm plots
#' }
#' @export
pool_posteriors <- function(bace_final_object, variable = NULL) {
  
  # Check inputs
  if (!inherits(bace_final_object, "bace_final")) {
    stop("Input must be an object of class 'bace_final' from bace_final_imp function")
  }
  
  all_models <- bace_final_object$all_models
  n_runs <- length(all_models)
  
  # Get variable names
  if (is.null(variable)) {
    variables <- names(all_models[[1]])
  } else {
    if (!variable %in% names(all_models[[1]])) {
      stop(paste("Variable", variable, "not found in models"))
    }
    variables <- variable
  }
  
  # Storage for pooled models
  pooled_models <- list()
  
  for (var in variables) {
    
    # Extract all models for this variable across runs
    var_models <- lapply(all_models, function(x) x[[var]])
    
    # Check if any are NULL (no missing data for this variable)
    if (all(sapply(var_models, is.null))) {
      warning(paste("Variable", var, "has no models (no missing data). Skipping."))
      next
    }
    
    # Remove NULL entries
    var_models <- var_models[!sapply(var_models, is.null)]
    n_valid <- length(var_models)
    
    if (n_valid == 0) {
      warning(paste("No valid models for variable", var))
      next
    }
    
    # Get dimensions from first model
    first_model <- var_models[[1]]
    n_iter_per_model <- nrow(first_model$Sol)
    n_fixed <- ncol(first_model$Sol)
    n_vcv <- ncol(first_model$VCV)
    
    # Check what components are present in the first model
    has_cp <- !is.null(first_model$CP) && length(first_model$CP) > 0
    has_liab <- !is.null(first_model$Liab) && length(first_model$Liab) > 0
    has_deviance <- !is.null(first_model$Deviance) && length(first_model$Deviance) > 0
    
    # Calculate total number of posterior samples
    total_samples <- n_iter_per_model * n_valid
    
    # Preallocate matrices for pooled posteriors
    pooled_sol <- matrix(NA, nrow = total_samples, ncol = n_fixed)
    colnames(pooled_sol) <- colnames(first_model$Sol)
    
    pooled_vcv <- matrix(NA, nrow = total_samples, ncol = n_vcv)
    colnames(pooled_vcv) <- colnames(first_model$VCV)
    
    if (has_cp) {
      n_cp <- ncol(first_model$CP)
      pooled_cp <- matrix(NA, nrow = total_samples, ncol = n_cp)
      colnames(pooled_cp) <- colnames(first_model$CP)
    }
    
    if (has_liab) {
      n_liab <- ncol(first_model$Liab)
      pooled_liab <- matrix(NA, nrow = total_samples, ncol = n_liab)
      colnames(pooled_liab) <- colnames(first_model$Liab)
    }
    
    if (has_deviance) {
      pooled_deviance <- numeric(total_samples)
    }
    
    # Combine posteriors by stacking
    for (i in 1:n_valid) {
      start_idx <- (i - 1) * n_iter_per_model + 1
      end_idx <- i * n_iter_per_model
      
      pooled_sol[start_idx:end_idx, ] <- var_models[[i]]$Sol
      pooled_vcv[start_idx:end_idx, ] <- var_models[[i]]$VCV
      
      if (has_cp) {
        pooled_cp[start_idx:end_idx, ] <- var_models[[i]]$CP
      }
      
      if (has_liab) {
        pooled_liab[start_idx:end_idx, ] <- var_models[[i]]$Liab
      }
      
      if (has_deviance) {
        pooled_deviance[start_idx:end_idx] <- var_models[[i]]$Deviance
      }
    }
    
    # Calculate mean DIC across imputations
    dics <- sapply(var_models, function(m) {
      if (!is.null(m$DIC)) m$DIC else NA
    })
    mean_dic <- mean(dics, na.rm = TRUE)
    
    # Convert pooled matrices to mcmc objects (required by MCMCglmm methods)
    # Set MCMC parameters: start, end, thin
    # Since we're pooling, we treat the concatenated samples as a single chain
    pooled_sol <- coda::as.mcmc(pooled_sol)
    pooled_vcv <- coda::as.mcmc(pooled_vcv)
    
    if (has_cp) {
      pooled_cp <- coda::as.mcmc(pooled_cp)
    }
    
    if (has_liab) {
      pooled_liab <- coda::as.mcmc(pooled_liab)
    }
    
    # Create pooled model object as a proper MCMCglmm object
    # Start with a copy of the first model's structure
    pooled_model <- first_model
    
    # Replace posterior samples with pooled versions (now as mcmc objects)
    pooled_model$Sol <- pooled_sol
    pooled_model$VCV <- pooled_vcv
    pooled_model$DIC <- if (is.finite(mean_dic)) mean_dic else NULL
    
    if (has_cp) {
      pooled_model$CP <- pooled_cp
    }
    
    if (has_liab) {
      pooled_model$Liab <- pooled_liab
    }
    
    if (has_deviance) {
      pooled_model$Deviance <- pooled_deviance
    }
    
    # Add metadata about the pooling process
    pooled_model$BACE_pooling <- list(
      n_imputations = n_valid,
      n_samples_per_imputation = n_iter_per_model,
      total_samples = total_samples,
      variable = var,
      pooled = TRUE
    )
    
    # Set class to bace_pooled_MCMCglmm with MCMCglmm as second class
    # Order matters for S3 dispatch - custom methods are checked first
    class(pooled_model) <- c("bace_pooled_MCMCglmm", "MCMCglmm")
    
    pooled_models[[var]] <- pooled_model
  }
  
  # Prepare output
  out <- list(
    models = pooled_models,
    n_imputations = n_runs,
    variables = names(pooled_models),
    call = match.call()
  )
  
  class(out) <- "bace_pooled"
  
  return(out)
}


#' @title print.bace_pooled
#' @description Print method for bace_pooled objects
#' @param x Object of class bace_pooled
#' @param ... Additional arguments
#' @export
print.bace_pooled <- function(x, ...) {
  cat("\n=== BACE Pooled Posterior Results ===\n\n")
  cat("Number of imputations pooled:", x$n_imputations, "\n")
  cat("Number of variables:", length(x$models), "\n")
  cat("Variables:", paste(x$variables, collapse = ", "), "\n\n")
  
  cat("Each pooled model is an MCMCglmm object accounting for imputation uncertainty.\n")
  cat("Use standard MCMCglmm methods:\n")
  cat("  - summary(object$models$variable_name)\n")
  cat("  - print(object$models$variable_name)\n")
  cat("  - plot(object$models$variable_name)\n\n")
  
  for (var in x$variables) {
    model <- x$models[[var]]
    cat(paste0("\n--- ", var, " ---\n"))
    if (!is.null(model$BACE_pooling)) {
      cat("Total posterior samples:", model$BACE_pooling$total_samples, "\n")
      cat("Samples per imputation:", model$BACE_pooling$n_samples_per_imputation, "\n")
      cat("Number of imputations:", model$BACE_pooling$n_imputations, "\n")
    }
    cat("Number of fixed effects:", ncol(model$Sol), "\n")
    if (!is.null(model$DIC)) {
      cat("Mean DIC:", round(model$DIC, 2), "\n")
    }
  }
  
  cat("\n")
  invisible(x)
}


#' @title print.bace_pooled_MCMCglmm
#' @description Print method for individual pooled MCMCglmm model objects
#' @param x Object of class bace_pooled_MCMCglmm
#' @param ... Additional arguments passed to MCMCglmm print method
#' @export
print.bace_pooled_MCMCglmm <- function(x, ...) {
  # Add header indicating this is a pooled model
  if (!is.null(x$BACE_pooling)) {
    cat("\n+--------------------------------------------------------------+\n")
    cat("|  BACE Pooled MCMCglmm Model (Imputation Uncertainty Included) |\n")
    cat("+--------------------------------------------------------------+\n\n")
    cat("Pooled from", x$BACE_pooling$n_imputations, "imputations\n")
    cat("Total posterior samples:", x$BACE_pooling$total_samples, "\n")
    cat("  (=", x$BACE_pooling$n_samples_per_imputation, "samples/imputation x", 
        x$BACE_pooling$n_imputations, "imputations)\n\n")
    cat("Note: Posterior distribution accounts for both estimation and imputation uncertainty.\n")
    cat("--------------------------------------------------------------\n\n")
  }
  
  # Call the MCMCglmm print method
  # Temporarily remove bace_pooled_MCMCglmm class to call MCMCglmm method
  class(x) <- "MCMCglmm"
  print(x, ...)
  invisible(x)
}


#' @title summary.bace_pooled_MCMCglmm
#' @description Summary method for individual pooled MCMCglmm model objects
#' @param object Object of class bace_pooled_MCMCglmm
#' @param ... Additional arguments passed to MCMCglmm summary method
#' @export
summary.bace_pooled_MCMCglmm <- function(object, ...) {
  # Add header indicating this is a pooled model
  if (!is.null(object$BACE_pooling)) {
    cat("\n+--------------------------------------------------------------+\n")
    cat("|  BACE Pooled MCMCglmm Summary (Imputation Uncertainty Included)|\n")
    cat("+--------------------------------------------------------------+\n\n")
    cat("Pooled from", object$BACE_pooling$n_imputations, "imputations\n")
    cat("Total posterior samples:", object$BACE_pooling$total_samples, "\n")
    cat("  (=", object$BACE_pooling$n_samples_per_imputation, "samples/imputation x", 
        object$BACE_pooling$n_imputations, "imputations)\n\n")
    cat("Note: Posterior summaries account for both estimation and imputation uncertainty.\n")
    cat("--------------------------------------------------------------\n\n")
  }
  
  # Call the MCMCglmm summary method
  # Temporarily remove bace_pooled_MCMCglmm class to call MCMCglmm method
  class(object) <- "MCMCglmm"
  result <- summary(object, ...)
  
  # Return invisibly so it can be captured if needed
  invisible(result)
}
