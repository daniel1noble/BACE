#' @title pool_posteriors
#' @description Combine per-imputation MCMC chains into a single set of posterior
#'   draws that approximate the marginal posterior integrating over the missing-data
#'   distribution.
#' @param bace_final_object An object of class 'bace_final' from bace_final_imp
#' @param variable Character string specifying which variable's model to pool.
#'   If NULL (default), pools all variables
#' @param sample_size Integer specifying how many posterior samples to draw from each
#'   imputation before pooling. If NULL (default), uses all posterior samples. Setting
#'   this to a smaller value (e.g., 1000) can greatly reduce the memory footprint of
#'   the pooled object. The total number of draws in the pooled posterior will be
#'   sample_size * n_imputations.
#' @return A list of class 'bace_pooled' containing pooled models for each variable.
#'   Each pooled model is an MCMCglmm object with class c("bace_pooled_MCMCglmm",
#'   "MCMCglmm") containing:
#'   - Sol: Stacked fixed-effect draws (combined across imputations)
#'   - VCV: Stacked variance-component draws
#'   - CP: Stacked cutpoint draws (for categorical/ordinal models, if present)
#'   - Liab: Stacked latent-variable draws (for categorical models, if present)
#'   - DIC: Mean DIC across imputations
#'   - Deviance: Stacked deviance samples
#'   - BACE_pooling: Metadata about the pooling (n_imputations,
#'     n_samples_per_imputation, etc.)
#'   - All other MCMCglmm components (Fixed, Random, etc.) carried over unchanged
#'     from the first imputation's fit.
#' @details
#'   \strong{What this function does.} It concatenates the MCMC draws from each
#'   per-imputation run produced by `bace_final_imp()` into a single mcmc object
#'   per parameter. It does \emph{not} apply Rubin's rules: no between- vs.
#'   within-imputation variance decomposition is computed, and no scalar point
#'   estimates are combined.
#'
#'   \strong{Why stacking is the right combiner for Bayesian MI.} Under a
#'   fully Bayesian imputation model, each per-imputation chain targets the
#'   conditional posterior \eqn{p(\theta \mid Y_{obs}, Y_{mis}^{(m)})}, where
#'   \eqn{Y_{mis}^{(m)}} is the m-th draw from the imputation distribution
#'   \eqn{p(Y_{mis} \mid Y_{obs})}. Concatenating these chains and treating the
#'   result as draws from a single distribution is a Monte Carlo approximation
#'   to the marginal posterior
#'   \deqn{p(\theta \mid Y_{obs}) = \int p(\theta \mid Y_{obs}, Y_{mis})\, p(Y_{mis} \mid Y_{obs})\, dY_{mis}.}
#'   Posterior summaries (means, quantiles, intervals) computed from the
#'   stacked draws therefore reflect uncertainty about both the model parameters
#'   and the imputed values in a single, internally consistent way. This is the
#'   standard Bayesian-MI combiner; see Zhou & Reiter (2010, The American
#'   Statistician 64:159-163) and Gelman et al. (2013, Bayesian Data Analysis,
#'   3rd ed., Ch. 18). It differs from Rubin's rules, which operate on scalar
#'   estimates and variances and are the standard combiner when per-imputation
#'   analyses produce frequentist point estimates rather than full posterior
#'   draws.
#'
#'   \strong{Assumptions you should check.} Validity of the stacked posterior
#'   rests on two conditions that `pool_posteriors()` does not verify:
#'   \enumerate{
#'     \item Each per-imputation MCMC chain has converged to and adequately
#'       mixed within its own conditional posterior. If individual MCMCglmm
#'       runs have not converged, the pooled draws inherit that bias. Run
#'       standard MCMCglmm diagnostics (trace plots, effective sample size,
#'       Geweke) on a handful of the per-imputation fits in
#'       `bace_final_object$all_models` before trusting pooled summaries.
#'     \item The chained-equations loop has reached its own stationary
#'       distribution over imputed values — i.e. the sequence of imputed
#'       datasets is no longer drifting. Use `assess_convergence()` on the
#'       initial `bace_imp()` output for this.
#'   }
#'
#'   \strong{Class and methods.} Pooled models inherit both
#'   "bace_pooled_MCMCglmm" and "MCMCglmm" classes, so all standard MCMCglmm
#'   methods (`summary`, `plot`, posterior extraction) continue to work on the
#'   concatenated chain.
#' @references
#'   Zhou, X. and Reiter, J.P. (2010) A note on Bayesian inference after
#'   multiple imputation. \emph{The American Statistician} 64(2):159-163.
#'
#'   Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A. and
#'   Rubin, D.B. (2013) \emph{Bayesian Data Analysis}, 3rd edition, Chapter 18.
#'   Chapman & Hall/CRC.
#' @importFrom coda as.mcmc
#' @examples \dontrun{
#' # After running bace_final_imp
#' final <- bace_final_imp(bace_obj, ...)
#'
#' # Stack all draws from every imputation (may create large objects)
#' pooled <- pool_posteriors(final)
#'
#' # Stack with per-imputation subsampling to reduce memory footprint
#' pooled <- pool_posteriors(final, sample_size = 1000)
#'
#' # Extract the pooled model for a specific variable — standard MCMCglmm
#' # methods operate on the concatenated chain.
#' pooled_y <- pooled$models$y
#' summary(pooled_y)
#' plot(pooled_y)
#'
#' # Before trusting pooled summaries, sanity-check a few of the
#' # per-imputation fits for MCMC convergence:
#' plot(final$all_models[[1]]$y)
#' coda::effectiveSize(final$all_models[[1]]$y$Sol)
#' }
#' @export
pool_posteriors <- function(bace_final_object, variable = NULL, sample_size = NULL) {
  
  # Check inputs
  if (!inherits(bace_final_object, "bace_final")) {
    stop("Input must be an object of class 'bace_final' from bace_final_imp function")
  }
  
  # Validate sample_size
  if (!is.null(sample_size)) {
    if (!is.numeric(sample_size) || length(sample_size) != 1) {
      stop("sample_size must be a single numeric value or NULL")
    }
    if (sample_size < 1) {
      warning("sample_size must be >= 1. Using all samples instead.")
      sample_size <- NULL
    } else if (sample_size != floor(sample_size)) {
      sample_size <- floor(sample_size)
      warning("sample_size must be an integer. Rounding down to ", sample_size)
    }
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
    
    # Determine sampling strategy
    use_sampling <- !is.null(sample_size) && sample_size < n_iter_per_model
    if (use_sampling) {
      n_samples_per_imp <- sample_size
    } else {
      n_samples_per_imp <- n_iter_per_model
    }
    
    # Calculate total number of posterior samples
    total_samples <- n_samples_per_imp * n_valid
    
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
    
    # Combine posteriors by stacking (with optional sampling)
    for (i in 1:n_valid) {
      start_idx <- (i - 1) * n_samples_per_imp + 1
      end_idx <- i * n_samples_per_imp
      
      # Sample indices if using sampling, otherwise use all
      if (use_sampling) {
        sample_idx <- sample(1:n_iter_per_model, n_samples_per_imp, replace = FALSE)
      } else {
        sample_idx <- 1:n_iter_per_model
      }
      
      pooled_sol[start_idx:end_idx, ] <- var_models[[i]]$Sol[sample_idx, , drop = FALSE]
      pooled_vcv[start_idx:end_idx, ] <- var_models[[i]]$VCV[sample_idx, , drop = FALSE]
      
      if (has_cp) {
        pooled_cp[start_idx:end_idx, ] <- var_models[[i]]$CP[sample_idx, , drop = FALSE]
      }
      
      if (has_liab) {
        pooled_liab[start_idx:end_idx, ] <- var_models[[i]]$Liab[sample_idx, , drop = FALSE]
      }
      
      if (has_deviance) {
        pooled_deviance[start_idx:end_idx] <- var_models[[i]]$Deviance[sample_idx]
      }
    }
    
    # Calculate mean DIC across imputations
    dics <- sapply(var_models, function(m) {
      if (!is.null(m$DIC)) m$DIC else NA
    })
    mean_dic <- mean(dics, na.rm = TRUE)
    
    # Convert pooled matrices to mcmc objects (required by MCMCglmm methods).
    # Note: the concatenated samples are NOT a single MCMC chain in the
    # stationarity sense -- they are Monte Carlo draws from the marginal
    # posterior p(theta | Y_obs) obtained by stacking draws from each
    # conditional posterior p(theta | Y_obs, Y_mis^(m)). Wrapping them in an
    # mcmc object is purely so that downstream MCMCglmm methods (summary,
    # plot, HPDinterval, etc.) accept the result; per-chain convergence
    # diagnostics (Geweke, autocorrelation) applied to the stacked object
    # are not meaningful and should be run on the per-imputation chains in
    # bace_final_object$all_models instead.
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
      n_samples_per_imputation = n_samples_per_imp,
      original_samples_per_imputation = n_iter_per_model,
      total_samples = total_samples,
      variable = var,
      pooled = TRUE,
      sampled = use_sampling
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
  cat("Number of imputations stacked:", x$n_imputations, "\n")
  cat("Number of variables:", length(x$models), "\n")
  cat("Variables:", paste(x$variables, collapse = ", "), "\n\n")

  cat("Each pooled model is an MCMCglmm-compatible object whose draws\n")
  cat("approximate the marginal posterior integrating over the\n")
  cat("imputation distribution (stacked per-imputation chains; see\n")
  cat("?pool_posteriors for details and references).\n\n")
  cat("Use standard MCMCglmm methods:\n")
  cat("  - summary(object$models$variable_name)\n")
  cat("  - print(object$models$variable_name)\n")
  cat("  - plot(object$models$variable_name)\n\n")
  
  for (var in x$variables) {
    model <- x$models[[var]]
    cat(paste0("\n--- ", var, " ---\n"))
    if (!is.null(model$BACE_pooling)) {
      cat("Total posterior samples:", model$BACE_pooling$total_samples, "\n")
      if (model$BACE_pooling$sampled) {
        cat("Samples per imputation:", model$BACE_pooling$n_samples_per_imputation, 
            "(sampled from", model$BACE_pooling$original_samples_per_imputation, ")\n")
      } else {
        cat("Samples per imputation:", model$BACE_pooling$n_samples_per_imputation, "\n")
      }
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
    cat("\n+----------------------------------------------------------------+\n")
    cat("|  BACE Pooled MCMCglmm (stacked per-imputation posterior draws)  |\n")
    cat("+----------------------------------------------------------------+\n\n")
    cat("Stacked from", x$BACE_pooling$n_imputations, "imputations\n")
    cat("Total posterior draws:", x$BACE_pooling$total_samples, "\n")
    if (x$BACE_pooling$sampled) {
      cat("  (=", x$BACE_pooling$n_samples_per_imputation, "sampled per imputation x",
          x$BACE_pooling$n_imputations, "imputations)\n")
      cat("  Original samples per imputation:", x$BACE_pooling$original_samples_per_imputation, "\n")
    } else {
      cat("  (=", x$BACE_pooling$n_samples_per_imputation, "samples/imputation x",
          x$BACE_pooling$n_imputations, "imputations)\n")
    }
    cat("\nDraws approximate p(theta | Y_obs) by Monte Carlo integration over\n")
    cat("the imputation distribution. This is NOT Rubin's rules -- see\n")
    cat("?pool_posteriors for the assumptions this combiner relies on.\n")
    cat("----------------------------------------------------------------\n\n")
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
    cat("\n+----------------------------------------------------------------+\n")
    cat("|  BACE Pooled MCMCglmm summary (stacked per-imputation draws)    |\n")
    cat("+----------------------------------------------------------------+\n\n")
    cat("Stacked from", object$BACE_pooling$n_imputations, "imputations\n")
    cat("Total posterior draws:", object$BACE_pooling$total_samples, "\n")
    if (object$BACE_pooling$sampled) {
      cat("  (=", object$BACE_pooling$n_samples_per_imputation, "sampled per imputation x",
          object$BACE_pooling$n_imputations, "imputations)\n")
      cat("  Original samples per imputation:", object$BACE_pooling$original_samples_per_imputation, "\n")
    } else {
      cat("  (=", object$BACE_pooling$n_samples_per_imputation, "samples/imputation x",
          object$BACE_pooling$n_imputations, "imputations)\n")
    }
    cat("\nPosterior summaries below are computed from stacked draws that\n")
    cat("approximate p(theta | Y_obs) (Monte Carlo integration over imputed\n")
    cat("values). Validity requires each per-imputation chain to have\n")
    cat("converged on its own; check a subset via\n")
    cat("plot(final$all_models[[1]]$<var>) and coda::effectiveSize().\n")
    cat("----------------------------------------------------------------\n\n")
  }
  
  # Call the MCMCglmm summary method
  # Temporarily remove bace_pooled_MCMCglmm class to call MCMCglmm method
  class(object) <- "MCMCglmm"
  result <- summary(object, ...)
  
  # Return invisibly so it can be captured if needed
  print(result)
}
