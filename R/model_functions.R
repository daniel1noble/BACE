#' @title .model_fit
#' @description Function takes incomplete data and a formula string and fits a model using MCMCglmm
#' @param data A dataframe containing missing data.
#' @param tree A phylogenetic tree of class 'phylo' from the ape package.
#' @param fixformula A string that specifies the fixed effect structure in the model.
#' @param randformula A string that specifies the random effect structure in the model.
#' @param type A string that specifies the type of model to fit. Options are "normal", "binary", "count", "categorical", "ordinal" and the appropriate family will be used in MCMCglmm.
#' @param nitt An integer specifying the number of iterations to run the MCMC algorithm. Default is 6000
#' @param thin An integer specifying the thinning rate for the MCMC algorithm. Default is 5.
#' @param burnin An integer specifying the number of iterations to discard as burnin. Default is 1000.
#' @param prior A list specifying the prior distributions for the MCMCglmm model.
#' @return A list of draws from the posterior distribution of the model parameters.
#' @export

.model_fit <- function(data, tree, fixformula, randformula, type, prior, nitt = 6000, thin = 5, burnin = 1000) {

	# Create sparse matrix of phylogeny
		   A  <- MCMCglmm::inverseA(tree, nodes = "TIPS")$Ainv
	
	# Name of the column in the data corresponding to the phylogeny
		name  <- all.vars(randformula)
		
	# Fit the model using MCMCglmm
  if(type != "categorical"){
  	model <- MCMCglmm::MCMCglmm(fixed = fixformula,
                               random = randformula,
                                 data = data,
                               family = type,
							               ginverse = setNames(list(A), name),
                              verbose = FALSE, pr = TRUE, pl = TRUE,
                                saveX = TRUE, saveZ = TRUE,
                                 nitt = nitt,
                                 thin = thin,
                               burnin = burnin, 
							    prior = prior, singular.ok=TRUE)	
        } else {

      # Categorical model needs special treatment. Append -1 to the right size of ~ formula to remove intercept
      fixformula_cat <- as.formula(paste0(as.character(fixformula)[2], "~ -1 + ", as.character(fixformula)[3]))
      ranformula_cat <- as.formula(paste0("~", "idh(trait):",as.character(randformula)[2]))

      model <- MCMCglmm::MCMCglmm(fixed = fixformula_cat,
                                 random = ranformula_cat,
                                 data = data,
                               family = type,
                               rcov = ~us(trait):units,
							               ginverse = setNames(list(A), name),
                              verbose = FALSE, pr = TRUE, pl = TRUE,
                                saveX = TRUE, saveZ = TRUE,
                                 nitt = nitt,
                                 thin = thin,
                               burnin = burnin, 
							    prior = prior, singular.ok=TRUE)	     
                  }								
  	return(model)
}

#' @title .list_of_G
#' @description Function creates a list of G priors for the MCMCglmm model
#' @param n_rand An integer specifying the number of random effects in the model.
#' @param nu A numeric specifying the nu parameter for the prior.
#' @param par_expand A logical indicating whether to use parameter expansion.
#' @param diag An integer specifying the dimension of the variance-covariance matrix.
#' @return A list of G priors for the MCMCglmm model.

.list_of_G <- function(n_rand, nu = NULL, par_expand = FALSE, diag = 1) {
  
  if(is.null(nu)) {nu <- 0.002}
  
  prior_G <- list()

  if (par_expand) {
    for (i in 1:n_rand) {
      prior_G[[paste0("G", i)]] <- list(
        V = diag(diag), nu = nu,
        alpha.mu = rep(0, diag), alpha.V = 1e4 * diag(diag)
      )
    }
  } else {
    for (i in 1:n_rand) {
      prior_G[[paste0("G", i)]] <- list(V = diag(diag), nu = nu)
    }
  }
  
  return(prior_G)
}

#' @title .make_prior
#' @description Function creates the prior for the MCMCglmm model
#' @param n_rand An integer specifying the number of random effects in the model.
#' @param type A string that specifies the type of model to fit.
#' @param nu A numeric specifying the nu parameter for the prior.
#' @param n_levels An integer specifying the number of levels for categorical or ordinal data.
#' @param par_expand A logical indicating whether to use parameter expansion.
#' @return A list of priors for the MCMCglmm model.
#' @export

.make_prior <- function(n_rand, type, nu = NULL, n_levels = NULL, par_expand = FALSE) {
  if (type == "gaussian") {
    if (is.null(nu)) {
      nu <- 0.002
    }

    prior_G <- .list_of_G(n_rand, nu, par_expand)
    prior <- list(
      R = list(V = 1, nu = nu),
      G = prior_G
    )
  }

  if (type == "poisson") {
    if (is.null(nu)) {
      nu <- 0.02
    }
    prior_G <- .list_of_G(n_rand, nu, par_expand)
    prior <- list(
      R = list(V = 1, nu = nu),
      G = prior_G
    )
  }

  if (type  == "categorical") {

    stopifnot(!is.null(n_levels))

    J <- n_levels - 1

    J_matrix <- array(1, dim = c(J, J)) # matrix of ones
    I_matrix <- diag(J) # identity matrix

    IJ <- (I_matrix + J_matrix) / n_levels # see Hadfield's Course notes p. 97

    prior_G <- .list_of_G(n_rand, nu = nu, par_expand, diag = J)
    prior <- list(
      R = list(V = IJ, fix = 1),
      G = prior_G
    )
  }

  if (type == "threshold" || type == "ordinal") {
    # stopifnot(!is.null(n_levels)) # this check not needed for ordinal prior

    if (is.null(nu)) {
      nu <- 0.002
    }

    # Fix residual variance R at 1
    # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm
    # priors from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q2/018194.html
    # The residual variance has to be fixed, because the latent variable is without scale

    # J <- n_levels # number of categories not needed for ordinal prior

    prior_G <- .list_of_G(n_rand, nu = nu, par_expand)
    prior <- list(
      R = list(V = 1, fix = 1),
      G = prior_G
    )
  }

  return(prior)
}

#' @title .predict_bace
#' @description Function creates a predcition from MCMCglmm model
#' @param model A MCMCglmm model object
#' @param dat_prep A list containing the prepared data frame and attributes for continuous variables.
#' @param response_var A string specifying the name of the response variable.
#' @param type A string that specifies the type of model to fit.
#' @param sample A logical indicating whether to sample from the distribution for categorical/ordinal variables.
#' @param ... Additional arguments (not used).
#' @return A vector of predicted values from the MCMCglmm model.
#' @export
.predict_bace <- function(model, dat_prep, response_var, type = NULL, sample = FALSE, ...) {				

				if(type == "gaussian"){
					# z-transformed so need to back-transform
					mean_val <- dat_prep[[2]][[response_var]]$mean
					sd_val   <- dat_prep[[2]][[response_var]]$sd
					
					# Predict from model and back-transform
				 		pred_values <- .pred_cont(model)[,1] * sd_val + mean_val
			     }

				 if(type == "poisson"){  
					# Predict from model and round to nearest integer to retain count data
						pred_values <- round(.pred_count(model)[,1], digits = 0)
				 }

				 if(type == "threshold" || type == "ordinal"){
					# Identify number of categories and their levels from the data
			      	     lv  <- dat_prep[[1]][[response_var]]
				  levels_var <- sort(unique(as.character(lv)))
				  
				   # Predicts probabilities for each category
				 pred_prob <- .pred_threshold(model, level_names = levels_var)

				   # For each observation, sample from the categorical distribution based on the predicted probabilities. TO DO: Note we could also just take the max probability for baseline level
             pred_values <- .impute_levels(pred_prob, levels_var, sample = sample)
          }
				 
				 
				 if(type == "categorical"){
					# Identify number of categories and their levels from the data
			      	     lv  <- dat_prep[[1]][[response_var]]
				  levels_var <- sort(unique(as.character(lv)))
					
					# Predict category probabilities
					pred_prob <- .pred_cat(model, baseline_name = levels_var[1])

					# For each observation, sample from the categorical distribution based on the predicted probabilities
				   pred_values <- .impute_levels(pred_prob, levels_var, sample = sample)
				 }

	return(pred_values)
}

#' @title .impute_levels
#' @description Function samples predicted levels for categorical or ordinal variables based on predicted probabilities
#' @param pred_prob A data frame of predicted probabilities for each category
#' @param levels_var A character vector specifying the names of the levels/categories
#' @param sample A logical indicating whether to sample from the distribution or take the maximum probability
#' @return A vector of predicted levels for each observation
.impute_levels <- function(pred_prob, levels_var, sample = FALSE) {
  if (sample) {
    pred_values <- apply(pred_prob, 1, function(probs) {
      sample(levels_var, size = 1, prob = probs)
    })
  } else {
    pred_values <- levels_var[apply(pred_prob, 1, which.max)]
  }
  return(pred_values)
}

#' @title .pred_cat
#' @description Function calculates predicted probabilities for each category from a categorical MCMCglmm model
#' @param model A MCMCglmm model object
#' @param baseline_name A string specifying the name of the baseline category
#' @return A data frame of predicted probabilities for each category
#' @export
.pred_cat <- function(model, baseline_name = "Baseline") {
  
  # 1. Basic Dimensions
     n_obs <- model$Residual$nrl
  n_traits <- (ncol(model$Sol) - model$Fixed$nll) / (model$Fixed$nfl / (ncol(model$Sol) / (ncol(model$Liab)/n_obs)))
  
  # Simpler way to get n_traits for categorical:
  n_traits <- ncol(model$Liab) / n_obs
  
  # 2. Calculate the c2 Scaling Factor
  # c = 16 * sqrt(3) / (15 * pi)
        c2 <- (16 * sqrt(3) / (15 * pi))^2
  
  # Get Total Variance (G + R) for each sample in the chain
  # We sum the Phylo and Units variances. 
  # Note: MCMCglmm categorical models usually have fixed residual variance (e.g. 1)
  
  # since residual variance is fixed - we can use directly the IJ matrix
  IJ <- IJ <- 1 / (n_traits+1) * (diag(n_traits) + matrix(1, nrow = n_traits, ncol = n_traits))
  scaling_factor <- sqrt(1 + c2 * diag(IJ))
  
  # 3. Extract and Scale Liabilities
  exp_liab_list <- list()
  
  # Initialize exp_sum with 1 (which is exp(0) for the baseline)
  exp_sum <- 1 
  
  for (i in 1:n_traits) {
    cols <- ((i - 1) * n_obs + 1):(i * n_obs)
    
    # Apply c2 scaling to the liabilities before exponentiating
    # We divide each sample's liabilities by that sample's scaling factor
    scaled_liab <- model$Liab[, cols] / scaling_factor[i]
    
    exp_liab_list[[i]] <- exp(scaled_liab)
    exp_sum <- exp_sum + exp_liab_list[[i]]
  }
  
  # 4. Calculate Mean Probabilities (%)
  prop_results <- list()
  for (i in 1:n_traits) {
    prop_results[[i]] <- colMeans(exp_liab_list[[i]] / exp_sum)
  }
  
  # Calculate Baseline level %
  prop_results[[n_traits + 1]] <- colMeans(1 / exp_sum) 
  
  # 5. Extract Names Generically
  # Look at Fixed effects to get trait/variable names
  raw_names <- colnames(model$Sol)[1:n_traits]
  
  # Removes "trait", the variable name, and the following dot
  clean_names <- gsub("^trait.*?\\.", "", raw_names)
  
  # 6. Final Data Frame Assembly
  df_final <- as.data.frame(do.call(cbind, prop_results))
  colnames(df_final) <- c(clean_names, baseline_name)
  
  # Rearrange columns alphabetically (Aquatic, Insessorial, Terrestrial)
  df_ordered <- df_final[, order(colnames(df_final))]
  rownames(df_ordered) <- paste0("Obs_", 1:n_obs)
  
  return(df_ordered)
}


#' @title .pred_threshold
#' @description Function calculates predicted probabilities for each category from a threshold MCMCglmm model
#' @param model A MCMCglmm model object
#' @param level_names A character vector specifying the names of the levels/categories	
#' @return A data frame of predicted probabilities for each category
#' @export
.pred_threshold <- function(model, level_names = NULL) {
  # 1. Dimensions and Data
           n_obs <- model$Residual$nrl
            liab <- model$Liab
         v_total <- rowSums(model$VCV)
  scaling_factor <- 1
  
  # 2. Thresholds (Cut-points)
  # MCMCglmm fixes the first threshold at 0.
  # Additional cut-points are in model$CP.
  if (!is.null(model$CP)) {
    # Ordinal case: J > 2 levels
    cp_samples <- cbind(0, model$CP) # First CP is 0
      n_levels <- ncol(cp_samples) + 1
  } else {
    # Binary case: 2 levels
    cp_samples <- matrix(0, nrow = nrow(liab), ncol = 1)
      n_levels <- 2
  }
  
  # 3. Calculate Probabilities
  # We iterate through iterations and calculate area under normal curve
  # Prob(category j) = Phi((CP_j - Liab) / scale) - Phi((CP_j-1 - Liab) / scale)
  all_probs <- list()
  for (j in 1:n_levels) {
    # Define upper and lower bounds for the current category
    # Lower bound (T_low)
    if (j == 1) {
      t_low <- -Inf
    } else {
      t_low <- cp_samples[, j - 1]
    }
    # Upper bound (T_high)
    if (j == n_levels) {
      t_high <- Inf
    } else {
      t_high <- cp_samples[, j]
    }
    
	# Calculate probability for this category across all MCMC samples
    # Probability = pnorm(Upper) - pnorm(Lower)
    # We use mapply or row-wise math to handle the scaling per iteration
    # Prob per sample: pnorm(t_high, mean=liab, sd=scaling_factor) - pnorm(t_low, ...)
             p_cat <- pnorm(t_high, mean = liab, sd = scaling_factor) -
                      pnorm(t_low, mean = liab, sd = scaling_factor)
    all_probs[[j]] <- colMeans(p_cat) 
  }
  
  # 4. Final Data Frame
  df_final <- as.data.frame(do.call(cbind, all_probs))
  
  # 5. Naming and Ordering
  if (is.null(level_names)) {
    colnames(df_final) <- paste0("Level_", 1:n_levels)
  } else {
    colnames(df_final) <- level_names
  }
  return(df_final)
}


#' @title .pred_cont
#' @description Posterior mean and posterior SD of fitted values for a Gaussian (identity-link) MCMCglmm model
#' @param model A MCMCglmm model object
#' @return A data frame with two columns: post_mean and post_sd (row per observation used in the fit)
#' @export
.pred_cont <- function(model) {
  
  # 1. Need X and Sol (and Z if you want conditional fitted values)
  if (is.null(model$X))   stop("model$X is missing: fit with saveX=TRUE.")
  if (is.null(model$Sol)) stop("model$Sol is missing.")
  
  X   <- as.matrix(model$X)
  Sol <- as.matrix(model$Sol)
  
  # Use Z if present (pr=TRUE typically required to have RE columns in Sol)
  if (!is.null(model$Z)) {
    W <- cbind(X, as.matrix(model$Z))
  } else {
    W <- X
  }
  
  # 2. Align columns (most robust)
  if (!is.null(colnames(W)) && !is.null(colnames(Sol))) {
    common <- intersect(colnames(W), colnames(Sol))
    if (length(common) == 0) {
      stop("No matching coefficient names between design matrix (X/Z) and Sol.")
    }
    W   <- W[, common, drop = FALSE]
    Sol <- Sol[, common, drop = FALSE]
  } else {
    # Fallback: assume ordering matches for first min columns
      p <- min(ncol(W), ncol(Sol))
    W   <- W[, seq_len(p), drop = FALSE]
    Sol <- Sol[, seq_len(p), drop = FALSE]
  }
  
  # 4. eta draws: [n_iter x n_rows_in_model_matrix]
  eta <- Sol %*% t(W)
  
  out <- data.frame(
    post_mean = as.numeric(colMeans(eta)),
    post_sd   = as.numeric(apply(eta, 2, sd))
  )
  
  rownames(out) <- paste0("Obs_", seq_len(nrow(out)))
  return(out)
}


#' @title .pred_count
#' @description Posterior mean and posterior SD of fitted values for a Poisson (log-link) MCMCglmm model
#' @param model A MCMCglmm model object
#' @return A data frame with two columns: post_mean and post_sd (row per observation)
#' @export
.pred_count <- function(model) {
  
  if (is.null(model$Liab)) stop("model$Liab is missing.")
  liab <- as.matrix(model$Liab)
  
  mu <- exp(liab)
  
  out <- data.frame(
    post_mean = colMeans(mu),
    post_sd   = apply(mu, 2, sd)
  )
  
  rownames(out) <- paste0("Obs_", seq_len(nrow(out)))
  return(round(out, 4))
}


#' @title .check_mcmc_diagnostics
#' @description Function checks MCMC convergence and mixing for all parameters in MCMCglmm models
#' @param bace_output Output from bace_imp function containing models_last_run
#' @return A list containing diagnostic statistics and a summary data frame
#' @details This function performs the following diagnostics:
#' \itemize{
#'   \item Effective Sample Size (ESS) for all parameters
#'   \item Geweke convergence diagnostic
#'   \item Autocorrelation at lag 1
#'   \item Summary statistics for fixed and random effects
#' }
#' @export
.check_mcmc_diagnostics <- function(bace_output) {
  
  if (!inherits(bace_output, "bace")) {
    stop("Input must be a 'bace' object from bace_imp()")
  }
  
  if (is.null(bace_output$models_last_run)) {
    stop("No models found in bace_output. Make sure models were saved during imputation.")
  }
  
  models <- bace_output$models_last_run
  diagnostics <- list()
  
  for (var_name in names(models)) {
    model <- models[[var_name]]
    
    # Extract ONLY fixed effects and variance components (exclude BLUPs)
    # Fixed effects: First nfl columns of Sol (number of fixed effects levels)
    n_fixed <- model$Fixed$nfl
    fixed_effects <- model$Sol[, 1:n_fixed, drop = FALSE]
    
    # Random effect variance components
    random_variances <- model$VCV
    
    # Combine fixed effects and variance components only
    all_params <- cbind(fixed_effects, random_variances)
    param_names <- colnames(all_params)
    
    # Initialize results for this model
    var_diagnostics <- data.frame(
      parameter = param_names,
      mean = numeric(length(param_names)),
      sd = numeric(length(param_names)),
      ess = numeric(length(param_names)),
      geweke_z = numeric(length(param_names)),
      geweke_pval = numeric(length(param_names)),
      autocorr_lag1 = numeric(length(param_names)),
      convergence = character(length(param_names)),
      stringsAsFactors = FALSE
    )
    
    # Calculate diagnostics for each parameter
    for (i in seq_along(param_names)) {
      param <- all_params[, i]
      
      # Basic statistics
      var_diagnostics$mean[i] <- mean(param)
      var_diagnostics$sd[i] <- sd(param)
      
      # Effective Sample Size
      var_diagnostics$ess[i] <- coda::effectiveSize(param)
      
      # Geweke diagnostic (tests for equality of means in first 10% and last 50% of chain)
      geweke_result <- coda::geweke.diag(param)
      var_diagnostics$geweke_z[i] <- geweke_result$z
      var_diagnostics$geweke_pval[i] <- 2 * pnorm(-abs(geweke_result$z))
      
      # Autocorrelation at lag 1
      acf_result <- stats::acf(param, lag.max = 1, plot = FALSE)
      var_diagnostics$autocorr_lag1[i] <- acf_result$acf[2]
      
      # Convergence assessment
      # ESS > 200 is generally considered adequate
      # Geweke p-value > 0.05 suggests convergence
      # Low autocorrelation (< 0.4) suggests good mixing
      
      ess_ok <- var_diagnostics$ess[i] > 200
      geweke_ok <- var_diagnostics$geweke_pval[i] > 0.05
      autocorr_ok <- abs(var_diagnostics$autocorr_lag1[i]) < 0.4
      
      if (ess_ok && geweke_ok && autocorr_ok) {
        var_diagnostics$convergence[i] <- "Good"
      } else if (ess_ok && geweke_ok) {
        var_diagnostics$convergence[i] <- "Acceptable"
      } else {
        var_diagnostics$convergence[i] <- "Poor"
      }
    }
    
    # Identify fixed vs random variance components
    # Fixed effects are the first n_fixed parameters
    var_diagnostics$effect_type <- ifelse(
      seq_along(param_names) <= n_fixed,
      "Fixed",
      "Variance"
    )
    
    diagnostics[[var_name]] <- var_diagnostics
  }
  
  # Create overall summary
  all_diagnostics <- do.call(rbind, lapply(names(diagnostics), function(var) {
    df <- diagnostics[[var]]
    df$response_variable <- var
    df
  }))
  
  # Summary statistics
  summary_stats <- list(
    n_models = length(models),
    n_parameters_total = nrow(all_diagnostics),
    n_fixed_effects = sum(all_diagnostics$effect_type == "Fixed"),
    n_variance_components = sum(all_diagnostics$effect_type == "Variance"),
    convergence_summary = table(all_diagnostics$convergence),
    mean_ess_fixed = mean(all_diagnostics$ess[all_diagnostics$effect_type == "Fixed"]),
    mean_ess_variance = mean(all_diagnostics$ess[all_diagnostics$effect_type == "Variance"]),
    prop_good_convergence = mean(all_diagnostics$convergence == "Good"),
    prop_poor_convergence = mean(all_diagnostics$convergence == "Poor")
  )
  
  result <- list(
    diagnostics_by_model = diagnostics,
    all_diagnostics = all_diagnostics,
    summary = summary_stats
  )
  
  class(result) <- c("bace_diagnostics", "list")
  return(result)
}


#' @title print.bace_diagnostics
#' @description Print method for MCMC diagnostics from bace_imp models
#' @param x A bace_diagnostics object from .check_mcmc_diagnostics()
#' @param ... Additional arguments (not used)
#' @export
print.bace_diagnostics <- function(x, ...) {
  
  cat("\n========================================\n")
  cat("BACE MCMC Diagnostics Summary\n")
  cat("========================================\n\n")
  
  # Overall summary
  cat("Number of models:", x$summary$n_models, "\n")
  cat("Total parameters:", x$summary$n_parameters_total, "\n")
  cat("  Fixed effects:", x$summary$n_fixed_effects, "\n")
  cat("  Variance components:", x$summary$n_variance_components, "\n\n")
  
  # Effective Sample Size
  cat("Mean Effective Sample Size:\n")
  cat(sprintf("  Fixed effects:      %.1f\n", x$summary$mean_ess_fixed))
  cat(sprintf("  Variance components: %.1f\n", x$summary$mean_ess_variance))
  cat("\n")
  
  # Convergence summary
  cat("Convergence Assessment:\n")
  conv_table <- x$summary$convergence_summary
  for (status in names(conv_table)) {
    prop <- conv_table[status] / sum(conv_table) * 100
    cat(sprintf("  %-12s: %3d parameters (%.1f%%)\n", 
                status, conv_table[status], prop))
  }
  cat("\n")
  
  # Overall assessment
  if (x$summary$prop_good_convergence >= 0.8) {
    cat("Overall: GOOD - Most parameters show good convergence\n")
  } else if (x$summary$prop_poor_convergence <= 0.2) {
    cat("Overall: ACCEPTABLE - Most parameters converged adequately\n")
  } else {
    cat("Overall: POOR - Many parameters show convergence issues\n")
    cat("  Consider increasing nitt, adjusting priors, or checking model specification\n")
  }
  cat("\n")
  
  # Model-by-model summary
  cat("Convergence by Response Variable:\n")
  cat("----------------------------------------\n")
  
  for (var_name in names(x$diagnostics_by_model)) {
    diag <- x$diagnostics_by_model[[var_name]]
    
    n_good <- sum(diag$convergence == "Good")
    n_poor <- sum(diag$convergence == "Poor")
    n_total <- nrow(diag)
    
    status_symbol <- if (n_good / n_total >= 0.8) "[OK]" else if (n_poor / n_total > 0.3) "[X]" else "[~]"
    
    cat(sprintf("%s %-15s: %2d/%2d parameters good\n", 
                status_symbol, var_name, n_good, n_total))
    
    # Highlight poorly converged parameters
    poor_params <- diag$parameter[diag$convergence == "Poor"]
    if (length(poor_params) > 0) {
      cat(sprintf("    Poor convergence: %s\n", 
                  paste(poor_params, collapse=", ")))
    }
  }
  
  cat("========================================\n\n")
  
  invisible(x)
}


#' @title summary.bace_diagnostics
#' @description Detailed summary method for MCMC diagnostics
#' @param object A bace_diagnostics object from .check_mcmc_diagnostics()
#' @param ... Additional arguments (not used)
#' @export
summary.bace_diagnostics <- function(object, ...) {
  
  cat("\n========================================\n")
  cat("Detailed MCMC Diagnostics\n")
  cat("========================================\n\n")
  
  for (var_name in names(object$diagnostics_by_model)) {
    cat("\nResponse Variable:", var_name, "\n")
    cat("----------------------------------------\n")
    
    diag <- object$diagnostics_by_model[[var_name]]
    
    # Print table
    print(diag[, c("parameter", "effect_type", "mean", "sd", "ess", 
                   "geweke_pval", "autocorr_lag1", "convergence")], 
          row.names = FALSE, digits = 3)
    
    cat("\n")
  }
  
  invisible(object)
}
