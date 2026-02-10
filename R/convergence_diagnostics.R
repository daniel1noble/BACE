#' @title assess_convergence
#' @description Assess convergence of BACE imputation using multiple diagnostic criteria
#' @param bace_object An object of class 'bace' returned from bace_imp function
#' @importFrom stats cor cor.test var approx
#' @param method Character string specifying the convergence method: "summary" (all summary tests), 
#'   "summary.acf" (ACF only), "summary.percentage" (% change only), "summary.trend" (trend only),
#'   "summary.geweke" (Geweke only), "energy", "wasserstein", or "all" for all methods
#' @param variables Character vector of variable names to assess. If NULL (default), 
#'   all variables with missing data are assessed. A warning is issued for variables 
#'   without missing data, which are ignored.
#' @param use_all_data Logical. If TRUE, use all data from each iteration for convergence assessment.
#'   If FALSE (default), use only the imputed values (originally missing data).
#' @param alpha Significance level for convergence decision. Default is 0.05
#' @param min_iterations Minimum number of iterations before testing convergence. Default is 3
#' @param lag Maximum lag for autocorrelation assessment. Default is 1
#' @param pct_change_threshold Threshold for percentage change between consecutive iterations. 
#'   Variables are considered converged if the mean absolute percentage change in the 
#'   last half of iterations is below this threshold. Default is 0.05 (5%).
#' @return A list containing:
#'   - converged: Logical indicating overall convergence status
#'   - method_results: Detailed results for each method
#'   - summary_stats: Data frame with summary statistics across iterations
#'   - diagnostics: Additional diagnostic information
#' @examples \dontrun{
#' # After running bace_imp
#' result <- bace_imp(fixformula = "y ~ x1 + x2", ...)
#' conv <- assess_convergence(result, method = "all")
#' print(conv$converged)
#' }
#' @export
assess_convergence <- function(bace_object, 
                               method = c("summary", "summary.acf", "summary.percentage", 
                                         "summary.trend", "summary.geweke",
                                         "energy", "wasserstein", "all"),
                               variables = NULL,
                               use_all_data = FALSE,
                               alpha = 0.05,
                               min_iterations = 3,
                               lag = 1,
                               pct_change_threshold = 0.05) {
  
  # Check inputs
  if (!inherits(bace_object, "bace")) {
    stop("Input must be an object of class 'bace' from bace_imp function")
  }
  
  
  method <- match.arg(method)
  
  # Extract imputed datasets
  imputed_data <- bace_object$data[-1]
  n_iter <- length(imputed_data)
  
  if (n_iter < min_iterations) {
    warning(paste("Only", n_iter, "iterations available. Minimum of", 
                  min_iterations, "recommended for convergence assessment"))
    return(list(converged = NA, 
                message = "Insufficient iterations for convergence assessment"))
  }
  
  # Get variable types and missing data info
  types <- bace_object$types
  miss_dat <- bace_object$miss_dat
  
  # Filter variables if specified
  vars_with_missing <- unique(miss_dat$colname)
  
  if (!is.null(variables)) {
    # Check which requested variables don't have missing data
    vars_without_missing <- setdiff(variables, vars_with_missing)
    
    if (length(vars_without_missing) > 0) {
      warning(paste("The following variables have no missing data and will be ignored:",
                    paste(vars_without_missing, collapse = ", ")))
    }
    
    # Filter to only variables with missing data
    variables <- intersect(variables, vars_with_missing)
    
    if (length(variables) == 0) {
      stop("None of the specified variables have missing data")
    }
  }
  
  # Initialize results
  results <- list(
    converged = FALSE,
    method_results = list(),
    summary_stats = NULL,
    diagnostics = list()
  )
  
  # Parse method to extract summary criterion if specified
  summary_criterion <- NULL
  if (grepl("^summary\\.", method)) {
    parts <- strsplit(method, "\\.")[[1]]
    summary_criterion <- parts[2]  # Extract criterion (acf, percentage, trend, geweke)
    base_method <- "summary"
  } else {
    base_method <- method
  }
  
  # Method 1: Summary statistics convergence
  if (base_method %in% c("summary", "all")) {
    summary_result <- .assess_summary_convergence(imputed_data, types, miss_dat, 
                                                   variables = variables,
                                                   use_all_data = use_all_data,
                                                   alpha = alpha, lag = lag,
                                                   pct_change_threshold = pct_change_threshold,
                                                   criterion = summary_criterion)
    results$method_results$summary <- summary_result
    results$summary_stats <- summary_result$stats_df
  }
  
  # Method 2: Energy distance
  if (base_method %in% c("energy", "all")) {
    energy_result <- .assess_energy_convergence(imputed_data, miss_dat, 
                                                 variables = variables,
                                                 use_all_data = use_all_data,
                                                 alpha = alpha)
    results$method_results$energy <- energy_result
  }
  
  # Method 3: Wasserstein distance
  if (base_method %in% c("wasserstein", "all")) {
    wasserstein_result <- .assess_wasserstein_convergence(imputed_data, miss_dat, 
                                                           types, 
                                                           variables = variables,
                                                           use_all_data = use_all_data,
                                                           alpha = alpha)
    results$method_results$wasserstein <- wasserstein_result
  }
  
  # Overall convergence decision
  if (base_method == "all") {
    # Converged if at least 2 out of 3 methods indicate convergence
    convergence_votes <- c(
      summary = results$method_results$summary$converged,
      energy = results$method_results$energy$converged,
      wasserstein = results$method_results$wasserstein$converged
    )
    results$converged <- sum(convergence_votes, na.rm = TRUE) >= 2
    results$diagnostics$votes <- convergence_votes
  } else if (base_method == "summary") {
    results$converged <- results$method_results$summary$converged
  } else if (base_method == "energy") {
    results$converged <- results$method_results$energy$converged
  } else if (base_method == "wasserstein") {
    results$converged <- results$method_results$wasserstein$converged
  }
  
  class(results) <- c("bace_convergence", "list")
  return(results)
}


#' @title .assess_summary_convergence
#' @description Internal function to assess convergence using summary statistics
#' @param imputed_data List of imputed datasets across iterations
#' @param types Named list of variable types
#' @param miss_dat Data frame with missing data information
#' @param variables Character vector of variables to check, or NULL for all
#' @param use_all_data Logical. Use all data (TRUE) or only imputed values (FALSE)
#' @param alpha Significance level
#' @param lag Maximum lag for autocorrelation
#' @param pct_change_threshold Threshold for mean absolute percentage change (default 0.05)
#' @param criterion Which specific criterion to use (NULL for all, "acf", "percentage", "trend", or "geweke")
#' @return List with convergence assessment
.assess_summary_convergence <- function(imputed_data, types, miss_dat, 
                                        variables = NULL,
                                        use_all_data = FALSE,
                                        alpha = 0.05, lag = 1,
                                        pct_change_threshold = 0.05,
                                        criterion = NULL) {
  
  n_iter <- length(imputed_data)
  
  # Get unique variables with missing data
  vars_to_check <- unique(miss_dat$colname)
  
  # Filter to specified variables if provided
  if (!is.null(variables)) {
    vars_to_check <- intersect(vars_to_check, variables)
  }
  
  # Initialize storage
  stats_list <- list()
  
  # Calculate summary statistics for each variable across iterations
  for (var in vars_to_check) {
    var_type <- types[[var]]
    var_stats <- numeric(n_iter)
    
    # Get indices with missing data for this variable
    miss_idx <- miss_dat[miss_dat$colname == var, "row"]
    
    for (i in 1:n_iter) {
      # Use all data or only imputed values
      if (use_all_data) {
        imputed_values <- imputed_data[[i]][[var]]
      } else {
        imputed_values <- imputed_data[[i]][[var]][miss_idx]
      }
      
      if (var_type == "gaussian" || var_type == "poisson") {
        # Use mean for continuous/count variables
        var_stats[i] <- mean(imputed_values, na.rm = TRUE)
      } else if (var_type %in% c("categorical", "threshold")) {
        # Use proportion of first category (alphabetically) for categorical
        if (is.factor(imputed_values)) {
          first_level <- levels(imputed_values)[1]
        } else {
          first_level <- sort(unique(imputed_values))[1]
        }
        var_stats[i] <- mean(imputed_values == first_level, na.rm = TRUE)
      }
    }
    
    stats_list[[var]] <- var_stats
  }
  
  # Create data frame with all statistics
  stats_matrix <- do.call(cbind, stats_list)
  stats_df <- data.frame(
    iteration = 1:n_iter,
    stats_matrix,
    check.names = FALSE
  )
  # Ensure column names are preserved from stats_list
  names(stats_df)[-1] <- names(stats_list)
  
  # Test for stationarity using multiple criteria
  stationarity_tests <- list()
  
  for (var in vars_to_check) {
    var_series <- stats_list[[var]]
    
    # 1. Geweke diagnostic (compare first 10% vs last 50%)
    geweke_result <- .geweke_test(var_series, alpha = alpha)
    
    # 2. Check autocorrelation at specified lag
    if (n_iter > lag + 1 && sum(is.finite(var_series)) >= lag + 2) {
      acf_val <- tryCatch({
        cor(var_series[1:(n_iter - lag)], 
            var_series[(1 + lag):n_iter], 
            use = "complete.obs")
      }, error = function(e) {
        NA
      })
      # Low autocorrelation indicates mixing/convergence
      acf_converged <- if (is.finite(acf_val)) abs(acf_val) < 0.3 else NA
    } else {
      acf_val <- NA
      acf_converged <- NA
    }
    
    # 3. Check if trend is stable (no significant linear trend)
    if (n_iter >= 3 && sum(is.finite(var_series)) >= 3) {
      trend_test <- tryCatch({
        cor.test(1:n_iter, var_series)
      }, error = function(e) {
        list(p.value = NA)
      })
      trend_converged <- if (!is.na(trend_test$p.value)) trend_test$p.value > alpha else NA
    } else {
      trend_converged <- NA
    }
    
    # 4. Check percentage change between consecutive iterations
    if (n_iter >= 3 && sum(is.finite(var_series)) >= 3) {
      # Calculate absolute percentage changes between consecutive iterations
      pct_changes <- abs(diff(var_series) / (var_series[-n_iter] + 1e-10)) * 100
      
      # Focus on last half of iterations
      last_half_idx <- ceiling(length(pct_changes) / 2)
      last_half_changes <- pct_changes[last_half_idx:length(pct_changes)]
      
      # Mean absolute percentage change
      mean_pct_change <- mean(last_half_changes, na.rm = TRUE)
      
      # Converged if mean change is below threshold
      pct_change_converged <- if (is.finite(mean_pct_change)) {
        mean_pct_change < (pct_change_threshold * 100)
      } else {
        NA
      }
    } else {
      mean_pct_change <- NA
      pct_change_converged <- NA
    }
    
    # TODO - Consider additional tests like KPSS for stationarity, but may be too strict for short series
    # Temporary turn-off Geweke stats
    geweke_result$converged <- NA
    geweke_result$z_score <- NA

    stationarity_tests[[var]] <- list(
      geweke_converged = geweke_result$converged,
      geweke_z = geweke_result$z_score,
      acf = acf_val,
      acf_converged = acf_converged,
      trend_converged = trend_converged,
      mean_pct_change = mean_pct_change,
      pct_change_converged = pct_change_converged
    )
  }
  
  # Overall convergence: all variables should show convergence
  # Build convergence matrix based on selected criterion
  if (!is.null(criterion)) {
    # Use only the specified criterion
    convergence_matrix <- switch(criterion,
      "acf" = sapply(stationarity_tests, function(x) x$acf_converged),
      "percentage" = sapply(stationarity_tests, function(x) x$pct_change_converged),
      "trend" = sapply(stationarity_tests, function(x) x$trend_converged),
      "geweke" = sapply(stationarity_tests, function(x) x$geweke_converged),
      stop(paste("Unknown criterion:", criterion))
    )
    # For single criterion, convergence_matrix is a vector
    var_converged <- convergence_matrix
    names(var_converged) <- names(stationarity_tests)
  } else {
    # Use all criteria
    convergence_matrix <- sapply(stationarity_tests, function(x) {
      c(x$geweke_converged, x$acf_converged, x$trend_converged, x$pct_change_converged)
    })
    # Converged if majority of tests pass for each variable
    var_converged <- colMeans(convergence_matrix, na.rm = TRUE) > 0.5
  }
  
  # Overall: converged if majority of variables converged
  overall_converged <- mean(var_converged, na.rm = TRUE) > 0.7
  
  return(list(
    converged = overall_converged,
    stats_df = stats_df,
    stationarity_tests = stationarity_tests,
    var_converged = var_converged
  ))
}


#' @title .geweke_test
#' @description Geweke convergence diagnostic comparing early vs late chain segments
#' @param x Numeric vector of summary statistics across iterations
#' @param alpha Significance level
#' @param frac1 Fraction of chain for first window. Default 0.1
#' @param frac2 Fraction of chain for second window. Default 0.5
#' @return List with convergence status and z-score
.geweke_test <- function(x, alpha = 0.05, frac1 = 0.1, frac2 = 0.5) {
  n <- length(x)
  
  if (n < 10 || sum(is.finite(x)) < 10) {
    return(list(converged = NA, z_score = NA, p_value = NA))
  }
  
  # First and last segments
  n1 <- floor(frac1 * n)
  n2_start <- ceiling((1 - frac2) * n)
  
  seg1 <- x[1:n1]
  seg2 <- x[n2_start:n]
  
  # Check for sufficient finite observations
  if (sum(is.finite(seg1)) < 2 || sum(is.finite(seg2)) < 2) {
    return(list(converged = NA, z_score = NA, p_value = NA))
  }
  
  # Calculate means and variances
  mean1 <- mean(seg1, na.rm = TRUE)
  mean2 <- mean(seg2, na.rm = TRUE)
  var1 <- var(seg1, na.rm = TRUE) / sum(is.finite(seg1))
  var2 <- var(seg2, na.rm = TRUE) / sum(is.finite(seg2))
  
  # Check for valid variance
  if (!is.finite(var1) || !is.finite(var2) || (var1 + var2) <= 0) {
    return(list(converged = NA, z_score = NA, p_value = NA))
  }
  
  # Z-score
  z_score <- (mean1 - mean2) / sqrt(var1 + var2)
  
  # Check if z_score is finite
  if (!is.finite(z_score)) {
    return(list(converged = NA, z_score = NA, p_value = NA))
  }
  
  # Two-tailed test
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  converged <- p_value > alpha
  
  return(list(
    converged = converged,
    z_score = z_score,
    p_value = p_value
  ))
}


#' @title .assess_energy_convergence
#' @description Assess convergence using energy distance between consecutive iterations
#' @param imputed_data List of imputed datasets
#' @param miss_dat Missing data information
#' @param variables Character vector of variables to check, or NULL for all
#' @param use_all_data Logical. Use all data (TRUE) or only imputed values (FALSE)
#' @param alpha Significance level
#' @return List with convergence assessment
.assess_energy_convergence <- function(imputed_data, miss_dat, 
                                       variables = NULL,
                                       use_all_data = FALSE, alpha = 0.05) {
  
  n_iter <- length(imputed_data)
  
  if (n_iter < 2) {
    return(list(converged = NA, message = "Need at least 2 iterations"))
  }
  
  # Calculate energy distance between consecutive iterations
  energy_distances <- numeric(n_iter - 1)
  
  for (i in 1:(n_iter - 1)) {
    # Get imputed values or all data
    data1_imputed <- .extract_imputed_values(imputed_data[[i]], miss_dat, variables, use_all_data)
    data2_imputed <- .extract_imputed_values(imputed_data[[i + 1]], miss_dat, variables, use_all_data)
    
    # Calculate energy distance
    energy_distances[i] <- .energy_distance(data1_imputed, data2_imputed)
  }
  
  # Test for decreasing trend and stabilization
  # Converged if:
  # 1. No significant increasing trend
  # 2. Low variation in last iterations
  
  if (n_iter >= 4 && sum(is.finite(energy_distances)) >= 3) {
    # Test trend in energy distances
    trend_test <- tryCatch({
      cor.test(1:length(energy_distances), energy_distances)
    }, error = function(e) {
      list(estimate = 0, p.value = 1)
    })
    no_increase <- trend_test$estimate <= 0 || trend_test$p.value > alpha
    
    # Check if stabilized (low CV in last half)
    last_half <- energy_distances[ceiling(length(energy_distances) / 2):length(energy_distances)]
    if (sum(is.finite(last_half)) >= 2 && mean(last_half, na.rm = TRUE) != 0) {
      cv <- sd(last_half, na.rm = TRUE) / mean(last_half, na.rm = TRUE)
      stabilized <- cv < 0.5  # Coefficient of variation < 50%
    } else {
      stabilized <- NA
    }
    
    converged <- if (!is.na(no_increase) && !is.na(stabilized)) no_increase && stabilized else NA
  } else {
    converged <- NA
    no_increase <- NA
    stabilized <- NA
  }
  
  return(list(
    converged = converged,
    energy_distances = energy_distances,
    no_increase = no_increase,
    stabilized = stabilized
  ))
}


#' @title .energy_distance
#' @description Calculate energy distance between two multivariate samples
#' @param x1 Matrix or data frame (first sample)
#' @param x2 Matrix or data frame (second sample)
#' @return Numeric energy distance
.energy_distance <- function(x1, x2) {
  
  # Convert to matrix if needed
  if (is.data.frame(x1)) x1 <- as.matrix(x1)
  if (is.data.frame(x2)) x2 <- as.matrix(x2)
  
  # Remove rows with any NA
  x1 <- x1[complete.cases(x1), , drop = FALSE]
  x2 <- x2[complete.cases(x2), , drop = FALSE]
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  # Check if we have enough data
  if (n1 == 0 || n2 == 0) {
    return(NA_real_)
  }
  
  # Subsample if large (for computational efficiency)
  max_sample <- 100
  if (n1 > max_sample) {
    idx1 <- sample(n1, max_sample)
    x1 <- x1[idx1, , drop = FALSE]
    n1 <- max_sample
  }
  if (n2 > max_sample) {
    idx2 <- sample(n2, max_sample)
    x2 <- x2[idx2, , drop = FALSE]
    n2 <- max_sample
  }
  
  # Calculate pairwise distances
  # E[||X-Y||] - 0.5*E[||X-X'||] - 0.5*E[||Y-Y'||]
  
  # Between samples
  dist_between <- 0
  for (i in 1:n1) {
    for (j in 1:n2) {
      dist_between <- dist_between + sqrt(sum((x1[i, ] - x2[j, ])^2, na.rm = TRUE))
    }
  }
  dist_between <- dist_between / (n1 * n2)
  
  # Within sample 1
  dist_within1 <- 0
  n_pairs1 <- 0
  if (n1 > 1) {
    for (i in 1:(n1 - 1)) {
      for (j in (i + 1):n1) {
        dist_within1 <- dist_within1 + sqrt(sum((x1[i, ] - x1[j, ])^2, na.rm = TRUE))
        n_pairs1 <- n_pairs1 + 1
      }
    }
    dist_within1 <- dist_within1 / n_pairs1
  }
  
  # Within sample 2
  dist_within2 <- 0
  n_pairs2 <- 0
  if (n2 > 1) {
    for (i in 1:(n2 - 1)) {
      for (j in (i + 1):n2) {
        dist_within2 <- dist_within2 + sqrt(sum((x2[i, ] - x2[j, ])^2, na.rm = TRUE))
        n_pairs2 <- n_pairs2 + 1
      }
    }
    dist_within2 <- dist_within2 / n_pairs2
  }
  
  energy_dist <- 2 * dist_between - dist_within1 - dist_within2
  
  # Check for valid result
  if (!is.finite(energy_dist)) {
    return(NA_real_)
  }
  
  return(energy_dist)
}


#' @title .assess_wasserstein_convergence
#' @description Assess convergence using Wasserstein distance for each variable
#' @param imputed_data List of imputed datasets
#' @param miss_dat Missing data information
#' @param types Variable types
#' @param variables Character vector of variables to check, or NULL for all
#' @param use_all_data Logical. Use all data (TRUE) or only imputed values (FALSE)
#' @param alpha Significance level
#' @return List with convergence assessment
.assess_wasserstein_convergence <- function(imputed_data, miss_dat, types, 
                                            variables = NULL,
                                            use_all_data = FALSE, alpha = 0.05) {
  
  n_iter <- length(imputed_data)
  vars_to_check <- unique(miss_dat$colname)
  
  # Filter to specified variables if provided
  if (!is.null(variables)) {
    vars_to_check <- intersect(vars_to_check, variables)
  }
  
  if (n_iter < 2) {
    return(list(converged = NA, message = "Need at least 2 iterations"))
  }
  
  # Calculate Wasserstein distance for each variable
  wasserstein_results <- list()
  
  for (var in vars_to_check) {
    var_type <- types[[var]]
    miss_idx <- miss_dat[miss_dat$colname == var, "row"]
    
    # Only for continuous variables
    if (var_type %in% c("gaussian", "poisson")) {
      distances <- numeric(n_iter - 1)
      
      for (i in 1:(n_iter - 1)) {
        # Use all data or only imputed values
        if (use_all_data) {
          vals1 <- imputed_data[[i]][[var]]
          vals2 <- imputed_data[[i + 1]][[var]]
        } else {
          vals1 <- imputed_data[[i]][[var]][miss_idx]
          vals2 <- imputed_data[[i + 1]][[var]][miss_idx]
        }
        
        # 1D Wasserstein distance (simple implementation)
        distances[i] <- .wasserstein_1d(vals1, vals2)
      }
      
      # Test for stabilization
      if (length(distances) >= 3) {
        # Use last half of iterations to assess stability
        last_half_idx <- ceiling(length(distances) / 2)
        last_half <- distances[last_half_idx:length(distances)]
        
        # Check multiple criteria for stabilization:
        # 1. Low coefficient of variation (CV < 0.3) in later iterations
        cv <- sd(last_half) / mean(last_half)
        low_variation <- !is.na(cv) && cv < 0.3
        
        # 2. No significant increasing trend in last half
        if (length(last_half) >= 3) {
          trend_test <- tryCatch({
            cor.test(seq_along(last_half), last_half)
          }, error = function(e) list(estimate = 0, p.value = 1))
          no_increase <- trend_test$estimate <= 0 || trend_test$p.value > alpha
        } else {
          no_increase <- TRUE
        }
        
        # 3. Later iterations are not substantially higher than earlier ones
        first_quarter <- distances[1:max(1, floor(length(distances) / 4))]
        not_increasing <- mean(last_half) <= mean(first_quarter) * 1.5
        
        # Stabilized if low variation AND (no increase OR not increasing overall)
        stabilized <- low_variation && (no_increase || not_increasing)
      } else {
        stabilized <- NA
      }
      
      wasserstein_results[[var]] <- list(
        distances = distances,
        stabilized = stabilized
      )
    }
  }
  
  # Overall convergence: most variables should be stabilized
  if (length(wasserstein_results) > 0) {
    stabilized_vars <- sapply(wasserstein_results, function(x) x$stabilized)
    converged <- mean(stabilized_vars, na.rm = TRUE) > 0.7
  } else {
    converged <- NA
  }
  
  return(list(
    converged = converged,
    wasserstein_results = wasserstein_results
  ))
}


#' @title .wasserstein_1d
#' @description Calculate 1-dimensional Wasserstein distance
#' @param x First sample
#' @param y Second sample
#' @return Numeric Wasserstein distance
.wasserstein_1d <- function(x, y) {
  x <- sort(x)
  y <- sort(y)
  
  n <- length(x)
  m <- length(y)
  
  # Make equal length by linear interpolation
  if (n != m) {
    if (n < m) {
      x <- approx(1:n, x, n = m)$y
    } else {
      y <- approx(1:m, y, n = n)$y
    }
  }
  
  # L1 distance between sorted values
  mean(abs(x - y))
}


#' @title .extract_imputed_values
#' @description Extract imputed values or all values from a dataset
#' @param data Data frame with imputed values
#' @param miss_dat Missing data information
#' @param variables Character vector of variables to extract, or NULL for all
#' @param use_all_data Logical. If TRUE, extract all data; if FALSE, extract only imputed values
#' @return Matrix of values
.extract_imputed_values <- function(data, miss_dat, variables = NULL, use_all_data = FALSE) {
  
  vars <- unique(miss_dat$colname)
  
  # Filter to specified variables if provided
  if (!is.null(variables)) {
    vars <- intersect(vars, variables)
  }
  
  imputed_list <- list()
  
  for (var in vars) {
    if (use_all_data) {
      values <- data[[var]]
    } else {
      miss_idx <- miss_dat[miss_dat$colname == var, "row"]
      values <- data[[var]][miss_idx]
    }
    
    # Convert factors to numeric for distance calculation
    if (is.factor(values)) {
      values <- as.numeric(values)
    }
    
    imputed_list[[var]] <- values
  }
  
  # Create matrix with equal length (pad with NA if needed)
  max_len <- max(sapply(imputed_list, length))
  imputed_matrix <- sapply(imputed_list, function(x) {
    c(x, rep(NA, max_len - length(x)))
  })
  
  # Remove rows with all NA
  imputed_matrix <- imputed_matrix[rowSums(is.na(imputed_matrix)) < ncol(imputed_matrix), , drop = FALSE]
  
  return(imputed_matrix)
}


#' @title print.bace_convergence
#' @description Print method for bace_convergence objects
#' @param x Object of class bace_convergence
#' @param ... Additional arguments
#' @export
print.bace_convergence <- function(x, ...) {
  cat("\n=== BACE Imputation Convergence Assessment ===\n\n")
  
  cat("Overall Convergence:", ifelse(x$converged, "YES", "NO"), "\n\n")
  
  if (!is.null(x$method_results$summary)) {
    cat("--- Summary Statistics Method ---\n")
    cat("Converged:", x$method_results$summary$converged, "\n")
    cat("Variables converged:", 
        sum(x$method_results$summary$var_converged, na.rm = TRUE), "/", 
        length(x$method_results$summary$var_converged), "\n\n")
  }
  
  if (!is.null(x$method_results$energy)) {
    cat("--- Energy Distance Method ---\n")
    cat("Converged:", x$method_results$energy$converged, "\n\n")
  }
  
  if (!is.null(x$method_results$wasserstein)) {
    cat("--- Wasserstein Distance Method ---\n")
    cat("Converged:", x$method_results$wasserstein$converged, "\n\n")
  }
  
  if (!is.null(x$diagnostics$votes)) {
    cat("--- Method Votes ---\n")
    print(x$diagnostics$votes)
  }
  
  cat("\nUse plot() to visualize convergence diagnostics\n")
}
