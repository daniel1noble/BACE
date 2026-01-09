#' @title assess_convergence
#' @description Assess convergence of BACE imputation using multiple diagnostic criteria
#' @param bace_object An object of class 'bace' returned from bace_imp function
#' @param method Character string specifying the convergence method: "summary" (default), 
#'   "energy", "wasserstein", or "all" for all methods
#' @param use_all_data Logical. If TRUE, use all data from each iteration for convergence assessment.
#'   If FALSE (default), use only the imputed values (originally missing data).
#' @param alpha Significance level for convergence decision. Default is 0.05
#' @param min_iterations Minimum number of iterations before testing convergence. Default is 3
#' @param lag Maximum lag for autocorrelation assessment. Default is 1
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
                               method = c("summary", "energy", "wasserstein", "all"),
                               use_all_data = FALSE,
                               alpha = 0.05,
                               min_iterations = 3,
                               lag = 1) {
  
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
  
  # Initialize results
  results <- list(
    converged = FALSE,
    method_results = list(),
    summary_stats = NULL,
    diagnostics = list()
  )
  
  # Method 1: Summary statistics convergence
  if (method %in% c("summary", "all")) {
    summary_result <- .assess_summary_convergence(imputed_data, types, miss_dat, 
                                                   use_all_data = use_all_data,
                                                   alpha = alpha, lag = lag)
    results$method_results$summary <- summary_result
    results$summary_stats <- summary_result$stats_df
  }
  
  # Method 2: Energy distance
  if (method %in% c("energy", "all")) {
    energy_result <- .assess_energy_convergence(imputed_data, miss_dat, 
                                                 use_all_data = use_all_data,
                                                 alpha = alpha)
    results$method_results$energy <- energy_result
  }
  
  # Method 3: Wasserstein distance
  if (method %in% c("wasserstein", "all")) {
    wasserstein_result <- .assess_wasserstein_convergence(imputed_data, miss_dat, 
                                                           types, 
                                                           use_all_data = use_all_data,
                                                           alpha = alpha)
    results$method_results$wasserstein <- wasserstein_result
  }
  
  # Overall convergence decision
  if (method == "all") {
    # Converged if at least 2 out of 3 methods indicate convergence
    convergence_votes <- c(
      summary = results$method_results$summary$converged,
      energy = results$method_results$energy$converged,
      wasserstein = results$method_results$wasserstein$converged
    )
    results$converged <- sum(convergence_votes, na.rm = TRUE) >= 2
    results$diagnostics$votes <- convergence_votes
  } else if (method == "summary") {
    results$converged <- results$method_results$summary$converged
  } else if (method == "energy") {
    results$converged <- results$method_results$energy$converged
  } else if (method == "wasserstein") {
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
#' @param use_all_data Logical. Use all data (TRUE) or only imputed values (FALSE)
#' @param alpha Significance level
#' @param lag Maximum lag for autocorrelation
#' @return List with convergence assessment
.assess_summary_convergence <- function(imputed_data, types, miss_dat, 
                                        use_all_data = FALSE,
                                        alpha = 0.05, lag = 1) {
  
  n_iter <- length(imputed_data)
  
  # Get unique variables with missing data
  vars_to_check <- unique(miss_dat$colname)
  
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
  stats_df <- data.frame(
    iteration = 1:n_iter,
    do.call(cbind, stats_list)
  )
  
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
    
    stationarity_tests[[var]] <- list(
      geweke_converged = geweke_result$converged,
      geweke_z = geweke_result$z_score,
      acf = acf_val,
      acf_converged = acf_converged,
      trend_converged = trend_converged
    )
  }
  
  # Overall convergence: all variables should show convergence
  convergence_matrix <- sapply(stationarity_tests, function(x) {
    c(x$geweke_converged, x$acf_converged, x$trend_converged)
  })
  
  # Converged if majority of tests pass for majority of variables
  var_converged <- colMeans(convergence_matrix, na.rm = TRUE) > 0.5
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
#' @param use_all_data Logical. Use all data (TRUE) or only imputed values (FALSE)
#' @param alpha Significance level
#' @return List with convergence assessment
.assess_energy_convergence <- function(imputed_data, miss_dat, 
                                       use_all_data = FALSE, alpha = 0.05) {
  
  n_iter <- length(imputed_data)
  
  if (n_iter < 2) {
    return(list(converged = NA, message = "Need at least 2 iterations"))
  }
  
  # Calculate energy distance between consecutive iterations
  energy_distances <- numeric(n_iter - 1)
  
  for (i in 1:(n_iter - 1)) {
    # Get imputed values or all data
    data1_imputed <- .extract_imputed_values(imputed_data[[i]], miss_dat, use_all_data)
    data2_imputed <- .extract_imputed_values(imputed_data[[i + 1]], miss_dat, use_all_data)
    
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
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  # Calculate pairwise distances
  # E[||X-Y||] - 0.5*E[||X-X'||] - 0.5*E[||Y-Y'||]
  
  # Between samples
  dist_between <- 0
  for (i in 1:min(n1, 100)) {  # Subsample if large
    for (j in 1:min(n2, 100)) {
      dist_between <- dist_between + sqrt(sum((x1[i, ] - x2[j, ])^2))
    }
  }
  dist_between <- dist_between / (min(n1, 100) * min(n2, 100))
  
  # Within sample 1
  dist_within1 <- 0
  if (n1 > 1) {
    for (i in 1:min(n1 - 1, 100)) {
      for (j in (i + 1):min(n1, 100)) {
        dist_within1 <- dist_within1 + sqrt(sum((x1[i, ] - x1[j, ])^2))
      }
    }
    n_pairs1 <- min(n1 - 1, 100) * (min(n1, 100) - min(n1 - 1, 100)) / 2
    dist_within1 <- dist_within1 / n_pairs1
  }
  
  # Within sample 2
  dist_within2 <- 0
  if (n2 > 1) {
    for (i in 1:min(n2 - 1, 100)) {
      for (j in (i + 1):min(n2, 100)) {
        dist_within2 <- dist_within2 + sqrt(sum((x2[i, ] - x2[j, ])^2))
      }
    }
    n_pairs2 <- min(n2 - 1, 100) * (min(n2, 100) - min(n2 - 1, 100)) / 2
    dist_within2 <- dist_within2 / n_pairs2
  }
  
  energy_dist <- 2 * dist_between - dist_within1 - dist_within2
  
  return(energy_dist)
}


#' @title .assess_wasserstein_convergence
#' @description Assess convergence using Wasserstein distance for each variable
#' @param imputed_data List of imputed datasets
#' @param miss_dat Missing data information
#' @param types Variable types
#' @param use_all_data Logical. Use all data (TRUE) or only imputed values (FALSE)
#' @param alpha Significance level
#' @return List with convergence assessment
.assess_wasserstein_convergence <- function(imputed_data, miss_dat, types, 
                                            use_all_data = FALSE, alpha = 0.05) {
  
  n_iter <- length(imputed_data)
  vars_to_check <- unique(miss_dat$colname)
  
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
        # Check if distances are decreasing and small in later iterations
        last_half <- distances[ceiling(length(distances) / 2):length(distances)]
        stabilized <- mean(last_half) < median(distances) * 0.5
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
#' @param use_all_data Logical. If TRUE, extract all data; if FALSE, extract only imputed values
#' @return Matrix of values
.extract_imputed_values <- function(data, miss_dat, use_all_data = FALSE) {
  
  vars <- unique(miss_dat$colname)
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
