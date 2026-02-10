#' @title plot.bace_convergence
#' @description Plot convergence diagnostics for BACE imputation
#' @param x Object of class bace_convergence from assess_convergence()
#' @importFrom graphics abline barplot boxplot grid layout legend lines mtext par plot.new points text
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom stats median na.pass
#' @param type Type of plot: "trace", "density", "acf", "pct_change", "energy", "all"
#' @param variables Character vector of variable names to plot. If NULL (default), plots all
#' @param ... Additional arguments passed to plotting functions
#' @return Plots convergence diagnostics
#' @examples \dontrun{
#' conv <- assess_convergence(bace_result)
#' plot(conv, type = "all")
#' plot(conv, type = "trace", variables = c("y", "x2"))
#' plot(conv, type = "pct_change")
#' }
#' @export
plot.bace_convergence <- function(x, 
                                   type = c("trace", "density", "acf", "pct_change", "energy", "all"),
                                   variables = NULL,
                                   ...) {
  
  type <- match.arg(type)
  
  if (type == "all") {
    # Plot each type sequentially (each sets its own layout internally)
    message("Plotting trace diagnostics...")
    plot_trace_convergence(x, variables = variables, ...)
    
    if (!is.null(x$summary_stats)) {
      readline(prompt = "Press [Enter] to see percentage change plots...")
      plot_pct_change_convergence(x, variables = variables, ...)
      
      readline(prompt = "Press [Enter] to see ACF plots...")
      plot_acf_convergence(x, variables = variables, ...)
    }
    
    if (!is.null(x$method_results$energy)) {
      readline(prompt = "Press [Enter] to see energy distance plot...")
      plot_energy_convergence(x, ...)
    }
    
  } else if (type == "trace") {
    plot_trace_convergence(x, variables = variables, ...)
  } else if (type == "density") {
    plot_density_convergence(x, variables = variables, ...)
  } else if (type == "pct_change") {
    plot_pct_change_convergence(x, variables = variables, ...)
  } else if (type == "acf") {
    plot_acf_convergence(x, variables = variables, ...)
  } else if (type == "energy") {
    plot_energy_convergence(x, ...)
  }
}


#' @title plot_trace_convergence
#' @description Plot trace plots of summary statistics across iterations
#' @param conv_object Object of class bace_convergence
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @export
plot_trace_convergence <- function(conv_object, variables = NULL, ...) {
  
  if (is.null(conv_object$summary_stats)) {
    stop("No summary statistics available. Run assess_convergence with method='summary' or 'all'")
  }
  
  stats_df <- conv_object$summary_stats
  
  # Select variables to plot
  var_cols <- setdiff(names(stats_df), "iteration")
  if (!is.null(variables)) {
    var_cols <- intersect(variables, var_cols)
  }
  
  if (length(var_cols) == 0) {
    stop("No valid variables to plot")
  }
  
  # Create trace plots
  n_vars <- length(var_cols)
  
  message(paste("Plotting", n_vars, "variables:", paste(var_cols, collapse = ", ")))
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Calculate proper layout
  n_rows <- ceiling(sqrt(n_vars))
  n_cols <- ceiling(n_vars / n_rows)
  
  message(paste("Layout:", n_rows, "rows x", n_cols, "cols"))
  
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  
  for (var in var_cols) {
    message(paste("Plotting variable:", var))
    plot(stats_df$iteration, stats_df[[var]], 
         type = "b", 
         pch = 19,
         col = "steelblue",
         xlab = "Iteration",
         ylab = "Summary Statistic",
         main = paste("Trace:", var),
         ...)
    
    # Add smoothed line
    if (nrow(stats_df) >= 3) {
      smooth_vals <- stats::lowess(stats_df$iteration, stats_df[[var]])
      lines(smooth_vals$x, smooth_vals$y, col = "red", lwd = 2)
    }
    
    # Add convergence indicator
    if (!is.null(conv_object$method_results$summary$var_converged)) {
      if (var %in% names(conv_object$method_results$summary$var_converged)) {
        converged <- conv_object$method_results$summary$var_converged[var]
        col_indicator <- ifelse(converged, "darkgreen", "darkred")
        
        # Get percentage change info if available
        pct_info <- ""
        if (!is.null(conv_object$method_results$summary$stationarity_tests[[var]]$mean_pct_change)) {
          pct_val <- conv_object$method_results$summary$stationarity_tests[[var]]$mean_pct_change
          if (is.finite(pct_val)) {
            pct_info <- sprintf(" (%.2f%% chg)", pct_val)
          }
        }
        
        status_text <- ifelse(converged, 
                             paste0("\u2713 Converged", pct_info), 
                             paste0("\u2717 Not Converged", pct_info))
        mtext(status_text, side = 3, line = 0.5, col = col_indicator, cex = 0.8)
      }
    }
    
    grid()
  }
}


#' @title plot_density_convergence
#' @description Plot density evolution across iterations
#' @param conv_object Object of class bace_convergence  
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @export
plot_density_convergence <- function(conv_object, variables = NULL, ...) {
  
  if (is.null(conv_object$summary_stats)) {
    stop("No summary statistics available")
  }
  
  stats_df <- conv_object$summary_stats
  var_cols <- setdiff(names(stats_df), "iteration")
  
  if (!is.null(variables)) {
    var_cols <- intersect(variables, var_cols)
  }
  
  n_vars <- length(var_cols)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(ceiling(sqrt(n_vars)), ceiling(n_vars / ceiling(sqrt(n_vars)))),
      mar = c(4, 4, 3, 1))
  
  n_iter <- nrow(stats_df)
  colors <- colorRampPalette(c("lightblue", "darkblue"))(n_iter)
  
  for (var in var_cols) {
    plot(NULL, 
         xlim = range(stats_df[[var]], na.rm = TRUE),
         ylim = c(0, n_iter + 1),
         xlab = "Value",
         ylab = "Iteration",
         main = paste("Density Evolution:", var),
         ...)
    
    # Plot density for each iteration
    for (i in 1:n_iter) {
      # Simple representation: vertical lines
      points(stats_df[[var]][i], i, pch = 19, col = colors[i], cex = 1.5)
    }
    
    # Add smoothed trend
    if (n_iter >= 3) {
      smooth_vals <- stats::lowess(stats_df[[var]], 1:n_iter)
      lines(smooth_vals$x, smooth_vals$y, col = "red", lwd = 2)
    }
    
    grid()
    
    # Add color legend
    legend("topright", 
           legend = c("Early", "Late"), 
           col = c(colors[1], colors[n_iter]),
           pch = 19,
           cex = 0.7)
  }
}


#' @title plot_pct_change_convergence
#' @description Plot percentage change between consecutive iterations
#' @param conv_object Object of class bace_convergence
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @export
plot_pct_change_convergence <- function(conv_object, variables = NULL, ...) {
  
  if (is.null(conv_object$summary_stats)) {
    stop("No summary statistics available. Run assess_convergence with method='summary' or 'all'")
  }
  
  stats_df <- conv_object$summary_stats
  
  # Select variables to plot
  var_cols <- setdiff(names(stats_df), "iteration")
  if (!is.null(variables)) {
    var_cols <- intersect(variables, var_cols)
  }
  
  if (length(var_cols) == 0) {
    stop("No valid variables to plot")
  }
  
  n_vars <- length(var_cols)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Calculate proper layout
  n_rows <- ceiling(sqrt(n_vars))
  n_cols <- ceiling(n_vars / n_rows)
  
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  
  for (var in var_cols) {
    var_series <- stats_df[[var]]
    n_iter <- length(var_series)
    
    # Calculate percentage changes
    if (n_iter >= 2) {
      pct_changes <- abs(diff(var_series) / (var_series[-n_iter] + 1e-10)) * 100
      
      # Plot percentage changes
      plot(2:n_iter, pct_changes,
           type = "b",
           pch = 19,
           col = "purple",
           xlab = "Iteration",
           ylab = "% Change from Previous",
           main = paste("% Change:", var),
           ylim = c(0, max(pct_changes, na.rm = TRUE) * 1.1),
           ...)
      
      # Add horizontal line at threshold if available
      if (!is.null(conv_object$method_results$summary$stationarity_tests[[var]])) {
        # Try to extract threshold from convergence test
        # Default to 5% if not available
        threshold <- 5  # default
        abline(h = threshold, col = "red", lty = 2, lwd = 2)
        text(n_iter * 0.7, threshold * 1.1, 
             paste0("Threshold: ", threshold, "%"), 
             col = "red", cex = 0.8)
      }
      
      # Add mean for last half
      last_half_idx <- ceiling(length(pct_changes) / 2)
      last_half_changes <- pct_changes[last_half_idx:length(pct_changes)]
      mean_change <- mean(last_half_changes, na.rm = TRUE)
      
      abline(h = mean_change, col = "blue", lty = 3, lwd = 1.5)
      
      # Add convergence indicator
      if (!is.null(conv_object$method_results$summary$stationarity_tests[[var]]$pct_change_converged)) {
        pct_converged <- conv_object$method_results$summary$stationarity_tests[[var]]$pct_change_converged
        mean_pct <- conv_object$method_results$summary$stationarity_tests[[var]]$mean_pct_change
        
        if (!is.na(pct_converged) && is.finite(mean_pct)) {
          col_indicator <- ifelse(pct_converged, "darkgreen", "darkred")
          status_text <- sprintf("%s Mean: %.2f%%", 
                                ifelse(pct_converged, "\u2713", "\u2717"),
                                mean_pct)
          mtext(status_text, side = 3, line = 0.5, col = col_indicator, cex = 0.8)
        }
      }
      
      grid()
    } else {
      plot.new()
      text(0.5, 0.5, "Insufficient iterations", cex = 1.2)
    }
  }
}


#' @title plot_acf_convergence
#' @description Plot autocorrelation diagnostics
#' @param conv_object Object of class bace_convergence
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @export
plot_acf_convergence <- function(conv_object, variables = NULL, ...) {
  
  if (is.null(conv_object$summary_stats)) {
    stop("No summary statistics available")
  }
  
  stats_df <- conv_object$summary_stats
  var_cols <- setdiff(names(stats_df), "iteration")
  
  if (!is.null(variables)) {
    var_cols <- intersect(variables, var_cols)
  }
  
  n_vars <- length(var_cols)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Calculate proper layout
  n_rows <- ceiling(sqrt(n_vars))
  n_cols <- ceiling(n_vars / n_rows)
  
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  
  for (var in var_cols) {
    if (length(stats_df[[var]]) > 3) {
      acf(stats_df[[var]], 
          main = paste("ACF:", var),
          na.action = na.pass,
          ...)
      
      # Add convergence indicator
      if (!is.null(conv_object$method_results$summary$stationarity_tests[[var]])) {
        acf_conv <- conv_object$method_results$summary$stationarity_tests[[var]]$acf_converged
        if (!is.na(acf_conv)) {
          col_indicator <- ifelse(acf_conv, "darkgreen", "darkred")
          mtext(ifelse(acf_conv, "\u2713 Low ACF", "\u2717 High ACF"), 
                side = 3, line = 0.5, col = col_indicator, cex = 0.8)
        }
      }
    }
  }
}


#' @title plot_energy_convergence
#' @description Plot energy distance evolution
#' @param conv_object Object of class bace_convergence
#' @param ... Additional plotting arguments
#' @export
plot_energy_convergence <- function(conv_object, ...) {
  
  if (is.null(conv_object$method_results$energy)) {
    stop("No energy distance results available. Run assess_convergence with method='energy' or 'all'")
  }
  
  energy_res <- conv_object$method_results$energy
  distances <- energy_res$energy_distances
  
  plot(1:length(distances), distances,
       type = "b",
       pch = 19,
       col = "darkgreen",
       xlab = "Iteration Pair",
       ylab = "Energy Distance",
       main = "Energy Distance Between Consecutive Iterations",
       ...)
  
  # Add smoothed line
  if (length(distances) >= 3) {
    smooth_vals <- stats::lowess(1:length(distances), distances)
    lines(smooth_vals$x, smooth_vals$y, col = "red", lwd = 2)
  }
  
  # Add horizontal line at median
  abline(h = median(distances), lty = 2, col = "gray50")
  
  # Add convergence status
  if (!is.na(energy_res$converged)) {
    col_indicator <- ifelse(energy_res$converged, "darkgreen", "darkred")
    status <- ifelse(energy_res$converged, 
                     "\u2713 Converged (stable energy distance)", 
                     "\u2717 Not Converged")
    mtext(status, side = 3, line = 0.5, col = col_indicator, cex = 0.9)
  }
  
  grid()
}


#' @title plot_wasserstein_convergence
#' @description Plot Wasserstein distance evolution for each variable
#' @param conv_object Object of class bace_convergence
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @export
plot_wasserstein_convergence <- function(conv_object, variables = NULL, ...) {
  
  if (is.null(conv_object$method_results$wasserstein)) {
    stop("No Wasserstein distance results available")
  }
  
  wass_res <- conv_object$method_results$wasserstein$wasserstein_results
  
  if (length(wass_res) == 0) {
    stop("No variables with Wasserstein distances available")
  }
  
  # Filter variables
  var_names <- names(wass_res)
  if (!is.null(variables)) {
    var_names <- intersect(variables, var_names)
  }
  
  n_vars <- length(var_names)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(ceiling(sqrt(n_vars)), ceiling(n_vars / ceiling(sqrt(n_vars)))),
      mar = c(4, 4, 3, 1))
  
  for (var in var_names) {
    distances <- wass_res[[var]]$distances
    stabilized <- wass_res[[var]]$stabilized
    
    plot(1:length(distances), distances,
         type = "b",
         pch = 19,
         col = "purple",
         xlab = "Iteration Pair",
         ylab = "Wasserstein Distance",
         main = paste("Wasserstein Distance:", var),
         ...)
    
    # Add smoothed line
    if (length(distances) >= 3) {
      smooth_vals <- stats::lowess(1:length(distances), distances)
      lines(smooth_vals$x, smooth_vals$y, col = "red", lwd = 2)
    }
    
    # Add convergence indicator
    if (!is.na(stabilized)) {
      col_indicator <- ifelse(stabilized, "darkgreen", "darkred")
      mtext(ifelse(stabilized, "\u2713 Stabilized", "\u2717 Not Stabilized"), 
            side = 3, line = 0.5, col = col_indicator, cex = 0.8)
    }
    
    grid()
  }
}


#' @title plot_bace_imputation_comparison
#' @description Compare imputed values across iterations
#' @param bace_object Object of class 'bace' from bace_imp
#' @param variable Variable name to visualize
#' @param n_iter Number of iterations to display. Default is all
#' @param ... Additional plotting arguments
#' @export
plot_bace_imputation_comparison <- function(bace_object, 
                                            variable,
                                            n_iter = NULL,
                                            ...) {
  
  if (!inherits(bace_object, "bace")) {
    stop("Input must be an object of class 'bace'")
  }
  
  imputed_data <- bace_object$data
  miss_dat <- bace_object$miss_dat
  types <- bace_object$types
  
  # Check variable exists
  if (!variable %in% names(types)) {
    stop(paste("Variable", variable, "not found"))
  }
  
  # Get iterations to plot
  total_iter <- length(imputed_data)
  if (is.null(n_iter) || n_iter > total_iter) {
    n_iter <- total_iter
  }
  
  # Select evenly spaced iterations
  iter_indices <- unique(round(seq(1, total_iter, length.out = n_iter)))
  
  # Get missing indices for this variable
  miss_idx <- miss_dat[miss_dat$colname == variable, "row"]
  
  if (length(miss_idx) == 0) {
    stop(paste("No missing data for variable", variable))
  }
  
  var_type <- types[[variable]]
  
  if (var_type %in% c("gaussian", "poisson")) {
    # Continuous variable: violin/boxplot
    plot_data <- data.frame(
      value = numeric(),
      iteration = character(),
      stringsAsFactors = FALSE
    )
    
    for (idx in iter_indices) {
      iter_name <- names(imputed_data)[idx]
      values <- imputed_data[[idx]][[variable]][miss_idx]
      plot_data <- rbind(plot_data, 
                         data.frame(value = values, 
                                   iteration = iter_name))
    }
    
    # Create boxplot
    boxplot(value ~ iteration, data = plot_data,
            col = colorRampPalette(c("lightblue", "darkblue"))(length(iter_indices)),
            main = paste("Imputed Values Across Iterations:", variable),
            xlab = "Iteration",
            ylab = "Value",
            las = 2,
            ...)
    
  } else {
    # Categorical variable: bar plot of proportions
    prop_data <- matrix(0, nrow = 0, ncol = 0)
    
    for (idx in iter_indices) {
      iter_name <- names(imputed_data)[idx]
      values <- imputed_data[[idx]][[variable]][miss_idx]
      
      # Get proportions
      tab <- table(values)
      props <- prop.table(tab)
      
      if (nrow(prop_data) == 0) {
        prop_data <- matrix(props, nrow = 1)
        colnames(prop_data) <- names(props)
        rownames(prop_data) <- iter_name
      } else {
        # Ensure same categories
        new_row <- numeric(ncol(prop_data))
        names(new_row) <- colnames(prop_data)
        for (cat in names(props)) {
          if (cat %in% colnames(prop_data)) {
            new_row[cat] <- props[cat]
          }
        }
        prop_data <- rbind(prop_data, new_row)
        rownames(prop_data)[nrow(prop_data)] <- iter_name
      }
    }
    
    # Create stacked bar plot
    barplot(t(prop_data),
            col = rainbow(ncol(prop_data)),
            main = paste("Category Proportions Across Iterations:", variable),
            xlab = "Iteration",
            ylab = "Proportion",
            legend = TRUE,
            las = 2,
            ...)
  }
}


#' @title plot_convergence_summary
#' @description Create a comprehensive summary plot of convergence diagnostics
#' @param conv_object Object of class bace_convergence
#' @param ... Additional plotting arguments
#' @export
plot_convergence_summary <- function(conv_object, ...) {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Create a 2x2 layout
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
         widths = c(1, 1),
         heights = c(1, 1))
  
  # Panel 1: Overall convergence status
  par(mar = c(2, 2, 3, 2))
  plot(1, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
       axes = FALSE, xlab = "", ylab = "", main = "Convergence Status")
  
  status_text <- ifelse(conv_object$converged, 
                       "CONVERGED", 
                       "NOT CONVERGED")
  status_col <- ifelse(conv_object$converged, "darkgreen", "darkred")
  
  text(0.5, 0.6, status_text, cex = 2.5, col = status_col, font = 2)
  
  # Show method votes if available
  if (!is.null(conv_object$diagnostics$votes)) {
    votes <- conv_object$diagnostics$votes
    n_pass <- sum(votes, na.rm = TRUE)
    text(0.5, 0.35, 
         paste(n_pass, "out of", length(votes), "methods agree"),
         cex = 1.2)
  }
  
  # Panel 2: Variable-level convergence
  if (!is.null(conv_object$method_results$summary$var_converged)) {
    par(mar = c(4, 8, 3, 2))
    
    var_conv <- conv_object$method_results$summary$var_converged
    barplot(as.numeric(var_conv),
            horiz = TRUE,
            names.arg = names(var_conv),
            col = ifelse(var_conv, "darkgreen", "darkred"),
            xlim = c(0, 1),
            main = "Variable Convergence",
            xlab = "Convergence Score",
            las = 1)
    abline(v = 0.5, lty = 2, col = "gray50")
  }
  
  # Panel 3: Energy distance (if available)
  if (!is.null(conv_object$method_results$energy$energy_distances)) {
    par(mar = c(4, 4, 3, 2))
    
    distances <- conv_object$method_results$energy$energy_distances
    plot(1:length(distances), distances,
         type = "b",
         pch = 19,
         col = "steelblue",
         main = "Energy Distance",
         xlab = "Iteration",
         ylab = "Distance")
    grid()
  }
  
  # Panel 4: Summary statistics trends
  if (!is.null(conv_object$summary_stats)) {
    par(mar = c(4, 4, 3, 2))
    
    stats_df <- conv_object$summary_stats
    var_cols <- setdiff(names(stats_df), "iteration")
    
    # Normalize and plot all variables
    plot(NULL,
         xlim = range(stats_df$iteration),
         ylim = c(0, 1),
         main = "Normalized Summary Statistics",
         xlab = "Iteration",
         ylab = "Normalized Value")
    
    colors <- rainbow(length(var_cols))
    for (i in seq_along(var_cols)) {
      var <- var_cols[i]
      vals <- stats_df[[var]]
      # Normalize to 0-1
      vals_norm <- (vals - min(vals, na.rm = TRUE)) / 
                   (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE))
      lines(stats_df$iteration, vals_norm, col = colors[i], lwd = 2)
    }
    
    legend("topright", legend = var_cols, col = colors, lwd = 2, cex = 0.6)
    grid()
  }
}
