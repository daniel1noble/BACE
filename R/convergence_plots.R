#' @title .conv_plot_vars
#' @description Internal: resolve which variables to draw, with the same
#'   message everywhere, and set up the panel grid. Returns NULL-safe vector.
#' @noRd
.conv_plot_vars <- function(stats_df, variables) {
  var_cols <- setdiff(names(stats_df), "iteration")
  if (!is.null(variables)) var_cols <- intersect(variables, var_cols)
  if (length(var_cols) == 0)
    stop("No valid variables to plot (none of the specified variables were assessed)")
  var_cols
}

#' @title .conv_finite_range
#' @description Internal: a finite plotting range even when a monitored
#'   series is constant or all non-finite (both are states the stability
#'   check deliberately tolerates, so the plots must not die on them).
#' @noRd
.conv_finite_range <- function(y) {
  y <- y[is.finite(y)]
  if (length(y) == 0) return(c(-1, 1))
  r <- range(y)
  if (diff(r) < 1e-12) r <- r + c(-1, 1) * max(abs(r), 1) * 0.05
  r
}

#' @title plot.bace_convergence
#' @description Plot the stability diagnostics for a BACE convergence chain.
#'   The trace panels are the primary diagnostic — the flag from
#'   \code{assess_convergence} is a pathology detector, not a certificate,
#'   so the traces deserve an eyeball (van Buuren 2018, ch. 6.5.2).
#' @param x Object of class bace_convergence from assess_convergence()
#' @importFrom graphics abline barplot boxplot grid layout legend lines mtext par plot.new points text rect title
#' @importFrom grDevices colorRampPalette rainbow adjustcolor
#' @importFrom stats median na.pass
#' @param type Type of plot: "trace" (the monitored mean per chained run, with
#'   each variable's verdict), "steps" (run-to-run movement in units of the
#'   variable's imputation spread, against its Monte-Carlo noise floor),
#'   "density", "acf", or "all". The pre-1.0.1 types "pct_change" and "energy"
#'   are defunct: the quantities they drew belonged to the defective pre-1.0.1
#'   criteria (see NEWS).
#' @param variables Character vector of variable names to plot. If NULL (default), plots all
#' @param ... Additional arguments passed to plotting functions
#' @return No return value (\code{NULL}, invisibly); called for its side
#'   effect of drawing the requested convergence diagnostic plot(s).
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- sim_bace(response_type = "gaussian", predictor_types = "gaussian",
#'                 phylo_signal = c(0.6, 0.6), missingness = c(0.2, 0),
#'                 n_cases = 40, n_species = 40)
#' fit <- bace_imp(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
#'                 phylo = sim$tree, data = sim$data,
#'                 runs = 3, nitt = 600, burnin = 100, thin = 2, verbose = FALSE)
#' conv <- assess_convergence(fit)
#' plot(conv, type = "trace")
#' }
#' @export
plot.bace_convergence <- function(x,
                                   type = c("trace", "steps", "density", "acf", "all"),
                                   variables = NULL,
                                   ...) {

  if (length(type) == 1 && type %in% c("pct_change", "energy")) {
    stop("plot type \"", type, "\" is defunct: it drew the pre-1.0.1 criteria, ",
         "which measured each variable's scale rather than its stability. ",
         "Use type = \"steps\" for the repaired movement panel. See NEWS.")
  }
  type <- match.arg(type)

  if (type == "all") {
    message("Plotting trace diagnostics...")
    plot_trace_convergence(x, variables = variables, ...)

    if (!is.null(x$summary_stats)) {
      if (interactive()) readline(prompt = "Press [Enter] to see step-size plots...")
      plot_steps_convergence(x, variables = variables, ...)

      if (interactive()) readline(prompt = "Press [Enter] to see ACF plots...")
      plot_acf_convergence(x, variables = variables, ...)
    }

  } else if (type == "trace") {
    plot_trace_convergence(x, variables = variables, ...)
  } else if (type == "density") {
    plot_density_convergence(x, variables = variables, ...)
  } else if (type == "steps") {
    plot_steps_convergence(x, variables = variables, ...)
  } else if (type == "acf") {
    plot_acf_convergence(x, variables = variables, ...)
  }
  invisible(NULL)
}


#' @title plot_trace_convergence
#' @description Trace of the monitored mean of each variable's imputed values
#'   across the chained runs, with the warm-up half shaded and each variable's
#'   verdict from the stability check in the panel header.
#' @param conv_object Object of class bace_convergence
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @return Invisibly returns \code{NULL}; called for its side effect of drawing
#'   trace panels (one per variable) on the current graphics device.
#' @export
plot_trace_convergence <- function(conv_object, variables = NULL, ...) {

  if (is.null(conv_object$summary_stats)) {
    stop("No summary statistics available")
  }

  stats_df <- conv_object$summary_stats
  var_cols <- .conv_plot_vars(stats_df, variables)

  n_vars <- length(var_cols)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  n_rows <- ceiling(sqrt(n_vars))
  n_cols <- ceiling(n_vars / n_rows)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

  warmup <- conv_object$diagnostics$warmup
  assessments <- conv_object$method_results$summary$assessments

  for (var in var_cols) {
    y <- stats_df[[var]]
    if (sum(is.finite(y)) == 0) {
      plot.new(); title(paste("Trace:", var))
      text(0.5, 0.5, "nothing to plot\n(no finite values)", cex = 1.1)
      next
    }
    ylim <- .conv_finite_range(y)

    plot(stats_df$iteration, y,
         type = "b",
         pch = 19,
         col = "steelblue",
         xlab = "Chained run",
         ylab = "Mean of imputed values",
         main = paste("Trace:", var),
         ylim = ylim,
         ...)

    # Shade the warm-up half that the stability check discards
    if (!is.null(warmup) && warmup >= 1) {
      rect(par("usr")[1], par("usr")[3], warmup + 0.5, par("usr")[4],
           col = adjustcolor("grey", 0.25), border = NA)
    }

    # Smoothed line over the retained window
    ok_i <- is.finite(y)
    if (sum(ok_i) >= 3) {
      smooth_vals <- stats::lowess(stats_df$iteration[ok_i], y[ok_i])
      lines(smooth_vals$x, smooth_vals$y, col = "red", lwd = 2)
    }

    # Verdict from the stability check
    if (!is.null(assessments) && var %in% names(assessments)) {
      a <- assessments[[var]]
      col_indicator <- ifelse(a$verdict == "ok", "darkgreen", "darkred")
      status_text <- if (a$verdict == "ok") "no pathology detected"
                     else paste0("[X] ", a$verdict)
      if (nzchar(a$note)) status_text <- paste0(status_text, " *")
      mtext(status_text, side = 3, line = 0.5, col = col_indicator, cex = 0.8)
    }

    grid()
  }
}


#' @title plot_density_convergence
#' @description Plot density evolution across iterations
#' @param conv_object Object of class bace_convergence  
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @return Invisibly returns \code{NULL}; called for its side effect of drawing
#'   density-evolution panels (one per variable) on the current graphics device.
#' @export
plot_density_convergence <- function(conv_object, variables = NULL, ...) {
  
  if (is.null(conv_object$summary_stats)) {
    stop("No summary statistics available")
  }
  
  stats_df <- conv_object$summary_stats
  var_cols <- .conv_plot_vars(stats_df, variables)

  n_vars <- length(var_cols)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(ceiling(sqrt(n_vars)), ceiling(n_vars / ceiling(sqrt(n_vars)))),
      mar = c(4, 4, 3, 1))
  
  n_iter <- nrow(stats_df)
  colors <- colorRampPalette(c("lightblue", "darkblue"))(n_iter)
  
  for (var in var_cols) {
    if (sum(is.finite(stats_df[[var]])) == 0) {
      plot.new(); title(paste("Density Evolution:", var))
      text(0.5, 0.5, "nothing to plot\n(no finite values)", cex = 1.1)
      next
    }
    plot(NULL, 
         xlim = .conv_finite_range(stats_df[[var]]),
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
    
    # Add smoothed trend (finite points only: lowess dies on NaN)
    ok_i <- is.finite(stats_df[[var]])
    if (sum(ok_i) >= 3) {
      smooth_vals <- stats::lowess(stats_df[[var]][ok_i], which(ok_i))
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


#' @title plot_steps_convergence
#' @description Run-to-run movement of each variable's monitored mean, in
#'   units of that variable's imputation spread — the statistic the stability
#'   check actually uses, drawn with its two reference lines: the Monte-Carlo
#'   noise floor sqrt(2/n_missing) (movement this size is expected by chance)
#'   and the loud-oscillation line at 3 noise floors. This panel replaces the
#'   pre-1.0.1 "% change" panel, which divided each step by the variable's own
#'   level and therefore exploded for any mean-centred variable (see NEWS).
#' @param conv_object Object of class bace_convergence
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @return Invisibly returns \code{NULL}; called for its side effect of drawing
#'   step-size panels (one per variable) on the current graphics device.
#' @export
plot_steps_convergence <- function(conv_object, variables = NULL, ...) {

  if (is.null(conv_object$summary_stats)) {
    stop("No summary statistics available")
  }

  stats_df <- conv_object$summary_stats
  sd_df    <- conv_object$sd_stats
  var_cols <- .conv_plot_vars(stats_df, variables)

  n_vars <- length(var_cols)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  n_rows <- ceiling(sqrt(n_vars))
  n_cols <- ceiling(n_vars / n_rows)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

  warmup <- conv_object$diagnostics$warmup
  assessments <- conv_object$method_results$summary$assessments

  for (var in var_cols) {
    m <- stats_df[[var]]
    n_iter <- length(m)

    if (n_iter < 2 || sum(is.finite(m)) < 2) {
      plot.new()
      text(0.5, 0.5, "Insufficient iterations", cex = 1.2)
      next
    }

    # the same yardstick the check uses: the late-chain imputation spread
    keep <- max(1, floor(n_iter / 2)):n_iter
    spread <- if (!is.null(sd_df)) mean(sd_df[[var]][keep]) else NA_real_
    if (!is.finite(spread) || spread < 1e-12) spread <- stats::sd(m[keep])
    if (!is.finite(spread) || spread < 1e-12) spread <- 1e-12

    steps <- abs(diff(m)) / spread
    if (sum(is.finite(steps)) == 0) {
      plot.new(); title(paste("Step size:", var))
      text(0.5, 0.5, "nothing to plot\n(no finite steps)", cex = 1.1)
      next
    }
    ref <- if (!is.null(assessments) && var %in% names(assessments)) {
      unname(assessments[[var]]$stats["ref"])
    } else {
      NA_real_
    }

    plot(2:n_iter, steps,
         type = "b",
         pch = 19,
         col = "purple",
         xlab = "Chained run",
         ylab = "|change| / imputation spread",
         main = paste("Step size:", var),
         ylim = c(0, max(c(steps[is.finite(steps)], ref * 3.5),
                         na.rm = TRUE) * 1.1),
         ...)

    # Shade the warm-up half that the stability check discards
    if (!is.null(warmup) && warmup >= 1) {
      rect(par("usr")[1], par("usr")[3], warmup + 0.5, par("usr")[4],
           col = adjustcolor("grey", 0.25), border = NA)
    }

    # Reference lines: the Monte-Carlo noise floor, and the loud line
    if (is.finite(ref)) {
      abline(h = ref, col = "darkgreen", lty = 2, lwd = 1.5)
      abline(h = 3 * ref, col = "red", lty = 2, lwd = 1.5)
      legend("topright",
             legend = c("noise floor", "3x floor (loud)"),
             col = c("darkgreen", "red"), lty = 2, lwd = 1.5, cex = 0.7)
    }

    # Verdict from the stability check
    if (!is.null(assessments) && var %in% names(assessments)) {
      a <- assessments[[var]]
      col_indicator <- ifelse(a$verdict == "ok", "darkgreen", "darkred")
      status_text <- if (a$verdict == "ok") {
        sprintf("mean step %.2f", unname(a$stats["step"]))
      } else {
        paste0("[X] ", a$verdict)
      }
      mtext(status_text, side = 3, line = 0.5, col = col_indicator, cex = 0.8)
    }

    grid()
  }
}


#' @title plot_pct_change_convergence
#' @description Defunct. The percentage-change panel drew the pre-1.0.1
#'   criterion, which divided each step by the variable's own level and so
#'   exploded for mean-centred variables. Use
#'   \code{\link{plot_steps_convergence}}.
#' @param conv_object Object of class bace_convergence
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @return Never returns; always throws an error pointing to the replacement.
#' @export
plot_pct_change_convergence <- function(conv_object, variables = NULL, ...) {
  stop("plot_pct_change_convergence() is defunct: it drew the pre-1.0.1 ",
       "%-change criterion, which measured each variable's scale rather ",
       "than its stability. Use plot_steps_convergence(). See NEWS.")
}


#' @title plot_acf_convergence
#' @description Autocorrelation panels for each variable's monitored mean
#'   series. Interpretation note: a chained sampler alternates a little by
#'   nature, so a negative lag-1 autocorrelation is not by itself a failure —
#'   the stability check flags oscillation only when it is strong AND the
#'   movement is well above the Monte-Carlo noise floor.
#' @param conv_object Object of class bace_convergence
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @return Invisibly returns \code{NULL}; called for its side effect of drawing
#'   autocorrelation (ACF) panels (one per variable) on the current graphics device.
#' @export
plot_acf_convergence <- function(conv_object, variables = NULL, ...) {
  
  if (is.null(conv_object$summary_stats)) {
    stop("No summary statistics available")
  }
  
  stats_df <- conv_object$summary_stats
  var_cols <- .conv_plot_vars(stats_df, variables)
  
  n_vars <- length(var_cols)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Calculate proper layout
  n_rows <- ceiling(sqrt(n_vars))
  n_cols <- ceiling(n_vars / n_rows)
  
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  
  assessments <- conv_object$method_results$summary$assessments

  for (var in var_cols) {
    y <- stats_df[[var]]
    ## acf() of a constant or all-non-finite series is all NaN, which makes
    ## plot.window fail on infinite ylim. Both are states the check tolerates.
    if (length(y) <= 3 || sum(is.finite(y)) < 3 ||
        stats::sd(y[is.finite(y)]) < 1e-12) {
      plot.new(); title(paste("ACF:", var))
      text(0.5, 0.5, "no autocorrelation to show\n(series constant or too short)",
           cex = 1.0)
      next
    }
    {
      acf(y, 
          main = paste("ACF:", var),
          na.action = na.pass,
          ...)
      
      # Report the kept-window lag-1 autocorrelation the check used
      if (!is.null(assessments) && var %in% names(assessments)) {
        a1 <- unname(assessments[[var]]$stats["acf1"])
        if (is.finite(a1)) {
          mtext(sprintf("lag-1 (post warm-up): %.2f", a1),
                side = 3, line = 0.5, col = "grey30", cex = 0.8)
        }
      }
    }
  }
}


#' @title plot_energy_convergence
#' @description Defunct. The energy-distance method's verdict logic passed
#'   steadily drifting chains and coin-flipped on settled ones, so the method
#'   was removed. See NEWS.
#' @param conv_object Object of class bace_convergence
#' @param ... Additional plotting arguments
#' @return Never returns; always throws an error pointing to the replacement.
#' @export
plot_energy_convergence <- function(conv_object, ...) {
  stop("plot_energy_convergence() is defunct: the energy method was removed ",
       "(its verdict passed drifting chains and coin-flipped on settled ",
       "ones). Use plot(x, type = \"trace\") and type = \"steps\". See NEWS.")
}


#' @title plot_wasserstein_convergence
#' @description Defunct. The Wasserstein method's verdict logic had the same
#'   divide-by-a-vanishing-mean defect as the pre-1.0.1 %-change criterion, so
#'   the method was removed. See NEWS.
#' @param conv_object Object of class bace_convergence
#' @param variables Variables to plot
#' @param ... Additional plotting arguments
#' @return Never returns; always throws an error pointing to the replacement.
#' @export
plot_wasserstein_convergence <- function(conv_object, variables = NULL, ...) {
  stop("plot_wasserstein_convergence() is defunct: the Wasserstein method ",
       "was removed (its verdict divided by a quantity that vanishes exactly ",
       "at convergence). Use plot(x, type = \"trace\") and type = \"steps\". ",
       "See NEWS.")
}


#' @title plot_bace_imputation_comparison
#' @description Compare imputed values across iterations
#' @param bace_object Object of class 'bace' from bace_imp
#' @param variable Variable name to visualize
#' @param n_iter Number of iterations to display. Default is all
#' @param ... Additional plotting arguments
#' @return Called for its side effect of drawing the plot. Invisibly returns the
#'   underlying graphics call's value: a \code{\link[graphics]{boxplot}} statistics
#'   list (gaussian/poisson) or \code{\link[graphics]{barplot}} bar midpoints (categorical).
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
            # barplot()'s formal is `legend.text`; `legend = TRUE` only worked
            # by partial argument matching, which would break silently if
            # graphics::barplot ever gained a `legend` argument of its own.
            legend.text = TRUE,
            las = 2,
            ...)
  }
}




#' @title plot_convergence_summary
#' @description One-page dashboard: the headline verdict, each variable's
#'   verdict from the stability check, and the normalised traces.
#' @param conv_object Object of class bace_convergence
#' @param ... Additional plotting arguments
#' @return Invisibly returns \code{NULL}; called for its side effect of
#'   drawing the summary dashboard on the current graphics device.
#' @export
plot_convergence_summary <- function(conv_object, ...) {

  if (!inherits(conv_object, "bace_convergence")) {
    stop("Input must be an object of class 'bace_convergence'")
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))

  # Panel 1: headline
  par(mar = c(1, 1, 3, 1))
  plot.new()
  title("Stability check")

  if (is.na(conv_object$converged)) {
    text(0.5, 0.6, "NOT ASSESSED", cex = 2, col = "grey40", font = 2)
  } else if (conv_object$converged) {
    text(0.5, 0.6, "no pathology detected", cex = 1.8, col = "darkgreen", font = 2)
  } else {
    text(0.5, 0.6, "PATHOLOGY DETECTED", cex = 1.8, col = "darkred", font = 2)
  }
  text(0.5, 0.3, "a pass is not a certificate - inspect the traces", cex = 0.9,
       col = "grey30")

  # Panel 2: per-variable verdicts
  verdicts <- conv_object$verdicts
  if (!is.null(verdicts) && length(verdicts)) {
    par(mar = c(4, 8, 3, 2))
    verdict_col <- c(ok = "darkgreen", drift = "darkorange",
                     oscillation = "darkred", stuck = "purple4",
                     diverging = "purple4", unassessable = "grey60")
    barplot(rep(1, length(verdicts)),
            horiz = TRUE,
            names.arg = paste0(names(verdicts), ": ", verdicts),
            col = verdict_col[verdicts],
            xlim = c(0, 1),
            main = "Per-variable verdict",
            xlab = "",
            axes = FALSE,
            las = 1)
  }

  # Panel 3: normalised traces
  if (!is.null(conv_object$summary_stats)) {
    par(mar = c(4, 4, 3, 2))

    stats_df <- conv_object$summary_stats
    var_cols <- setdiff(names(stats_df), "iteration")

    plot(NULL,
         xlim = range(stats_df$iteration),
         ylim = c(0, 1),
         main = "Normalised traces (grey band = warm-up)",
         xlab = "Chained run",
         ylab = "Normalised mean of imputed values")

    warmup <- conv_object$diagnostics$warmup
    if (!is.null(warmup) && warmup >= 1) {
      rect(par("usr")[1], par("usr")[3], warmup + 0.5, par("usr")[4],
           col = adjustcolor("grey", 0.25), border = NA)
    }

    colors <- rainbow(length(var_cols))
    for (i in seq_along(var_cols)) {
      var  <- var_cols[i]
      vals <- stats_df[[var]]
      fin  <- vals[is.finite(vals)]
      if (length(fin) == 0) next          # nothing finite to draw
      rng <- max(fin) - min(fin)
      vals_norm <- if (rng < 1e-12) rep(0.5, length(vals))
                   else (vals - min(fin)) / rng
      lines(stats_df$iteration, vals_norm, col = colors[i], lwd = 2)
    }

    legend("topright", legend = var_cols, col = colors, lwd = 2, cex = 0.6)
    grid()
  }
  invisible(NULL)
}
