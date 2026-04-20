
#' @title bace
#' @description Wrapper function to perform complete Bayesian imputation using BACE with 
#'   convergence checking and posterior pooling across multiple imputations.
#' @param fixformula A character string specifying the fixed effects formula that is of major 
#'   interest or a list of formulas that one wishes to estimate. This should be of the form: 
#'   y ~ x. Providing a list gives users complex flexibility in the types of models fit for 
#'   the different variables. Note that users can also include interaction terms using the 
#'   standard R formula syntax (e.g., x1*x2). 
#' @param ran_phylo_form A character string specifying the random effects and phylogenetic 
#'   structure formula used in the model.
#' @param phylo A phylogenetic tree of class 'phylo' from the ape package.
#' @param data A data frame containing the dataset with missing values to be imputed.
#' @param runs An integer specifying the number of initial imputation iterations for convergence 
#'   checking. Default is 10.
#' @param nitt An integer or list specifying the number of iterations to run the MCMC algorithm. 
#'   Can be a single value (applied to all models) or a list of values (one per formula). 
#'   Default is 6000.
#' @param thin An integer or list specifying the thinning rate for the MCMC algorithm. Can be a 
#'   single value (applied to all models) or a list of values (one per formula). Default is 5.
#' @param burnin An integer or list specifying the number of iterations to discard as burnin. 
#'   Can be a single value (applied to all models) or a list of values (one per formula). 
#'   Default is 1000.
#' @param n_final An integer specifying the number of final imputation runs to perform after 
#'   convergence is achieved. Default is 10. These runs are used for posterior pooling.
#' @param species A logical indicating whether to decompose phylogenetic and non-phylogenetic 
#'   species effects. Default is FALSE. When TRUE, the random effects structure is modified to 
#'   include both a phylogenetic effect and a non-phylogenetic species effect (with identity 
#'   matrix in ginverse). This requires sufficient replicated observations per species.
#' @param verbose A logical indicating whether to print progress messages. Default is TRUE.
#' @param plot A logical indicating whether to plot convergence diagnostics. Default is FALSE.
#' @param max_attempts Maximum number of attempts to achieve convergence by increasing runs. 
#'   Default is 3.
#' @param skip_conv A logical indicating whether to skip convergence retry logic. When TRUE, 
#'   convergence is still assessed, but if it fails, the function proceeds directly to final 
#'   imputation instead of retrying with more runs. Default is FALSE.
#' @param sample_size Integer specifying how many posterior samples to draw from each
#'   imputation before pooling. If NULL (default), uses all posterior samples. Setting this
#'   to a smaller value (e.g., 1000) can greatly reduce memory usage of the final pooled
#'   models while still properly accounting for imputation and parameter uncertainty.
#' @param nitt_cat_mult Integer multiplier applied to nitt and burnin for categorical and
#'   threshold/ordinal variables. Default is 1 (no change). Set to 2 or 3 to give harder-to-mix
#'   categorical models proportionally longer chains.
#' @param ovr_categorical Logical. If TRUE, categorical variables are modelled using
#'   one-vs-rest binary threshold MCMCglmm models (J models per variable, one per level)
#'   instead of a single multinomial probit. Binary threshold models mix more reliably
#'   and are the recommended default. Default is TRUE.
#' @param phylo_signal Logical. If TRUE, run \code{\link{phylo_signal_summary}} for every
#'   variable in the formula instead of performing the full BACE imputation pipeline. The
#'   function returns the signal-only preview object (class \code{c("phylo_signal",
#'   "bace_preview")}); no call to \code{bace_imp()}, \code{bace_final_imp()}, or
#'   \code{pool_posteriors()} is made. Use this to check whether the data carry enough
#'   phylogenetic signal to justify a full BACE run. Default is FALSE.
#' @param n_cores Integer specifying the number of parallel cores to use for the final
#'   imputation runs. Default is 1 (serial). Values > 1 use \code{parallel::mclapply}.
#'   Note: parallel execution may be unstable on macOS with multithreaded BLAS; the
#'   function falls back to serial automatically if worker failures are detected.
#' @param ... Additional arguments to be passed to the underlying modeling functions.
#' @return A list of class 'bace_complete' containing:
#'   - pooled_models: Pooled posterior distributions accounting for imputation uncertainty
#'   - final_results: Results from final imputation runs (bace_final object)
#'   - imputed_datasets: List of n_final imputed datasets from final runs (for convenient access)
#'   - initial_results: Results from initial convergence runs (bace object)
#'   - convergence: Convergence assessment results
#'   - converged: Logical indicating if convergence was achieved
#'   - n_attempts: Number of attempts needed to achieve convergence
#'   - call: The function call
#' @examples \dontrun{
#' # Complete BACE analysis with convergence checking and posterior pooling
#' result <- bace(
#'   fixformula = "y ~ x1 + x2",
#'   ran_phylo_form = "~1|Species",
#'   phylo = phylo_tree,
#'   data = my_data,
#'   runs = 10,
#'   n_final = 10
#' )
#' 
#' # With posterior sampling to reduce memory usage
#' result <- bace(
#'   fixformula = "y ~ x1 + x2",
#'   ran_phylo_form = "~1|Species",
#'   phylo = phylo_tree,
#'   data = my_data,
#'   runs = 10,
#'   n_final = 10,
#'   sample_size = 1000  # Sample 1000 draws from each imputation
#' )
#' 
#' # Access pooled results
#' summary(result$pooled_models$models$y)
#' 
#' # Access final imputed datasets
#' imputed_data1 <- result$imputed_datasets[[1]]
#' 
#' # Check convergence
#' print(result$convergence)
#' }
#' @export
bace <- function(fixformula, ran_phylo_form, phylo, data, nitt = 6000, thin = 5,
                burnin = 1000, runs = 10, n_final = 10, species = FALSE,
                verbose = TRUE, plot = FALSE, max_attempts = 3, skip_conv = FALSE,
                sample_size = NULL, n_cores = 1L,
                nitt_cat_mult = 1L, ovr_categorical = TRUE,
                phylo_signal = FALSE, ...) {

##-----------------------##
# phylo_signal = TRUE: short-circuit to signal-only preview
##-----------------------##
if (isTRUE(phylo_signal)) {
  if (verbose) {
    cat("\n=======================================================\n")
    cat("BACE: phylo_signal = TRUE (signal-only preview)\n")
    cat("=======================================================\n\n")
  }
  # Collect all variables across fix formula(s)
  if (is.list(fixformula) && !inherits(fixformula, "formula")) {
    sig_vars <- unique(unlist(lapply(fixformula,
                                     function(f) .get_variables(f, fix = TRUE)$fix)))
  } else {
    sig_vars <- .get_variables(fixformula, fix = TRUE)$fix
  }
  # Extract species column from ran_phylo_form
  sig_cluster <- .get_variables(ran_phylo_form, fix = FALSE)$cluster
  sig_species_col <- if (length(sig_cluster) >= 1) sig_cluster[[1]] else "Species"

  preview <- phylo_signal_summary(
    data            = data,
    tree            = phylo,
    species_col     = sig_species_col,
    variables       = sig_vars,
    species         = species,
    ovr_categorical = ovr_categorical,
    verbose         = verbose
  )
  if (verbose) print(preview)
  class(preview) <- c("phylo_signal", "bace_preview")
  message("phylo_signal = TRUE: returning signal-only preview. ",
          "Rerun with phylo_signal = FALSE to proceed with imputation.")
  return(invisible(preview))
}

##-----------------------##
# Run bace_imp first
##-----------------------##
if (verbose) {
  cat("\n=======================================================\n")
  cat("BACE: Bayesian Analysis with Chained Equations\n")
  cat("=======================================================\n\n")
  cat("Step 1: Initial imputation for convergence assessment\n")
  cat("-------------------------------------------------------\n")
}

start <- bace_imp(fixformula = fixformula, ran_phylo_form = ran_phylo_form,
                 phylo = phylo, data = data, nitt = nitt, thin = thin,
                 burnin = burnin, runs = runs, species = species,
                 verbose = verbose, nitt_cat_mult = nitt_cat_mult,
                 ovr_categorical = ovr_categorical, ...)

##-----------------------## 
# Check convergence 
##-----------------------##
if (verbose) {
  cat("\n\nStep 2: Convergence assessment\n")
  cat("-------------------------------------------------------\n")
}

converge <- assess_convergence(start, method = "summary")

if (plot) { 
  plot(converge) 
}

# Track convergence attempts
n_attempts <- 1
converged <- converge$converged

# Treat NA as not converged
if (is.na(converged)) converged <- FALSE

##-----------------------## 
# If not converged, try again with more runs (unless skip_conv is TRUE)
##-----------------------##
while (!converged && n_attempts < max_attempts && !skip_conv) {
  
  if (verbose) {
    cat("\n\nConvergence not achieved. Running additional iterations...\n")
    cat("Attempt", n_attempts + 1, "of", max_attempts, "\n")
    cat("-------------------------------------------------------\n")
  }
  
  # Increase runs by 50% each attempt
  runs_new <- round(runs * (1.5 ^ n_attempts))
  
  start <- bace_imp(fixformula = fixformula, ran_phylo_form = ran_phylo_form,
                   phylo = phylo, data = data, nitt = nitt, thin = thin,
                   burnin = burnin, runs = runs_new, species = species,
                   verbose = verbose, nitt_cat_mult = nitt_cat_mult,
                   ovr_categorical = ovr_categorical, ...)
  
  converge <- assess_convergence(start, method = "summary")
  
  if (plot) { 
    plot(converge) 
  }
  
  converged <- converge$converged
  # Treat NA as not converged
  if (is.na(converged)) converged <- FALSE
  n_attempts <- n_attempts + 1
}

##-----------------------## 
# Starting from last dataset, do final runs and pool posteriors
##-----------------------##
if (converged) {
  
  if (verbose) {
    cat("\n\nStep 3: Final imputation runs for posterior pooling\n")
    cat("-------------------------------------------------------\n")
  }
  
  # Run final imputations starting from converged data
  final_results <- bace_final_imp(
    bace_object = start,
    fixformula = fixformula,
    ran_phylo_form = ran_phylo_form,
    phylo = phylo,
    nitt = nitt,
    thin = thin,
    burnin = burnin,
    n_final = n_final,
    species = species,
    verbose = verbose,
    n_cores = n_cores,
    nitt_cat_mult = nitt_cat_mult,
    ovr_categorical = ovr_categorical,
    ...
  )

  # Pool posteriors across imputations
  if (verbose) {
    cat("\n\nStep 4: Pooling posteriors across imputations\n")
    cat("-------------------------------------------------------\n")
    if (!is.null(sample_size)) {
      cat("Using posterior sampling:", sample_size, "samples per imputation\n")
    }
  }

  pooled <- pool_posteriors(final_results, sample_size = sample_size)

  if (verbose) {
    cat("\nPosterior pooling complete!\n")
    cat("Pooled", length(pooled$models), "variable model(s) across",
        n_final, "imputations\n")
  }

} else {

  if (verbose) {
    if (!skip_conv) {
      cat("\n\nWARNING: Convergence not achieved after", max_attempts, "attempts\n")
      cat("Consider increasing 'runs' or 'max_attempts' parameters\n")
    } else {
      cat("\n\nSkipping convergence retry (skip_conv = TRUE)\n")
    }
    cat("Proceeding with final imputation...\n")
    cat("-------------------------------------------------------\n")
  }

  # Still proceed with final runs, but warn user (only if not explicitly skipping)
  final_results <- bace_final_imp(
    bace_object = start,
    fixformula = fixformula,
    ran_phylo_form = ran_phylo_form,
    phylo = phylo,
    nitt = nitt,
    thin = thin,
    burnin = burnin,
    n_final = n_final,
    species = species,
    verbose = verbose,
    n_cores = n_cores,
    nitt_cat_mult = nitt_cat_mult,
    ovr_categorical = ovr_categorical,
    ...
  )
  
  if (verbose && !is.null(sample_size)) {
    cat("\nPooling posteriors with sampling:", sample_size, "samples per imputation\n")
  }
  
  pooled <- pool_posteriors(final_results, sample_size = sample_size)
  
  # Only warn if convergence was attempted but not achieved (not when explicitly skipped)
  if (!skip_conv) {
    warning("BACE completed but convergence was not achieved. Results should be interpreted with caution.")
  }
}

##-----------------------## 
# Prepare and return output
##-----------------------##
out <- list(
  pooled_models = pooled,
  final_results = final_results,
  imputed_datasets = final_results$all_datasets,   # For convenient access
  prob_preds       = final_results$all_prob_preds, # Pooled probability matrices (cat/threshold)
  initial_results = start,
  convergence = converge,
  converged = converged,
  n_attempts = n_attempts,
  call = match.call()
)

class(out) <- "bace_complete"

if (verbose) {
  cat("\n=======================================================\n")
  cat("BACE Analysis Complete!\n")
  cat("=======================================================\n\n")
  cat("Convergence achieved:", converged, "\n")
  cat("Number of attempts:", n_attempts, "\n")
  cat("Final imputations:", n_final, "\n")
  cat("\nAccess results via:\n")
  cat("  - $pooled_models: Pooled posterior distributions\n")
  cat("  - $imputed_datasets: List of", n_final, "imputed datasets\n")
  cat("  - $convergence: Convergence diagnostics\n")
  cat("  - $final_results: Full final imputation results\n\n")
}

return(out)
}


#' @title print.bace_complete
#' @description Print method for bace_complete objects
#' @param x Object of class bace_complete
#' @param ... Additional arguments
#' @export
print.bace_complete <- function(x, ...) {
  cat("\n=== BACE Complete Analysis Results ===\n\n")
  cat("Convergence Status:", ifelse(x$converged, "ACHIEVED", "NOT ACHIEVED"), "\n")
  cat("Number of attempts:", x$n_attempts, "\n")
  cat("Final imputations:", x$pooled_models$n_imputations, "\n")
  cat("Variables analyzed:", length(x$pooled_models$models), "\n")
  cat("  -", paste(x$pooled_models$variables, collapse = ", "), "\n\n")
  
  if (!x$converged) {
    cat("WARNING: Convergence was not achieved. Results should be interpreted with caution.\n\n")
  }
  
  cat("Components:\n")
  cat("  $pooled_models - Contains MCMCglmm objects with imputation uncertainty\n")
  cat("  $imputed_datasets - List of", length(x$imputed_datasets), "imputed datasets for additional analyses\n")
  cat("  $convergence - Convergence assessment details\n")
  cat("  $final_results - Full final imputation results\n")
  cat("  $initial_results - Initial convergence runs\n\n")
  
  cat("Access pooled MCMCglmm models:\n")
  cat("  result$pooled_models$models$variable_name\n\n")
  cat("Access imputed datasets:\n")
  cat("  result$imputed_datasets[[1]]  # First imputed dataset\n\n")
  cat("Use standard MCMCglmm methods:\n")
  cat("  summary(result$pooled_models$models$y)\n")
  cat("  print(result$pooled_models$models$y)\n")
  cat("  plot(result$pooled_models$models$y)\n")
}