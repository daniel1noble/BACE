
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
#' @param ... Additional arguments to be passed to the underlying modeling functions.
#' @return A list of class 'bace_complete' containing:
#'   - pooled_models: Pooled posterior distributions accounting for imputation uncertainty
#'   - final_results: Results from final imputation runs (bace_final object)
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
#' # Access pooled results
#' summary(result$pooled_models$models$y)
#' 
#' # Check convergence
#' print(result$convergence)
#' }
#' @export
bace <- function(fixformula, ran_phylo_form, phylo, data, nitt = 6000, thin = 5, 
                burnin = 1000, runs = 10, n_final = 10, species = FALSE, 
                verbose = TRUE, plot = FALSE, max_attempts = 3, ...) {

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
                 verbose = verbose, ...)

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
# If not converged, try again with more runs
##-----------------------##
while (!converged && n_attempts < max_attempts) {
  
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
                   verbose = verbose, ...)
  
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
    ...
  )
  
  # Pool posteriors across imputations
  if (verbose) {
    cat("\n\nStep 4: Pooling posteriors across imputations\n")
    cat("-------------------------------------------------------\n")
  }
  
  pooled <- pool_posteriors(final_results)
  
  if (verbose) {
    cat("\nPosterior pooling complete!\n")
    cat("Pooled", length(pooled$models), "variable model(s) across", 
        n_final, "imputations\n")
  }
  
} else {
  
  if (verbose) {
    cat("\n\nWARNING: Convergence not achieved after", max_attempts, "attempts\n")
    cat("Consider increasing 'runs' or 'max_attempts' parameters\n")
    cat("Proceeding with final imputation despite lack of convergence...\n")
    cat("-------------------------------------------------------\n")
  }
  
  # Still proceed with final runs, but warn user
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
    ...
  )
  
  pooled <- pool_posteriors(final_results)
  
  warning("BACE completed but convergence was not achieved. Results should be interpreted with caution.")
}

##-----------------------## 
# Prepare and return output
##-----------------------##
out <- list(
  pooled_models = pooled,
  final_results = final_results,
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
  cat("  - $convergence: Convergence diagnostics\n")
  cat("  - $final_results: Final imputation datasets and models\n\n")
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
  cat("  $convergence - Convergence assessment details\n")
  cat("  $final_results - Final imputation runs\n")
  cat("  $initial_results - Initial convergence runs\n\n")
  
  cat("Access pooled MCMCglmm models:\n")
  cat("  result$pooled_models$models$variable_name\n\n")
  cat("Use standard MCMCglmm methods:\n")
  cat("  summary(result$pooled_models$models$y)\n")
  cat("  print(result$pooled_models$models$y)\n")
  cat("  plot(result$pooled_models$models$y)\n")
}