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
#' @param n_cores Integer specifying the number of parallel cores to use. Default is 1
#'   (serial). Values > 1 use \code{parallel::mclapply}. The function falls back to serial
#'   automatically if any parallel worker returns an error.
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
                           n_final = 10, species = FALSE, verbose = TRUE,
                           n_cores = 1, ...) {

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

  # Clamp n_cores
  n_cores <- max(1L, min(as.integer(n_cores), n_final))

  if (verbose) {
    if (n_cores > 1L) {
      cat("\n=== Running", n_final, "final imputation iterations",
          "(parallel:", n_cores, "cores) ===\n\n")
    } else {
      cat("\n=== Running", n_final, "final imputation iterations ===\n\n")
    }
  }

  # Inner function: one complete final imputation starting from data_start.
  # Each run is independent (all start from the converged last_data), giving
  # truly independent posterior draws suitable for Rubin's rules pooling.
  .one_final_run <- function(run, data_start, worker_verbose) {
    data_current <- data_start
    models_this_run <- list()

    for (i in seq_along(formulas)) {

      response_var <- all.vars(formulas[[i]][[2]])

      dat_prep <- .data_prep(
        data = data_current,
        formula = formulas[[i]],
        types = types,
        ran_cluster = phylo_ran[["cluster"]]
      )

      data_i <- dat_prep[[1]]

      if (species) {
        species_col_name <- phylo_ran[["cluster"]]
        data_i[[paste0(species_col_name, "2")]] <- data_i[[species_col_name]]
      }

      has_missing <- response_var %in% miss_dat$colname

      if (has_missing) {
        fixform    <- formulas[[i]]
        # Compute number of factor levels from the working data (dat_prep$levels is not set)
        levels <- if (types[[response_var]] %in% c("categorical", "threshold")) {
          nlevels(data_current[[response_var]])
        } else {
          NULL
        }
        n_rand_eff <- if (species) 2L else 1L

        gelman_val <- if (types[[response_var]] == "categorical") 2L else 0L
        prior_i <- .make_prior(
          n_rand = n_rand_eff,
          n_levels = levels,
          type = types[[response_var]],
          fixform = fixform,
          data = data_i,
          gelman = gelman_val
        )

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

        models_this_run[[response_var]] <- model

        predictions <- .predict_bace(
          model,
          dat_prep,
          response_var = response_var,
          type         = types[[response_var]],
          sample       = TRUE,               # draw from posterior predictive each run
          formula      = formulas[[i]],
          data_full    = data_i,
          cluster_col  = phylo_ran[["cluster"]]
        )

        id     <- miss_dat[miss_dat$colname == response_var, "row"]
        data_id <- which(colnames(data_current) == response_var)
        data_current[id, data_id] <- predictions[["pred_values"]][id]

        if (worker_verbose) {
          cat(paste0("Run ", run, "/", n_final,
                     " - Imputed variable: ", response_var, "\n"))
        }
      } else {
        if (worker_verbose) {
          cat(paste0("Run ", run, "/", n_final,
                     " - Variable ", response_var,
                     " has no missing data (skipped imputation)\n"))
        }
      }
    }

    list(models = models_this_run, dataset = data_current)
  }

  # Run serially or in parallel
  if (n_cores > 1L) {
    # Suppress per-worker output to avoid garbled console; progress printed below
    results_list <- parallel::mclapply(
      seq_len(n_final),
      function(run) .one_final_run(run, last_data, worker_verbose = FALSE),
      mc.cores    = n_cores,
      mc.set.seed = TRUE
    )

    # Check for worker failures (mclapply returns try-error objects on failure)
    n_failed <- sum(vapply(results_list, inherits, logical(1L), "try-error"))
    if (n_failed > 0L) {
      if (verbose) {
        cat(sprintf(
          "\nWARNING: %d/%d parallel workers failed -- falling back to serial execution.\n\n",
          n_failed, n_final
        ))
      }
      results_list <- lapply(
        seq_len(n_final),
        function(run) .one_final_run(run, last_data, worker_verbose = verbose)
      )
    } else if (verbose) {
      cat(paste0("Completed ", n_final, " parallel imputation runs\n"))
    }
  } else {
    results_list <- lapply(
      seq_len(n_final),
      function(run) .one_final_run(run, last_data, worker_verbose = verbose)
    )
  }

  all_models   <- lapply(results_list, `[[`, "models")
  all_datasets <- lapply(results_list, `[[`, "dataset")
  
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
