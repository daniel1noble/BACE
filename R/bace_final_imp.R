#' @title bace_final_imp
#' @description Run final imputation iterations starting from converged data, saving all models
#' @param bace_object An object of class 'bace' from bace_imp
#' @param fixformula A character string or list specifying the fixed effects formulas
#' @param ran_phylo_form A character string specifying the random effects and phylogenetic structure formula
#' @param phylo A phylogenetic tree of class 'phylo' from the ape package
#' @param nitt An integer or list specifying the number of MCMC iterations. Default is 6000
#' @param thin An integer or list specifying the thinning rate. Default is 5
#' @param burnin An integer or list specifying the burn-in iterations. Default is 1000
#' @param n_final An integer specifying the number of final imputation runs.
#'   Default is 50 -- large enough that the empirical 2.5 and 97.5 percent
#'   quantiles of per-cell imputed values give a reasonable estimate of the
#'   95 percent posterior predictive interval. Smaller values undercover:
#'   n_final = 5 gives roughly 67 percent effective coverage, n_final = 20
#'   gives roughly 88 percent, n_final = 50 gives roughly 92 percent.
#' @param species A logical indicating whether to decompose phylogenetic and non-phylogenetic species effects. Default is FALSE
#' @param verbose A logical indicating whether to print progress messages. Default is TRUE
#' @param n_cores Integer specifying the number of parallel cores to use. Default is 1
#'   (serial). Values > 1 use \code{parallel::mclapply}. The function falls back to serial
#'   automatically if any parallel worker returns an error.
#' @param nitt_cat_mult Integer multiplier applied to nitt and burnin for categorical and
#'   threshold/ordinal variables. Default is 1 (no change). Set to 2 or 3 to give harder-to-mix
#'   categorical models proportionally longer chains.
#' @param ovr_categorical Logical. If TRUE, categorical variables are modelled using
#'   one-vs-rest binary threshold MCMCglmm models (J models per variable, one per level)
#'   instead of a single multinomial probit. Default is FALSE: in the common
#'   one-observation-per-species regime the OVR binary models are prone to
#'   quasi-separation (overconfident, jointly incoherent normalized
#'   probabilities), while the single multinomial probit is both better
#'   calibrated and faster. Set to TRUE only if a specific multinomial model
#'   mixes too poorly to use.
#' @param priors Optional named list of MCMCglmm prior overrides, one element
#'   per imputed variable; same structure and rules as in \code{bace_imp()}
#'   (blocks \code{G}/\code{R}/\code{B} replaced wholesale; no \code{R}
#'   override for threshold/ordinal/categorical variables). Use the same
#'   \code{priors} in both phases so the convergence chain and the final
#'   imputations target the same models.
#' @param sweeps Integer (default 3): number of stochastic passes over the
#'   variables within each final run. Variables are imputed one at a time, so
#'   with a single pass the first-visited variable would condition on the SAME
#'   chained-phase point estimates of the other variables in every run,
#'   understating between-imputation variance and making results depend on
#'   formula order. Repeating the pass lets every variable condition on fresh
#'   posterior-predictive draws of the others (a Gibbs sweep over the missing
#'   data); only the last pass's values and models are kept. Warm-up passes
#'   use a shortened post-burn-in tail (full burn-in, 50 stored draws), so the
#'   default costs roughly 1.4x a single pass, not 3x. Automatically reduced
#'   to 1 when at most one variable has missing data (extra passes are then
#'   statistically inert). \code{sweeps = 1} reproduces the pre-0.1.0
#'   single-pass behaviour.
#' @param seed Optional integer. If supplied, each final run gets its own
#'   L'Ecuyer-CMRG random-number stream: results are reproducible and
#'   identical whether runs execute serially or in parallel (any
#'   \code{n_cores}), and the serial retry of a failed parallel worker reuses
#'   exactly the stream that worker should have consumed. The caller's RNG
#'   state is restored on exit. Default NULL (previous behaviour).
#' @param ... Additional arguments passed to modeling functions
#' @return A list of class 'bace_final' containing:
#'   - all_models: List of length n_final, each containing models for all variables
#'   - all_ovr_models: List of length n_final; for variables imputed via the
#'     one-vs-rest path, the J per-level binary MCMCglmm models (empty otherwise)
#'   - all_datasets: List of length n_final with imputed datasets
#'   - all_prob_preds: List of length n_final; per run, for each categorical/
#'     threshold/binary variable, the probability matrix at that variable's
#'     originally missing rows (rownames = species), with category
#'     probabilities averaged over the model's posterior draws (NOT evaluated
#'     at posterior-mean coefficients, which would be overconfident). The
#'     across-run average gives MI-pooled probabilities.
#'   - formulas: The formulas used
#'   - types: Variable types
#'   - phylo_ran: Phylogenetic random effects structure
#'   - call: The function call
#' @examples
#' \donttest{
#' # Simulate a small phylogenetic dataset with missing values
#' set.seed(1)
#' sim <- sim_bace(response_type = "gaussian", predictor_types = "gaussian",
#'                 phylo_signal = c(0.6, 0.6), missingness = c(0.2, 0),
#'                 n_cases = 40, n_species = 40)
#'
#' # Convergence chain (deterministic point estimates)
#' fit <- bace_imp(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
#'                 phylo = sim$tree, data = sim$data,
#'                 runs = 3, nitt = 600, burnin = 100, thin = 2, verbose = FALSE)
#' conv <- assess_convergence(fit, method = "summary")
#'
#' # Final independent posterior-predictive imputations
#' final <- bace_final_imp(fit, fixformula = "y ~ x1",
#'                         ran_phylo_form = "~ 1 | Species",
#'                         phylo = sim$tree, n_final = 2, nitt = 600,
#'                         burnin = 100, thin = 2, verbose = FALSE)
#' length(final$all_datasets)
#' }
#' @export
bace_final_imp <- function(bace_object, fixformula, ran_phylo_form, phylo,
                           nitt = 6000, thin = 5, burnin = 1000,
                           n_final = 50, species = FALSE, verbose = TRUE,
                           n_cores = 1, nitt_cat_mult = 1L, ovr_categorical = FALSE,
                           priors = NULL, sweeps = 3L, seed = NULL, ...) {

  # Check inputs
  if (!inherits(bace_object, "bace")) {
    stop("bace_object must be an object of class 'bace' from bace_imp function")
  }

  # Extract components from bace_object
  last_data <- bace_object$data[[length(bace_object$data)]]
  miss_dat <- bace_object$miss_dat
  types <- bace_object$types
  phylo_ran <- bace_object$phylo_ran

  # Validate any user-supplied priors BEFORE fitting anything (fail fast).
  .validate_priors(priors, types, miss_vars = unique(miss_dat$colname))

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

  # Apply nitt/burnin multiplier for categorical/threshold (harder to mix)
  if (nitt_cat_mult != 1L) {
    for (i_mult in seq_along(formulas)) {
      rv_mult <- all.vars(formulas[[i_mult]][[2]])
      if (!is.null(types[[rv_mult]]) && types[[rv_mult]] %in% c("categorical", "threshold", "ordinal")) {
        nitt_list[[i_mult]]   <- round(nitt_list[[i_mult]]   * nitt_cat_mult)
        burnin_list[[i_mult]] <- round(burnin_list[[i_mult]] * nitt_cat_mult)
      }
    }
  }

  # Variables to impute
  fix <- names(types)

  # Clamp n_cores
  n_cores <- max(1L, min(as.integer(n_cores), n_final))

  # --- Stochastic sweeps (proper MI under multi-variable missingness) -------
  # The variables are imputed one at a time, so whichever is visited first
  # would otherwise condition on the SAME chained-phase point fills of the
  # others in every run (between-imputation variance collapses for it, and
  # formula order changes results). Cycling the variables `sweeps` times with
  # sample = TRUE throughout is a Gibbs sweep over the missing data: by pass
  # 2+ every variable conditions on fresh random draws of the others, and the
  # dataset kept from the LAST pass is (approximately) one draw from the
  # joint posterior predictive distribution.
  sweeps <- as.integer(sweeps)
  if (length(sweeps) != 1L || is.na(sweeps) || sweeps < 1L) {
    stop("sweeps must be a single integer >= 1")
  }
  # With a single missing variable, one pass already IS an exact draw from
  # the joint predictive (there is nothing else missing to condition on), so
  # extra passes change nothing and are skipped.
  n_miss_vars <- length(unique(miss_dat$colname))
  eff_sweeps  <- if (n_miss_vars <= 1L) 1L else sweeps

  # --- Reproducibility: one L'Ecuyer-CMRG stream per run --------------------
  # With `seed` set, serial and parallel execution give identical results,
  # reruns reproduce, and the serial fallback reuses exactly the stream a
  # failed worker should have consumed. The caller's RNG state is restored
  # on exit. With seed = NULL, behaviour is unchanged (mc.set.seed = TRUE).
  streams <- NULL
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    old_kind <- RNGkind()
    on.exit({
      do.call(RNGkind, as.list(old_kind))
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    streams <- vector("list", n_final)
    s_cur <- get(".Random.seed", envir = globalenv())
    for (r in seq_len(n_final)) {
      streams[[r]] <- s_cur
      s_cur <- parallel::nextRNGStream(s_cur)
    }
  }

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
  .one_final_run <- function(run, data_start, worker_verbose, rng_state = NULL) {
    if (!is.null(rng_state)) {
      assign(".Random.seed", rng_state, envir = globalenv())
    }
    data_current <- data_start
    models_this_run <- list()
    ovr_models_this_run <- list()   # per-level binary models for OVR categorical vars
    prob_preds_this_run <- list()   # probability matrices for cat/threshold vars

    for (s in seq_len(eff_sweeps)) {
    last_sweep <- (s == eff_sweeps)

    for (i in seq_along(formulas)) {

      response_var <- all.vars(formulas[[i]][[2]])

      # Warm-up sweeps exist only to randomise the conditioning values; they
      # need ONE approximately-stationary posterior draw, not the long
      # post-burn-in tail that buys ESS for pooled inference. So keep the
      # full burn-in (never shorten it: a draw from an unburned chain
      # re-anchors to the inits, recreating the very bias the sweeps remove)
      # and a short 50-draw tail. Only the last sweep's full-length models
      # are kept and pooled.
      nitt_run <- if (last_sweep) nitt_list[[i]] else
        burnin_list[[i]] + 50L * thin_list[[i]]

      # Reinstate NAs on the current response (mirroring bace_imp) so that
      # MCMCglmm performs data augmentation on the originally missing cells
      # rather than conditioning on the convergence chain's point estimates
      # as if they were observed data. Without this the final "posterior-
      # predictive" draws are anchored to the deterministic fills: the
      # residual variance is deflated, the BLUPs are anchored, and the
      # between-imputation variance collapses (improper MI; Rubin 1987 Ch.4).
      # Predictors keep their most recent imputed values.
      data_for_model <- data_current
      resp_na_rows <- miss_dat[miss_dat$colname == response_var, "row"]
      if (length(resp_na_rows) > 0) {
        data_for_model[resp_na_rows, response_var] <- NA
      }

      dat_prep <- .data_prep(
        data = data_for_model,
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

        # Use the same BACE.gelman option as bace_imp (this phase used to
        # hard-code gelman = 2, silently disagreeing with the chained phase).
        # A user-supplied B block replaces the default anyway, so skip the
        # Gelman construction in that case.
        cat_gelman <- if (is.null(getOption("BACE.gelman"))) {
          bace_option_defaults()$gelman
        } else {
          getOption("BACE.gelman")
        }
        gelman_val <- if (types[[response_var]] != "categorical") 0L
          else if (!is.null(priors[[response_var]]$B)) 0L
          else cat_gelman
        prior_i <- .resolve_prior(
          .make_prior(
            n_rand = n_rand_eff,
            n_levels = levels,
            type = types[[response_var]],
            fixform = fixform,
            data = data_i,
            gelman = gelman_val
          ),
          override = priors[[response_var]],
          type = types[[response_var]], var = response_var
        )

        if (ovr_categorical && types[[response_var]] == "categorical") {
          # One-vs-rest: J binary threshold models, aggregate probabilities
          predictions <- .fit_predict_ovr(
            data_i       = data_i,
            response_var = response_var,
            fixformula   = formulas[[i]],
            phylo        = phylo,
            ran_phylo_form = ran_phylo_form,
            nitt         = nitt_run,
            thin         = thin_list[[i]],
            burnin       = burnin_list[[i]],
            n_rand_eff   = n_rand_eff,
            cluster_col  = phylo_ran[["cluster"]],
            dat_prep     = dat_prep,
            sample       = TRUE,   # draw from posterior predictive each final run
            prior_override = priors[[response_var]]
          )
          # Keep the per-variable slot NULL so pool_posteriors() (which
          # expects one MCMCglmm object per variable) skips OVR variables,
          # but retain the J per-level binary models for diagnostics.
          if (last_sweep) {
            models_this_run[[response_var]] <- NULL
            ovr_models_this_run[[response_var]] <- predictions[["models"]]
          }

        } else {
          model <- .model_fit(
            data = data_i,
            tree = phylo,
            fixformula = formulas[[i]],
            randformula = ran_phylo_form,
            type = types[[response_var]],
            prior = prior_i,
            nitt = nitt_run,
            thin = thin_list[[i]],
            burnin = burnin_list[[i]]
          )

          if (last_sweep) models_this_run[[response_var]] <- model

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
        }

        id      <- miss_dat[miss_dat$colname == response_var, "row"]
        data_id <- which(colnames(data_current) == response_var)
        data_current[id, data_id] <- predictions[["pred_values"]][id]

        # Store probability matrix for categorical/threshold (used for pooled prediction)
        # Attach animal names as rownames so downstream pooling joins by name, not position.
        # make.unique() handles the multi-row-per-species case (e.g. within-species
        # replicates in simulated data); AVONET-like single-row-per-species data is
        # unaffected because make.unique is a no-op on unique inputs.
        if (last_sweep &&
            types[[response_var]] %in% c("categorical", "threshold") &&
            !is.null(predictions[["full_prediction"]])) {
          prob_mat <- predictions[["full_prediction"]][id, , drop = FALSE]
          rownames(prob_mat) <- make.unique(
            as.character(data_current[[phylo_ran[["cluster"]]]][id]))
          prob_preds_this_run[[response_var]] <- prob_mat
        }

        if (worker_verbose) {
          sweep_tag <- if (eff_sweeps > 1L)
            paste0(" (sweep ", s, "/", eff_sweeps, ")") else ""
          cat(paste0("Run ", run, "/", n_final, sweep_tag,
                     " - Imputed variable: ", response_var, "\n"))
        }
      } else {
        if (worker_verbose && s == 1L) {
          cat(paste0("Run ", run, "/", n_final,
                     " - Variable ", response_var,
                     " has no missing data (skipped imputation)\n"))
        }
      }
    }
    }

    # NOTE (intentional, do not "fix"): within a run, the kept model for
    # variable k was fit while variables AFTER k in the visit order still
    # held their previous-sweep draws, whereas the returned dataset holds the
    # final sweep's draws for every variable. Both pooling pathways are
    # consistent (pool_posteriors uses only the kept models; with_imputations
    # / pool_mi use only the datasets), so this one-generation offset is
    # harmless - it is the standard Gibbs state at the time of each update.
    list(models = models_this_run, dataset = data_current,
         prob_preds = prob_preds_this_run, ovr_models = ovr_models_this_run)
  }

  # Run serially or in parallel
  run_one <- function(run, worker_verbose) {
    .one_final_run(run, last_data, worker_verbose = worker_verbose,
                   rng_state = if (is.null(streams)) NULL else streams[[run]])
  }

  if (n_cores > 1L) {
    # Suppress per-worker output to avoid garbled console; progress printed below
    results_list <- parallel::mclapply(
      seq_len(n_final),
      function(run) run_one(run, worker_verbose = FALSE),
      mc.cores    = n_cores,
      # With per-run streams the workers set their own RNG state; without a
      # seed keep mclapply's own per-child seeding (previous behaviour).
      mc.set.seed = is.null(seed)
    )

    # Re-run ONLY the failed runs serially (mclapply returns try-error
    # objects on failure). Successful draws are kept - redoing all n_final
    # runs would waste sweeps x fits of MCMC per healthy run.
    failed <- which(vapply(results_list, inherits, logical(1L), "try-error"))
    if (length(failed) > 0L) {
      if (verbose) {
        cat(sprintf(
          "\nWARNING: %d/%d parallel workers failed -- re-running those runs serially.\n\n",
          length(failed), n_final
        ))
      }
      for (r in failed) {
        results_list[[r]] <- run_one(r, worker_verbose = verbose)
      }
    } else if (verbose) {
      cat(paste0("Completed ", n_final, " parallel imputation runs\n"))
    }
  } else {
    results_list <- lapply(
      seq_len(n_final),
      function(run) run_one(run, worker_verbose = verbose)
    )
  }

  all_models     <- lapply(results_list, `[[`, "models")
  all_datasets   <- lapply(results_list, `[[`, "dataset")
  all_prob_preds <- lapply(results_list, `[[`, "prob_preds")
  all_ovr_models <- lapply(results_list, `[[`, "ovr_models")

  # Prepare output
  out <- list(
    all_models     = all_models,
    all_datasets   = all_datasets,
    all_prob_preds = all_prob_preds,   # prob matrices per run, for pooled prediction
    all_ovr_models = all_ovr_models,   # per-level binary models per run (OVR path only)
    formulas = formulas,
    types = types,
    phylo_ran = phylo_ran,
    miss_dat = miss_dat,
    n_final = n_final,
    sweeps = eff_sweeps,            # sweeps actually performed per run
    sweeps_requested = sweeps,      # as asked (collapses to 1 if <=1 missing var)
    seed = seed,
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
#' @return Returns \code{x} invisibly; called for its side effect of printing
#'   a summary of the \code{bace_final} object to the console.
#' @export
print.bace_final <- function(x, ...) {
  cat("\n=== BACE Final Imputation Results ===\n\n")
  cat("Number of imputation runs:", x$n_final, "\n")
  cat("Number of variables with models:", 
      length(x$all_models[[1]]), "\n")
  cat("Variables:", paste(names(x$all_models[[1]]), collapse = ", "), "\n\n")
  cat("Use pool_posteriors() to combine model results across imputations\n")
  invisible(x)
}
