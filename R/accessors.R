#' @title get_pooled_model
#' @description Extract a pooled MCMCglmm model for a specific variable from a 
#'   \code{bace_complete} or \code{bace_pooled} object. The returned model has 
#'   posteriors pooled across all final imputation runs so that inference properly 
#'   accounts for imputation uncertainty.
#' @param object A \code{bace_complete} object (returned by \code{bace()}) or a 
#'   \code{bace_pooled} object (returned by \code{pool_posteriors()}).
#' @param variable A character string giving the name of the response variable 
#'   whose pooled model should be extracted. If \code{NULL} (the default), a named 
#'   list of all pooled models is returned.
#' @return When \code{variable} is specified, a single MCMCglmm model object of 
#'   class \code{c("MCMCglmm", "bace_pooled_MCMCglmm")}. When \code{variable} is 
#'   \code{NULL}, a named list of all such models. Each model contains pooled 
#'   posteriors (Sol, VCV, etc.) and a \code{BACE_pooling} element with metadata 
#'   about the pooling process.
#' @examples \donttest{
#' # Tiny toy analysis; use far larger nitt/burnin/runs/n_final in practice
#' set.seed(1)
#' sim <- sim_bace(response_type = "gaussian", predictor_types = "gaussian",
#'                 phylo_signal = c(0.6, 0.6), missingness = c(0.2, 0),
#'                 n_cases = 40, n_species = 40)
#' result <- bace(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
#'                phylo = sim$tree, data = sim$data,
#'                runs = 3, nitt = 600, burnin = 100, thin = 2,
#'                n_final = 2, skip_conv = TRUE, verbose = FALSE)
#'
#' # Extract the pooled model for the response variable
#' y_model <- get_pooled_model(result, variable = "y")
#' summary(y_model)
#'
#' # Extract all pooled models as a named list
#' all_models <- get_pooled_model(result)
#' names(all_models)
#' }
#' @export
get_pooled_model <- function(object, variable = NULL) {

  # --- Resolve to the bace_pooled layer ---
  if (inherits(object, "bace_complete")) {
    pooled <- object$pooled_models
  } else if (inherits(object, "bace_pooled")) {
    pooled <- object
  } else {
    stop("'object' must be a 'bace_complete' object (from bace()) ",
         "or a 'bace_pooled' object (from pool_posteriors()).")
  }

  models <- pooled$models

  # --- Return all models when variable is NULL ---
  if (is.null(variable)) {
    return(models)
  }

  # --- Validate variable name ---
  if (!is.character(variable) || length(variable) != 1L) {
    stop("'variable' must be a single character string.")
  }

  available <- names(models)
  if (!variable %in% available) {
    stop("Variable '", variable, "' not found. Available variables: ",
         paste(available, collapse = ", "))
  }

  models[[variable]]
}


#' @title get_imputed_data
#' @description Extract the complete (imputed) datasets from a \code{bace_complete} 
#'   or \code{bace_final} object. These are the datasets produced during the final 
#'   imputation runs after convergence, with all missing values filled in.
#' @param object A \code{bace_complete} object (returned by \code{bace()}) or a 
#'   \code{bace_final} object (returned by \code{bace_final_imp()}).
#' @param format Character string specifying the output format. One of 
#'   \code{"list"} (default) or \code{"data.frame"}. When \code{"list"}, a list of 
#'   data frames is returned (one per imputation run). When \code{"data.frame"}, all 
#'   imputed datasets are row-bound into a single data frame with an additional 
#'   \code{.imputation} column indicating the imputation run number.
#' @return Depending on \code{format}:
#'   \describe{
#'     \item{\code{"list"}}{A list of data frames, one per final imputation run.}
#'     \item{\code{"data.frame"}}{A single data frame with all imputed datasets 
#'       stacked and an \code{.imputation} column (integer) identifying the run.}
#'   }
#' @examples \donttest{
#' # Tiny toy analysis; use far larger nitt/burnin/runs/n_final in practice
#' set.seed(1)
#' sim <- sim_bace(response_type = "gaussian", predictor_types = "gaussian",
#'                 phylo_signal = c(0.6, 0.6), missingness = c(0.2, 0),
#'                 n_cases = 40, n_species = 40)
#' result <- bace(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
#'                phylo = sim$tree, data = sim$data,
#'                runs = 3, nitt = 600, burnin = 100, thin = 2,
#'                n_final = 2, skip_conv = TRUE, verbose = FALSE)
#'
#' # Get imputed datasets as a list
#' imp_list <- get_imputed_data(result)
#' length(imp_list)       # number of imputations
#' head(imp_list[[1]])    # first imputed dataset
#'
#' # Get as a single stacked data frame
#' imp_df <- get_imputed_data(result, format = "data.frame")
#' table(imp_df$.imputation)
#' }
#' @export
get_imputed_data <- function(object, format = c("list", "data.frame")) {

  format <- match.arg(format)

  # --- Resolve to the list of datasets ---
  if (inherits(object, "bace_complete")) {
    datasets <- object$imputed_datasets
  } else if (inherits(object, "bace_final")) {
    datasets <- object$all_datasets
  } else {
    stop("'object' must be a 'bace_complete' object (from bace()) ",
         "or a 'bace_final' object (from bace_final_imp()).")
  }

  if (is.null(datasets) || length(datasets) == 0L) {
    stop("No imputed datasets found in the supplied object.")
  }

  # --- Return in requested format ---
  if (format == "list") {
    return(datasets)
  }

  # Stack into a single data.frame with an .imputation identifier
  for (i in seq_along(datasets)) {
    datasets[[i]][[".imputation"]] <- i
  }
  do.call(rbind, datasets)
}


#' Scaling constants used to standardise each variable
#'
#' @description
#' BACE z-scores continuous variables before fitting (subtract the mean, divide
#' by the SD) because MCMC mixes far better on tidy numbers. Imputed values are
#' returned on your original scale automatically, but the **fitted models** --
#' \code{$pooled_models}, and the per-run models in \code{$final_results} -- are
#' on the standardised scale. Their coefficients are in SD units, not data
#' units.
#'
#' This returns the constants that were actually used, so a coefficient can be
#' put back on the data scale:
#'
#' \deqn{\beta_{data} = \beta_{z} \times sd(y) / sd(x)}
#'
#' @details
#' The constants are recorded **per imputation run and per variable**, because
#' that is how they are computed: every run re-standardises using its own
#' completed data, so the divisor differs from run to run. There is no single
#' constant for the whole fit, and guessing one from the returned datasets is
#' wrong -- measured at -0.6 per cent to +3.3 per cent depending on which guess
#' is made. Non-gaussian variables are never z-scored and are absent.
#'
#' Before 1.0.2 these were discarded once fitting finished, so no exact
#' back-transform was possible from BACE's output at all.
#'
#' @param object A \code{bace_complete} (from \code{\link{bace}}) or a
#'   \code{bace_final} (from \code{\link{bace_final_imp}}).
#' @param variable Optional variable name. When given, returns a data.frame with
#'   one row per run (\code{run}, \code{mean}, \code{sd}); otherwise the full
#'   nested list, keyed run -> model -> variable.
#' @param model Optional: the model (i.e. the response) whose fit you want
#'   \code{variable}'s constants from. Defaults to \code{variable}'s own model.
#'   A coefficient is in units of BOTH sides, so converting the slope of
#'   \code{y} on \code{x1} needs \code{get_scaling(fit, "x1", model = "y")} --
#'   as a predictor, \code{x1} is standardised by its current IMPUTED values,
#'   which differ from its observed-only spread.
#'
#' @return A data.frame when \code{variable} is supplied, else a list of length
#'   \code{n_final}, each a named list of \code{list(mean =, sd =)}.
#' @seealso \code{\link{get_pooled_model}}, \code{\link{get_imputed_data}}
#' @examples
#' \dontrun{
#' fit <- bace("y ~ x1", "~ 1 | Species", phylo = tree, data = d)
#' get_scaling(fit, "y")          # one row per imputation
#'
#' # put the pooled slope back on the data scale
#' sy <- mean(get_scaling(fit, "y")$sd)
#' sx <- mean(get_scaling(fit, "x1")$sd)
#' beta_z <- mean(as.matrix(get_pooled_model(fit, "y")$Sol)[, "x1"])
#' beta_z * sy / sx
#' }
#' @export
get_scaling <- function(object, variable = NULL, model = NULL) {
  sc <- if (inherits(object, "bace_complete")) {
    object[["scaling", exact = TRUE]] %||%
      object$final_results[["scaling", exact = TRUE]]
  } else if (inherits(object, "bace_final")) {
    object[["scaling", exact = TRUE]]
  } else {
    stop("`object` must be a 'bace_complete' or a 'bace_final'.", call. = FALSE)
  }

  if (is.null(sc)) {
    stop("This fit carries no scaling constants. They are recorded from BACE ",
         "1.0.2; a fit made with an earlier version cannot be back-transformed ",
         "exactly, because the divisor each run used was not kept.",
         call. = FALSE)
  }
  if (is.null(variable)) return(sc)

  # sc is run -> model (by response) -> variable -> list(mean, sd).
  # Default to the variable's constants in ITS OWN model; `model` asks for the
  # constants it had as a PREDICTOR in another variable's model, which is what
  # a cross-variable coefficient needs.
  mdl <- if (is.null(model)) variable else model
  models_present <- unique(unlist(lapply(sc, names)))
  if (!mdl %in% models_present) {
    stop("No model recorded for '", mdl, "'. This fit has models for: ",
         paste(models_present, collapse = ", "), ".", call. = FALSE)
  }
  vars_present <- unique(unlist(lapply(sc, function(r) names(r[[mdl]]))))
  if (!variable %in% vars_present) {
    stop("No scaling recorded for '", variable, "' in the model of '", mdl,
         "'. Only z-scored (gaussian) variables have any; that model has: ",
         paste(vars_present, collapse = ", "), ".", call. = FALSE)
  }
  rows <- lapply(seq_along(sc), function(m) {
    e <- sc[[m]][[mdl]][[variable]]
    data.frame(run  = m,
               mean = if (is.null(e$mean)) NA_real_ else e$mean,
               sd   = if (is.null(e$sd))   NA_real_ else e$sd)
  })
  out <- do.call(rbind, rows)
  out[!is.na(out$mean) | !is.na(out$sd), , drop = FALSE]
}
