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
#' @examples \dontrun{
#' # Run a full bace analysis
#' result <- bace(fixformula = "y ~ x1 + x2",
#'                ran_phylo_form = "~1|species",
#'                phylo = tree, data = data)
#'
#' # Extract pooled model for the response variable
#' y_model <- get_pooled_model(result, variable = "y")
#' summary(y_model)
#' plot(y_model)
#'
#' # Extract all pooled models as a named list
#' all_models <- get_pooled_model(result)
#' names(all_models)
#'
#' # Also works with bace_pooled objects directly
#' pooled <- pool_posteriors(final_results)
#' y_model <- get_pooled_model(pooled, variable = "y")
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
#' @examples \dontrun{
#' # Run a full bace analysis
#' result <- bace(fixformula = "y ~ x1 + x2",
#'                ran_phylo_form = "~1|species",
#'                phylo = tree, data = data)
#'
#' # Get imputed datasets as a list
#' imp_list <- get_imputed_data(result)
#' length(imp_list)       # number of imputations
#' head(imp_list[[1]])    # first imputed dataset
#'
#' # Get as a single stacked data frame
#' imp_df <- get_imputed_data(result, format = "data.frame")
#' table(imp_df$.imputation)
#'
#' # Also works with bace_final objects directly
#' final <- bace_final_imp(bace_obj, ...)
#' imp_list <- get_imputed_data(final)
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
