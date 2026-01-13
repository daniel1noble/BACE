#' Get or set BACE options
#'
#' @param ... Named options to set. If empty, returns current options.
#' @param .reset Logical; if TRUE, reset all mypkg.* options (in this session) to defaults.
#' @return A named list of effective option values (after defaults + overrides).
#' @examples
#' \dontrun{
#' # get current options
#' bace_options()
#' #' # set verbose to TRUE
#' bace_options(verbose = TRUE)
#' bace_options()
#' # reset all options to defaults
#' bace_options(.reset = TRUE)
#' # user overriden option
#' getOption("BACE.verbose")
#' }
#' @export
bace_options <- function(..., .reset = FALSE) {
  defaults <- BACE:::bace_option_defaults()

  if (.reset) {
    # remove any user overrides; do NOT set defaults into options()
    to_drop <- paste0("BACE.", names(defaults))
    op <- options()
    existing <- intersect(names(op), to_drop)
    if (length(existing)) options(structure(rep(list(NULL), length(existing)), names = existing))
    return(bace_options())
  }

  dots <- list(...)
  if (length(dots) == 0L) {
    # effective options = user overrides in options() + defaults
    user <- BACE:::bace_option_get_raw(names(defaults))
    eff <- utils::modifyList(defaults, user)
    return(eff)
  }

  # set path
  BACE:::bace_option_set(dots, defaults = defaults)
  invisible(bace_options())
}


#' default option values
#' @description Internal function returning default option values. Currently:
#' \itemize{
#'   \item{verbose: FALSE}
#'   \item{digits: 3L}
#'   \item{gelman: 1}
#' }
#' @return named list of default option values
#' @export
bace_option_defaults <- function() {
  list(
    verbose = TRUE,
    digits = 3L,
    gelman = 1
  )
}

#' internal: get raw user-set options from options()
#' @param keys character vector of option names (without "BACE." prefix)
#' @return named list of option values set by user (may be empty)
#' @noRd
bace_option_get_raw <- function(keys) {
  pref <- paste0("BACE.", keys)
  vals <- lapply(pref, getOption, default = NULL)
  names(vals) <- keys
  # drop NULLs (not set by user)
  vals[!vapply(vals, is.null, logical(1))]
}

#' internal: set options after validation
#' @param dots named list of options to set
#' @param defaults named list of default option values
#' @noRd
bace_option_set <- function(dots, defaults) {
  if (is.null(names(dots)) || any(names(dots) == "")) {
    stop("All options must be named, e.g. bace_options(verbose = TRUE).", call. = FALSE)
  }

  # allow short names (verbose) but store as bace_verbose
  known <- names(defaults)
  unknown <- setdiff(names(dots), known)
  if (length(unknown)) {
    stop("Unknown BACE option(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  validated <- Map(function(key, val) BACE:::bace_option_validate(key, val),
                   names(dots), dots)
  opt_names <- paste0("BACE.", names(validated))
  options(structure(validated, names = opt_names))
  invisible(NULL)
}

#' internal: validate option values
#' @param key option name (without "BACE." prefix)
#' @param val option value
#' @return validated option value
#' @noRd
bace_option_validate <- function(key, val) {
  switch(key,
    verbose = {
      if (!is.logical(val) || length(val) != 1) stop("verbose must be a single TRUE/FALSE.", call. = FALSE)
      val
    },
    digits = {
      if (!is.numeric(val) || length(val) != 1) stop("digits must be length-1 numeric.", call. = FALSE)
      as.integer(val)
    },
    gelman = {
      if (!is.numeric(val) || length(val) != 1) stop("gelman must be a length-1 number.", call. = FALSE)
      allowed <- c(0, 1, 2)
      if (!val %in% allowed) stop("gelman must be one of: ", paste(allowed, collapse = ", "), call. = FALSE)
      val
    },
    stop("No validator defined for option '", key, "'.", call. = FALSE)
  )
}
