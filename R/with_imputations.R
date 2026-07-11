#' Fit a downstream model on every imputed dataset
#'
#' Apply a user-supplied model-fitting function to each of the imputed
#' datasets produced by \code{\link{bace}} / \code{\link{bace_final_imp}} and
#' return the list of fits. This is the middle step of the multiple-imputation
#' workflow \code{bace()} -> \code{with_imputations()} -> \code{\link{pool_mi}}
#' (Rubin's rules) or \code{\link{pool_posteriors}} (Bayesian stacking).
#'
#' @details
#' \code{with_imputations()} is model-agnostic. \code{.f} may fit any model:
#' frequentist (\code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#' \code{nlme::gls}, \code{phylolm::phylolm}, ...) for pooling with
#' \code{\link{pool_mi}}, or \code{MCMCglmm::MCMCglmm} for either
#' \code{\link{pool_mi}} (Rubin on posterior summaries) or
#' \code{\link{pool_posteriors}} (exact stacked posterior). For frequentist
#' fits, \code{\link[stats]{coef}} and \code{\link[stats]{vcov}} must work on
#' the return value.
#'
#' If \code{.f} declares a \code{tree} argument and \code{tree} is supplied
#' here, the phylogeny is passed through to every fit (e.g. for a phylogenetic
#' GLS). This also leaves room for a future posterior-tree workflow without an
#' API change.
#'
#' @param object A \code{bace_complete} (from \code{\link{bace}}), a
#'   \code{bace_final} (from \code{\link{bace_final_imp}}), or a plain list of
#'   completed data.frames.
#' @param .f A function \code{function(dataset, ...)} that fits a model to one
#'   completed data.frame and returns the fitted model.
#' @param ... Additional arguments passed to \code{.f} for every imputation.
#' @param tree Optional phylogeny (class \code{phylo}) passed to \code{.f} when
#'   \code{.f} declares a \code{tree} argument. Default \code{NULL}.
#' @param .progress Logical; show a progress line (default: interactive only).
#' @param .on_error One of \code{"continue"} (default) or \code{"stop"}. With
#'   \code{"continue"}, a fit that errors is captured and the loop proceeds;
#'   \code{\link{pool_mi}} drops the failures.
#'
#' @return A list of length \code{M} with class \code{"bace_mi_fits"}.
#' @seealso \code{\link{pool_mi}}, \code{\link{pool_posteriors}}
#' @examples \dontrun{
#' res  <- bace("y ~ x1 + x2", "~1|Species", tree, dat, n_final = 50)
#' fits <- with_imputations(res, function(d) lm(y ~ x1 + x2, data = d))
#' pool_mi(fits)
#' }
#' @export
with_imputations <- function(object, .f, ..., tree = NULL,
                             .progress = interactive(),
                             .on_error = c("continue", "stop")) {

  .on_error <- match.arg(.on_error)
  if (!is.function(.f)) {
    stop("`.f` must be a function of one argument (a data.frame).", call. = FALSE)
  }

  datasets <- .extract_imputed_datasets(object)
  M <- length(datasets)
  if (M < 2L) {
    stop("Need at least 2 imputed datasets to pool later; got ", M, ".",
         call. = FALSE)
  }

  dots      <- list(...)
  f_formals <- names(formals(.f))
  fits      <- vector("list", M)
  failures  <- integer(0)

  for (i in seq_len(M)) {
    if (isTRUE(.progress)) {
      message(sprintf("\r  fitting imputation %d/%d", i, M), appendLF = FALSE)
    }
    call_args <- c(list(datasets[[i]]), dots)
    if (!is.null(tree) && "tree" %in% f_formals && !("tree" %in% names(dots))) {
      call_args[["tree"]] <- tree
    }
    res <- tryCatch(
      do.call(.f, call_args),
      error = function(e) structure(
        list(index = i, message = conditionMessage(e)),
        class = c("bace_mi_error", "condition"))
    )
    if (inherits(res, "bace_mi_error")) {
      failures <- c(failures, i)
      if (.on_error == "stop") {
        if (isTRUE(.progress)) message("")
        stop("`.f` failed on imputation ", i, ": ", res$message, call. = FALSE)
      }
    }
    fits[[i]] <- res
  }
  if (isTRUE(.progress)) message("")

  if (length(failures) > 0L) {
    warning(sprintf(
      "%d of %d fits failed (indices: %s). They will be dropped by pool_mi().",
      length(failures), M, paste(failures, collapse = ", ")), call. = FALSE)
  }

  first_ok <- NULL
  for (i in seq_along(fits)) {
    if (!inherits(fits[[i]], "bace_mi_error")) { first_ok <- fits[[i]]; break }
  }
  attr(fits, "n_fits")      <- M
  attr(fits, "n_failed")    <- length(failures)
  attr(fits, "failed")      <- failures
  attr(fits, "is_mcmcglmm") <- !is.null(first_ok) && inherits(first_ok, "MCMCglmm")
  class(fits) <- c("bace_mi_fits", "list")
  fits
}

#' Extract the list of completed datasets from a BACE object
#' @keywords internal
.extract_imputed_datasets <- function(object) {
  if (inherits(object, "bace_complete")) return(object$imputed_datasets)
  if (inherits(object, "bace_final"))    return(object$all_datasets)
  if (is.list(object) && length(object) > 0L &&
      all(vapply(object, is.data.frame, logical(1)))) {
    return(object)
  }
  stop("`object` must be a 'bace_complete', a 'bace_final', or a list of ",
       "data.frames.", call. = FALSE)
}

#' @export
print.bace_mi_fits <- function(x, ...) {
  n  <- attr(x, "n_fits") %||% length(x)
  nf <- attr(x, "n_failed") %||% 0L
  cat(sprintf("BACE multiple-imputation fits: %d/%d successful\n", n - nf, n))
  if (nf > 0L) {
    cat(sprintf("  failed at indices: %s\n",
                paste(attr(x, "failed"), collapse = ", ")))
  }
  first_ok <- NULL
  for (i in seq_along(x)) {
    if (!inherits(x[[i]], "bace_mi_error")) { first_ok <- x[[i]]; break }
  }
  if (!is.null(first_ok)) {
    cat(sprintf("  model class: %s\n", paste(class(first_ok), collapse = " / ")))
  }
  if (isTRUE(attr(x, "is_mcmcglmm"))) {
    cat("\n  Pool with: pool_mi(fits)  (Rubin, gives fmi)  OR\n")
    cat("             pool_posteriors() on the bace_final for the exact stacked posterior\n")
  } else {
    cat("\n  Pool with: pool_mi(fits)\n")
  }
  invisible(x)
}

# Local null-coalesce (avoids a hard rlang dependency at load).
`%||%` <- function(a, b) if (is.null(a)) b else a
