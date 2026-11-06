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
#' here, the phylogeny is passed through to the fits (e.g. for a phylogenetic
#' GLS). A single \code{phylo} is passed to every fit. A \code{multiPhylo} (or
#' plain list of \code{phylo}) assigns tree \code{m} to imputation \code{m},
#' recycled in balanced rotation when there are fewer trees than imputations --
#' the same deterministic rule \code{\link{bace_final_imp}} uses, so a
#' downstream phylogenetic model can condition on the same tree that generated
#' each imputation (tree-per-imputation workflow). If \code{.f}'s own
#' \code{tree} argument expects a whole tree set, wrap it instead:
#' \code{.f = function(d) g(d, my_trees)}.
#'
#' @param object A \code{bace_complete} (from \code{\link{bace}}), a
#'   \code{bace_final} (from \code{\link{bace_final_imp}}), or a plain list of
#'   completed data.frames.
#' @param .f A function \code{function(dataset, ...)} that fits a model to one
#'   completed data.frame and returns the fitted model.
#' @param ... Additional arguments passed to \code{.f} for every imputation.
#' @param tree Optional phylogeny passed to \code{.f} when \code{.f} declares a
#'   \code{tree} argument: a single \code{phylo} (same tree for every fit) or a
#'   \code{multiPhylo}/list of \code{phylo} (tree \code{m} for imputation
#'   \code{m}, recycled if shorter). Default \code{NULL}.
#' @param .progress Logical; show a progress line (default: interactive only).
#' @param .on_error One of \code{"continue"} (default) or \code{"stop"}. With
#'   \code{"continue"}, a fit that errors is captured and the loop proceeds;
#'   \code{\link{pool_mi}} drops the failures.
#'
#' @return A list of length \code{M} with class \code{"bace_mi_fits"}.
#' @seealso \code{\link{pool_mi}}, \code{\link{pool_posteriors}}
#' @examples
#' # with_imputations() also accepts a plain list of completed data.frames,
#' # so the fit-and-pool steps can be shown without running MCMC. In the
#' # full workflow, `object` is the result of bace() or bace_final_imp().
#' set.seed(1)
#' imps <- lapply(1:3, function(i) {
#'   d <- data.frame(x1 = rnorm(30))
#'   d$y <- 1 + 0.5 * d$x1 + rnorm(30, sd = 0.2)
#'   d
#' })
#' fits <- with_imputations(imps, function(d) lm(y ~ x1, data = d),
#'                          .progress = FALSE)
#' fits
#' pool_mi(fits)
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

  # Per-imputation trees: a multiPhylo or bare list of phylo assigns tree m to
  # imputation m (balanced recycling, mirroring bace_final_imp's rule); a
  # single phylo keeps the same-tree-for-all behaviour.
  tree_list <- NULL
  if (!is.null(tree) && !inherits(tree, "phylo")) {
    tree_list <- .as_tree_list(tree, context = "with_imputations")
    gen_hashes <- if (inherits(object, "bace_complete")) {
      object$final_results$tree_hashes
    } else if (inherits(object, "bace_final")) {
      object$tree_hashes
    } else NULL
    if (!is.null(gen_hashes) && length(gen_hashes) > 1L &&
        length(tree_list) != length(gen_hashes)) {
      message("with_imputations: ", length(tree_list), " trees supplied, but ",
              "the imputations were generated under ", length(gen_hashes),
              " trees (see $tree_index on the bace_final object). Pass the ",
              "same tree set in the same order for each fit to condition on ",
              "the tree that generated its imputation.")
    }
  }

  for (i in seq_len(M)) {
    if (isTRUE(.progress)) {
      message(sprintf("\r  fitting imputation %d/%d", i, M), appendLF = FALSE)
    }
    call_args <- c(list(datasets[[i]]), dots)
    if (!is.null(tree) && "tree" %in% f_formals && !("tree" %in% names(dots))) {
      call_args[["tree"]] <- if (is.null(tree_list)) tree else
        tree_list[[((i - 1L) %% length(tree_list)) + 1L]]
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
#' @noRd
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
