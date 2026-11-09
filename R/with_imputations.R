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
#' plain list of \code{phylo}) gives one tree per imputation
#' (tree-per-imputation workflow).
#'
#' Which tree goes to which fit is \emph{read from the imputation object}: a
#' \code{bace_final} / \code{bace_complete} records the tree each run was
#' imputed under in \code{$tree_index}, and when the trees supplied here are
#' the same set in the same order (checked by content hash, not by count) that
#' recorded pairing is used exactly, so every downstream model conditions on
#' the tree that generated its own imputation. Otherwise BACE falls back to
#' \code{bace_final_imp}'s own assignment rule and says so with a message.
#' Supplying a plain list of data.frames records nothing, so the fallback rule
#' applies.
#'
#' If \code{.f}'s own \code{tree} argument expects a whole tree set, wrap it
#' instead: \code{.f = function(d) g(d, my_trees)}.
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

  # A tree that cannot be delivered must never be dropped in silence. Before
  # 1.0.1 the tree was forwarded only when `.f` declared a formal literally
  # named `tree`; a user who wrote function(d, phy = my_tree) got their
  # DEFAULT tree in every fit, so every model conditioned on the wrong
  # phylogeny with no error and no warning. Measured on one such pair: the
  # reported coefficient was 0.4407 with the tree dropped against 0.9785 with
  # it passed. Nothing downstream looked wrong. (GitHub issue #18.)
  if (!is.null(tree) && !("tree" %in% f_formals) && !("..." %in% f_formals) &&
      !("tree" %in% names(dots))) {
    looks_like_tree <- intersect(f_formals, c("phy", "phylo", "phy_tree",
                                              "phylogeny", "treeobj"))
    stop("`tree` was supplied but `.f` cannot receive it: `.f` has no `tree` ",
         "argument and no `...`.\n",
         "  `.f` takes: ", paste(f_formals, collapse = ", "), "\n",
         if (length(looks_like_tree))
           paste0("  Rename `", looks_like_tree[1], "` to `tree`, or ")
         else "  Name the argument `tree`, or ",
         "wrap it: .f = function(d) your_fun(d, my_tree).\n",
         "  Passing it anyway would fit every imputation against the wrong ",
         "phylogeny without warning.",
         call. = FALSE)
  }

  # Per-imputation trees: a multiPhylo or bare list of phylo assigns one tree
  # per imputation; a single phylo keeps the same-tree-for-all behaviour.
  #
  # The pairing is READ from the imputation object whenever the object records
  # it ($tree_index), never re-derived here. Before 1.0.1 it was re-derived as
  # plain modular recycling, which disagrees with the rule bace_final_imp
  # actually uses (.select_tree_indices) whenever MORE trees were supplied than
  # there are imputations: with 100 trees and 50 imputations the generator used
  # trees 1, 3, 5, ..., 99 while the fits used 1, 2, 3, ..., 50, so 49 of the 50
  # models conditioned on a tree that had not generated their data. The old
  # guard compared only the tree COUNTS, so passing the same 100 trees -- the
  # correct thing to do -- matched on count and said nothing at all.
  tree_list    <- NULL
  tree_pairing <- NULL
  if (!is.null(tree) && !inherits(tree, "phylo")) {
    tree_list <- .as_tree_list(tree, context = "with_imputations")
    gen <- if (inherits(object, "bace_complete")) object$final_results
           else if (inherits(object, "bace_final")) object
           else NULL
    gen_hashes <- if (is.list(gen)) gen[["tree_hashes", exact = TRUE]] else NULL
    gen_index  <- if (is.list(gen)) gen[["tree_index",  exact = TRUE]] else NULL

    # Same trees, same order? Compare contents, not counts.
    same_set <- !is.null(gen_hashes) &&
      length(tree_list) == length(gen_hashes) &&
      identical(as.character(vapply(tree_list, .tree_cache_key, character(1L))),
                as.character(gen_hashes))

    if (same_set && is.numeric(gen_index) && length(gen_index) == M &&
        all(gen_index >= 1L & gen_index <= length(tree_list))) {
      tree_pairing <- as.integer(gen_index)        # exact recorded pairing
    } else {
      tree_pairing <- .select_tree_indices(length(tree_list), M)
      if (!is.null(gen_hashes) && length(gen_hashes) > 1L) {
        message("with_imputations: the ", length(tree_list), " tree(s) supplied ",
                "are not the ", length(gen_hashes), " tree(s) the imputations ",
                "were generated under (same order, same trees), so no fit can ",
                "be matched to the tree that generated its data. Trees are ",
                "being assigned by bace_final_imp's own rule instead. Pass the ",
                "same tree set, in the same order, to pair them exactly (see ",
                "$tree_index on the bace_final object).")
      }
    }
  }

  for (i in seq_len(M)) {
    if (isTRUE(.progress)) {
      message(sprintf("\r  fitting imputation %d/%d", i, M), appendLF = FALSE)
    }
    call_args <- c(list(datasets[[i]]), dots)
    if (!is.null(tree) &&
        ("tree" %in% f_formals || "..." %in% f_formals) &&
        !("tree" %in% names(dots))) {
      call_args[["tree"]] <- if (is.null(tree_list)) tree else
        tree_list[[tree_pairing[i]]]
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
