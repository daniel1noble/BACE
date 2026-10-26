#' @title .validate_priors
#' @description Shallow structural validation of the user-facing \code{priors}
#'   argument (named list: one prior specification per imputed variable). Run
#'   once per phase BEFORE any MCMC so misspecifications fail fast. Each
#'   element may contain any of the MCMCglmm prior blocks \code{G}, \code{R},
#'   \code{B}; a supplied block replaces the corresponding default block
#'   wholesale (no element-wise merging - MCMCglmm prior components are
#'   internally coupled, and partial merges produce undiagnosable errors).
#'   Overriding \code{R} is refused for threshold/ordinal/categorical
#'   variables: their residual variance is fixed on the latent probit scale
#'   and BACE's category-probability computations hard-code that scale
#'   (unit residual for threshold; the (I+J)/K rescaling for categorical),
#'   so a custom R would silently corrupt every probability.
#' @param priors NULL, or a named list keyed by variable name.
#' @param types Named list of detected variable types (all modelled variables).
#' @param miss_vars Character vector of variables that actually have missing
#'   data (only these get imputation models); priors for other variables
#'   trigger a warning, not an error.
#' @return invisible(NULL); called for its side effect (errors/warnings).
#' @noRd
.validate_priors <- function(priors, types, miss_vars = NULL) {
  if (is.null(priors)) return(invisible(NULL))

  if (!is.list(priors) || is.null(names(priors)) || any(names(priors) == "")) {
    stop("'priors' must be a named list with one element per variable, ",
         "e.g. priors = list(y = list(G = ...), x = list(G = ..., B = ...)).",
         call. = FALSE)
  }
  if (anyDuplicated(names(priors))) {
    stop("'priors' has duplicated names: ",
         paste(unique(names(priors)[duplicated(names(priors))]), collapse = ", "),
         call. = FALSE)
  }

  unknown <- setdiff(names(priors), names(types))
  if (length(unknown) > 0) {
    stop("'priors' names variable(s) not in the model formulas: ",
         paste(unknown, collapse = ", "),
         ". Modelled variables are: ", paste(names(types), collapse = ", "), ".",
         call. = FALSE)
  }

  if (!is.null(miss_vars)) {
    unused <- setdiff(names(priors), miss_vars)
    if (length(unused) > 0) {
      warning("'priors' supplied for variable(s) with no missing data (no ",
              "imputation model is fitted for them): ",
              paste(unused, collapse = ", "), ". These priors are ignored.",
              call. = FALSE)
    }
  }

  for (v in names(priors)) {
    el <- priors[[v]]
    if (!is.list(el) || length(el) == 0) {
      stop("priors$", v, " must be a non-empty list containing one or more of ",
           "the MCMCglmm prior blocks 'G', 'R', 'B'.", call. = FALSE)
    }
    bad <- setdiff(names(el), c("G", "R", "B"))
    if (length(bad) > 0 || is.null(names(el)) || any(names(el) == "")) {
      stop("priors$", v, " may only contain named blocks 'G', 'R', 'B' ",
           "(found: ", paste(if (length(bad)) bad else "<unnamed>",
                             collapse = ", "), "). ",
           "Each supplied block replaces the default block wholesale.",
           call. = FALSE)
    }
    if ("R" %in% names(el) &&
        !is.null(types[[v]]) &&
        types[[v]] %in% c("threshold", "ordinal", "categorical")) {
      stop("priors$", v, ": the residual (R) prior cannot be changed for ",
           types[[v]], " variables. Their residual variance is fixed at 1 on ",
           "the latent probit scale, and BACE's category-probability ",
           "computations hard-code that scale - a custom R prior would ",
           "silently change every predicted probability. Customise the G ",
           "(random-effect) and/or B (fixed-effect) blocks instead.",
           call. = FALSE)
    }
  }
  invisible(NULL)
}


#' @title .resolve_prior
#' @description Merge a user prior override into the default prior built by
#'   \code{.make_prior} for one variable: each supplied top-level block
#'   (\code{G}, \code{R}, \code{B}) replaces the default block wholesale,
#'   after dimension checks against the default (which encodes what MCMCglmm
#'   will accept for this variable's family, number of random effects and,
#'   for categorical, J = n_levels - 1).
#' @param default Prior list from .make_prior.
#' @param override NULL (returns default) or one element of the validated
#'   user 'priors' list.
#' @param type Variable type (for error messages / R-ban belt-and-braces).
#' @param var Variable name (for error messages).
#' @return The merged prior list.
#' @noRd
.resolve_prior <- function(default, override, type, var) {
  if (is.null(override)) return(default)

  out <- default

  if (!is.null(override$R)) {
    if (type %in% c("threshold", "ordinal", "categorical")) {
      stop("priors$", var, ": R prior override is not allowed for ", type,
           " variables (fixed latent scale).", call. = FALSE)
    }
    if (!is.list(override$R) || is.null(override$R$V)) {
      stop("priors$", var, "$R must be a list with at least element V ",
           "(e.g. list(V = 1, nu = 0.002)).", call. = FALSE)
    }
    out$R <- override$R
  }

  if (!is.null(override$G)) {
    G <- override$G
    # Accept a single G element supplied bare (list(V = ..., nu = ...))
    # by wrapping it, but only when one random effect is expected.
    if (!is.null(G$V) && length(default$G) == 1L) G <- list(G)
    if (!is.list(G) || length(G) != length(default$G)) {
      stop("priors$", var, "$G must be a list of ", length(default$G),
           " element(s) (one per random effect: phylogenetic",
           if (length(default$G) == 2L) " + species" else "",
           "), each a list like list(V = ..., nu = ...).", call. = FALSE)
    }
    exp_dim <- nrow(as.matrix(default$G[[1]]$V))
    for (g in seq_along(G)) {
      if (!is.list(G[[g]]) || is.null(G[[g]]$V)) {
        stop("priors$", var, "$G[[", g, "]] must be a list with at least ",
             "element V.", call. = FALSE)
      }
      got_dim <- nrow(as.matrix(G[[g]]$V))
      if (got_dim != exp_dim) {
        stop("priors$", var, "$G[[", g, "]]$V has dimension ", got_dim,
             " but this ", type, " variable requires dimension ", exp_dim,
             if (type == "categorical")
               " (J = n_levels - 1 latent traits)" else "",
             ".", call. = FALSE)
      }
    }
    out$G <- G
  }

  if (!is.null(override$B)) {
    B <- override$B
    if (!is.list(B) || is.null(B$mu) || is.null(B$V)) {
      stop("priors$", var, "$B must be a list with elements mu and V.",
           call. = FALSE)
    }
    nB <- length(B$mu)
    VB <- as.matrix(B$V)
    if (nrow(VB) != nB || ncol(VB) != nB) {
      stop("priors$", var, "$B: length(mu) = ", nB, " but V is ",
           nrow(VB), " x ", ncol(VB), "; they must agree.", call. = FALSE)
    }
    if (!is.null(default$B)) {
      nExp <- length(default$B$mu)
      if (nB != nExp) {
        stop("priors$", var, "$B has ", nB, " fixed-effect entries but this ",
             "model expects ", nExp,
             " (categorical models expand to J latent traits; count ",
             "intercept + slopes per trait).", call. = FALSE)
      }
    }
    out$B <- B
  }

  out
}
