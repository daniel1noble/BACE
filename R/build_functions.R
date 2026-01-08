
#' @title .build_formula
#' @description Function takes a string specifying a formula and converts this to a formula object to be used in the models.
#' @param x A character string specifying the formula used in the model. This should be of the form: y ~ x.
#' @return A formula object
#' #' @examples \dontrun{
#' .build_formula("y ~ x1 + x2")
#' }
#' @export

.build_formula <- function(x) {

	# Some checks. Formula must have a '`~`' in it
	if(!grepl("~", as.character(x))) { stop("Make sure you specify a formula string. This involves a structure of the form: y ~ x") } 

  return(stats::as.formula(x))
}

#' @title .build_formula_string
#' @description Function takes a string specifying a formula and creates all combinations of models using the variables.
#' @param x A character string specifying the formula used in the model. This should be of the form: y ~ x.
#' @return A list of formulas
#' #' @examples \dontrun{
#' .build_formula_string("y ~ x1 + x2 + x3")
#' .build_formula_string("y ~ x1 + x2*x3")
#' .build_formula_string("y ~ x1 + x2:x3")
#' .build_formula_string("y ~ x1 + x2 + x3 + x1:x2 + x1:x3 + x2:x3 + x1:x2:x3")
#' }
#' @export
.build_formula_string <- function(x) {

  f <- if (inherits(x, "formula")) x else stats::as.formula(x)

  if (length(f) < 3L) {
    stop("Make sure you specify a formula of the form: y ~ ...")
  }

  lhs <- f[[2]]
  rhs <- f[[3]]

  # all variables appearing anywhere in the formula
  vars <- unique(all.vars(f))

  # keep original as the first element
  out <- vector("list", length(vars))
  names(out) <- vars

  for (i in seq_along(vars)) {
    v <- vars[[i]]

    if (identical(v, as.character(lhs))) {
      # original formula unchanged
      out[[i]] <- f
    } else {
      # swap v with original lhs inside the RHS expression
      rhs2 <- .swap_symbols(rhs, a = v, b = as.character(lhs))

      # new formula: v ~ (modified RHS)
      out[[i]] <- stats::as.formula(call("~", as.name(v), rhs2), env = environment(f))
    }
  }

  out
}

#' @title .replace_symbols
#' @description Recursively replace symbol names in an expression
#' @param expr An expression
#' @param map A named list specifying the replacements to be made
#' @return An expression with the symbol names replaced
.replace_symbols <- function(expr, map) {
  if (is.null(expr)) return(expr)

  # symbol / name
  if (is.name(expr)) {
    nm <- as.character(expr)
    if (nm %in% names(map)) return(as.name(map[[nm]]))
    return(expr)
  }

  # atomic (numeric/character/etc.)
  if (!is.language(expr)) return(expr)

  # call / language: recurse over arguments
  as.call(lapply(as.list(expr), .replace_symbols, map = map))
}

#' @title .swap_symbols
#' @description Swap two symbol names in an expression
#' @param expr An expression
#' @param a The first symbol name to swap
#' @param b The second symbol name to swap
#' @param tmp A temporary symbol name used during the swap
#' @return An expression with the two symbol names swapped
.swap_symbols <- function(expr, a, b, tmp = ".TMP_SWAP_SYMBOL") {
  if (identical(a, b)) return(expr)

  # a -> tmp, b -> a, tmp -> b
  expr <- .replace_symbols(expr, setNames(tmp, a))
  expr <- .replace_symbols(expr, setNames(a, b))
  expr <- .replace_symbols(expr, setNames(b, tmp))
  expr
}

#' @title .build_formula_string_random
#' @description Function takes a string specifying a random effects formula and converts this to a formula object to be used in the models.
#' @param ran_phylo_form A character string specifying the random effects and phylogenetic structure formula used in the model.
#' @return A formula	 object
#' #' @examples \dontrun{
#' .build_formula_string_random("~ 1 | Species")
#' }
#' @export	
.build_formula_string_random <- function(ran_phylo_form) {

  f <- if (inherits(ran_phylo_form, "formula")) {
    ran_phylo_form
  } else {
    stats::as.formula(ran_phylo_form)
  }

  # Extract RHS as character
  rhs <- as.character(f)[length(as.character(f))]

  # Remove parentheses
  rhs <- gsub("[()]", "", rhs)

  # Split on pipe with optional whitespace
  parts <- strsplit(rhs, "\\s*\\|\\s*", perl = TRUE)[[1]]

  if (length(parts) < 2L) {
    stop("No random-effects structure detected (missing '|').")
  }

  cluster <- trimws(parts[2])

  # Return formula of the form: ~ Species
  stats::as.formula(paste("~", cluster))
}

