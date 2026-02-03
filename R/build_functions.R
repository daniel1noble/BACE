
#' @title .build_formula
#' @description Function takes a string specifying a formula and converts this to a formula object to be used in the models.
#' @param x A character string specifying the formula used in the model. This should be of the form: y ~ x.
#' @return A formula object
#' @examples
#' \dontrun{
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
#' @examples
#' \dontrun{
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
#' @description Function takes a string specifying a random effects formula and converts this to a formula object to be used in MCMCglmm models. Supports both lme4-style (~ 1 | cluster) and MCMCglmm-style (~us(1):cluster) formulas, including multiple random effects.
#' @param ran_phylo_form A character string or formula specifying the random effects and phylogenetic structure formula used in the model.
#' @param return_list A logical indicating whether to return a detailed list (TRUE) or just the formula object (FALSE, default). For backward compatibility, the default is FALSE.
#' @return If return_list = FALSE (default): A formula object for use in MCMCglmm.
#'   If return_list = TRUE: A list containing:
#'   \itemize{
#'     \item formula: The random effects formula for MCMCglmm
#'     \item terms: A list of parsed random effect terms, each containing structure type, variables, and cluster
#'     \item clusters: Character vector of cluster/grouping variables
#'   }
#' @examples
#' \dontrun{
#' # lme4-style formula (converted to MCMCglmm format)
#' .build_formula_string_random("~ 1 | Species")
#' 
#' # MCMCglmm-style formulas
#' .build_formula_string_random("~ us(1):Species")
#' .build_formula_string_random("~ us(1):Species + us(1):Genus")
#' .build_formula_string_random("~ idh(trait):Species")
#' .build_formula_string_random("~ us(trait1 + trait2):Species")
#' 
#' # Get detailed information about the random effects structure
#' .build_formula_string_random("~ us(1):Species + us(1):Genus", return_list = TRUE)
#' }
#' @export	
.build_formula_string_random <- function(ran_phylo_form, return_list = FALSE) {

  f <- if (inherits(ran_phylo_form, "formula")) {
    ran_phylo_form
  } else {
    stats::as.formula(ran_phylo_form)
  }

  # Extract RHS as character
  rhs <- as.character(f)[length(as.character(f))]
  rhs <- trimws(rhs)

  # Check if this is lme4-style (contains '|') or MCMCglmm-style (contains ':')
  is_lme4_style <- grepl("\\|", rhs) && !grepl(":", gsub("\\|.*", "", rhs))

  if (is_lme4_style) {
    # Convert lme4-style to MCMCglmm-style
    # Handle multiple random effects: (1|Species) + (1|Genus) or 1|Species + 1|Genus
    # First split on '+' to get individual terms
    terms_str <- strsplit(rhs, "\\s*\\+\\s*")[[1]]
    
    parsed_terms <- lapply(terms_str, function(term) {
      # Remove parentheses
      term_clean <- gsub("[()]", "", trimws(term))
      
      # Split on pipe with optional whitespace
      parts <- strsplit(term_clean, "\\s*\\|\\s*", perl = TRUE)[[1]]
      
      if (length(parts) < 2L) {
        stop(paste0("No random-effects structure detected in term '", term, "' (missing '|')."))
      }
      
      cluster <- trimws(parts[2])
      vars <- trimws(parts[1])
      
      # Convert to MCMCglmm format: us(1):cluster or us(vars):cluster
      if (vars == "1") {
        mcmc_rhs <- paste0("us(1):", cluster)
      } else {
        mcmc_rhs <- paste0("us(", vars, "):", cluster)
      }
      
      list(
        structure = "us",
        variables = if(vars == "1") "1" else strsplit(gsub("\\s*\\+\\s*|\\s*\\*\\s*", "+", vars), "\\+")[[1]],
        cluster = cluster,
        original = mcmc_rhs
      )
    })
    
    # Combine all terms into a single MCMCglmm formula
    mcmc_rhs_all <- sapply(parsed_terms, function(x) x$original)
    mcmc_formula <- stats::as.formula(paste("~", paste(mcmc_rhs_all, collapse = " + ")))
    
  } else {
    # MCMCglmm-style formula
    # Split on '+' to get individual random effects terms
    terms_str <- strsplit(rhs, "\\s*\\+\\s*")[[1]]
    
    parsed_terms <- lapply(terms_str, function(term) {
      term <- trimws(term)
      
      # Match pattern: structure(variables):cluster
      # e.g., us(1):Species, idh(trait1+trait2):Species, us(1):animal
      pattern <- "^([a-z]+)\\(([^)]+)\\):([a-zA-Z0-9_\\.]+)$"
      
      if (grepl(pattern, term)) {
        matches <- regmatches(term, regexec(pattern, term))[[1]]
        structure_type <- matches[2]
        vars_str <- matches[3]
        cluster_var <- matches[4]
        
        # Parse variables (may be "1", or "trait1 + trait2", etc.)
        if (vars_str == "1") {
          variables <- "1"
        } else {
          variables <- trimws(strsplit(vars_str, "\\s*\\+\\s*")[[1]])
        }
        
        list(
          structure = structure_type,
          variables = variables,
          cluster = cluster_var,
          original = term
        )
      } else {
        # Fallback: try to extract cluster after colon
        if (grepl(":", term)) {
          parts <- strsplit(term, ":")[[1]]
          cluster_var <- trimws(parts[length(parts)])
          
          list(
            structure = "us",  # default
            variables = "1",
            cluster = cluster_var,
            original = term
          )
        } else {
          stop(paste("Unable to parse random effects term:", term, 
                     "\nExpected format: structure(variables):cluster (e.g., us(1):Species)"))
        }
      }
    })
    
    mcmc_formula <- f
  }
  
  # Extract all cluster variables
  clusters <- sapply(parsed_terms, function(x) x$cluster)
  
  # Return either just the formula (default, for backward compatibility) or the full list
  if (return_list) {
    return(list(
      formula = mcmc_formula,
      terms = parsed_terms,
      clusters = unique(clusters)
    ))
  } else {
    return(mcmc_formula)
  }
}

#' @title .list_of_G
#' @description Helper function to create a named list of inverse matrices for use in MCMCglmm's ginverse argument. This is useful when you have multiple random effects that require inverse relationship matrices (e.g., phylogenetic effects for different taxonomic levels).
#' @param ... Named inverse matrices. Each argument should be named with the cluster variable name and contain the inverse matrix for that cluster.
#' @return A named list suitable for use as the ginverse argument in MCMCglmm::MCMCglmm()
#' @examples
#' \dontrun{
#' library(ape)
#' library(MCMCglmm)
#' 
#' # Create phylogenetic tree
#' tree <- rtree(10)
#' A <- inverseA(tree, nodes = "ALL")$Ainv
#' 
#' # Single random effect (phylogeny)
#' ginv <- .list_of_G(Species = A)
#' 
#' # Multiple random effects (e.g., Species and Genus phylogenies)
#' # You would need separate trees/matrices for each level
#' genus_tree <- rtree(5)
#' A_genus <- inverseA(genus_tree, nodes = "ALL")$Ainv
#' ginv <- .list_of_G(Species = A, Genus = A_genus)
#' }
#' @export
.list_of_G <- function(...) {
  args <- list(...)
  
  if (length(args) == 0) {
    stop("At least one inverse matrix must be provided")
  }
  
  if (is.null(names(args)) || any(names(args) == "")) {
    stop("All arguments must be named with the cluster variable name")
  }
  
  return(args)
}

