#' @title get_variables
#' @description Function takes a formula string and identifies the variables in the formula that relate to the data
#' @param x A string
#' @param fix A logical value. If TRUE, the function returns the fixed effects variables. If FALSE, the function returns the random effects variables.
#' @return A vector of variable names. This vector is used to subset out the 'fixed effects' and 'random effect' columns within a given dataframe. 
#' @examples {
#' # All should return the same three variables: y, x1, x2
#' form <- "y ~ x1 + x2"
#' get_variables(form)
#' form <- "y ~ x1 * x2"
#' get_variables(form)
#' form <- "y ~ x1 + x2 + x1 : x2"
#' get_variables(form)
#' form <- "y ~ x1 + x2 + x1:x2"
#' get_variables(form)
#' form <- "y ~ x1 + x2 + x1*x2"
#' get_variables(form)
#' form <- "~ 1 + x1|Species"
#' get_variables(form, fix = FALSE)
#' form <- "~ 1 + x1| Species"
#' get_variables(form, fix = FALSE)
#' form <- "~ 1 |Species"
#' get_variables(form, fix = FALSE)
#' form <- "~ 1 | Species"
#' get_variables(form, fix = FALSE)}
#' @export
get_variables <- function(x, fix = TRUE) {

  # helper: split into "tokens" (variable-ish strings), dropping empties
  tokens <- function(s) {
    out <- unlist(strsplit(s, "\\W+"))
    out <- out[out != ""]
    unique(out)
  }

  if (fix) {
    # For fixed effects, just return the tokens
    return(list(fix = tokens(x)))
  }

  # For random effects, allow multiple terms:
  # e.g. "(1 | g) + (x1 + x2 | h)" or "1 | g"
  # Split on '+' at top-level-ish (good enough for typical formulas)
  terms <- trimws(unlist(strsplit(x, "\\s*\\+\\s*")))

  ran_vars <- character()
  clusters <- character()

  for (term in terms) {
    # remove surrounding parentheses if present
    term2 <- trimws(gsub("^\\(|\\)$", "", term))

    # split on pipe with optional whitespace around it
    parts <- strsplit(term2, "\\s*\\|\\s*", perl = TRUE)[[1]]
    if (length(parts) < 2) next  # not a random effect term

    lhs <- trimws(parts[1])  # e.g. "1 + x"
    rhs <- trimws(parts[2])  # e.g. "group"

    clusters <- c(clusters, rhs)

    # random slope vars from LHS (drop intercept & empties)
    lhs_vars <- tokens(lhs)
    lhs_vars <- lhs_vars[!lhs_vars %in% c("1")]

    ran_vars <- c(ran_vars, lhs_vars)
  }

  clusters <- unique(clusters)
  ran_vars <- unique(ran_vars)

  # Return intercept + slopes; intercept always included for random term
  list(ran = c(1, ran_vars), cluster = clusters)
}

#' @title get_type
#' @description Function takes a dataframe and identifies the class of variables (type)
#' @param fixed_form The fixed effect formula
#' @param cluster_name The name of the cluster variable
#' @param random_form The random effect formula
#' @return A list of variable names and their corresponding class. 
#' @examples \dontrun{
#' data <- data.frame(y = c(1,2,3,NA,5), x1 = factor(c("A","B","A","B","A")), x2 = c(10,20,30,40,50))
#' get_type(data)
#' }
#' @export
get_type <- function(data) {
  types  <- lapply(data, function(x) check_type(x))
  return(types)

}

#' @title check_type
#' @description Function takes a variable in a dataframe and checks / classifies what type of variable it is (binary, continuous, count, categorical, ordered_categorical) so that it can be modeled appropriately.
#' @param x A vector of data
#' @return Returns a character string specifying the class of the variable. 
#' @examples \dontrun{
#' x1 <- as.factor(c(1,0,1,1,0,NA))
#' check_type(x1) # returns "binary"
#' x1.1 <- c(1,0,1,1,0,NA)
#' check_type(x1.1) # returns "warning" about binary numeric coding
#' x2 <- c(1.5, 2.3, 3.1, NA, 4.0)
#' check_type(x2) # returns "continuous"
#' x3 <- c(0,1,2,3,4,NA)	
#' check_type(x3) # returns "count"
#' x4 <- as.factor(c("A","B","C","A","B",NA))
#' check_type(x4) # returns "categorical"
#' x5 <- as.ordered(c("low","medium","high","medium",NA))
#' check_type(x5) # returns "ordered_categorical"
#' }
#' @export
check_type <- function(x) {
	
	# Use only complete cases to fit the models
		x_complete <- x[!is.na(x)]

	# Add a warning to check binary variables dummy coded as 0/1 numeric
		if(length(unique(x_complete)) == 2 && is.numeric(x_complete)){
			warning("Check that binary / categorical variables are coded as factors / ordered factors, not numeric (e.g., 0/1 dummy coding).")
		}

	# Classification of variable type
		type <- NULL 
        
		# Binary
		if(length(unique(x_complete)) == 2 && is.factor(x_complete)){type <- "ordinal"}

		# Continuous or count
		if(is.numeric(x_complete) && all(x_complete >= 0)){
			
			fit_norm <- suppressWarnings(stats::glm(x_complete ~ 1, family = "gaussian"))
			fit_pois <- suppressWarnings(stats::glm(x_complete ~ 1, family = "poisson"))
			
			if((stats::AIC(fit_norm) - stats::AIC(fit_pois) >= 2)){
				type <- "poisson" } else {type <- "gaussian"}			
		} 
    
		# Multinomial categorical
		if(is.factor(x_complete) && length(levels(x_complete)) > 2){type <- "categorical"}
      
		# Ordered categorical
	    if(is.ordered(x_complete)){type <- "ordinal"}
  
  	return(type)
}

#' @title data_prep
#' @description Function prepares the data for imputation by handling missing values and standardizing continuous variables.
#' @param formula A formula specifying the response and predictor variables.
#' @param data A data frame containing the dataset to be prepared.
#' @param types A list specifying the type of each variable in the dataset.
#' @return A list containing the prepared data frame and attributes for continuous variables.
#' @examples \dontrun{
#' data <- data.frame(y = c(1,2,3,NA,5), x1 = factor(c("A","B","A","B","A")), x2 = c(10,20,30,NA,50))
#' formula <- as.formula("y ~ x1 + x2")
#' types <- list(y = "gaussian", x1 = "categorical", x2 = "gaussian")
#' data_prep(formula, data, types)
#' }	
#' @export
data_prep  <- function(formula, data, types) {
				# Identify response variable in formula
			response_var <- all.vars(formula[[2]])
			
			# Identify predictors in formula
			  predictors <- all.vars(formula[[3]])
			
			# Check missing data in predictors and impute values for the model for. For count and gaussian we use the mean for now, for categorical we sample from the empirical distribution.
			predictor_data <- data[, predictors, drop = FALSE]
			for(p in 1:ncol(predictor_data)){
				if(any(is.na(predictor_data[, p]))){
					if(types[[predictors[p]]] %in% c("gaussian")){
						predictor_data[is.na(predictor_data[, p]), p] <- mean(predictor_data[, p], na.rm = TRUE)
					}
					
					if(types[[predictors[p]]] %in% c("poisson")){
						predictor_data[is.na(predictor_data[, p]), p] <- round(mean(predictor_data[, p], na.rm = TRUE))
					}

					if(types[[predictors[p]]] %in% c("categorical", "threshold")){
						obs_values <- predictor_data[!is.na(predictor_data[, p]), p]
						predictor_data[is.na(predictor_data[, p]), p] <- sample(obs_values, sum(is.na(predictor_data[, p])), replace = TRUE)
					}
				}
			}
			
			# Create the data frame for the model fitting
			  data_i <- data.frame(data[, response_var],
			                     predictor_data,
			                   data[, phylo_ran[["cluster"]], drop = FALSE])
			colnames(data_i)[1] <- response_var

			# z-transform all gaussian variables for better mixing and store attributes to revert later name the slots with variable names

			data_i_attrs <- extract_gaussian_attrs(data_i, types)

			data_i <- data_i %>%
			          dplyr::mutate(dplyr::across(.cols = dplyr::where(is.numeric) & dplyr::all_of(names(types)[types == "gaussian"]),
			                                      .fns = ~ as.numeric(scale(.x))))

return(list = (list(data_i,
				   data_i_attrs)))
}


#' @title extract_gaussian_attrs
#' @description Function extracts the mean and standard deviation of gaussian variables in a dataframe for later use in back-transformation.
#' @param data A data frame containing the dataset.
#' @param types A list specifying the type of each variable in the dataset.
#' @return A list containing the mean and standard deviation of each gaussian variable.
#' @examples \dontrun{
#' data <- data.frame(y = c(1,2,3,NA,5), x1 = factor(c("A","B","A","B","A")), x2 = c(10,20,30,NA,50))
#' types <- list(y = "gaussian", x1 = "categorical", x2 = "gaussian")
#' extract_gaussian_attrs(data, types)
#' }
#' @export
extract_gaussian_attrs <- function(data, types) {
  out <- lapply(names(types), function(v) {
    if (!is.null(types[[v]]) && types[[v]] == "gaussian") {
      x <- data[[v]]
      list(
        mean = mean(x, na.rm = TRUE),
        sd   = sd(x, na.rm = TRUE)
      )
    } else {
      NULL
    }
  })
  
  names(out) <- names(types)
  out
}

#' @title summarise_var_types
#' @description Function summarizes the types and characteristics of variables in a dataframe.
#' @param df A data frame containing the dataset.
#' @return A data frame summarizing the types and characteristics of each variable.
#' @examples \dontrun{
#' data <- data.frame(y = c(1,2,3,NA,5), x1 = factor(c("A","B","A","B","A")), x2 = c(10,20,30,NA,50))
#' summarise_var_types(data)
#' }
#' @export
summarise_var_types <- function(df, store_levels = TRUE, max_levels_store = 200) {
  stopifnot(is.data.frame(df))
  
  # helper to store "levels" for factors, and optionally for low-cardinality character columns
  get_levels <- function(x) {
    if (!store_levels) return(NULL)
    if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) return(NULL)
    
    if (is.factor(x)) {
      lv <- levels(x)
      return(if (length(lv) <= max_levels_store) lv else NULL)
    }
    if (is.character(x)) {
      ux <- sort(unique(x[!is.na(x)]))
      return(if (length(ux) <= max_levels_store) ux else NULL)
    }
    NULL
  }
  
  out <- data.frame(
    variable = names(df),
    class = vapply(df, function(x) paste(class(x), collapse = ", "), character(1)),
    typeof = vapply(df, typeof, character(1)),
    is_numeric = vapply(df, is.numeric, logical(1)),
    is_integer = vapply(df, is.integer, logical(1)),
    is_factor = vapply(df, is.factor, logical(1)),
    is_character = vapply(df, is.character, logical(1)),
    is_logical = vapply(df, is.logical, logical(1)),
    n_levels = vapply(df, function(x) if (is.factor(x)) nlevels(x) else NA_integer_, integer(1)),
    n_unique = vapply(df, function(x) length(unique(x[!is.na(x)])), integer(1)),
    has_na = vapply(df, function(x) any(is.na(x)), logical(1)),
    stringsAsFactors = FALSE
  )
  
  # list-column of levels
  out$levels <- lapply(df, get_levels)
  out
}

