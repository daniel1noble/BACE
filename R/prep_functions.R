#' @importFrom stats pnorm predict runif sd setNames
#' @importFrom utils combn head
#' @keywords internal
"_PACKAGE"

#' @title .get_variables
#' @description Function takes a formula string and identifies the variables in the formula that relate to the data
#' @param x A string
#' @param fix A logical value. If TRUE, the function returns the fixed effects variables. If FALSE, the function returns the random effects variables.
#' @return A vector of variable names. This vector is used to subset out the 'fixed effects' and 'random effect' columns within a given dataframe. 
#' @examples {
#' # All should return the same three variables: y, x1, x2
#' form <- "y ~ x1 + x2"
#' .get_variables(form)
#' form <- "y ~ x1 * x2"
#' .get_variables(form)
#' form <- "y ~ x1 + x2 + x1 : x2"
#' .get_variables(form)
#' form <- "y ~ x1 + x2 + x1:x2"
#' .get_variables(form)
#' form <- "y ~ x1 + x2 + x1*x2"
#' .get_variables(form)
#' form <- "~ 1 + x1|Species"
#' .get_variables(form, fix = FALSE)
#' form <- "~ 1 + x1| Species"
#' .get_variables(form, fix = FALSE)
#' form <- "~ 1 |Species"
#' .get_variables(form, fix = FALSE)
#' form <- "~ 1 | Species"
#' .get_variables(form, fix = FALSE)}
#' @export
.get_variables <- function(x, fix = TRUE) {

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

#' @title .get_type
#' @description Function takes a variable in a dataframe and checks / classifies what type of variable it is (binary, continuous, count, categorical, ordered_categorical) so that it can be modeled appropriately.
#' @param x A summary of a variables in a dataframe
#' @param var The name of the variable
#' @param data The dataframe containing the variables
#' @return Returns a character string specifying the class of the variable. 
#' @keywords internal
#' @export
.get_type <- function(x, var, data) {
	
      x_complete <- data[[var]][!is.na(data[[var]])]

	# Classification of variable type
		   type <- NULL 
    sum_var <- x[x$variable == var, ]
		
    # Binary
		if(sum_var$is_factor == TRUE && sum_var$is_ordered == TRUE && sum_var$n_levels == 2){type <- "threshold"}

    if(sum_var$is_factor == TRUE && sum_var$is_ordered == FALSE && sum_var$n_levels == 2){type <- "threshold"}

    # Gaussian / Continuous
    if(sum_var$is_numeric == TRUE & sum_var$is_integer == FALSE){type <- "gaussian"}

		# Poisson?
		if(sum_var$is_numeric == TRUE & sum_var$is_integer == TRUE){
			
			fit_norm <- suppressWarnings(stats::glm(x_complete ~ 1, family = "gaussian"))
			fit_pois <- suppressWarnings(stats::glm(x_complete ~ 1, family = "poisson"))
			
			if((stats::AIC(fit_norm) - stats::AIC(fit_pois) >= 2)){
				type <- "poisson" } else {type <- "gaussian"}			
		} 
    
		# Multinomial categorical
		if(sum_var$is_factor == TRUE && sum_var$is_ordered == FALSE && sum_var$n_levels > 2){type <- "categorical"}
    if(sum_var$is_factor == TRUE && sum_var$is_ordered == TRUE && sum_var$n_levels > 2) {type <- "threshold"}
      
  	return(type)
}

#' @title data_prep
#' @description Function prepares the data for imputation by handling missing values and standardizing continuous variables.
#' @param formula A formula specifying the response and predictor variables.
#' @param data A data frame containing the dataset to be prepared.
#' @param types A list specifying the type of each variable in the dataset.
#' @param ran_cluster A string specifying the random effect cluster variable.
#' @return A list containing the prepared data frame and attributes for continuous variables.
#' @examples \dontrun{
#' data <- data.frame(y = c(1,2,3,NA,5), x1 = factor(c("A","B","A","B","A")), x2 = c(10,20,30,NA,50))
#' formula <- as.formula("y ~ x1 + x2")
#' types <- list(y = "gaussian", x1 = "categorical", x2 = "gaussian")
#' data_prep(formula, data, types)
#' }	
#' @export
.data_prep  <- function(formula, data, types, ran_cluster) {
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
			                       data[, ran_cluster, drop = FALSE])
			colnames(data_i)[1] <- response_var

			# z-transform all gaussian variables for better mixing and store attributes to revert later name the slots with variable names

			data_i_attrs <- .extract_gaussian_attrs(data_i, types)

		# Get gaussian variables that actually exist in data_i
		gaussian_vars_in_data <- intersect(names(types)[types == "gaussian"], colnames(data_i))

		data_i <- data_i %>%
		          dplyr::mutate(dplyr::across(.cols = dplyr::all_of(gaussian_vars_in_data),
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
.extract_gaussian_attrs <- function(data, types) {
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
#' @param store_levels A logical indicating whether to store levels for factor and character variables. Default is TRUE.
#' @param max_levels_store An integer specifying the maximum number of levels to store. Default is 200.
#' @return A data frame summarizing the types and characteristics of each variable.
#' @examples \dontrun{
#' data <- data.frame(y = c(1,2,3,NA,5), x1 = factor(c("A","B","A","B","A")), x2 = c(10,20,30,NA,50))
#' summarise_var_types(data)
#' }
#' @export
.summarise_var_types <- function(df, store_levels = TRUE, max_levels_store = 200) {
  
  stopifnot(is.data.frame(df))

  out <- data.frame(
        variable = names(df),
           class = vapply(df, function(x) paste(class(x), collapse = ", "), character(1)),
          typeof = vapply(df, typeof, character(1)),
      is_numeric = vapply(df, is.numeric, logical(1)),
      is_integer = vapply(df, is.integer, logical(1)),
       is_factor = vapply(df, is.factor, logical(1)),
      is_ordered = vapply(df, is.ordered, logical(1)),
    is_character = vapply(df, is.character, logical(1)),
      is_logical = vapply(df, is.logical, logical(1)),
        n_levels = vapply(df, function(x) if (is.factor(x)) nlevels(x) else NA_integer_, integer(1)),
        n_unique = vapply(df, function(x) length(unique(x[!is.na(x)])), integer(1)),
          has_na = vapply(df, function(x) any(is.na(x)), logical(1)),
    stringsAsFactors = FALSE
  )
  
  # list-column of levels
  out$levels <- lapply(df, function(x) .get_levels(x, store_levels, max_levels_store))
  out
}

 
  #' @title .get_levels
  #' @description Helper function to get levels of a variable if applicable
  #' @param x A variable
  #' @param store_levels A logical indicating whether to store levels
  #' @param max_levels_store An integer specifying the maximum number of levels to store
  #' @return A vector of levels or NULL if not applicable
  .get_levels <- function(x, store_levels, max_levels_store) {
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
  

  #' Pipe operator
#'
#' See magrittr::%>% for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
NULL