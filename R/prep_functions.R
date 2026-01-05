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
#' form <- "~ 1 |Species"
#' get_variables(form, fix = FALSE)}
#' @export
get_variables <- function(x, fix = TRUE) {
	
	if(fix){
		vars <- unique(unlist(strsplit(x, "\\W+")))
		return(list(fix = vars))
	} else{
		cluster <- unique(unlist(strsplit(x, "\\|")))[2]
		vars <- unique(unlist(strsplit(x, "\\W+")))
		vars <- vars[!vars %in% c(cluster,"1", "")]
		return(list(ran = c(1, vars), cluster = cluster))
	}
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
		if(length(unique(x_complete)) == 2 && is.factor(x_complete)){type <- "binary"}

		# Continuous or count
		if(is.numeric(x_complete) && all(x_complete >= 0)){
			
			fit_norm <- suppressWarnings(stats::glm(x_complete ~ 1, family = "gaussian"))
			fit_pois <- suppressWarnings(stats::glm(x_complete ~ 1, family = "poisson"))
			
			if((stats::AIC(fit_norm) - stats::AIC(fit_pois) >= 2)){
				type <- "count" } else {type <- "continuous"}			
		} 
    
		# Multinomial categorical
		if(is.factor(x_complete) && length(levels(x_complete)) > 2){type <- "categorical"}
      
		# Ordered categorical
	    if(is.ordered(x_complete)){type <- "ordered_categorical"}
  
  	return(type)
}

