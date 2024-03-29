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
		return(list(ran = vars, cluster = cluster))
	}
}

#' @title get_type
#' @description Function takes a dataframe and identifies the class of variables (type)
#' @param fixed_form The fixed effect formula
#' @param cluster_name The name of the cluster variable
#' @param random_form The random effect formula
#' @return A list of variable names and their correpsonding class. 
#' @examples {
#' data(iris)}
#' @export
#' 

get_type <- function(data) {
 
  types  <- lapply(data, function(x) check_type(x))

  return(types)

}

#' @title check_type
#' @description Function takes a variable in a dataframe and checks/ classifies what tope
#' @param x A vector of data
#' @return Returns a character string specifying the class of the variable. 
#' @export
 
check_type <- function(x) {

        if(length(unique(x[!is.na(x)])) == 2){type <- "binary"}
		
		if(is.numeric(x)){type <- "continuous"}

		if(is.numeric(x) && all(x >= 0)){
			fit_norm <- stats::glm(x ~ 1, family = "gaussian")
			fit_pois <- suppressWarnings(stats::glm(x ~ 1, family = "poisson"))
			
			if((stats::AIC(fit_norm) - stats::AIC(fit_pois) >= 2)){type <- "count"} else{type <- "continuous"}			
		} # Need to work on identifying counts from continuous
    
		if(is.factor(x) && length(levels(x)) > 2){type <- "categorical"}
      
	    if(is.ordered(x)){type <- "ordered_categorical"}
  
  	return(type)
}

