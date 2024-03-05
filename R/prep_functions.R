
#' @title get_class
#' @description Function identifies the class of all variables in the dataframe
#' @param x A dataframe containing missing data.
#' @return A list identifying the class of all columns in the dataframe.
#' @examples {
#' data(iris)
#' get_class(iris)}
#' @export
get_class <- function(x) {
  return(lapply(x, class))	
}

#' @title get_variables
#' @description Function takes a formula string and identifies the variables in the formula that relate to the data
#' @param x A string
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
#' get_variables(form)}
#' @export
get_variables <- function(x) {
	vars <- unique(unlist(strsplit(x, "\\W+")))
  return(vars)
}