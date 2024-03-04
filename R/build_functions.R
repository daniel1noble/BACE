
#' @title build_formula
#' @description Function takes a string specifying a formula and converts this to a formula object to be used in the models.
#' @param x A character string specifying the formula used in the model. This should be of the form: y ~ x.
#' @return A formula object
#' @export

build_formula <- function(x) {

	# Some checks. Formula must have a '`~`' in it
	if(!grepl("~", x)) { stop("Make sure you specify a formula string. This involves a structure of the form: y ~ x") } 

  return(stats::as.formula(x))
}




