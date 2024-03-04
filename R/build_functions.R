
#' @title build_formula
#' @description Function takes a string specifying a formula and converts this to a formula object to be used in the models.
#' @param x A character strong specifying the formula one wishes to use in the model
#' @return A formula object
#' @export

build_formula <- function(x) {

	# Some checks. Formula must have a '`~`' in it
	if(!grepl("~", x)) { stop("Make sure you specify a formula string. This involves a structure of the form: y ~ x") } 

  return(as.formula(x))
}



