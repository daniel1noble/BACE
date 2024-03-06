
#' @title build_formula
#' @description Function takes a string specifying a formula and converts this to a formula object to be used in the models.
#' @param x A character string specifying the formula used in the model. This should be of the form: y ~ x.
#' @return A formula object
#' @export

build_formula <- function(x) {

	# Some checks. Formula must have a '`~`' in it
	if(!grepl("~", as.character(x))) { stop("Make sure you specify a formula string. This involves a structure of the form: y ~ x") } 

  return(stats::as.formula(x))
}

#' @title build_formula_string
#' @description Function takes a string specifying a formula and creates all combinations of models using the variables.
#' @param x A character string specifying the formula used in the model. This should be of the form: y ~ x.
#' @return A list of formulas
#' @export

build_formula_string <- function(x) {

	# Some checks. Formula must have a '`~`' in it
	if(!grepl("~", as.character(x))) { stop("Make sure you specify a formula string. This involves a structure of the form: y ~ x") } 

	vars <- get_variables(x) 

	formulas <- list()
	# Create all combinations of the variables
	for(i in 1:length(vars)) {

	 	if(i == 1) { formulas[[i]] <- x}
		else{
		formulas[[i]] <-  paste(vars[i], "~", paste(vars[-i], collapse = " + "))
	}}

  return(lapply(formulas, function(x) build_formula(x)))
}






