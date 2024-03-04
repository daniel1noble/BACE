
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

