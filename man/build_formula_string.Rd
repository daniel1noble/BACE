% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_functions.R
\name{build_formula_string}
\alias{build_formula_string}
\title{build_formula_string}
\usage{
build_formula_string(x)
}
\arguments{
\item{x}{A character string specifying the formula used in the model. This should be of the form: y ~ x.}
}
\value{
A list of formulas
}
\description{
Function takes a string specifying a formula and creates all combinations of models using the variables.
}
