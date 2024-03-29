% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_functions.R
\name{model_fit}
\alias{model_fit}
\title{model_fit}
\usage{
model_fit(data, fixformula, randformula, type, nitt, thin, burnin, n_rand)
}
\arguments{
\item{data}{A dataframe containing missing data.}

\item{fixformula}{A string that specifies the fixed effects in the model.}

\item{randformula}{A string that specifies the random effects in the model.}

\item{type}{A string that specifies the type of model to fit.}

\item{nitt}{An integer specifying the number of iterations to run the MCMC algorithm.}

\item{thin}{An integer specifying the thinning rate for the MCMC algorithm.}

\item{burnin}{An integer specifying the number of iterations to discard as burnin.}

\item{n.par.rand}{An integer specifying the number of random effects in the model.}
}
\value{
A list of draws from the posterior distribution of the model parameters.
}
\description{
Function takes incomplete data and a formula string and fits a model using MCMCglmm
}
