% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods-McmcParams.R
\docType{methods}
\name{burnin}
\alias{burnin}
\alias{burnin,MixtureModel-method}
\alias{burnin,MultiBatch-method}
\alias{burnin,MultiBatchList-method}
\alias{burnin<-}
\alias{burnin<-,MixtureModel,ANY-method}
\alias{burnin<-,MultiBatch,numeric-method}
\alias{burnin<-,MultiBatchList,numeric-method}
\alias{burnin,McmcParams-method}
\alias{burnin<-,McmcParams,ANY-method}
\alias{burnin<-,McmcParams-method}
\title{Number of burnin iterations.}
\usage{
burnin(object)

burnin(object) <- value

\S4method{burnin}{McmcParams}(object)

\S4method{burnin}{McmcParams,ANY}(object) <- value
}
\arguments{
\item{object}{A MultiBatch or MultiBatchList object}

\item{value}{new number of burnin iterations}
}
\value{
The number of burnin simulations.
}
\description{
This function retrieves the number of burnin simulations to be discarded.

This function changes the number of burnin simulations to be discarded.
}
\examples{
burnin(SingleBatchModelExample)
mp <- mcmcParams(SingleBatchModelExample)
burnin(mp)
}
