% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-McmcChains.R
\name{zFreq,McmcChains-method}
\alias{zFreq,McmcChains-method}
\alias{zfreq,McmcChains-method}
\title{Accessor for frequency of mixture component label in McmcChains class}
\usage{
\S4method{zFreq}{McmcChains}(object)
}
\arguments{
\item{object}{a McmcChains instance}
}
\description{
Accessor for frequency of mixture component label in McmcChains class
}
\details{
A running tally of the mixture component assignments is obtained by the `zFreq` accessor.  The running tally can also be extracted from a `MultiBatch` instance by \code{zFreq(chains(mb))}, where `mb` is a MultiBatch object.
}
