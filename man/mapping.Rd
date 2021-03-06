% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{mapping}
\alias{mapping}
\alias{mapping,MultiBatch-method}
\alias{mapping<-}
\alias{mapping<-,MultiBatch-method}
\title{Map mixture components to copy number states}
\usage{
mapping(object)

mapping(object) <- value
}
\arguments{
\item{object}{a `MultiBatch` model}

\item{value}{a k-length numeric vector with values in {1, 2, ..., k}, where k is the number of mixture components}
}
\value{
a character vector of length k (number of components) indicating the inferred copy number of each component.  The mapping is not neccessarily one-to-one as multiple mixture components can be mapped to a single copy number state.
}
\description{
Map mixture components to copy number states
}
