% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiBatch.R, R/MultiBatchList.R,
%   R/methods-McmcChains.R, R/methods-MixtureModel.R,
%   R/methods-MultiBatchModel.R
\docType{methods}
\name{[[,CnList,numeric,ANY-method}
\alias{[[,CnList,numeric,ANY-method}
\alias{[,MultiBatch,numeric,ANY,ANY-method}
\alias{[,MultiBatch,logical,ANY,ANY-method}
\alias{[[,MultiBatchList,numeric,ANY-method}
\alias{"[[",MultiBatchList,numeric,ANY}
\alias{[[,MultiBatchList,character,ANY-method}
\alias{[[,MultiBatchList,character,ANY}
\alias{[,MultiBatchList,numeric,ANY,ANY-method}
\alias{[,MultiBatchList,numeric,ANY,ANY}
\alias{[,MultiBatchList,character,ANY,ANY-method}
\alias{[,MultiBatchList,character,ANY,ANY}
\alias{[,MultiBatchList,logical,ANY,ANY-method}
\alias{[,MultiBatchList,logical,ANY,ANY}
\alias{[[<-,MultiBatchList,ANY,ANY,ANY-method}
\alias{[[<-,MultiBatchList,ANY,ANY}
\alias{[[<-,MultiBatchList,ANY,ANY-method}
\alias{[<-,MultiBatchList,ANY,ANY,ANY-method}
\alias{[,MultiBatchList,ANY,ANY,ANY}
\alias{[,McmcChains,ANY,ANY,ANY-method}
\alias{[,McmcChains-method}
\alias{[,McmcChains,ANY-method}
\alias{[,MixtureModel,logical,ANY,ANY-method}
\alias{"[",MixtureModel,logical,ANY,ANY}
\alias{[,MixtureModel,numeric,ANY,ANY-method}
\alias{"[",MixtureModel,numeric,ANY,ANY}
\alias{[,MultiBatchModel,ANY,ANY,ANY-method}
\alias{"[",MultiBatchModel,ANY,ANY,ANY}
\title{Subsetting methods for CNPBayes objects}
\usage{
\S4method{[[}{CnList,numeric,ANY}(x, i)

\S4method{[}{MultiBatch,numeric,ANY,ANY}(x, i, j, ..., drop = FALSE)

\S4method{[}{MultiBatch,logical,ANY,ANY}(x, i, j, ..., drop = FALSE)

\S4method{[[}{MultiBatchList,numeric,ANY}(x, i)

\S4method{[[}{MultiBatchList,character,ANY}(x, i)

\S4method{[}{MultiBatchList,numeric,ANY,ANY}(x, i, j, ..., drop = TRUE)

\S4method{[}{MultiBatchList,character,ANY,ANY}(x, i, j, ..., drop = TRUE)

\S4method{[}{MultiBatchList,logical,ANY,ANY}(x, i, j, ..., drop = TRUE)

\S4method{[[}{MultiBatchList,ANY,ANY,ANY}(x, i) <- value

\S4method{[}{MultiBatchList,ANY,ANY,ANY}(x, i, j) <- value

\S4method{[}{McmcChains,ANY,ANY,ANY}(x, i, j, ..., drop = FALSE)

\S4method{[}{MixtureModel,logical,ANY,ANY}(x, i, j, ..., drop = TRUE)

\S4method{[}{MixtureModel,numeric,ANY,ANY}(x, i, j, ..., drop = TRUE)

\S4method{[}{MultiBatchModel,ANY,ANY,ANY}(x, i, j, ..., drop = FALSE)
}
\arguments{
\item{x}{a MultiBatch instance}

\item{i}{elements to select}

\item{j}{ignored}

\item{...}{ignored}

\item{drop}{ignored}

\item{value}{a MultiBatch instance}
}
\value{
An object of class 'McmcChains'
}
\description{
Many of the classes defined in CNPBayes can be subset using the "[" operator.
}
