% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{McmcChains-class}
\alias{McmcChains-class}
\title{An object to hold estimated paraeters.}
\description{
An object of this class holds estimates of each parameter at each iteration of the MCMC simulation.
}
\section{Slots}{

\describe{
\item{\code{theta}}{means of each batch and component}

\item{\code{sigma2}}{variances of each batch and component}

\item{\code{pi}}{mixture probabilities}

\item{\code{mu}}{overall mean in a marginal. In batch model, averaged across batches}

\item{\code{tau2}}{overall variance in a marginal model. In a batch model, weighted average by precision across batches.}

\item{\code{nu.0}}{shape parameter for sigma.2 distribution}

\item{\code{sigma2.0}}{rate parameter for sigma.2 distribution}

\item{\code{logprior}}{log likelihood of prior.}

\item{\code{loglik}}{log likelihood.}

\item{\code{zfreq}}{table of z.}

\item{\code{predictive}}{posterior predictive distribution}

\item{\code{zstar}}{needed for plotting posterior predictive distribution}

\item{\code{k}}{integer specifying number of components}

\item{\code{iter}}{integer specifying number of MCMC simulations}

\item{\code{B}}{integer specifying number of batches}
}}

