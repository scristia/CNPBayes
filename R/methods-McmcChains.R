.initializeMcmc <- function(object){
  ## add 1 for starting values (either the last run from the burnin,
  ## or default values if no burnin
  mcmc.params <- mcmcParams(object)
  nr <- iter(mcmc.params)
  ns <- length(y(object))
  K <- k(object)
  mati <- matrix(as.integer(NA), nr, K)
  vec <- numeric(nr)
  new("McmcChains",
      theta=matrix(NA, nr, K),
      sigma2=matrix(NA, nr, K),
      pi=matrix(NA, nr, K),
      mu=numeric(nr),
      tau2=numeric(nr),
      nu.0=numeric(nr),
      sigma2.0=numeric(nr),
      logprior=numeric(nr),
      loglik=numeric(nr),
      zfreq=mati,
      z=matrix(NA, nr, ns))
}

.initializeMcmcPooledVar <- function(object){
  ## add 1 for starting values (either the last run from the burnin,
  ## or default values if no burnin
  mcmc.params <- mcmcParams(object)
  nr <- iter(mcmc.params)
  ns <- length(y(object))
  K <- k(object)
  mati <- matrix(as.integer(NA), nr, K)
  vec <- numeric(nr)
  new("McmcChains",
      theta=matrix(NA, nr, K),
      sigma2=matrix(NA, nr, 1), ## this is the only difference from non-pooled
      pi=matrix(NA, nr, K),
      mu=numeric(nr),
      tau2=numeric(nr),
      nu.0=numeric(nr),
      sigma2.0=numeric(nr),
      logprior=numeric(nr),
      loglik=numeric(nr),
      zfreq=mati,
      z=matrix(NA, nr, ns))
}

setMethod("McmcChains", "missing", function(object){
  new("McmcChains", theta=matrix(), sigma2=matrix(),
      pi=matrix(), mu=numeric(), tau2=numeric(),
      nu.0=numeric(), sigma2.0=numeric(),
      zfreq=matrix(),
      logprior=numeric(),
      loglik=numeric(),
      z=matrix())
})




setMethod("McmcChains", "MixtureModel", function(object){
  .initializeMcmc(object)
})

setMethod("McmcChains", "SingleBatchPooledVar", function(object){
  .initializeMcmcPooledVar(object)
})

.initializeMcmcBatch <- function(object){
  mcmc.params <- mcmcParams(object)
  nr <- iter(mcmc.params)[1]
  ns <- length(y(object))
  K <- k(object)
  B <- nBatch(object)
  mati <- matrix(as.integer(NA), nr, K)
  new("McmcChains",
      theta=matrix(NA, nr, K*B),
      sigma2=matrix(NA, nr, K*B),
      pi=matrix(NA, nr, K),
      mu=matrix(NA, nr, K),
      tau2=matrix(NA, nr, K),
      nu.0=numeric(nr),
      sigma2.0=numeric(nr),
      logprior=numeric(nr),
      loglik=numeric(nr),
      zfreq=mati,
      z=matrix(NA, nr, ns))
}

setMethod("McmcChains", "BatchModel", function(object){
  .initializeMcmcBatch(object)
})

#' @rdname mu-method
#' @aliases mu,McmcChains-method
setMethod("mu", "McmcChains", function(object) object@mu)

#' @rdname tau2-method
#' @aliases tau2,McmcChains-method
setMethod("tau2", "McmcChains", function(object) object@tau2)

#' @rdname theta-method
#' @aliases theta,McmcChains-method
setMethod("theta", "McmcChains", function(object) object@theta)

#' @rdname sigma2-method
#' @aliases sigma2,missing-method
setMethod("sigma2", "McmcChains", function(object) object@sigma2)

setMethod("show", "McmcChains", function(object){
  cat("An object of class 'McmcChains'\n")
  cat("    chain dim:", nrow(theta(object)), "x", ncol(theta(object)), "\n")
  cat("    see theta(), sigma2(), p(), ...\n")
})

#' extract estimated parameters at particular iteration of simulation.
#' @aliases [,McmcChains-method [,McmcChains,ANY-method
#' @return An object of class 'McmcChains'
#' @docType methods
#' @rdname extract-methods
setMethod("[", "McmcChains", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@theta <- x@theta[i, , drop=FALSE]
    x@sigma2 <- x@sigma2[i, , drop=FALSE]
    x@pi <- x@pi[i, , drop=FALSE]
    if(is.matrix(x@mu)){
      x@mu <- x@mu[i, , drop=FALSE]
    } else x@mu <- x@mu[i]
    if(is.matrix(x@tau2)){
      x@tau2 <- x@tau2[i, , drop=FALSE]
    } else  x@tau2 <- x@tau2[i]
    x@nu.0 <- x@nu.0[i]
    x@sigma2.0 <- x@sigma2.0[i]
    x@logprior <- x@logprior[i]
    x@loglik <- x@loglik[i]
    x@zfreq <- x@zfreq[i, , drop=FALSE]
    x@z <- x@z[i, , drop=FALSE]
  }
  x
})

#' @rdname nu.0-method
#' @aliases nu.0,McmcChains-method
setMethod("nu.0", "McmcChains", function(object) object@nu.0)

#' @rdname sigma2.0-method
#' @aliases sigma2.0,McmcChains-method
setMethod("sigma2.0", "McmcChains", function(object) object@sigma2.0)

setReplaceMethod("theta", "McmcChains", function(object, value){
  object@theta <- value
  object
})

setReplaceMethod("sigma2", "McmcChains", function(object, value){
  object@sigma2 <- value
  object
})

setReplaceMethod("p", "McmcChains", function(object, value){
  object@pi <- value
  object
})

setReplaceMethod("mu", "McmcChains", function(object, value){
  object@mu <- value
  object
})

setReplaceMethod("tau2", "McmcChains", function(object, value){
  object@tau2 <- value
  object
})

setReplaceMethod("nu.0", "McmcChains", function(object, value){
  object@nu.0 <- value
  object
})

setReplaceMethod("sigma2.0", "McmcChains", function(object, value){
  object@sigma2.0 <- value
  object
})

setReplaceMethod("log_lik", "McmcChains", function(object, value){
  object@loglik <- value
  object
})

#' @rdname log_lik-method
#' @aliases log_lik,McmcChains-method
setMethod("log_lik", "McmcChains", function(object){
  object@loglik
})

#' Retrieve the names of the parameters estimated in the MCMC chain.
#' @param x an object of class 'McmcChains'
#' @return A vector of strings containing the names of each parameter
#' @aliases names,McmcChains-method
#' @docType methods
#' @rdname names-methods
setMethod("names", "McmcChains", function(x) slotNames(x))

#' @rdname zfreq-method
#' @aliases zfreq,McmcChains-method
setMethod("zFreq", "McmcChains", function(object) object@zfreq )

#' @rdname z-method
#' @aliases z,McmcChains-method
setMethod("z", "McmcChains", function(object) object@z )

#' @rdname logPrior-method
#' @aliases logPrior,McmcChains-method
setMethod("logPrior", "McmcChains", function(object) object@logprior)

setReplaceMethod("logPrior", "McmcChains", function(object, value) {
  object@logprior <- value
  object
})

setReplaceMethod("zFreq", "McmcChains", function(object, value){
  object@zfreq <- value
  object
})
