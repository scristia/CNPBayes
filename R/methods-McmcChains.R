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
      pi_chd=matrix(NA, nr, K),
      mu=numeric(nr),
      tau2=numeric(nr),
      nu.0=numeric(nr),
      sigma2.0=numeric(nr),
      logprior=numeric(nr),
      loglik=numeric(nr),
      zfreq=mati,
      zfreq_parents=mati)
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
      zfreq=mati)
}

setMethod("McmcChains", "missing", function(object){
  new("McmcChains", theta=matrix(0), sigma2=matrix(0),
      pi=matrix(0), mu=numeric(),  tau2=numeric(), 
      nu.0=numeric(), sigma2.0=numeric(),
      zfreq=matrix(0),
      logprior=numeric(),
      loglik=numeric())
})

setMethod("McmcChains", "MixtureModel", function(object){
  .initializeMcmc(object)
})

setMethod("McmcChains", "SingleBatchPooled", function(object){
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
      zfreq=mati)
}


setMethod("McmcChains", "MultiBatchModel", function(object){
  .initializeMcmcBatch(object)
})

.initializeMcmcBatch2 <- function(object){
  mcmc.params <- mcmcParams(object)
  nr <- iter(mcmc.params)[1]
  ns <- length(y(object))
  K <- k(object)
  B <- nBatch(object)
  mati <- matrix(as.integer(NA), nr, K)
  new("McmcChains",
      theta=matrix(NA, nr, K*B),
      theta_chd=matrix(NA, nr, K*B),
      sigma2=matrix(NA, nr, K*B),
      sigma2_chd=matrix(NA, nr, K*B),
      pi=matrix(NA, nr, K),
      pi_chd=matrix(NA, nr, K),
      mu=matrix(NA, nr, K),
      mu_chd=matrix(NA, nr, K),
      tau2=matrix(NA, nr, K),
      tau2_chd=matrix(NA, nr, K),
      nu.0=numeric(nr),
      nu.0_chd=numeric(nr),
      sigma2.0=numeric(nr),
      sigma2.0_chd=numeric(nr),
      logprior=numeric(nr),
      loglik=numeric(nr),
      zfreq=mati,
      zfreq_parents=mati,
      zfreq_chd=mati)
}

setMethod("McmcChains", "TrioBatchModel", function(object){
  .initializeMcmcBatch2(object)
})

chains_mb <- function(object){
  mcmc.params <- mcmcParams(object)
  nr <- iter(mcmc.params)[1]
  ns <- length(y(object))
  K <- k(object)
  B <- nBatch(object)
  mati <- matrix(as.integer(NA), nr, K)
  new("McmcChains",
      theta=matrix(NA, nr, K*B),
      sigma2=matrix(NA, nr, B),
      pi=matrix(NA, nr, K),
      mu=matrix(NA, nr, K),
      tau2=matrix(NA, nr, K),
      nu.0=numeric(nr),
      sigma2.0=numeric(nr),
      logprior=numeric(nr),
      loglik=numeric(nr),
      zfreq=mati)
}

setMethod("McmcChains", "MultiBatchPooled", function(object){
  chains_mb(object)
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
    x@pi_chd <- x@pi_chd[i, , drop=FALSE]
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
    x@zfreq_parents <- x@zfreq_parents[i, , drop=FALSE]
  }
  x
})

#' @rdname nu.0-method
#' @aliases nu.0,McmcChains-method
setMethod("nu.0", "McmcChains", function(object) object@nu.0)

#' @rdname sigma2.0-method
#' @aliases sigma2.0,McmcChains-method
setMethod("sigma2.0", "McmcChains", function(object) object@sigma2.0)

setReplaceMethod("pp", "McmcChains", function(object, value){
  object@pi_chd <- value
  object
})

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

#' @rdname zfreqpar-method
#' @aliases zfreqpar,McmcChains-method
setMethod("zFreqPar", "McmcChains", function(object) object@zfreq_parents )

#' @rdname zfreqchd-method
#' @aliases zfreqchd,McmcChains-method
setMethod("zFreqChd", "McmcChains", function(object) object@zfreq_chd )

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

setReplaceMethod("zFreqPar", "McmcChains", function(object, value){
  object@zfreq_parents <- value
  object
})

setReplaceMethod("zFreqChd", "McmcChains", function(object, value){
  object@zfreq_chd <- value
  object
})

#' @rdname thetachd-method
#' @aliases thetachd,McmcChains-method
setMethod("thetachd", "McmcChains", function(object) object@theta_chd )

setReplaceMethod("thetachd", "McmcChains", function(object, value){
  object@theta_chd <- value
  object
})

#' @rdname sigma2chd-method
#' @aliases sigma2chd,McmcChains-method
setMethod("sigma2chd", "McmcChains", function(object) object@sigma2_chd )

setReplaceMethod("sigma2chd", "McmcChains", function(object, value){
  object@sigma2_chd <- value
  object
})

#' @rdname muchd-method
#' @aliases muchd,McmcChains-method
setMethod("muchd", "McmcChains", function(object) object@mu_chd )

setReplaceMethod("muchd", "McmcChains", function(object, value){
  object@mu_chd <- value
  object
})

#' @rdname muchd-method
#' @aliases muchd,McmcChains-method
setMethod("muchd", "McmcChains", function(object) object@mu_chd )

setReplaceMethod("muchd", "McmcChains", function(object, value){
  object@mu_chd <- value
  object
})

#' @rdname tau2chd-method
#' @aliases tau2chd,McmcChains-method
setMethod("tau2chd", "McmcChains", function(object) object@tau2_chd )

setReplaceMethod("tau2chd", "McmcChains", function(object, value){
  object@tau2_chd <- value
  object
})

#' @rdname pp-method
#' @aliases pp,McmcChains-method
setMethod("pp", "McmcChains", function(object) object@pi_chd )

setReplaceMethod("pp", "McmcChains", function(object, value){
  object@pi_chd <- value
  object
})

#' @rdname sigma2.0chd-method
#' @aliases sigma2.0chd,McmcChains-method
setMethod("sigma2.0chd", "McmcChains", function(object) object@sigma2.0_chd )

setReplaceMethod("sigma2.0chd", "McmcChains", function(object, value){
  object@sigma2.0_chd <- value
  object
})

#' @rdname nuchd-method
#' @aliases nuchd,McmcChains-method
setMethod("nuchd", "McmcChains", function(object) object@nu.0_chd )

setReplaceMethod("nuchd", "McmcChains", function(object, value){
  object@nu.0_chd <- value
  object
})


#' @rdname probzchd-method
#' @aliases probzchd,McmcChains-method
setMethod("probzchd", "McmcChains", function(object) object@probz_chd )

setReplaceMethod("probzchd", "McmcChains", function(object, value){
  object@probz_chd <- value
  object
})

#' @rdname probzpar-method
#' @aliases probzpar,McmcChains-method
setMethod("probzpar", "McmcChains", function(object) object@probz_par )

setReplaceMethod("probzpar", "McmcChains", function(object, value){
  object@probz_par <- value
  object
})
