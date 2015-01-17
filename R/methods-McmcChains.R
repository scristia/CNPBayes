##McmcChains <- function(mcmc.params, hyper.param){
##  if(missing(mcmc.params)){
##    return(new("McmcChains", theta=matrix(), sigma2=matrix(),
##               pi=matrix(), mu=numeric(), tau2=numeric(),
##               nu.0=numeric(), sigma2.0=numeric(),
##               logpotential=numeric()))
##  }
##  nr <- iter(mcmc.params)/thin(mcmc.params)
##  K <- k(hyper.param)
##  mat <- matrix(NA, nr, K)
##  vec <- numeric(nr)
##  new("McmcChains", theta=mat,
##      sigma2=mat,
##      pi=mat,
##      mu=vec,
##      tau2=vec,
##      nu.0=vec,
##      sigma2.0=vec,
##      logpotential=vec)
##}

.initializeMcmc <- function(object, mcmc.params){
  ## add 1 for starting values (either the last run from the burnin, or default values if no burnin
  nr <- iter(mcmc.params)/thin(mcmc.params) + 1
  K <- k(object)
  mat <- matrix(NA, nr, K)
  vec <- numeric(nr)
  new("McmcChains", theta=mat,
      sigma2=mat,
      pi=mat,
      mu=vec,
      tau2=vec,
      nu.0=vec,
      sigma2.0=vec,
      logpotential=vec)
}

setMethod("McmcChains", "Hyperparameters", function(object, mcmc.params){
  .initializeMcmc(object, mcmc.params)
})

setMethod("McmcChains", "missing", function(object, mcmc.params){
  new("McmcChains", theta=matrix(), sigma2=matrix(),
      pi=matrix(), mu=numeric(), tau2=numeric(),
      nu.0=numeric(), sigma2.0=numeric(),
      logpotential=numeric())
})

setMethod("McmcChains", "MixtureModel", function(object, mcmc.params){
  .initializeMcmc(object, mcmc.params)
})

.initializeMcmcBatch <- function(object, mcmc.params){
  nr <- iter(mcmc.params[1])/thin(mcmc.params)[1] + 1
  K <- k(object)
  B <- nBatch(object)
  mat_batch <- matrix(NA, nr, K*B)
  mat <- matrix(NA, nr, K)
  vec <- numeric(nr)
  new("McmcChains",
      theta=mat_batch,
      sigma2=mat_batch,
      pi=mat,
      mu=mat,
      tau2=mat,
      nu.0=vec,
      sigma2.0=vec,
      logpotential=vec)
}

setMethod("McmcChains", c("BatchModel", "missing"), function(object, mcmc.params){
  .initializeMcmcBatch(object, mcmc.params)
})

setMethod("McmcChains", c("BatchModel", "McmcParams"), function(object, mcmc.params){
  .initializeMcmcBatch(object, mcmc.params)
})

setMethod("mu", "McmcChains", function(object) object@mu)
setMethod("tau2", "McmcChains", function(object) object@tau2)
setMethod("theta", "McmcChains", function(object) object@theta)
setMethod("sigma2", "McmcChains", function(object) object@sigma2)


##posteriorMean <- function(object) object@mean
##posteriorPrec <- function(object) object@prec
##posteriorProb <- function(object) object@prob

setMethod("show", "McmcChains", function(object){
  cat("An object of class 'McmcChains'\n")
  cat("    chain dim:", nrow(theta(object)), "x", ncol(theta(object)), "\n")
  cat("    see theta(), sigma2(), p(), ...\n")
})

setMethod("[", "McmcChains", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@theta <- x@theta[i, , drop=FALSE]
    x@sigma2 <- x@sigma2[i, , drop=FALSE]
    x@pi <- x@pi[i, , drop=FALSE]
    x@mu <- x@mu[i]
    x@tau2 <- x@tau2[i]
    x@nu.0 <- x@nu.0[i]
    x@sigma2.0 <- x@sigma2.0[i]
    x@logpotential <- x@logpotential[i]
  }
  x
})

setMethod("nu.0", "McmcChains", function(object) object@nu.0)

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

setReplaceMethod("logpotential", "McmcChains", function(object, value){
  object@logpotential <- value
  object
})

setMethod("names", "McmcChains", function(x) slotNames(x))
