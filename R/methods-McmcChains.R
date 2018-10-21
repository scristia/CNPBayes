## number comonents (K)
## number MCMC (S)
## number batches (B)
initialize_mcmc <- function(K, S, B){
  new("McmcChains",
      theta=matrix(numeric(), S, K*B),
      sigma2=matrix(numeric(), S, K*B),
      pi=matrix(numeric(), S, K),
      mu=matrix(numeric(), S, K),
      tau2=matrix(numeric(), S, K),
      nu.0=numeric(S),
      sigma2.0=numeric(S),
      logprior=numeric(S),
      loglik=numeric(S),
      zfreq=matrix(as.integer(NA), S, K))
}

## pooled
initialize_mcmcP <- function(K, S, B){
  new("McmcChains",
      theta=matrix(numeric(), S, K*B),
      sigma2=matrix(numeric(), S, B),
      pi=matrix(numeric(), S, K),
      mu=matrix(numeric(), S, K),
      tau2=matrix(numeric(), S, K),
      nu.0=numeric(S),
      sigma2.0=numeric(S),
      logprior=numeric(S),
      loglik=numeric(S),
      zfreq=matrix(as.integer(NA), S, K))
}

.initializeMcmc <- function(object){
  ## add 1 for starting values (either the last run from the burnin,
  ## or default values if no burnin
  mcmc.params <- mcmcParams(object)
  nr <- iter(mcmc.params)
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
      zfreq=mati)
}

.initializeMcmcPooledVar <- function(object){
  ## add 1 for starting values (either the last run from the burnin,
  ## or default values if no burnin
  mcmc.params <- mcmcParams(object)
  nr <- iter(mcmc.params)
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
  new("McmcChains", theta=matrix(), sigma2=matrix(),
      pi=matrix(), mu=numeric(), tau2=numeric(),
      nu.0=numeric(), sigma2.0=numeric(),
      zfreq=matrix(),
      logprior=numeric(),
      loglik=numeric())
})


setMethod("nu.0", "McmcChains", function(object) object@nu.0)

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

setMethod("p", "McmcChains", function(object){
  object@pi
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

McmcChains2 <- function(mc, iter, k, batch){
  new("McmcChains2",
      iter=iter,
      k=k,
      batch=batch,
      theta=theta(mc),
      sigma2=sigma2(mc),
      pi=p(mc),
      mu=mu(mc),
      tau2=tau2(mc),
      nu.0=nu.0(mc),
      sigma2.0=sigma2.0(mc),
      logprior=logPrior(mc),
      loglik=log_lik(mc),
      zfreq=zFreq(mc))
}

longFormatKB <- function(x, K, B){
  col_names <- rep(seq_len(B), B) %>%
    paste(rep(seq_len(K), each=K), sep=",")
  x <- x %>%
    as.tibble %>%
    set_colnames(col_names) %>%
    mutate(s=seq_len(nrow(.))) %>%
    gather("bk", "value", -s) %>%
    mutate(b=sapply(strsplit(bk, ","), "[", 1),
           k=sapply(strsplit(bk, ","), "[", 2)) %>%
    mutate(b=factor(paste("batch", b)),
           k=factor(paste("k", k))) %>%
    select(-bk)
  x
}

longFormatK <- function(x, K){
  col_names <- seq_len(K) %>%
    as.character
  x <- x %>%
    as.tibble %>%
    set_colnames(col_names) %>%
    mutate(s=seq_len(nrow(.))) %>%
    gather("k", "value", -s) %>%
    mutate(k=factor(paste("k ", k)))
  x
}

setAs("McmcChains2", "list", function(from){
  K <- from@k
  B <- from@batch
  S <- from@iter
  theta <- longFormatKB(theta(from), K, B)
  sigma2 <- longFormatKB(sigma2(from), K, B)
  p <- longFormatK(p(from), K)
  mu <- longFormatK(mu(from), K)
  tau2 <- longFormatK(tau2(from), K)
  zfreq <- longFormatK(zFreq(from), K)
  params <- tibble(s=seq_len(S),
                   nu.0=nu.0(from),
                   sigma2.0=sigma2.0(from),
                   logprior=logPrior(from),
                   loglik=log_lik(from)) %>%
    gather("parameter", "value", -s)
  list(theta=theta,
       sigma2=sigma2,
       p=p,
       mu=mu,
       tau2=tau2,
       zfreq=zfreq,
       scalars=params)
})

setMethod("listChains", "MultiBatch", function(object){
  ch.list <- McmcChains2(mc=chains(object),
                     iter=iter(object),
                     k=k(object),
                     batch=nBatch(object)) %>%
    as("list")
  ch.list
})
