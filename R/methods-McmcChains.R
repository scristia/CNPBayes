## number comonents (K)
## number MCMC (S)
## number batches (B)
initialize_mcmc <- function(K, S, B){
  K <- as.integer(K)
  B <- as.integer(B)
  S <- as.integer(S)
  new("McmcChains",
      theta=matrix(numeric(), S, K*B),
      sigma2=matrix(numeric(), S, K*B),
      pi=matrix(numeric(), S, K*B),
      mu=matrix(numeric(), S, K),
      tau2=matrix(numeric(), S, K),
      nu.0=numeric(S),
      sigma2.0=numeric(S),
      logprior=numeric(S),
      loglik=numeric(S),
      zfreq=matrix(as.integer(NA), S, K),
      predictive=matrix(as.numeric(NA), S, K*B),
      zstar=matrix(as.integer(NA), S, K*B),
      iter=S,
      k=K,
      B=B)
}

## pooled
initialize_mcmcP <- function(K, S, B){
  K <- as.integer(K)
  B <- as.integer(B)
  S <- as.integer(S)
  new("McmcChains",
      theta=matrix(numeric(), S, K*B),
      sigma2=matrix(numeric(), S, B),
      pi=matrix(numeric(), S, K*B),
      mu=matrix(numeric(), S, K),
      tau2=matrix(numeric(), S, K),
      nu.0=numeric(S),
      sigma2.0=numeric(S),
      logprior=numeric(S),
      loglik=numeric(S),
      zfreq=matrix(as.integer(NA), S, K),
      predictive=matrix(as.numeric(NA), S, K*B),
      zstar=matrix(as.integer(NA), S, K*B),
      iter=S,
      k=K,
      B=B)
}

# trios
initialize_mcmcT <- function(K, S, B, T){
  K <- as.integer(K)
  B <- as.integer(B)
  S <- as.integer(S)
  new("McmcChainsTrios",
      theta=matrix(numeric(), S, K*B),
      sigma2=matrix(numeric(), S, K*B),
      pi=matrix(numeric(), S, K),
      pi_parents=matrix(numeric(), S, K),
      mu=matrix(numeric(), S, K),
      tau2=matrix(numeric(), S, K),
      nu.0=numeric(S),
      sigma2.0=numeric(S),
      logprior=numeric(S),
      loglik=numeric(S),
      zfreq=matrix(as.integer(NA), S, K),
      zfreq_parents=matrix(as.integer(NA), S, K),
      predictive=matrix(as.numeric(NA), S, K*B),
      zstar=matrix(as.integer(NA), S, K*B),
      is_mendelian=integer(T),
      iter=S,
      k=K,
      B=B)
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
      zfreq=mati,
      predictive=matrix(as.numeric(NA), nr, K*1),
      zstar=matrix(as.integer(NA), nr, K*1),
      iter=iter(mcmc.params),
      k=k(object),
      B=nBatch(object))
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
      zfreq=mati,
      predictive=matrix(as.numeric(NA), nr, K),
      zstar=matrix(as.integer(NA), nr, K),
      iter=iter(mcmc.params),
      k=k(object),
      B=nBatch(object))
}

setMethod("McmcChains", "missing", function(object){
  new("McmcChains",
      theta=matrix(),
      sigma2=matrix(),
      pi=matrix(),
      mu=numeric(),
      tau2=numeric(),
      nu.0=numeric(),
      sigma2.0=numeric(),
      zfreq=matrix(),
      logprior=numeric(),
      loglik=numeric(),
      predictive=matrix(),
      zstar=matrix(),
      iter=integer(),
      k=integer(),
      B=integer())
})

setMethod("McmcChainsTrios", "missing", function(object){
  initialize_mcmcT(0, 0, 0, 0)
})

setValidity("McmcChains", function(object){
  msg <- TRUE
  if(length(iter(object)) > 0){
    if(iter(object) != nrow(predictive(object))){
      msg <- "predictive slot has incorrect dimension"
      return(msg)
    }
  }
  msg
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
      pi=matrix(NA, nr, K*B),
      mu=matrix(NA, nr, K),
      tau2=matrix(NA, nr, K),
      nu.0=numeric(nr),
      sigma2.0=numeric(nr),
      logprior=numeric(nr),
      loglik=numeric(nr),
      zfreq=mati,
      predictive=matrix(as.numeric(NA), nr, K*B),
      zstar=matrix(as.integer(NA), nr, K*B),
      iter=iter(object),
      k=k(object),
      B=nBatch(object))

}


setMethod("McmcChains", "MultiBatchModel", function(object){
  .initializeMcmcBatch(object)
})

.initializeMcmcTrios <- function(object){
  mcmc.params <- mcmcParams(object)
  S <- iter(mcmc.params)[1]
  ns <- length(y(object))
  K <- k(object)
  B <- nBatch(object)
  dat <- triodata(object)
  T <- length(unique(dat$id))
  initialize_mcmcT(K, S, B, T)
}

setMethod("McmcChainsTrios", "TrioBatchModel", function(object){
  K <- k(object)
  S <- iter(object)
  B <- numBatch(object)
  T <- nTrios(object)
  initialize_mcmcT(K, S, B, T)
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
      pi=matrix(NA, nr, K*B),
      mu=matrix(NA, nr, K),
      tau2=matrix(NA, nr, K),
      nu.0=numeric(nr),
      sigma2.0=numeric(nr),
      logprior=numeric(nr),
      loglik=numeric(nr),
      zfreq=mati,
      predictive=matrix(as.numeric(NA), nr, K*B),
      zstar=matrix(as.integer(NA), nr, K*B),
      iter=iter(object),
      k=k(object),
      B=nBatch(object))
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

#' @aliases sigma2,missing-method
setMethod("sigma_", "McmcChains", function(object) sqrt(object@sigma2))

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
    x@predictive <- x@predictive[i, , drop=FALSE]
    x@zstar <- x@zstar[i, , drop=FALSE]
  }
  x@iter <- nrow(x@theta)
  x
})

#' extract estimated parameters at particular iteration of simulation.
#' @aliases [,McmcChainsTrios-method [,McmcChainsTrios,ANY-method
#' @return An object of class 'McmcChainsTrios'
#' @docType methods
#' @rdname extract-methods
setMethod("[", "McmcChainsTrios", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@theta <- x@theta[i, , drop=FALSE]
    x@sigma2 <- x@sigma2[i, , drop=FALSE]
    x@pi <- x@pi[i, , drop=FALSE]
    x@pi_parents <- x@pi_parents[i, , drop=FALSE]
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
    x@predictive <- x@predictive[i, , drop=FALSE]
    x@zstar <- x@zstar[i, , drop=FALSE]
  }
  x@iter <- nrow(x@theta)
  x
})

#' @rdname nu.0-method
#' @aliases nu.0,McmcChains-method
setMethod("nu.0", "McmcChains", function(object) object@nu.0)

#' @rdname sigma2.0-method
#' @aliases sigma2.0,McmcChains-method
setMethod("sigma2.0", "McmcChains", function(object) object@sigma2.0)

setReplaceMethod("pp", "McmcChains", function(object, value){
  object@pi_parents <- value
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

#' @rdname zfreqpar-method
#' @aliases zfreqpar,McmcChains-method
setMethod("zFreqPar", "McmcChains", function(object) object@zfreq_parents )

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

longFormatKB <- function(x, K, B){
  col_names <- expand.grid(seq_len(B), seq_len(K)) %>%
    mutate(col_names=paste(Var1, Var2, sep=",")) %$%
    col_names
  ##col_names <- col_names[ !duplicated(col_names) ]
  x <- x %>%
    as_tibble %>%
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

longFormatKB2 <- function(x, K, B){
  col_names <- rep(seq_len(B), B) %>%
    paste(rep(seq_len(K), each=K), sep=",")
  col_names <- col_names[ !duplicated(col_names) ]
  x <- x %>%
    as_tibble %>%
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
    as_tibble %>%
    set_colnames(col_names) %>%
    mutate(s=seq_len(nrow(.))) %>%
    gather("k", "value", -s) %>%
    mutate(k=factor(paste("k ", k)))
  x
}

setMethod("k", "McmcChains",  function(object) object@k)
setMethod("iter", "McmcChains",  function(object) object@iter)
setReplaceMethod("iter", "McmcChains",  function(object, value){
  object@iter <- value
  object
})
setReplaceMethod("k", "McmcChains",  function(object, value){
  object@k <- value
  object
})
setMethod("numBatch", "McmcChains",  function(object) object@B)
setReplaceMethod("numBatch", "McmcChains",  function(object, value){
  object@B <- value
})

setAs("McmcChains", "list", function(from){
  K <- k(from)
  B <- from@B
  S <- iter(from)
  theta <- longFormatKB(theta(from), K, B)
  if(ncol(sigma2(from)) == K*B){
    sigma2 <- longFormatKB(sigma2(from), K, B)
  } else {
    sigma2 <- longFormatK(sigma2(from), K)
  }
  p <- longFormatKB(p(from), K, B)
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


setMethod("predictive", "McmcChains", function(object) object@predictive)
setMethod("zstar", "McmcChains", function(object) object@zstar)
setMethod("predictive", "McmcChainsTrios", function(object) object@predictive)
setMethod("zstar", "McmcChainsTrios", function(object) object@zstar)
setMethod("predictive", "MultiBatchModel", function(object) predictive(chains(object)))
setMethod("zstar", "MultiBatchModel", function(object) zstar(chains(object)))
setMethod("predictive", "TrioBatchModel", function(object) predictive(chains(object)))
setMethod("zstar", "TrioBatchModel", function(object) zstar(chains(object)))
setMethod("predictive", "MultiBatchPooled", function(object) predictive(chains(object)))
setMethod("zstar", "MultiBatchPooled", function(object) zstar(chains(object)))


setReplaceMethod("predictive", c("McmcChains", "matrix"), function(object, value) {
  object@predictive <- value
  object
})

setReplaceMethod("predictive", c("McmcChainsTrios", "matrix"), function(object, value) {
  object@predictive <- value
  object
})


setMethod("updateObject", "McmcChains",
          function(object, verbose=FALSE){
            object <- callNextMethod(object)
            K <- k(object)
            S <- iter(object)
            B <- numBatch(object)
            predictive(object) <- matrix(as.numeric(NA),
                                         nrow=S,
                                         ncol=K*B)
            object@zstar <- matrix(as.integer(NA),
                                   nrow=S,
                                   ncol=K*B)
            object@k <- K
            object@iter <- S
            object@B <- B
            object
})

setMethod("isMendelian", "McmcChains", function(object) object@is_mendelian)
