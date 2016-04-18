#' Create an object for running hierarchical MCMC simulations.
#' @examples
#'      model <- BatchModel(rnorm(10), k=1, batch=rep(1:2, each=5))
#' @param data the data for the simulation.
#' @param k An integer value specifying the number of latent classes.
#' @param batch a vector of the different batch numbers (must be sorted)
#' @param hypp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mcmc.params An object of class 'McmcParams'
#' @return An object of class `BatchModel`
#' @export
BatchModel <- function(data=numeric(), k=2L, batch, hypp, mcmc.params){
  if(missing(batch)) batch <- as.integer(factor(rep("a", length(data))))
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  mcmc.chains <- McmcChains()
  bf <- factor(batch)
  batch <- as.integer(bf)
  ub <- unique(batch)
  ##ix <- order(batch)
  ix <- seq_along(batch)
  nbatch <- setNames(as.integer(table(batch)), levels(bf))
  B <- length(ub)
  if(B==1 && length(data) > 0){
    if(missing(hypp)) hypp <- HyperparametersMarginal(k=k)
    zz <- as.integer(factor(numeric(k)))
    zfreq <- as.integer(table(zz))
    obj <- MarginalModel(data, k, hypp, mcmc.params)
    return(obj)
  }
  if(k == 1) {
    if(missing(hypp)) hypp <- HyperparametersBatch(k=1)
    obj <- UnivariateBatchModel(data, k, batch, hypp, mcmc.params)
    return(obj)
  }
  if(missing(hypp)) hypp <- HyperparametersBatch(k=k)
  zz <- integer(length(data))
  zfreq <- as.integer(table(zz))
  if(length(data) != length(batch)) {
    stop("batch vector must be the same length as data")
  }
  obj <- new("BatchModel",
             k=as.integer(k),
             hyperparams=hypp,
             theta=matrix(NA, B, k),
             sigma2=matrix(NA, B, k),
             mu=numeric(k),
             tau2=numeric(k),
             nu.0=numeric(1),
             sigma2.0=numeric(1),
             pi=numeric(k),
             data=data[ix],
             data.mean=matrix(NA, B, k),
             data.prec=matrix(NA, B, k),
             z=zz,
             zfreq=zfreq,
             probz=matrix(0, length(data), k),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=mcmc.chains,
             mcmc.params=mcmc.params,
             batch=batch[ix],
             batchElements=nbatch,
             .internal.constraint=5e-4,
             .internal.counter=0L)
  obj <- startingValues(obj)
  obj <- ensureAllComponentsObserved(obj)
  obj
}

## BatchModel2 <- function(data=numeric(), k=2L, batch, hypp, mcmc.params){
##   if(missing(batch)) batch <- as.integer(factor(rep("a", length(data))))
##   if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
##   mcmc.chains <- McmcChains()
##   bf <- factor(batch)
##   batch <- as.integer(bf)
##   ub <- unique(batch)
##   ##ix <- order(batch)
##   ix <- seq_along(batch)
##   nbatch <- setNames(as.integer(table(batch)), levels(bf))
##   B <- length(ub)
##   if(B==1 && length(data) > 0){
##     if(missing(hypp)) hypp <- HyperparametersMarginal(k=k)
##     zz <- as.integer(factor(numeric(k)))
##     zfreq <- as.integer(table(zz))
##     obj <- MarginalModel(data, k, hypp, mcmc.params)
##     return(obj)
##   }
##   if(k == 1) {
##     if(missing(hypp)) hypp <- HyperparametersBatch(k=1)
##     obj <- UnivariateBatchModel(data, k, batch, hypp, mcmc.params)
##     return(obj)
##   }
##   if(missing(hypp)) hypp <- HyperparametersBatch(k=k)
##   zz <- integer(length(data))
##   zfreq <- as.integer(table(zz))
##   if(length(data) != length(batch)) {
##     stop("batch vector must be the same length as data")
##   }
##   obj <- new("BatchModel",
##              k=as.integer(k),
##              hyperparams=hypp,
##              theta=matrix(NA, B, k),
##              sigma2=matrix(NA, B, k),
##              mu=numeric(k),
##              tau2=numeric(k),
##              nu.0=numeric(1),
##              sigma2.0=numeric(1),
##              pi=numeric(k),
##              data=data[ix],
##              data.mean=matrix(NA, B, k),
##              data.prec=matrix(NA, B, k),
##              z=zz,
##              zfreq=zfreq,
##              probz=matrix(0, length(data), k),
##              logprior=numeric(1),
##              loglik=numeric(1),
##              mcmc.chains=mcmc.chains,
##              mcmc.params=mcmc.params,
##              batch=batch[ix],
##              batchElements=nbatch,
##              .internal.constraint=5e-4,
##              .internal.counter=0L)
##   obj <- startingValues(obj)
##   obj <- ensureAllComponentsObserved(obj)
##   obj
## }

ensureAllComponentsObserved <- function(obj){
  zz <- table(batch(obj), z(obj))
  K <- seq_len(k(obj))
  if(any(zz<=1)){
    index <- which(rowSums(zz<=1) > 0)
    for(i in seq_along(index)){
      j <- index[i]
      zup <- z(obj)[batch(obj) == j]
      zfact <- factor(zup, levels=K)
      minz <- as.integer(names(table(zfact))[table(zfact) <= 1])
      ##missingk <- K[!K %in% unique(zup)]
      maxk <- names(table(zfact))[which.max(table(zfact))]
      nreplace <- length(minz)*2
      zup[sample(which(zup == maxk), nreplace)] <- as.integer(minz)
      obj@z[batch(obj) == j] <- as.integer(zup)
    }
  }
  obj
}

##
## Multiple batches, but only 1 component
##
UnivariateBatchModel <- function(data, k=1, batch, hypp, mcmc.params){
  mcmc.chains <- McmcChains()
  zz <- integer(length(data))
  zfreq <- as.integer(table(zz))
  bf <- factor(batch)
  batch <- as.integer(bf)
  ix <- order(bf)
  ##B <- length(levels(batch))
  nbatch <- elementNROWS(split(batch, batch))
  B <- length(nbatch)
  if(missing(hypp)) hypp <- HyperparametersBatch(k=1)
  obj <- new("UnivariateBatchModel",
             k=as.integer(k),
             hyperparams=hypp,
             theta=matrix(NA, B, 1),
             sigma2=matrix(NA, B, 1),
             mu=numeric(k),
             tau2=numeric(k),
             nu.0=numeric(1),
             sigma2.0=numeric(1),
             pi=numeric(k),
             ##pi=matrix(NA, B, 1),
             data=data[ix],
             data.mean=matrix(NA, B, 1),
             data.prec=matrix(NA, B, 1),
             z=integer(length(data)),
             zfreq=zfreq,
             probz=matrix(1, length(data), 1),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=mcmc.chains,
             mcmc.params=mcmc.params,
             batch=batch[ix],
             batchElements=nbatch,
             .internal.constraint=5e-4)
  obj <- startingValues(obj)
  obj
}

setValidity("BatchModel", function(object){
  msg <- TRUE
  ztab <- table(batch(object), z(object))
  if(any(ztab < 1)){
    msg <- "All components in each batch must have 1 or more observations"
    return(msg)
  }
  msg
})

#' extract data, latent variable, and batch for given observation
#' @name extract
#' @param x An object of class BatchModel, McmcChains, or McmcParams
#' @param i An element of the instance to be extracted.
#' @param j Not used.
#' @param ... Not used.
#' @param drop Not used.
#' @return An object of class 'BatchModel'
#' @aliases [,BatchModel-method [,BatchModel,ANY-method [,BatchModel,ANY,ANY-method [,BatchModel,ANY,ANY,ANY-method
#' @docType methods
#' @rdname extract-methods
setMethod("[", "BatchModel", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    y(x) <- y(x)[i]
    z(x) <- z(x)[i]
    batch(x) <- batch(x)[i]
  }
  x
})

#' @rdname bic-method
#' @aliases bic,BatchModel-method
setMethod("bic", "BatchModel", function(object){
  object <- useModes(object)
  ## K: number of free parameters to be estimated
  ##   - component and batch-specific parameters:  theta, sigma2  ( k(model) * nBatch(model))
  ##   - mixing probabilities: (k-1)*nBatch
  ##   - component-specific parameters: mu, tau2                 2 x k(model)
  ##   - length-one parameters: sigma2.0, nu.0                   +2
  K <- 2*k(object)*nBatch(object) + (k(object)-1) + 2*k(object) + 2
  n <- length(y(object))
  bicstat <- -2*(log_lik(object) + logPrior(object)) + K*(log(n) - log(2*pi))
  bicstat
})

#' @rdname collapseBatch-method
#' @aliases collapseBatch,BatchModel-method
setMethod("collapseBatch", "BatchModel", function(object){
  collapseBatch(y(object), as.character(batch(object)))
})


batchLik <- function(x, p, mean, sd)  p*dnorm(x, mean, sd)


setMethod("computeMeans", "BatchModel", function(object){
  compute_means_batch(object)

})

setMethod("computePrec", "BatchModel", function(object){
  compute_prec_batch(object)
})

setMethod("computePrior", "BatchModel", function(object){
  compute_logprior_batch(object)
})

.computeModesBatch <- function(object){
  i <- argMax(object)
  mc <- chains(object)
  B <- nBatch(object)
  K <- k(object)
  thetamax <- matrix(theta(mc)[i, ], B, K)
  sigma2max <- matrix(sigma2(mc)[i, ], B, K)
  pmax <- p(mc)[i, ]
  mumax <- mu(mc)[i, ]
  tau2max <- tau2(mc)[i,]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mumax,
                tau2=tau2max,
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i],
                zfreq=zFreq(mc)[i, ],
                loglik=log_lik(mc)[i],
                logprior=logPrior(mc)[i])
  modes
}

setMethod("computeModes", "BatchModel", function(object){
  .computeModesBatch(object)
})

componentVariances <- function(y, z)  v <- sapply(split(y, z), var)

setMethod("computeVars", "BatchModel", function(object){
  compute_vars_batch(object)
})

#' @rdname mu-method
#' @aliases mu,BatchModel-method
setMethod("mu", "BatchModel", function(object) object@mu)

setReplaceMethod("mu", "BatchModel", function(object, value){
  object@mu <- value
  object
})

nBatch <- function(object) length(uniqueBatch(object))

batchElements <- function(object) object@batchElements

setReplaceMethod("p", "BatchModel", function(object, value){
  object@pi <- value
  object
})

setMethod("pMean", "BatchModel", function(object) {
  mns <- colMeans(pic(object))
  mns
})

setMethod("showMeans", "BatchModel", function(object){
  thetas <- round(theta(object), 2)
  mns <- c("\n", paste0(t(cbind(thetas, "\n")), collapse="\t"))
  mns <- paste0("\t", mns[2])
  mns <- paste0("\n", mns[1])
  mns
})

setMethod("showSigmas", "BatchModel", function(object){
  sigmas <- round(sqrt(sigma2(object)), 2)
  sigmas <- c("\n", paste0(t(cbind(sigmas, "\n")), collapse="\t"))
  sigmas <- paste0("\t", sigmas[2])
  sigmas <- paste0("\n", sigmas[1])
  sigmas
})


setReplaceMethod("sigma2", "BatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})

#' @rdname sigma2-method
#' @aliases sigma2,BatchModel-method
setMethod("sigma2", "BatchModel", function(object) {
  s2 <- object@sigma2
  s2 <- matrix(s2, nBatch(object), k(object))
  rownames(s2) <- uniqueBatch(object)
  s2
})

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})

setMethod("sigmaMean", "BatchModel", function(object) {
  mns <- colMeans(sigmac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

#' @rdname tau2-method
#' @aliases tau2,BatchModel-method
setMethod("tau2", "BatchModel", function(object) object@tau2)

setReplaceMethod("tau2", "BatchModel", function(object, value){
  object@tau2 <- value
  object
})

#' @rdname theta-method
#' @aliases theta,BatchModel-method
setMethod("theta", "BatchModel", function(object) {
  b <- object@theta
  b <- matrix(b, nBatch(object), k(object))
  rownames(b) <- uniqueBatch(object)
  b
})

setReplaceMethod("theta", "BatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@theta <- value
  object
})



setMethod("thetaMean", "BatchModel", function(object) {
  mns <- colMeans(thetac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})


setMethod("show", "BatchModel", function(object){
  ##callNextMethod()
  cls <- class(object)
  cat(paste0("An object of class ", cls), "\n")
  cat("     n. obs      :", length(y(object)), "\n")
  cat("     n. batches  :", nBatch(object), "\n")
  cat("     k           :", k(object), "\n")
  cat("     nobs/batch  :", table(batch(object)), "\n")
  cat("     loglik (s)  :", round(log_lik(object), 1), "\n")
  cat("     logprior (s):", round(logPrior(object), 1), "\n")
})

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab <- tab[uniqueBatch(object), , drop=FALSE]
  tab
})

uniqueBatch <- function(object) unique(batch(object))
