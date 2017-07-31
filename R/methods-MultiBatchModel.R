setValidity("BatchModel", function(object){
  msg <- TRUE
  if(length(p(object)) != k(object)){
    msg <- "Mixture probability vector must be the same length as k"
    return(msg)
  }
  if(ncol(theta(object)) != k(object)){
    msg <- "theta matrix must have k columns"
    return(msg)
  }
  if(ncol(sigma(object)) != k(object)){
    msg <- "sigma matrix must have k columns"
    return(msg)
  }
  if(length(mu(object)) != k(object)){
    msg <- "mu vector must be length k "
    return(msg)
  }
  if(length(tau(object)) != k(object)){
    msg <- "tau vector must be length k "
    return(msg)
  }
  if(k(object) != k(hyperParams(object))){
    msg <- "k must be the same in the hyperparameters and in the model object"
    return(msg)
  }
  msg
})

#' Constructor for list of batch models
#'
#' An object of class BatchModel is constructed for each k, creating a list of
#' BatchModels.
#'
#' @param data numeric vector of average log R ratios
#' @param batch vector of batch labels
#' @param mcmc.params a \code{McmcParams} object
#' @param k numeric vector indicating the number of mixture components for each model
#' @param ... additional arguments to \code{HyperparametersBatch}
#' @return a list. Each element of the list is a \code{BatchModel}
#' @seealso \code{\link{BatchModel}}.  For single-batch data, use \code{\link{MarginalModelList}}.
#' @examples
#' mlist <- BatchModelList(data=y(BatchModelExample), k=1:4, batch=batch(BatchModelExample))
#' mcmcParams(mlist) <- McmcParams(iter=1, burnin=1, nStarts=0)
#' mlist2 <- posteriorSimulation(mlist)
#' @export
BatchModelList <- function(data=numeric(),
                           k=numeric(),
                           batch,
                           mcmc.params=McmcParams(),
                           ...){
  .Deprecated("See MultiBatchModelList")
  model.list <- vector("list", length(k))
  for(i in seq_along(k)){
    hypp <- HyperparametersBatch(k=k[i], ...)
    model.list[[i]] <- BatchModel(data=data, k=k[i], batch=batch,
                                  mcmc.params=mcmc.params,
                                  hypp=hypp)
  }
  model.list
}

#' Constructor for list of batch models
#'
#' An object of class MultiBatchModel is constructed for each k, creating a list of
#' BatchModels.
#'
#' @param data numeric vector of average log R ratios
#' @param batch vector of batch labels
#' @param mcmc.params a \code{McmcParams} object
#' @param k numeric vector indicating the number of mixture components for each model
#' @param ... additional arguments to \code{HyperparametersBatch}
#' @return a list. Each element of the list is a \code{BatchModel}
#' @seealso \code{\link{BatchModel}}.  For single-batch data, use \code{\link{MarginalModelList}}.
#' @examples
#' mlist <- BatchModelList(data=y(MultiBatchModelExample), k=1:4, batch=batch(MultiBatchModelExample))
#' mcmcParams(mlist) <- McmcParams(iter=1, burnin=1, nStarts=0)
#' mlist2 <- posteriorSimulation(mlist)
#' @export
MultiBatchModelList <- function(data=numeric(),
                           k=numeric(),
                           batch,
                           mcmc.params=McmcParams(),
                           ...){
  model.list <- vector("list", length(k))
  for(i in seq_along(k)){
    hypp <- HyperparametersBatch(k=k[i], ...)
    model.list[[i]] <- MultiBatchModel(data=data, k=k[i], batch=batch,
                                  mcmc.params=mcmc.params,
                                  hypp=hypp)
  }
  model.list
}

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
BatchModel <- function(data=numeric(),
                       k=3,
                       batch,
                       hypp,
                       mcmc.params){
  if(missing(batch)) batch <- as.integer(factor(rep("a", length(data))))
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  if(missing(hypp)) hypp <- HyperparametersBatch(k=k)
  if(missing(k) & !missing(hypp)){
    k <- k(hypp)
  }
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
    obj <- SingleBatchModel(data, k, hypp, mcmc.params)
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
             label_switch=FALSE,
             .internal.constraint=5e-4,
             .internal.counter=0L)
  obj <- startingValues(obj)
  obj
}

.empty_batch_model <- function(hp){
  K <- k(hp)
  B <- 0
  N <- 0
  obj <- new("BatchModel",
             k=as.integer(K),
             hyperparams=hp,
             theta=matrix(NA, 0, K),
             sigma2=matrix(NA, 0, K),
             mu=numeric(K),
             tau2=numeric(K),
             nu.0=numeric(1),
             sigma2.0=numeric(1),
             pi=numeric(K),
             data=numeric(0),
             data.mean=matrix(NA, B, K),
             data.prec=matrix(NA, B, K),
             z=integer(0),
             zfreq=integer(K),
             probz=matrix(0, N, K),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=McmcChains(),
             mcmc.params=mp,
             batch=integer(0),
             batchElements=integer(0),
             label_switch=FALSE,
             marginal_lik=as.numeric(NA),
             .internal.constraint=5e-4,
             .internal.counter=0L)
  chains(obj) <- McmcChains(obj)
  obj
}



#' Create an object for running hierarchical MCMC simulations.
#' @examples
#'      model <- MultiBatchModel(rnorm(10), k=1, batch=rep(1:2, each=5))
#' @param data the data for the simulation.
#' @param k An integer value specifying the number of latent classes.
#' @param batch a vector of the different batch numbers (must be sorted)
#' @param hypp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mcmc.params An object of class 'McmcParams'
#' @return An object of class `MultiBatchModel`
#' @export
MultiBatchModel <- function(data=numeric(),
                       k=3,
                       batch,
                       hypp,
                       mcmc.params){
  if(missing(batch)) batch <- as.integer(factor(rep("a", length(data))))
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  if(missing(hypp)) hypp <- HyperparametersBatch(k=k)
  if(missing(k) & !missing(hypp)){
    k <- k(hypp)
  }
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
    obj <- SingleBatchModel(data, k, hypp, mcmc.params)
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
  obj <- new("MultiBatchModel",
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
             label_switch=FALSE,
             .internal.constraint=5e-4,
             .internal.counter=0L)
  obj <- startingValues(obj)
  obj
}


MultiBatchModel2 <- function(dat=numeric(),
                            hp=HyperparametersBatch(),
                            mp=McmcParams(),
                            batches=integer()){
  if(length(dat) == 0){
    return(.empty_batch_model(hp))
  }
  ub <- unique(batches)
  nbatch <- setNames(as.integer(table(batches)), ub)
  B <- length(ub)
  N <- length(dat)
  ## move to setValidity
  if(length(dat) != length(batches)) {
    stop("batch vector must be the same length as data")
  }
  K <- k(hp)
  ## mu_k is the average across batches of the thetas for component k
  ## tau_k is the sd of the batch means for component k
  mu <- sort(rnorm(k(hp), mu.0(hp), sqrt(tau2.0(hp))))
  tau2 <- 1/rgamma(k(hp), 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
  p <- rdirichlet(1, alpha(hp))[1, ]
  sim_theta <- function(mu, tau, B) sort(rnorm(B, mu, tau))
  ##library(magrittr)
  thetas <- map2(mu, sqrt(tau2), sim_theta, B) %>%
    do.call(cbind, .) %>%
    apply(., 1, sort) %>%
    t
  if(K == 1) thetas <- t(thetas)
  nu.0 <- 3.5
  sigma2.0 <- 0.25
  sigma2s <- 1/rgamma(k(hp) * B, 0.5 * nu.0, 0.5 * nu.0 * sigma2.0) %>%
    matrix(B, k(hp))
  obj <- new("BatchModel",
             k=as.integer(K),
             hyperparams=hp,
             theta=thetas,
             sigma2=sigma2s,
             mu=mu,
             tau2=tau2,
             nu.0=nu.0,
             sigma2.0=sigma2.0,
             pi=p,
             data=dat,
             data.mean=matrix(NA, B, K),
             data.prec=matrix(NA, B, K),
             z=integer(N),
             zfreq=integer(K),
             probz=matrix(0, N, K),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=McmcChains(),
             mcmc.params=mp,
             batch=batches,
             batchElements=nbatch,
             label_switch=FALSE,
             marginal_lik=as.numeric(NA),
             .internal.constraint=5e-4,
             .internal.counter=0L)
  chains(obj) <- McmcChains(obj)
  obj
}


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
             probz=matrix(0, length(data), 1),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=mcmc.chains,
             mcmc.params=mcmc.params,
             batch=batch[ix],
             batchElements=nbatch,
             label_switch=FALSE,
             .internal.constraint=5e-4)
  obj <- startingValues(obj)
  if(!is.null(obj)){
    obj@probz[, 1] <- 1
  }
  obj
}

##
## This is the problem with calling new("BatchModel", ) with placeholders for values
## setValidity("BatchModel", function(object){
##   msg <- TRUE
##   ztab <- table(batch(object), z(object))
## ##  if(any(ztab < 1)){
## ##    msg <- "All components in each batch must have 1 or more observations"
## ##    return(msg)
## ##  }
## ##  if(ncol(ztab) != nrow(ztab)){
## ##    msg <- "All batches much have at least one observation from each component"
## ##    return(msg)
## ##  }
##   ## A valid batch model should have data ordered by batch
## ##  deltas <- diff(batch(object))
## ##  if(!all(deltas > 0)){
## ##    msg <- "Constructor for BatchModel should return data and batch assignment in batch-order"
## ##    return(msg)
##   ##  }
## ##  if(length(y(object)) > 0){
## ##    pz <- probz(object)
## ##    maxprob <- max(pz)
## ##    if(maxprob > 1 ){
## ##      msg <- "Posterior probabilities exceed 1"
## ##      return(msg)
## ##    }
## ##  }
##   msg
## })

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

#' extract data, latent variable, and batch for given observation
#' @name extract
#' @param x An object of class MultiBatchModel, McmcChains, or McmcParams
#' @param i An element of the instance to be extracted.
#' @param j Not used.
#' @param ... Not used.
#' @param drop Not used.
#' @return An object of class 'MultiBatchModel'
#' @aliases [,MultiBatchModel-method [,MultiBatchModel,ANY-method [,MultiBatchModel,ANY,ANY-method [,MultiBatchModel,ANY,ANY,ANY-method
#' @docType methods
#' @rdname extract-methods
setMethod("[", "MultiBatchModel", function(x, i, j, ..., drop=FALSE){
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

#' @rdname bic-method
#' @aliases bic,MultiBatchModel-method
setMethod("bic", "MultiBatchModel", function(object){
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

#' @rdname collapseBatch-method
#' @aliases collapseBatch,MultiBatchModel-method
setMethod("collapseBatch", "MultiBatchModel", function(object){
  collapseBatch(y(object), as.character(batch(object)))
})


batchLik <- function(x, p, mean, sd)  p*dnorm(x, mean, sd)


setMethod("computeMeans", "BatchModel", function(object){
  compute_means_batch(object)

})

setMethod("computeMeans", "MultiBatchModel", function(object){
  compute_means_batch(object)

})

setMethod("computePrec", "BatchModel", function(object){
  compute_prec_batch(object)
})

setMethod("computePrec", "MultiBatchModel", function(object){
  compute_prec_batch(object)
})

setMethod("computePrior", "BatchModel", function(object){
  compute_logprior_batch(object)
})

setMethod("computePrior", "MultiBatchModel", function(object){
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

setMethod("computeModes", "MultiBatchModel", function(object){
  .computeModesBatch(object)
})

componentVariances <- function(y, z)  v <- sapply(split(y, z), var)

setMethod("computeVars", "BatchModel", function(object){
  compute_vars_batch(object)
})

setMethod("computeVars", "MultiBatchModel", function(object){
  compute_vars_batch(object)
})

#' @rdname mu-method
#' @aliases mu,BatchModel-method
setMethod("mu", "BatchModel", function(object) object@mu)

#' @rdname mu-method
#' @aliases mu,MultiBatchModel-method
setMethod("mu", "MultiBatchModel", function(object) object@mu)

setReplaceMethod("mu", "BatchModel", function(object, value){
  object@mu <- value
  object
})

setReplaceMethod("mu", "MultiBatchModel", function(object, value){
  object@mu <- value
  object
})

nBatch <- function(object) length(uniqueBatch(object))

batchElements <- function(object) object@batchElements

setReplaceMethod("p", "BatchModel", function(object, value){
  object@pi <- value
  object
})

setReplaceMethod("p", "MultiBatchModel", function(object, value){
  object@pi <- value
  object
})

setMethod("pMean", "BatchModel", function(object) {
  mns <- colMeans(pic(object))
  mns
})

setMethod("pMean", "MultiBatchModel", function(object) {
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

setMethod("showMeans", "MultiBatchModel", function(object){
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

setMethod("showSigmas", "MultiBatchModel", function(object){
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

setReplaceMethod("sigma2", "MultiBatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})

#' @rdname sigma2-method
#' @aliases sigma2,BatchModel-method
setMethod("sigma2", "BatchModel", function(object) {
  s2 <- object@sigma2
  ##s2 <- matrix(s2, nBatch(object), k(object))
  rownames(s2) <- uniqueBatch(object)
  s2
})

#' @rdname sigma2-method
#' @aliases sigma2,MultiBatchModel-method
setMethod("sigma2", "MultiBatchModel", function(object) {
  s2 <- object@sigma2
  ##s2 <- matrix(s2, nBatch(object), k(object))
  rownames(s2) <- uniqueBatch(object)
  s2
})

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})

setMethod("tablez", "MultiBatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})

setMethod("sigmaMean", "BatchModel", function(object) {
  mns <- colMeans(sigmac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

setMethod("sigmaMean", "MultiBatchModel", function(object) {
  mns <- colMeans(sigmac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})


#' @rdname tau2-method
#' @aliases tau2,BatchModel-method
setMethod("tau2", "BatchModel", function(object) object@tau2)

#' @rdname tau2-method
#' @aliases tau2,BatchModel-method
setMethod("tau2", "MultiBatchModel", function(object) object@tau2)


setReplaceMethod("tau2", "BatchModel", function(object, value){
  object@tau2 <- value
  object
})

setReplaceMethod("tau2", "MultiBatchModel", function(object, value){
  object@tau2 <- value
  object
})

#' @rdname theta-method
#' @aliases theta,BatchModel-method
setMethod("theta", "BatchModel", function(object) {
  b <- object@theta
  ##b <- matrix(b, nBatch(object), k(object))
  rownames(b) <- uniqueBatch(object)
  b
})

#' @rdname theta-method
#' @aliases theta,MultiBatchModel-method
setMethod("theta", "MultiBatchModel", function(object) {
  b <- object@theta
  ##b <- matrix(b, nBatch(object), k(object))
  rownames(b) <- uniqueBatch(object)
  b
})


setReplaceMethod("theta", "BatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@theta <- value
  object
})


setReplaceMethod("theta", "MultiBatchModel", function(object, value){
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

setMethod("thetaMean", "MultiBatchModel", function(object) {
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

setMethod("show", "MultiBatchModel", function(object){
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

setMethod("tablez", "MultiBatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab <- tab[uniqueBatch(object), , drop=FALSE]
  tab
})


uniqueBatch <- function(object) unique(batch(object))


#' Create a data.frame of the component densities for each batch
#'
#' @param object an object of class \code{BatchModel}
#' @return a \code{{data.frame}}
#' @export
#' @examples
#'    nbatch <- 3
#'    k <- 3
#'    means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
#'        0, 0.05), nbatch, k, byrow = FALSE)
#'    sds <- matrix(0.15, nbatch, k)
#'    sds[, 1] <- 0.3
#'    N <- 1000
#'    truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
#'                                                  length.out = N),
#'                               p = c(1/10, 1/5, 1 - 0.1 - 0.2), theta = means,
#'                               sds = sds)
#'    mcmcp <- McmcParams(iter = 1000, burnin = 500, thin = 1,
#'                        nStarts = 10)
#'
#'    ## this parameter setting for m2.0 allows a lot of varation of the thetas
#'    ## between batch
#'    hypp <- CNPBayes:::HyperparametersBatch(m2.0 = 1/60, eta.0 = 1800,
#'                                            k = 3, a = 1/6, b = 180)
#'    model <- BatchModel(data = y(truth), batch = batch(truth),
#'                        k = 3, mcmc.params = mcmcp, hypp = hypp)
#'    model <- posteriorSimulation(model)
#'    df <- multiBatchDensities(model)
#'    df.observed <- data.frame(y=observed(model), batch=batch(model))
#'    library(ggplot2)
#'    ggplot(df, aes(x, d)) +
#'    geom_histogram(data=df.observed,
#'                   aes(y, ..density..),
#'                   bins=300, inherit.aes=FALSE) +
#'    geom_area(stat="identity", aes(color=name, fill=name),
#'              alpha=0.4) +
#'    xlab("quantiles") + ylab("density") +
#'    scale_color_manual(values=colors) +
#'    scale_fill_manual(values=colors) +
#'    guides(fill=guide_legend(""), color=guide_legend("")) +
#'    facet_wrap(~batch, nrow=2)
multiBatchDensities <- function(object){
  probs <- p(object)
  thetas <- theta(object)
  sigmas <- sigma(object)
  P <- matrix(probs, nrow(thetas), ncol(thetas), byrow=TRUE)
  rownames(P) <- uniqueBatch(object)
  avglrrs <- observed(object)
  quantiles <- seq(min(avglrrs), max(avglrrs), length.out=500)
  batchPr <- table(batch(object))/length(y(object))
  dens.list <- batchDensities(quantiles, uniqueBatch(object), 
                              thetas, sigmas, P, batchPr)
  ##component <- lapply(dens.list, rowSums)
  ##overall <- rowSums(do.call(cbind, component))
  ix <- order(thetas[1, ])
  d <- do.call(rbind, dens.list[ix])
  K <- ncol(thetas)
  NB <- nBatch(object)
  over <- Reduce("+", dens.list)
  batches.overall <- rep(1:2, each=nrow(over))
  quantile.overall <- rep(quantiles, 2)
  overall <- as.numeric(over)

  d.vec <- as.numeric(d, overall)
  d.vec <- c(d.vec, overall)
  batches <- c(rep(uniqueBatch(object), each=nrow(d)),
               batches.overall)
  K <- seq_len(ncol(thetas))
  name <- paste0("cn", K-1)
  name <- rep(rep(name, elementNROWS(dens.list)), 2)
  name <- c(name, rep("overall", length(overall)))
  x <- rep(rep(quantiles, length(dens.list)), 2)
  x <- c(x, quantile.overall)
  df <- data.frame(x=x, d=d.vec, name=name, batch=batches)
  df$batch <- factor(df$batch, uniqueBatch(object))
  df$name <- factor(df$name, levels=c("overall", paste0("cn", K-1)))
  df
}
