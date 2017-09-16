### From AllClasses.R

#' An object for running MCMC simulations.
#'
#' Run hierarchical MCMC for batch model.
#' @slot k An integer value specifying the number of latent classes.
#' @slot hyperparams An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @slot theta the means of each component and batch
#' @slot sigma2 the variances of each component and batch
#' @slot nu.0 the shape parameter for sigma2
#' @slot sigma2.0 the rate parameter for sigma2
#' @slot pi mixture probabilities which are assumed to be the same for all batches
#' @slot mu means from batches, averaged across batches
#' @slot tau2 variances from batches,  weighted by precisions
#' @slot data the data for the simulation.
#' @slot data.mean the empirical means of the components
#' @slot data.prec the empirical precisions
#' @slot z latent variables
#' @slot zfreq table of latent variables
#' @slot probz n x k matrix of probabilities
#' @slot logprior log likelihood of prior: log(p(sigma2.0)p(nu.0)p(mu))
#' @slot loglik log likelihood: \eqn{\sum p_k \Phi(\theta_k, \sigma_k)}
#' @slot mcmc.chains an object of class 'McmcChains' to store MCMC samples
#' @slot batch a vector of the different batch numbers
#' @slot batchElements a vector labeling from which batch each observation came from
#' @slot modes the values of parameters from the iteration which maximizes log likelihood and log prior
#' @slot label_switch length-one logical vector indicating whether label-switching occurs (possibly an overfit model)
#' @slot mcmc.params An object of class 'McmcParams'
#' @slot .internal.constraint Constraint on parameters. For internal use only.
setClass("BatchModel", contains="MixtureModel")

#' The 'MarginalModel' class
#'
#' Run marginal MCMC simulation
#' @slot k An integer value specifying the number of latent classes.
#' @slot hyperparams An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @slot theta the means of each component and batch
#' @slot sigma2 the variances of each component and batch
#' @slot nu.0 the shape parameter for sigma2
#' @slot sigma2.0 the rate parameter for sigma2
#' @slot pi mixture probabilities which are assumed to be the same for all batches
#' @slot mu overall mean
#' @slot tau2 overall variance
#' @slot data the data for the simulation.
#' @slot data.mean the empirical means of the components
#' @slot data.prec the empirical precisions
#' @slot z latent variables
#' @slot zfreq table of latent variables
#' @slot probz n x k matrix of probabilities
#' @slot logprior log likelihood of prior: log(p(sigma2.0)p(nu.0)p(mu))
#' @slot loglik log likelihood: \eqn{\sum p_k \Phi(\theta_k, \sigma_k)}
#' @slot mcmc.chains an object of class 'McmcChains' to store MCMC samples
#' @slot batch a vector of the different batch numbers
#' @slot batchElements a vector labeling from which batch each observation came from
#' @slot modes the values of parameters from the iteration which maximizes log likelihood and log prior
#' @slot mcmc.params An object of class 'McmcParams'
#' @slot label_switch length-one logical vector indicating whether label-switching occurs (possibly an overfit model)
#' @slot .internal.constraint Constraint on parameters. For internal use only.
setClass("MarginalModel", contains="MixtureModel")


#' Constructor for list of single-batch models
#'
#' An object of class MarginalModel is constructed for each k, creating a list of
#' MarginalModels.
#'
#' @param data numeric vector of average log R ratios
#' @param k numeric vector indicating the number of mixture components for each model
#' @param mcmc.params an object of class \code{McmcParams}
#' @param ... additional arguments passed to \code{Hyperparameters}
#' @seealso \code{\link{MarginalModel}} \code{\link{BatchModelList}}
#' @return a list. Each element of the list is a \code{BatchModel}
#' @examples
#' mlist <- MarginalModelList(data=y(MarginalModelExample), k=1:4)
#' mcmcParams(mlist) <- McmcParams(iter=1, burnin=1, nStarts=0)
#' mlist2 <- posteriorSimulation(mlist)
#' @export
MarginalModelList <- function(data=numeric(), k=numeric(),
                              mcmc.params=McmcParams(),
                              ...){
  .Deprecated("See SingleBatchModelList")
  model.list <- vector("list", length(k))
  for(i in seq_along(k)){
    hypp <- Hyperparameters(k=k[i], ...)
    model.list[[i]] <- MarginalModel(data=data, k=k[i], mcmc.params=mcmc.params,
                                     hypp=hypp)
  }
  model.list
}

#' Create an object for running marginal MCMC simulations.
#' @examples
#'      model <- MarginalModel(data=rnorm(10), k=1)
#' @param data the data for the simulation.
#' @param k An integer value specifying the number of latent classes.
#' @param hypp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mcmc.params An object of class 'McmcParams'
#' @return An object of class 'MarginalModel'
#' @export
MarginalModel <- function(data=numeric(), k=3, hypp, mcmc.params){
  .Deprecated("See SingleBatchModel")
  if(missing(hypp)) hypp <- Hyperparameters(k=3)
  batch <- rep(1L, length(data))
  if(missing(k)){
    k <- k(hypp)
  } else{
    k(hypp) <- k
  }
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  if(missing(hypp)) hypp <- HyperparametersMarginal(k=k)
  nbatch <- setNames(as.integer(table(batch)), levels(batch))
  zz <- sample(seq_len(k), length(data), replace=TRUE)
  zfreq <- as.integer(table(zz))
  object <- new("MarginalModel",
                k=as.integer(k),
                hyperparams=hypp,
                theta=numeric(k),
                sigma2=numeric(k),
                mu=numeric(1),
                tau2=numeric(1),
                nu.0=1L,
                sigma2.0=1L,
                pi=rep(1/k, k),
                data=data,
                data.mean=numeric(k),
                data.prec=numeric(k),
                z=zz,
                zfreq=zfreq,
                probz=matrix(0, length(data), k),
                logprior=numeric(1),
                loglik=numeric(1),
                mcmc.chains=McmcChains(),
                batch=batch,
                batchElements=nbatch,
                modes=list(),
                mcmc.params=mcmc.params,
                label_switch=FALSE,
                marginal_lik=as.numeric(NA),
                .internal.constraint=5e-4,
                .internal.counter=0L)
  object <- startingValues(object)
}

#' @rdname mu-method
#' @aliases mu,MarginalModel-method
#' @export
setMethod("mu", "MarginalModel", function(object){
  .Deprecated("see SingleBatchModel")
  object@mu
})

#' @rdname tau2-method
#' @aliases tau2,MarginalModel-method
setMethod("tau2", "MarginalModel", function(object) object@tau2)

setMethod("show", "MarginalModel", function(object) callNextMethod())

setMethod("computeMeans", "MarginalModel", function(object){
  .Deprecated("See SingleBatchModel")
  compute_means(object)
})

setMethod("computeVars", "MarginalModel", function(object){
  .Deprecated("See SingleBatchModel")
  compute_vars(object)
})

setReplaceMethod("tau2", "MarginalModel", function(object, value){
  .Deprecated("See SingleBatchModel")
  object@tau2 <- value
  object
})

setReplaceMethod("mu", "MarginalModel", function(object, value){
  .Deprecated("See SingleBatchModel")
  object@mu <- value
  object
})


#' @rdname bic-method
#' @aliases bic,MarginalModel-method
setMethod("bic", "MarginalModel", function(object){
  .Deprecated("See SingleBatchModel")
  object <- useModes(object)
  ## K: number of free parameters to be estimated
  ##   - component-specific parameters:  theta, sigma2   (3 x k(model))
  ##   - mixing probabilities:  k-1
  ##   - length-one parameters: mu, tau2, sigma2.0, nu.0             +4
  K <- 2*k(object) + (k(object)-1) + 4
  n <- length(y(object))
  -2*(log_lik(object) + logPrior(object)) + K*(log(n) - log(2*pi))
})

#' @rdname theta-method
#' @aliases theta,MarginalModel-method
setMethod("theta", "MarginalModel", function(object) object@theta)

#' @rdname sigma2-method
#' @aliases sigma2,MarginalModel-method
setMethod("sigma2", "MarginalModel", function(object) object@sigma2)

setMethod("relabel", "MarginalModel", function(object, zindex){
  object <- newMarginalModel(object)
  if(identical(zindex, seq_len(k(object)))) return(object)
  ##
  ## Permute the labels for the components
  ##
  zz <- factor(z(object), levels=zindex)
  zz <- as.integer(zz)
  z(object) <- zz
  zFreq(object) <- as.integer(table(zz))
  dataMean(object) <- dataMean(object)[zindex]
  dataPrec(object) <- dataPrec(object)[zindex]
  object
})


setMethod("relabel", "BatchModel", function(object, zindex){
  object <- newBatchModel(object)
  if(identical(zindex, seq_len(k(object)))) return(object)
  ##
  ## Permute only the latent variables
  ##
  zz <- factor(z(object), levels=zindex)
  zz <- as.integer(zz)
  z(object) <- zz
  zFreq(object) <- as.integer(table(zz))
  dataMean(object) <- dataMean(object)[, zindex, drop=FALSE]
  dataPrec(object) <- dataPrec(object)[, zindex, drop=FALSE]
  object
})

setMethod("computeModes", "MarginalModel", function(object){
  .computeModesMarginal(object)
})


setMethod("showMeans", "MarginalModel", function(object){
  paste(round(theta(object), 2), collapse=", ")
})

setMethod("showSigmas", "MarginalModel", function(object){
  paste(round(sqrt(sigma2(object)), 2), collapse=", ")
})

setMethod("tablez", "MarginalModel", function(object) table(z(object)))

setMethod("updateMultinomialProb", "MarginalModel", function(object){
  update_multinomialPr(object)
})

setMethod("updateMultinomialProb", "BatchModel", function(object){
  update_multinomialPr_batch(object)
})

setMethod("computeLoglik", "BatchModel", function(object){
  compute_loglik_batch(object)
})

setMethod("computeLoglik", "MarginalModel", function(object){
  loglik(object)
})
