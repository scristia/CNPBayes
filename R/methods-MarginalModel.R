#' Create an object for running marginal MCMC simulations.
#' @examples
#'      model <- MarginalModel(data=rnorm(10), k=1)
#' @param data the data for the simulation.
#' @param k An integer value specifying the number of latent classes.
#' @param hypp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mcmc.params An object of class 'McmcParams'
#' @return An object of class 'MarginalModel'
#' @export
MarginalModel <- function(data=numeric(), k=2, hypp, mcmc.params){
  batch <- rep(1L, length(data))
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
                .internal.constraint=5e-4)
  object <- startingValues(object)
}

SingleBatchPooledVar <- function(data=numeric(), k=2, hypp, mcmc.params){
  obj <- MarginalModel(data, k, hypp, mcmc.params)
  obj@sigma2 <- mean(sigma2(obj))
  ch <- chains(obj)
  ch@sigma2 <- matrix(NA, iter(obj), 1)
  obj@mcmc.chains <- ch
  as(obj, "SingleBatchPooledVar")
}

setValidity("SingleBatchPooledVar", function(object){
  s2 <- sigma2(object)
  if(length(s2) != 1){
    return("sigma2 slot should be length-one numeric vector")
  }
  TRUE
})

getK <- function(object){
  hypp <- hyperParams(object)
  getK(hypp)
}

#' @rdname mu-method
#' @aliases mu,MarginalModel-method
#' @export
setMethod("mu", "MarginalModel", function(object) object@mu)

#' @rdname tau2-method
#' @aliases tau2,MarginalModel-method
setMethod("tau2", "MarginalModel", function(object) object@tau2)

setMethod("show", "MarginalModel", function(object) callNextMethod())

setMethod("computeMeans", "MarginalModel", function(object){
  compute_means(object)
})

setMethod("computeVars", "MarginalModel", function(object){
  compute_vars(object)
})

setReplaceMethod("tau2", "MarginalModel", function(object, value){
  object@tau2 <- value
  object
})

setReplaceMethod("mu", "MarginalModel", function(object, value){
  object@mu <- value
  object
})

#' @rdname bic-method
#' @aliases bic,MarginalModel-method
setMethod("bic", "MarginalModel", function(object){
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

newMarginalModel <- function(object){
  mp <- mcmcParams(object)
  object2 <- MarginalModel(y(object), k=k(object), mcmc.params=mp,
                           hypp=hyperParams(object))
  theta(object2) <- theta(object)
  sigma2(object2) <- sigma2(object)
  p(object2) <- p(object)
  z(object2) <- z(object)
  nu.0(object2) <- nu.0(object)
  mu(object2) <- mu(object)
  tau2(object2) <- tau2(object)
  zFreq(object2) <- zFreq(object)
  probz(object2) <- probz(object)
  sigma2.0(object2) <- sigma2.0(object)
  dataMean(object2) <- dataMean(object)
  dataPrec(object2) <- dataPrec(object)
  log_lik(object2) <- log_lik(object)
  logPrior(object2) <- logPrior(object)
  modes(object2) <- modes(object)
  object2
}

newBatchModel <- function(object){
  mp <- mcmcParams(object)
  object2 <- BatchModel(y(object), batch=batch(object),
                        k=k(object), mcmc.params=mp,
                        hypp=hyperParams(object))
  theta(object2) <- theta(object)
  sigma2(object2) <- sigma2(object)
  p(object2) <- p(object)
  z(object2) <- z(object)
  nu.0(object2) <- nu.0(object)
  mu(object2) <- mu(object)
  tau2(object2) <- tau2(object)
  zFreq(object2) <- zFreq(object)
  probz(object2) <- probz(object)
  sigma2.0(object2) <- sigma2.0(object)
  dataMean(object2) <- dataMean(object)
  dataPrec(object2) <- dataPrec(object)
  log_lik(object2) <- log_lik(object)
  logPrior(object2) <- logPrior(object)
  modes(object2) <- modes(object)
  object2
}

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

.computeModesMarginal <- function(object){
  i <- argMax(object)
  mc <- chains(object)
  thetamax <- theta(mc)[i, ]
  sigma2max <- sigma2(mc)[i,]
  pmax <- p(mc)[i, ]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mu(mc)[i],
                tau2=tau2(mc)[i],
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i],
                zfreq=zFreq(mc)[i, ],
                loglik=log_lik(mc)[i],
                logprior=logPrior(mc)[i])
  modes
}

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

modelOtherModes <- function(model, maxperm=5){
  kperm <- permnK(k(model), maxperm)
  model.list <- vector("list", nrow(kperm))
  for(i in seq_along(model.list)){
    model.list[[i]] <- relabel(model, kperm[i, ])
  }
  model.list
}

#' Compute the log bayes factor between models.
#'
#' Models of varying component sizes are compared. The log bayes factor is 
#' calculated comparing the two models with the highest marginal likelihood, 
#' as computed by \code{marginalLikelihood}.
#' @param x the result of a call to \code{computeMarginalLik}.
#' @return Log Bayes factor comparing the two models with highest likelihood.
#' @export
logBayesFactor <- function(x){
    top.two <- sort(x, decreasing=TRUE)[1:2]
    log.bf <- diff(-top.two)
    log.bf
}

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
