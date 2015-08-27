#' Create an object for running marginal MCMC simulations.
#' @param data the data for the simulation.
#' @param k An integer value specifying the number of latent classes.
#' @param hypp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mcmc.params An object of class 'McmcParams'
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
                hwe=numeric(),
                modes=list(),
                mcmc.params=mcmc.params,
                .internal.constraint=5e-4)
  object <- startingValues(object)
}

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



setMethod("initializeSigma2.0", "MarginalModel", function(object){
  hypp <- hyperParams(object)
  sum(alpha(hypp)*sigma2(object))/sum(alpha(hypp))
})


## just choose a big number
setMethod("initializeTau2", "MarginalModel", function(object)  1000)

setMethod("posteriorMultinomial", "MarginalModel", function(object){
  update_multinomialPr(object)
})

setMethod("show", "MarginalModel", function(object) callNextMethod())

setMethod("computeMeans", "MarginalModel", function(object){
  compute_means(object)
})

setMethod("computeVars", "MarginalModel", function(object){
  compute_vars(object)
})

setMethod("simulateY", "MarginalModel", function(object){
  zz <- z(object)
  yy <- rnorm(length(zz), mean=theta(object)[zz], sd=sigma(object)[zz])
})

setMethod("updateThetaCpp", "MarginalModel", function(object, constrain) {
  update_theta(object, constrain=constrain)
})

setMethod("updateTheta", "MarginalModel", function(object) {
  update_theta(object)
})

setMethod("updateSigma2", "MarginalModel", function(object) {
  update_sigma2(object)
})

setMethod("updateSigma2.0", "MarginalModel", function(object){
  update_sigma2_0(object)
})

setReplaceMethod("tau2", "MarginalModel", function(object, value){
  object@tau2 <- value
  object
})

setReplaceMethod("mu", "MarginalModel", function(object, value){
  object@mu <- value
  object
})

##
## For the marginal model, mu has already been initialized in the hyperparameters
##
setMethod("initializeMu", "MarginalModel", function(object)   mu(object))

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

setMethod("reorderComponents", "MarginalModel", function(object){
  ##
  ## First, update the model so that the components are ordered by theta
  ##
  ix <- order(theta(object))
  theta(object) <- sort(theta(object))
  sigma2(object) <- sigma2(object)[ix]
  p(object) <- p(object)[ix]
  zz <- z(object)
  ##
  ## Set the labels of the latent variable such that 1=first
  ## components, 2= second component, ...
  ##
  zz <- factor(as.integer(factor(zz, levels=ix)), levels=seq_len(k(object)))
  z(object) <- zz
  dataMean(object) <- dataMean(object)[ix]
  dataPrec(object) <- dataPrec(object)[ix]
  object
})

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

setMethod("sort", "MarginalModel", function(x, decreasing=FALSE, ...){
  mc <- chains(x)
  pot <- logpotential(mc)
  index <- which.max(pot)
  thetas <- theta(mc)[index, ]
  if(identical(thetas, sort(thetas))){
    ## nothing to do
    return(x)
  }
  cn <- order(thetas)
  theta(mc) <- theta(mc)[, cn]
  theta(x) <- theta(x)[cn]

  sigma2(mc) <- sigma2(mc)[, cn]
  sigma2(x) <- sigma2(x)[cn]

  p(mc) <- p(mc)[, cn]
  p(x) <- p(x)[cn]

  mu(x) <- mu(x)[cn]
  tau2(x) <- tau2(x)[cn]

  probz(x) <- probz(x)[, cn]

  zz <- as.integer(z(x))
  z(x) <- factor(as.integer(factor(zz, levels=cn)), levels=sort(unique(zz)))
  dataMean(x) <- dataMean(x)[cn]
  dataPrec(x) <- dataPrec(x)[cn]
  chains(x) <- mc
  x
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

sumSquares <- function(object){
  B <- batch(object)
  thetas <- theta(object)
  yy <- y(object)
  K <- k(object)
  zz <- z(object)
  ss <- matrix(NA, nBatch(object), k(object))
  rownames(ss) <- uniqueBatch(object)
  batch.index <- split(seq_len(length(yy)), B)
  zz <- z(object)
  for(b in uniqueBatch(object)){
    k <- batch.index[[b]]
    y <- yy[k]
    cn <- zz[k]
    m <- thetas[b, ]
    ## This could be tricky in C.  It works in R because of the factor to an integer:
    ##  as.integer(factor(c(1, 3), levels=c("1", "2", "3"))) evaluates to 1,3
    m <- m[as.integer(cn)]
    squares <- (y - m)^2
    ss[b, ] <- sapply(split(squares, factor(cn, levels=seq_len(K))), sum)
  }
  ss
}

setMethod("showMeans", "MarginalModel", function(object){
  paste(round(theta(object), 2), collapse=", ")
})

setMethod("showSigmas", "MarginalModel", function(object){
  paste(round(sqrt(sigma2(object)), 2), collapse=", ")
})

setMethod("tablez", "MarginalModel", function(object) table(z(object)))

permuteZ <- function(object, ix){
  zz <- zChain(object)
  zz2 <- zz
  for(i in seq_along(ix)){
    zz2[zz == ix[i]] <- i
  }
  zz2
}


setMethod("permuteModes", "MarginalModel", function(object, ix){
  modes(object)[["theta"]] <- modes(object)[["theta"]][ix]
  modes(object)[["sigma2"]] <- modes(object)[["sigma2"]][ix]
  modes(object)[["mixprob"]] <- modes(object)[["mixprob"]][ix]
  object
})

permuteZandModes <- function(kmod, ix){
  zChain(kmod) <- permuteZ(kmod, ix)
  kmod <- permuteModes(kmod, ix)
  kmod
}

## same for marginal and batch models

setMethod("pMixProb", "MixtureModel", function(object) {
  exp(p_pmix_reduced(object))
})

modelOtherModes <- function(model, maxperm=5){
  kperm <- permnK(k(model), maxperm)
  model.list <- vector("list", nrow(kperm))
  for(i in seq_along(model.list)){
    model.list[[i]] <- relabel(model, kperm[i, ])
  }
  model.list
}

#' Reorder models of varying component sizes.
#'
#' Models are ordered according to marginal likelihood. The marginal
#' likelihood is computed for each chain of each component size model
#' separately. The mean is taken by model, and ordering by this mean
#' marginal is performed.
#' @param x the result of a call to \code{computeMarginalLik}.
#' @export
orderModels <- function(x){
  models <- x$models
  ##K <- k(models)
  K <- names(models)
  marginal.est.list <- x$marginal
  m <- sapply(marginal.est.list, function(x) x["marginal"])
  K <- K[order(m, decreasing=TRUE)]
  ix <- match(K, names(models))
  models <- models[ix]
  return(models)
}

#' Compute the log bayes factor between models.
#'
#' Models of varying component sizes are compared. The log bayes factor is calculated comparing the two models with the highest marginal likelihood, as computed by \code{computeMarginalLik}.
#' @param x the result of a call to \code{computeMarginalLik}.
#' @export
logBayesFactor <- function(x){
  models <- orderModels(x)
  if(length(models) <= 1) {
    return(NA)
  }
  K <- k(models)
  ##nms <- names(orderModels(x))
  nms <- names(models)
  ml <- sapply(x$marginal[nms], function(x) x["marginal"])
  names(ml) <- substr(names(ml), 1, 2)
  bf <- setNames(ml[1]-ml[2], paste0(names(ml[1:2]), collapse="-"))
  bf
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

setMethod("computeLogLikxPrior", "MixtureModel", function(object){
  compute_llxprior(object)
})
