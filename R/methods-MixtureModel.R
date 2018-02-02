#' @include AllGenerics.R
NULL

setValidity("MixtureModel", function(object){
  msg <- TRUE
  if(length(p(object)) != k(object)){
    msg <- "Mixture probability vector must be the same length as k"
    return(msg)
  }
  if(k(object)!=k(hyperParams(object))){
    msg <- "disagreement of k in hyperparams and model"
    return(msg)
  }
  ## maximum value of nu0 is currently hard-coded in C as 100
  ##if(nu.0(object) > 100){
    ##return("nu.0 can not exceed 100")
## }
  msg
})

#' @rdname hyperParams-method
#' @aliases hyperParams,MixtureModel-method
setMethod("hyperParams", "MixtureModel", function(object) object@hyperparams)

#' @rdname hyperParams-method
#' @aliases hyperParams<-,MixtureModel-method
setReplaceMethod("hyperParams", c("MixtureModel", "Hyperparameters"),
  function(object, value) {
    object@hyperparams <- value
    object
})

setReplaceMethod("batch", "MixtureModel", function(object, value){
  object@batch <- value
  object
})

observed <- function(object) object@data

#' Retrieve standard deviations of each component/batch mean.
#'
#' @examples
#'      sigma(SingleBatchModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of length K, or a matrix of size B x K, where
#' K is the number of components and B is the number of batches
#' @export
sigma <- function(object) sqrt(sigma2(object))

#' Retrieve overall standard deviation.
#'
#' @examples
#'      tau(SingleBatchModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of standard deviations
#' @export
tau <- function(object) sqrt(tau2(object))

#' @rdname nu.0-method
#' @aliases nu.0,MixtureModel-method
setMethod("nu.0", "MixtureModel", function(object) object@nu.0)

#' @rdname sigma2.0-method
#' @aliases sigma2.0,MixtureModel-method
setMethod("sigma2.0", "MixtureModel", function(object) object@sigma2.0)

sigma.0 <- function(object) sqrt(sigma2.0(object))



#' Retrieve mixture proportions.
#'
#' @examples
#'      p(MarginalModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of length the number of components
#' @export
p <- function(object) object@pi

nComp <- function(object) length(p(object))
dataMean <- function(object) object@data.mean
dataPrec <- function(object) object@data.prec
dataSd <- function(object) sqrt(1/dataPrec(object))

#' @rdname k-method
#' @aliases k,MixtureModel-method
setMethod("k", "MixtureModel", function(object) object@k)

#' @rdname k-method
#' @aliases k<-,MixtureModel-method
setReplaceMethod("k", "MixtureModel",
    function(object, value) {
        k <- as.integer(value)
        hypp <- hyperParams(object)
        hypp@k <- k
        hypp@alpha <- rep(1, k)
        hyperParams(object) <- hypp
        object@k <- k
        object@pi <- rep(1/k, k)
        object@probz <- matrix(0, length(y(object)), k)
        object <- startingValues(object)
        object
    }
)

setReplaceMethod("z", "MixtureModel", function(object, value){
  ##object@z <- factor(value, levels=seq_len(k(object)))
  object@z <- value
  object
})

setReplaceMethod("theta", "MixtureModel", function(object, value){
  object@theta <- value
  object
})

setReplaceMethod("sigma2", "MixtureModel", function(object, value){
  object@sigma2 <- value
  object
})

setReplaceMethod("p", "MixtureModel", function(object, value){
  object@pi <- value
  object
})



setReplaceMethod("nu.0", "MixtureModel", function(object, value){
  object@nu.0 <- value
  object
})

setReplaceMethod("sigma2.0", "MixtureModel", function(object, value){
  object@sigma2.0 <- value
  object
})

#' @rdname chains-method
#' @aliases chains,MixtureModel-method
setMethod("chains", "MixtureModel", function(object) object@mcmc.chains)

setReplaceMethod("chains", "MixtureModel", function(object, value){
  object@mcmc.chains <- value
  object
})

setMethod("show", "MixtureModel", function(object){
  cat("An object of class '", class(object), "'\n")
  cat("  data: \n")
  cat("     n           :", length(y(object)), "\n")
  cat("     k           :", nComp(object), "\n")
  cat("     table(z)    :", paste(tablez(object), collapse=", "), "\n")
  cat("     mix prob (s):", paste(round(p(object), 2), collapse=", "), "\n")
  sds <- showSigmas(object)
  mns <- showMeans(object)
  cat("     sigma (s)   :", sds, "\n")
  cat("     theta (s)   :", mns, "\n")
  cat("     sigma2.0 (s):", round(sigma2.0(object), 2), "\n")
  cat("     nu.0 (s)    :", nu.0(object), "\n")
  cat("     logprior(s):", round(logPrior(object), 2), "\n")
  cat("     loglik (s)  :", round(log_lik(object), 2), "\n")
})

setMethod("alpha", "MixtureModel", function(object) alpha(hyperParams(object)))

#' @rdname y-method
#' @aliases y,MixtureModel-method
setMethod("y", "MixtureModel", function(object) object@data)

setReplaceMethod("y", "MixtureModel", function(object, value){
  object@data <- value
  object
})

#' @rdname oned-method
#' @aliases oned,MixtureModel-method
setMethod("oned", "MixtureModel", function(object) object@data)

#' @rdname batch-method
#' @aliases batch,MixtureModel-method
setMethod("batch", "MixtureModel", function(object) object@batch)

#' @rdname z-method
#' @aliases z,MixtureModel-method
setMethod("z", "MixtureModel", function(object) object@z)

setMethod("computePrec", "MarginalModel", function(object){
  compute_prec(object)
})

setMethod("computePrior", "SingleBatchModel", function(object){
  compute_logprior(object)
})


.computeProbZ <- function(object){
  pZ <- probz(object)
  zz <- as.integer(z(object))
  for(j in seq_len(k(object))){
    pZ[, j] <- pZ[, j] + as.integer(zz==j)
  }
  pZ
}

#' @aliases numberObs,MixtureModel-method
#' @rdname numberObs-method
setMethod("numberObs", "MixtureModel", function(model) length(y(model)))


setReplaceMethod("probz", "MixtureModel", function(object, value){
  object@probz <- value
  object
})

## TODO Dangerous to have accessor do something more than return the value of it
## slot.  Further
## probz(object) <- probz(object)
## will not behave as expected
#' @rdname probz-method
#' @aliases probz,MixtureModel-method
setMethod("probz", "MixtureModel", function(object) {
  ## because first iteration not saved
  object@probz/(iter(object)-1)
})

setMethod("runBurnin", "SingleBatchModel", function(object){
  mcmc_marginal_burnin(object, mcmcParams(object))
})

setMethod("runBurnin", "SingleBatchPooled", function(object){
  burnin_singlebatch_pooled(object, mcmcParams(object))
})

setMethod("runBurnin", "MultiBatchModel", function(object){
  mcmc_batch_burnin(object, mcmcParams(object))
})

setMethod("runBurnin", "MultiBatchPooled", function(object){
  burnin_multibatch_pvar(object, mcmcParams(object))
})

setMethod("runMcmc", "SingleBatchModel", function(object){
  mcmc_marginal(object, mcmcParams(object))
})

setMethod("runMcmc", "SingleBatchPooled", function(object){
  mcmc_singlebatch_pooled(object, mcmcParams(object))
})

setMethod("runMcmc", "MultiBatchModel", function(object){
  mcmc_batch(object, mcmcParams(object))
})

setMethod("runMcmc", "MultiBatchPooled", function(object){
  mcmc_multibatch_pvar(object, mcmcParams(object))
})

multipleStarts <- function(object){
  if(k(object)==1) return(object)
  mcmcp <- mcmcParams(object)
  mmod <- replicate(nStarts(mcmcp), SingleBatchModel(y(object), mcmc.params=mcmcp,
                                                  hypp=hyperParams(object), k=k(object)))
  models <- suppressMessages(lapply(mmod, runBurnin))
  lp <- sapply(models, log_lik)
  select <- which.max(lp)
  if(length(select) == 0) stop("No model selected")
  model <- models[[select]]
  if(isSB(object)) return(model)
  ##
  ##  initialize batch model
  ##
  bmodel <- MultiBatchModel(data=y(model), batch=batch(object), k=k(object), hypp=hyperParams(object))
  mcmcParams(bmodel, force=TRUE) <- mcmcParams(object)
  theta(bmodel) <- matrix(theta(model), nBatch(object), k(object), byrow=TRUE)
  mu(bmodel) <- theta(model)
  z(bmodel) <- z(model)
  bmodel <- ensureAllComponentsObserved(bmodel)
  zFreq(bmodel) <- as.integer(table(z(bmodel)))
  dataMean(bmodel) <- computeMeans(bmodel)
  dataPrec(bmodel) <- computePrec(bmodel)
  bmodel
}

selectByLogLik <- function(model.list){
  lp <- sapply(model.list, log_lik)
  isfin <- is.finite(lp)
  isnonzero <- lp != 0
  is.candidate <- isfin & isnonzero
  if(!any(is.candidate)){
    stop("Bad starting values.")
  }
  lp <- lp[is.candidate]
  model.list <- model.list[is.candidate]
  select <- which.max(lp)
  if(length(select) == 0 ) stop("No model selected")
  model <- model.list[[select]]
  model
}

multipleStarts2 <- function(object){
  if(k(object)==1) return(object)
  mcmcp <- mcmcParams(object)
  ##
  ##
  if(is(object, "MultiBatchModel")){
    model.list <- replicate(nStarts(mcmcp),
                            MultiBatchModel(y(object), mcmc.params=mcmcp,
                                       hypp=hyperParams(object), k=k(object),
                                       batch=batch(object)))
  }
  if(is(object, "SingleBatchModel")){
    model.list <- replicate(nStarts(mcmcp),
                            SingleBatchModel(y(object), mcmc.params=mcmcp,
                                          hypp=hyperParams(object), k=k(object)))
  }
  model <- selectByLogLik(model.list)
  model
}


psParams <- function(warnings=TRUE,
                     returnNULLonWarnings=FALSE){
  list(warnings=warnings)
}




.ordered_thetas_multibatch<- function(model){
  thetas <- theta(model)
  checkOrder <- function(theta) identical(order(theta), seq_along(theta))
  is_ordered <- apply(thetas, 1, checkOrder)
  all(is_ordered)
}

reorderMultiBatch <- function(model){
  is_ordered <- .ordered_thetas_multibatch(model)
  if(is_ordered) return(model)
  ## thetas are not all ordered
  thetas <- theta(model)
  s2s <- sigma2(model)
  K <- k(model)
  ix <- order(thetas[1, ])
  B <- nBatch(model)
  zlist <- split(z(model), batch(model))
  for(i in seq_len(B)){
    ix.next <- order(thetas[i, ])
    thetas[i, ] <- thetas[i, ix.next]
    s2s[i, ] <- s2s[i, ix]
    zlist[[i]] <- as.integer(factor(zlist[[i]], levels=ix.next))
  }
  zs <- unlist(zlist)
  ps <- p(model)[ix]
  mu(model) <- mu(model)[ix]
  tau2(model) <- tau2(model)[ix]
  sigma2(model) <- s2s
  theta(model) <- thetas
  p(model) <- ps
  z(model) <- zs
  dataMean(model) <- computeMeans(model)
  dataPrec(model) <- computePrec(model)
  log_lik(model) <- computeLoglik(model)
  model
}



reorderSingleBatch <- function(model){
  thetas <- theta(model)
  K <- k(model)
  ix <- order(thetas)
  if(identical(ix, seq_len(K))) return(model)
  thetas <- thetas[ix]
  s2s <- sigma2(model)[ix]
  zs <- as.integer(factor(z(model), levels=ix))
  ps <- p(model)[ix]
  sigma2(model) <- s2s
  theta(model) <- thetas
  p(model) <- ps
  z(model) <- zs
  dataPrec(model) <- 1/computeVars(model)
  dataMean(model) <- computeMeans(model)
  model
}


setMethod("sortComponentLabels", "SingleBatchModel", function(model){
  reorderSingleBatch(model)
})

setMethod("sortComponentLabels", "MultiBatchModel", function(model){
  reorderMultiBatch(model)
})

setMethod("isOrdered", "MixtureModel", function(object){
  identical(order(theta(object)), seq_along(theta(object)))
})

setMethod("isOrdered", "MultiBatchModel", function(object){
  .ordered_thetas_multibatch(object)
})

setReplaceMethod("dataMean", "MixtureModel", function(object, value){
  object@data.mean <- value
  object
})

setReplaceMethod("dataPrec", "MixtureModel", function(object, value){
  object@data.prec <- value
  object
})

setMethod("mu.0", "MixtureModel", function(object) mu.0(hyperParams(object)))
setMethod("mu.0", "Hyperparameters", function(object) object@mu.0)

setMethod("tau2.0", "MixtureModel", function(object) tau2.0(hyperParams(object)))
setMethod("tau2.0", "Hyperparameters", function(object) object@tau2.0)

tau.0 <- function(object) sqrt(tau2.0(object))

#' @rdname eta.0-method
#' @aliases eta.0,MixtureModel-method
setMethod("eta.0", "MixtureModel", function(object) eta.0(hyperParams(object)))

#' @rdname eta.0-method
#' @aliases eta.0,Hyperparameters-method
setMethod("eta.0", "Hyperparameters", function(object) object@eta.0)

#' @rdname m2.0-method
#' @aliases m2.0,MixtureModel-method
setMethod("m2.0", "MixtureModel", function(object) m2.0(hyperParams(object)))

#' @rdname m2.0-method
#' @aliases m2.0,Hyperparameters-method
setMethod("m2.0", "Hyperparameters", function(object) object@m2.0)

setReplaceMethod("eta.0", "MixtureModel", function(object, value){
  eta.0(hyperParams(object)) <- value
  object
})

setReplaceMethod("m2.0", "MixtureModel", function(object, value){
  m2.0(hyperParams(object)) <- value
  object
})

setReplaceMethod("eta.0", "Hyperparameters", function(object, value){
  object@eta.0 <- value
  object
})

setReplaceMethod("m2.0", "Hyperparameters", function(object, value){
  object@m2.0 <- value
  object
})




##
## the batch names tend to be much too long
##
makeUnique <- function(x){
  ub <- unique(x)
  ##names(ub) <- ub
  maxchar <- pmin(nchar(ub), 8)
  abbrv <- setNames(make.unique(substr(ub, 1, maxchar)), ub)
  as.character(abbrv[x])
}

#' Calculate the maximum a posteriori estimate of latent variable assignment.
#'
#' @examples
#'      map_z(SingleBatchModelExample)
#' @param object an object of class MixtureModel.
#' @return map estimate of latent variable assignment for each observation
#' @export
map_z <- function(object) {
  max.col(probz(object))
}

setMethod("thetac", "MixtureModel", function(object) theta(chains(object)))

setMethod("thetaMean", "MixtureModel", function(object) colMeans(thetac(object)))

setMethod("sigmaMean", "MixtureModel", function(object) colMeans(sigmac(object)))

logLikc <- function(object) log_lik(chains(object))


#' Retrieve standard deviation of each component/batch mean at each iteration of the MCMC.
#'
#' @examples
#'      sigmac(SingleBatchModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A matrix of size N x K where N is the number of observations
#' and K is the number of components
#' @export
sigmac <- function(object) sigma(chains(object))

#' Retrieve mixture proportions at each iteration of the MCMC.
#'
#' @examples
#'      pic(SingleBatchModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A matrix of size MCMC iterations x Number of components
#' @export
pic <- function(object) p(chains(object))

setMethod("pMean", "MixtureModel", function(object){
  colMeans(pic(object))
})

#' Retrieve overall mean at each iteration of the MCMC.
#'
#' @examples
#'      muc(SingleBatchModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of length N or matrix of size N x B, where N is the 
#' number of observations and B is the number of unique batches.
#' @export
muc <- function(object) mu(chains(object))

#' Retrieve overall mean averaged across MCMC simulations.
#'
#' @examples
#'      muMean(SingleBatchModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of size 1 or number of batches
#' @export
muMean <- function(object) {
  x <- muc(object)
  if(is(object, "MultiBatchModel")){
    return(colMeans(x))
  }
  mean(x)
}

#' Retrieve overall standard deviation at each iteration of the MCMC.
#'
#' @examples
#'      tauc(SingleBatchModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of length N or matrix of size N x B, where N is the 
#' number of observations and B is the number of unique batches.
#' @export
tauc <- function(object) sqrt(tau2(chains(object)))

#' Retrieve overall standard deviation averaged across MCMC simulations.
#'
#' @examples
#'      tauMean(SingleBatchModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of size 1 or number of batches
#' @export
tauMean <- function(object){
  x <- tauc(object)
  if(is(object, "MultiBatchModel")){
    return(colMeans(x))
  }
  mean(x)
}

#' @rdname modes-method
#' @aliases modes,MixtureModel-method
setMethod("modes", "MixtureModel", function(object) object@modes)

#' @rdname modes-method
#' @aliases modes<-,MixtureModel-method
setReplaceMethod("modes", "MixtureModel", function(object, value) {
  object@modes <- value
  object
})

#' @rdname log_lik-method
#' @aliases log_lik,MixtureModel-method
setMethod("log_lik", "MixtureModel", function(object){
  object@loglik
})

setReplaceMethod("log_lik", "MixtureModel", function(object, value){
  object@loglik <- value
  object
})

argMax <- function(object){
  ll <- log_lik(chains(object))
  if(length(ll) == 1) return(1)
  lp <- logPrior(chains(object))
  p <- ll+lp
  p <- p[is.finite(p)]
  if(length(p) == 0){
    return(1)
  }
  maxp <- max(p)
  which(p == maxp)
}

setMethod("isSB", "MarginalModel", function(object) TRUE)
setMethod("isSB", "SingleBatchModel", function(object) TRUE)
setMethod("isSB", "BatchModel", function(object) FALSE)
setMethod("isSB", "MultiBatchModel", function(object) FALSE)

startAtTrueValues <- function(model, truth){
  theta(model) <- theta(truth)
  sigma2(model) <- sigma2(truth)
  p(model) <- p(truth)
  z(model) <- z(truth)
  mu(model) <- mu(truth)
  tau2(model) <- tau2(truth)
  dataMean(model) <- computeMeans(model)
  dataPrec(model) <- 1/computeVars(model)
  zFreq(model) <- as.integer(table(z(model)))
  log_lik(model) <- computeLoglik(model)
  model
}

restartAtChainIndex <- function(model, index){
  ch <- chains(model)
  if(!isSB(model) ){
    B <- nBatch(model)
    K <- k(model)
    theta(model) <- matrix(theta(ch)[index, ], B, K)
    sigma2(model) <- matrix(sigma2(ch)[index, ], B, K)
    p(model) <- p(ch)[index, ]
    z(model) <- z(ch)[index, ]
    mu(model) <- mu(ch)[index, ]
    tau2(model) <- tau2(ch)[index, ]
    sigma2.0(model) <- sigma2.0(ch)[index]
    nu.0(model) <- nu.0(ch)[index]
    zFreq(model) <- as.integer(table(z(model)))
    dataMean(model) <- computeMeans(model)
    dataPrec(model) <- 1/computeVars(model)
    return(model)
  }
  theta(model) <- theta(ch)[index, ]
  sigma2(model) <- sigma2(ch)[index, ]
  p(model) <- p(ch)[index, ]
  z(model) <- z(ch)[index, ]
  mu(model) <- mu(ch)[index]
  tau2(model) <- tau2(ch)[index]
  sigma2.0(model) <- sigma2.0(ch)[index]
  nu.0(model) <- nu.0(ch)[index]
  zFreq(model) <- as.integer(table(z(model)))
  dataMean(model) <- computeMeans(model)
  dataPrec(model) <- 1/computeVars(model)
  model
}

#' @rdname zfreq-method
#' @aliases zfreq,MixtureModel-method
setMethod("zFreq", "MixtureModel", function(object) object@zfreq)

setReplaceMethod("zFreq", "MixtureModel", function(object, value){
  object@zfreq <- value
  object
})

#' @rdname mcmcParams-method
#' @aliases mcmcParams,MixtureModel-method
setMethod("mcmcParams", "MixtureModel", function(object) object@mcmc.params )

#' @rdname iter-method
#' @aliases iter,MixtureModel-method
setMethod("iter", "MixtureModel", function(object) iter(mcmcParams(object)))

#' @rdname nStarts-method
#' @aliases nStarts,MixtureModel-method
setMethod("nStarts", "MixtureModel", function(object) nStarts(mcmcParams(object)))

#' @rdname nStarts-method
#' @aliases nStarts<-,MixtureModel-method
setReplaceMethod("nStarts", "MixtureModel", function(object, value){
  mcmcParams(object)@nstarts <- as.integer(value)
  object
})

#' @rdname thin-method
#' @aliases thin,MixtureModel-method
setMethod("thin", "MixtureModel", function(object) thin(mcmcParams(object)))

#' @rdname burnin-method
#' @aliases burnin,MixtureModel-method
setMethod("burnin", "MixtureModel", function(object) burnin(mcmcParams(object)))

#' @rdname burnin-method
#' @aliases burnin<-,MixtureModel-method
setReplaceMethod("burnin", "MixtureModel", function(object, value){
  burnin(mcmcParams(object)) <- value
  object
})

#' @rdname iter-method
#' @aliases iter<-,MixtureModel-method
setReplaceMethod("iter", "MixtureModel", function(object, force=FALSE, value){
  mp <- mcmcParams(object)
  iter(mp) <- value
  mcmcParams(object, force=force) <- mp
  object
})

#' @rdname logPrior-method
#' @aliases logPrior,MixtureModel-method
setMethod("logPrior", "MixtureModel", function(object) object@logprior)

setReplaceMethod("logPrior", "MixtureModel", function(object, value) {
  object@logprior <- value
  object
})

setMethod("paramUpdates", "MixtureModel", function(x) paramUpdates(mcmcParams(x)))
setReplaceMethod("paramUpdates", "MixtureModel", function(x, value){
  paramUpdates(mcmcParams(x)) <- value
  x
})

nu0c <- function(object) nu.0(chains(object))
sigma20c <- function(object) sigma2.0(chains(object))

#' @rdname mcmcParams-method
#' @aliases mcmcParams,MixtureModel-method
setReplaceMethod("mcmcParams", "MixtureModel", function(object, force=TRUE, value){
  it <- iter(object)
  if(it != iter(value)){
    if(!force){
      msg <- "Replacement will change the size of the elements in mcmc.chains slot."
      msg2 <- "Force=TRUE will allow the replacement"
      stop(paste(msg, msg2, sep="\n"))
    } else {
      ## force is TRUE
      if(iter(value) > iter(object)){
        object@mcmc.params <- value
        ## create a new chain
        mcmc_chains <- McmcChains(object)
      } else {
        object@mcmc.params <- value
        index <- seq_len(iter(value))
        mcmc_chains <- chains(object)[index, ]
      }
      chains(object) <- mcmc_chains
      return(object)
    }
  }
  ## if we've got to this point, it must be safe to update mcmc.params
  ## (i.e., size of chains is not effected)
  object@mcmc.params <- value
  object
})

#' @rdname mcmcParams-method
#' @aliases mcmcParams,list-method
setReplaceMethod("mcmcParams", "list",
                 function(object, force=TRUE, value){
                   for(i in seq_along(object)){
                     mcmcParams(object[[i]], force=force) <- value
                   }
                   object
                 })

setMethod("zChain", "MixtureModel", function(object) chains(object)@z)

mapModel <- function(model){
  model2 <- restartAtChainIndex(model, argMax(model))
  model2
}

#' Probabiliistic copy number assigments.
#'
#' Calculate probabilistic copy number assignments using Bayes Rule applied at the MAP estimates of the cluster mean, variance, and class proportion parameters
#' @param model An object of class MixtureModel.
#' @return A matrix of size N x K where N is number of observations and K is the number of components.
#' @export
mapCnProbability <- function(model){
  ## Cardin et al. : calculate probabilistic copy number assignments
  ## using Bayes Rule applied at the MAP estimates of the cluster
  ## mean, variance, and class proportion parameters
  map_model <- mapModel(model)
  p <- updateMultinomialProb(map_model)
  if(isSB(model)){
    p <- p[, order(theta(map_model))]
  } else {
    p <- p[, order(mu(map_model))]
  }
  return(p)
}


#' @param value a length-one numeric vector indicating how often to save MCMC iterations to the chain.  For example, a thin of 10 means that every 10th MCMC simulation is saved to the chain.
#' @rdname thin-method
#' @aliases thin<-,MixtureModel,numeric-method
setReplaceMethod("thin", c("MixtureModel", "numeric"), function(object, value){
  mp <- mcmcParams(object)
  mp@thin <- as.integer(value)
  mcmcParams(object) <- mp
  object
})

#' @rdname mcmcParams-method
#' @aliases mcmcParams,list-method
setMethod("mcmcParams", "list", function(object){
  mcmcParams(object[[1]])
})

#' @aliases label_switch,MixtureModel-method
#' @rdname label_switch
setMethod("label_switch", "MixtureModel", function(object) object@label_switch)

setReplaceMethod("label_switch", c("MixtureModel", "logical"),
                 function(object, value){
                   object@label_switch <- value
                   object
                 })

#' @aliases marginal_lik,MixtureModel-method
#' @rdname marginal_lik
setMethod("marginal_lik", "MixtureModel", function(object){
  object@marginal_lik
})

#' @aliases marginal_lik<-,MixtureModel,numeric-method
#' @rdname marginal_lik
setReplaceMethod("marginal_lik", c("MixtureModel", "numeric"),
                 function(object, value){
                   object@marginal_lik <- value
                   object
                 })


#' @rdname tile-functions
#' @aliases upSample,MultiBatchModel-method
setMethod("upSample", "MultiBatchModel", function(model, tiles){
  tile.sum <- tileSummaries(tiles)
  ##stopifnot(all.equal(y(model), tile.sum$avgLRR))
  key <- match(tiles$tile, tile.sum$tile)
  model2 <- model
  y(model2) <- tiles$logratio
  probs <- probz(model)
  probs2 <- probs[key, , drop=FALSE]
  probz(model2) <- probs2 * (iter(model) - 1)
  z(model2) <- z(model)[key]
  batch(model2) <- tiles$batch
  model2
})

#' @rdname tile-functions
#' @aliases upSample,MixtureModel-method
setMethod("upSample", "MixtureModel", function(model, tiles){
  tile.sum <- tileSummaries(tiles)
  ##stopifnot(all.equal(y(model), tile.sum$avgLRR, tolerance=0.1, scale=1))
  key <- match(tiles$tile, tile.sum$tile)
  model2 <- model
  y(model2) <- tiles$logratio
  probs <- probz(model)
  probs2 <- probs[key, , drop=FALSE]
  probz(model2) <- probs2 * (iter(model) - 1)
  z(model2) <- z(model)[key]
  model2
})
