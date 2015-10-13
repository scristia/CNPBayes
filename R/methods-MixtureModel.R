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
  if(nu.0(object) > 100){
    return("nu.0 can not exceed 100")
  }
  msg
})

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
#'      sigma(MarginalModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of length K, or a matrix of size B x K, where
#' K is the number of components and B is the number of batches
#' @export
sigma <- function(object) sqrt(sigma2(object))

#' Retrieve overall standard deviation.
#'
#' @examples
#'      tau(MarginalModelExample)
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

setMethod("computePrior", "MarginalModel", function(object){
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

setReplaceMethod("probz", "MixtureModel", function(object, value){
  object@probz <- value
  object
})

#' @rdname probz-method
#' @aliases probz,MixtureModel-method
setMethod("probz", "MixtureModel", function(object) {
  object@probz/iter(object)
})


setMethod("runBurnin", "MarginalModel", function(object){
  mcmc_marginal_burnin(object, mcmcParams(object))

})

setMethod("runBurnin", "SingleBatchPooledVar", function(object){
  burnin_singlebatch_pooled(object, mcmcParams(object))
})

setMethod("runBurnin", "BatchModel", function(object){
  mcmc_batch_burnin(object, mcmcParams(object))
})

setMethod("runMcmc", "MarginalModel", function(object){
  mcmc_marginal(object, mcmcParams(object))
})

setMethod("runMcmc", "SingleBatchPooledVar", function(object){
  mcmc_singlebatch_pooled(object, mcmcParams(object))
})

setMethod("runMcmc", "BatchModel", function(object){
  mcmc_batch(object, mcmcParams(object))
})

multipleStarts <- function(object){
  if(k(object)==1) return(object)
  mcmcp <- mcmcParams(object)
  mmod <- replicate(nStarts(mcmcp), MarginalModel(y(object), mcmc.params=mcmcp,
                                                  hypp=hyperParams(object), k=k(object)))
  models <- suppressMessages(lapply(mmod, runBurnin))
  lp <- sapply(models, log_lik)
  model <- models[[which.max(lp)]]
  if(isMarginalModel(object)) return(model)
  ##
  ##  initialize batch model
  ##
  bmodel <- BatchModel(data=y(model), batch=batch(object), k=k(object), hypp=hyperParams(object))
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

#' @rdname posteriorSimulation-method
#' @aliases posteriorSimulation,MixtureModel-method
setMethod("posteriorSimulation", "MixtureModel", function(object){
  .posteriorSimulation(object)
})

#' @rdname posteriorSimulation-method
#' @aliases posteriorSimulation,MixtureModel-method
setMethod("posteriorSimulation", c("MixtureModel", "integer"),
    function(object, k) {
        if (length(k) > 1) {
            mlist <- vector("list", length(k))
            for (i in seq_along(k)) {
                k(object) <- k[i]
                mlist[[i]] <- .posteriorSimulation(object)
            }
            mlist
        } else {
            k(object) <- k
            .posteriorSimulation(object)
        }
    }
)

#' @rdname posteriorSimulation-method
#' @aliases posteriorSimulation,MixtureModel-method
setMethod("posteriorSimulation", c("MixtureModel", "numeric"),
    function(object, k) {
        posteriorSimulation(object, as.integer(k))
    }
)

.posteriorSimulation <- function(post){
  if(nStarts(post) > 1){
    post <- multipleStarts(post)
  }
  post <- runBurnin(post)
  if( iter(post)==0 ) return(post)
  post <- runMcmc(post)
  modes(post) <- computeModes(post)
  post
}

posteriorSimulationPooled <- function(object, iter=1000,
                                      burnin=1000,
                                      thin=10, param_updates){
  if(missing(param_updates)){
    param_updates <- paramUpdates(object)
  }
  mp <- McmcParams(iter=iter, burnin=burnin, thin=thin, param_updates=param_updates)
  mcmcParams(object, force=TRUE) <- mp
  object <- runBurnin(object)
  if(!iter(object) > 0) return(object)
  object <- runMcmc(object)
  modes(object) <- computeModes(object)
  object
}

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
#'      map(MarginalModelExample)
#' @param object an object of class MixtureModel.
#' @return map estimate of latent variable assignment for each observation
#' @export
map <- function(object) {
  estimates <- apply(probz(object), 1, which.max)
  names(estimates) <- names(y(object))
  estimates
}

setMethod("thetac", "MixtureModel", function(object) theta(chains(object)))

setMethod("thetaMean", "MixtureModel", function(object) colMeans(thetac(object)))

setMethod("sigmaMean", "MixtureModel", function(object) colMeans(sigmac(object)))

logLikc <- function(object) log_lik(chains(object))


#' Retrieve standard deviation of each component/batch mean at each iteration of the MCMC.
#'
#' @examples
#'      sigmac(MarginalModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A matrix of size N x K where N is the number of observations
#' and K is the number of components
#' @export
sigmac <- function(object) sigma(chains(object))

#' Retrieve mixture proportions at each iteration of the MCMC.
#'
#' @examples
#'      pic(MarginalModelExample)
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
#'      muc(MarginalModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of length N or matrix of size N x B, where N is the 
#' number of observations and B is the number of unique batches.
#' @export
muc <- function(object) mu(chains(object))

#' Retrieve overall mean averaged across MCMC simulations.
#'
#' @examples
#'      muMean(MarginalModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of size 1 or number of batches
#' @export
muMean <- function(object) {
  x <- muc(object)
  if(is(object, "BatchModel")){
    return(colMeans(x))
  }
  mean(x)
}

#' Retrieve overall standard deviation at each iteration of the MCMC.
#'
#' @examples
#'      tauc(MarginalModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of length N or matrix of size N x B, where N is the 
#' number of observations and B is the number of unique batches.
#' @export
tauc <- function(object) sqrt(tau2(chains(object)))

#' Retrieve overall standard deviation averaged across MCMC simulations.
#'
#' @examples
#'      tauMean(MarginalModelExample)
#' @param object an object of class MarginalModel or BatchModel
#' @return A vector of size 1 or number of batches
#' @export
tauMean <- function(object){
  x <- tauc(object)
  if(is(object, "BatchModel")){
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
  lp <- logPrior(chains(object))
  which.max(ll + lp)
}

setMethod("isMarginalModel", "MarginalModel", function(object) TRUE)
setMethod("isMarginalModel", "BatchModel", function(object) FALSE)

startAtTrueValues <- function(model, truth){
  theta(model) <- theta(truth)
  sigma2(model) <- sigma2(truth)
  p(model) <- p(truth)
  z(model) <- z(truth)
  if( isMarginalModel(truth) ){
    mu(model) <- mean(theta(truth))
    tau2(model) <- 1000
  } else {
    mu(model) <- colMeans(theta(truth))
    tau2(model) <- rep(1000, k(truth))
  }
  dataMean(model) <- computeMeans(model)
  dataPrec(model) <- 1/computeVars(model)
  zFreq(model) <- as.integer(table(z(model)))
  model
}

restartAtChainIndex <- function(model, index){
  ch <- chains(model)
  if(!isMarginalModel(model) ){
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
setReplaceMethod("mcmcParams", "MixtureModel", function(object, force=FALSE, value){
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


setMethod("zChain", "MixtureModel", function(object) chains(object)@z)

useModes <- function(object){
  m2 <- object
  theta(m2) <- modes(object)[["theta"]]
  sigma2(m2) <- modes(object)[["sigma2"]]
  tau2(m2) <- modes(object)[["tau2"]]
  nu.0(m2) <- modes(object)[["nu0"]]
  sigma2.0(m2) <- modes(object)[["sigma2.0"]]
  p(m2) <- modes(object)[["mixprob"]]
  zFreq(m2) <- as.integer(modes(object)[["zfreq"]])
  log_lik(m2) <- modes(object)[["loglik"]]
  logPrior(m2) <- modes(object)[["logprior"]]
  ##
  ## update z using the modal values from above
  ##
  if(is(object, "MarginalModel")){
    z(m2) <- update_z(m2)
  } else {
    z(m2) <- update_z_batch(m2)
  }
  m2
}



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
  if(isMarginalModel(model)){
    p <- p[, order(theta(map_model))]
  } else {
    p <- p[, order(mu(map_model))]
  }
  return(p)
}
