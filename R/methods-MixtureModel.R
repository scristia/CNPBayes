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

setMethod("hyperParams", "MixtureModel", function(object) object@hyperparams)

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



#' @export
sigma <- function(object) sqrt(sigma2(object))

#' @export
tau <- function(object) sqrt(tau2(object))

setMethod("nu.0", "MixtureModel", function(object) object@nu.0)

setMethod("sigma2.0", "MixtureModel", function(object) object@sigma2.0)

sigma.0 <- function(object) sqrt(sigma2.0(object))



#' @export
p <- function(object) object@pi

nComp <- function(object) length(p(object))
dataMean <- function(object) object@data.mean
dataPrec <- function(object) object@data.prec
dataSd <- function(object) sqrt(1/dataPrec(object))

#' @export
logpotential <- function(object) object@logpotential

setMethod("k", "MixtureModel", function(object) object@k)

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

setReplaceMethod("logpotential", "MixtureModel", function(object, value){
  object@logpotential <- value
  object
})

setMethod("mcmcChains", "MixtureModel", function(object) object@mcmc.chains)

setMethod("chains", "MixtureModel", function(object) object@mcmc.chains)

setMethod("dat", "MixtureModel", function(object) object@data)

setReplaceMethod("dat", "MixtureModel", function(object, value) {
  object@data <- value
  object
})

setReplaceMethod("mcmcChains", "MixtureModel", function(object, value){
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
  ##post.sd <- paste(round(sqrt(sigma2(object)), 2), collapse=", ")
  sds <- showSigmas(object)
  mns <- showMeans(object)
  cat("     sigma (s)   :", sds, "\n")
  cat("     theta (s)   :", mns, "\n")
  cat("     sigma2.0 (s):", round(sigma2.0(object), 2), "\n")
  cat("     nu.0 (s)    :", nu.0(object), "\n")
  ##cat("     tau2 (s)    :", round(tau2(object), 2),  "\n")
  ##cat("     mu (s)      :", round(mu(object), 2), "\n")
  cat("     logprior(s):", round(logPrior(object), 2), "\n")
  cat("     loglik (s)  :", round(logLik(object), 2), "\n")
})

setMethod("updateMixProbs", "MixtureModel", function(object){
  ##alpha.n <- updateAlpha(object)
  ##as.numeric(rdirichlet(1, alpha.n))
  .Call("update_p", object)
})

.updateMix <- function(object){
  alpha.n <- updateAlpha(object)
  as.numeric(rdirichlet(1, alpha.n))
}


setMethod("alpha", "MixtureModel", function(object) alpha(hyperParams(object)))

setMethod("updateAlpha", "MixtureModel", function(object){
  alpha(object) + table(z(object))
})

setMethod("y", "MixtureModel", function(object) object@data)
setMethod("batch", "MixtureModel", function(object) object@batch)
setMethod("z", "MixtureModel", function(object) object@z)

updateAll <- function(post, is_burnin=FALSE){
  K <- k(post)
  up <- paramUpdates(mcmcParams(post))
  if(up["theta"] > 0)
    theta(post) <- updateTheta(post)
  if(up["sigma2"] > 0)
    sigma2(post) <- updateSigma2(post)
  if(up ["mu"] > 0)
    mu(post) <- updateMu(post)
  if(up["tau2"] > 0 )
    tau2(post) <- updateTau2(post)
  if(up["sigma2.0"] > 0)
    sigma2.0(post) <- updateSigma2.0(post)
  if(up["nu.0"] > 0)
    nu.0(post) <- updateNu.0(post)
  if(up["p"] > 0)
    p(post) <- updateMixProbs(post)
  if(up["z"] > 0){
    z(post) <- updateZ(post)
    zFreq(post) <- as.integer(table(z(post)))
  }
  dataMean(post) <- computeMeans(post)
  ##dataPrec(post) <- 1/computeVars(post)
  dataPrec(post) <- computePrec(post)
  if(!is_burnin){
    ## faster if we omit these steps during burnin
    ##logpotential(post) <- computeLogLikxPrior(post)
    logPrior(post) <- computePrior(post)
    logLik(post) <- computeLoglik(post)
  }
  post
}


setMethod("computePrec", "MarginalModel", function(object){
  .Call("compute_prec", object)
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

setMethod("probz", "MixtureModel", function(object) {
  object@probz/iter(object)
})


setMethod("runBurnin", "MarginalModel", function(object){
  .Call("mcmc_marginal_burnin", object, mcmcParams(object))

})

setMethod("runBurnin", "BatchModel", function(object){
  .Call("mcmc_batch_burnin", object, mcmcParams(object))
  ##.runBurnin(object)
})

setMethod("runMcmc", "MarginalModel", function(object){
  .Call("mcmc_marginal", object, mcmcParams(object))
})

setMethod("runMcmc", "BatchModel", function(object){
  .Call("mcmc_batch", object, mcmcParams(object))
  ##.runMcmc(object, mcmcParams(object))
})

.runBurnin <- function(object){
  message("Running burnin...")
  for(b in seq_len(burnin(object))){
    object <- updateAll(object, is_burnin=TRUE)
  }
  object
}

.runMcmc <- function(object){
  if(burnin(object) > 0) message("Burnin finished. Running additional MCMC simulations...")
  S <- 2:iter(object)
  do_thin <- thin(object) > 1
  T <- seq_len(thin(object))
  object <- moveChain(object, 1)
  for(s in S){
    object <- updateAll(object, FALSE)
    object <- moveChain(object, s)
    ## update without moving chain
    if(do_thin){
      for(t in T) object <- updateAll(object, TRUE)
    }
  }
  object
}

multipleStarts <- function(object){
  if(k(object)==1) return(object)
  mcmcp <- mcmcParams(object)
  message("Running ", nStarts(mcmcp), " chains for model k=", k(object))
  ##  if(is(object, "BatchModel")){
  ##    mmod <- replicate(nStarts(mcmcp), BatchModel(y(object), mcmc.params=mcmcp,
  ##                                                 hypp=hyperParams(object), k=k(object),
  ##                                                 batch=batch(object)))
  ##  } else {
  mmod <- replicate(nStarts(mcmcp), MarginalModel(y(object), mcmc.params=mcmcp,
                                                  hypp=hyperParams(object), k=k(object)))
##  browser()
##  model1 <- mmod[[1]]
##  K <- k(object)
##  if(K==2){
##    x <- MarginalModel(y(object), k=2)
##    z(model1)[y(model1) < -0.4] <- 1L
##    z(model1)[y(model1) >= -0.4] <- 2L
##    zFreq(model1) <- as.integer(table(z(model1)))
##    dataMean(model1) <- computeMeans(model1)
##    dataPrec(model1) <- computePrec(model1)
##    mmod[[1]] <- model1
##  }
##  if(K >= 3){
##    if(isMarginalModel(object)){
##      hypp <- hyperParams(object)
##    } else{
##      ## use defaults for marginal model
##      hypp <- Hyperparameters("marginal", k=K)
##    }
##    x <- MarginalModel(y(object), k=K, hypp=hypp)
##    z(x)[y(x) < -1] <- 1L
##    z(x)[y(x) > -1 & y(x) <=-0.2] <- 2L
##    z(x)[y(x) > -0.2] <- 3L
##    if(K == 4){
##      z(x)[y(x) > 0.2] <- 4L
##    }
##    zFreq(x) <- as.integer(table(z(x)))
##    theta(x) <- dataMean(x) <- computeMeans(x)
##    dataPrec(x) <- 1/computeVars(x)
##    sigma2(x) <- 1/dataPrec(x)
##    computeLoglik(x)
##    mcmcParams(x, force=TRUE) <- mcmcParams(object)
##    mmod[[1]] <- x
##  }
  ##}
  ##
  ## TODO: Remove nStartIter slot
  ##
  models <- suppressMessages(lapply(mmod, runBurnin))
  lp <- sapply(models, logLik)
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

computeLogLikForRandomStarts <- function(seed, params, hypp, 
                                         return.model=FALSE) {
  set.seed(seed)
  model <- initializeModel(params, hypp)
  if(return.model) return(model)
  logLikData(model)
}



setMethod("posteriorSimulation", "MixtureModel", function(object){
  .posteriorSimulation(object)
})

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


.class_constructor <- function(class, ...){
  new(class, ...)
}

homozygousComponent <- function(y){
  ##mean(y < -1.5) > 0.005
  sum(y < -1.5) >= 3
}

setMethod("computeMeans", "UnivariateMarginalModel", function(object){
  median(y(object), na.rm=TRUE)
})

setMethod("computeVars", "UnivariateMarginalModel", function(object){
  var(y(object), na.rm=TRUE)
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

setMethod("eta.0", "MixtureModel", function(object) eta.0(hyperParams(object)))
setMethod("eta.0", "Hyperparameters", function(object) object@eta.0)
setMethod("m2.0", "MixtureModel", function(object) m2.0(hyperParams(object)))
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



HardyWeinberg <- function(object){
  warning("requires cn inference...")
  tab <- table(map(object))
  ##freq <- setNames(table(zz), c("AA", "AB", "BB"))
  if(length(tab) == 1) return(NA)
  if(length(tab)==2){
    tab <- setNames(c(sort(tab, decreasing=TRUE), 0), c("AA", "AB", "BB"))
  }
  if(length(tab)==3){
    tab <- setNames(tab, c("AA", "AB", "BB"))
  }
  if(length(tab) > 3){
    tab <- setNames(tab[1:3], c("AA", "AB", "BB"))
  }
  HWChisq(tab)
}

setMethod("hwe", "MixtureModel", function(object) object@hwe)

#' @export
map <- function(object) apply(probz(object), 1, which.max)

setMethod("fitMixtureModels", "numeric", function(object, mcmcp, K=1:5){
  message("Fitting ", length(K), " models")
  fit <- vector("list", length(K))
  for(j in seq_along(K)){
    cat(".")
    mp <- mcmcp[j]
    kk <- K[j]
    params <- ModelParams("marginal", y=object, k=kk, mcmc.params=mp)
    model <- posteriorSimulation(params, mp)
    ##m.y(model) <- marginalY(model, mp)
    fit[[j]] <- model
  }
  fit
})

setMethod("thetac", "MixtureModel", function(object) theta(mcmcChains(object)))

setMethod("thetaMean", "MixtureModel", function(object) colMeans(thetac(object)))

setMethod("sigmaMean", "MixtureModel", function(object) colMeans(sigmac(object)))

logpotentialc <- function(object) logpotential(mcmcChains(object))

logLikc <- function(object) logLik(mcmcChains(object))


#' @export
sigmac <- function(object) sigma(mcmcChains(object))

#' @export
pic <- function(object) p(mcmcChains(object))

setMethod("pMean", "MixtureModel", function(object){
  colMeans(pic(object))
})

#' @export
muc <- function(object) mu(mcmcChains(object))

#' @export
muMean <- function(object) colMeans(muc(object))

#' @export
tauc <- function(object) sqrt(tau2(mcmcChains(object)))

#' @export
tauMean <- function(object) colMeans(tauc(object))

setMethod("modes", "MixtureModel", function(object) object@modes)

setReplaceMethod("modes", "MixtureModel", function(object, value) {
  object@modes <- value
  object
})

setMethod("logLik", "MixtureModel", function(object){
  object@loglik
})

setReplaceMethod("logLik", "MixtureModel", function(object, value){
  object@loglik <- value
  object
})


setMethod("m.y", "MixtureModel", function(object) object@m.y)
setReplaceMethod("m.y", "MixtureModel", function(object,value){
  object@m.y <- value
  object
})


setMethod("orderTheta", "MixtureModel", function(object) object@theta_order)

setReplaceMethod("orderTheta", "MixtureModel", function(object, value) {
  object@theta_order <- value
  object
})

modalLogLik <- function(object){
  x <- mcmcChains(object)
  max(logLik(x))
}

argmaxLogLik <- function(object){
  x <- mcmcChains(object)
  which.max(logLik(x))
}

argMax <- function(object){
  ll <- logLik(chains(object))
  lp <- logPrior(chains(object))
  which.max(ll + lp)
}



setMethod("isMarginalModel", "MarginalModel", function(object) TRUE)
setMethod("isMarginalModel", "BatchModel", function(object) FALSE)

#' @export
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


setMethod("zFreq", "MixtureModel", function(object) object@zfreq)

setReplaceMethod("zFreq", "MixtureModel", function(object, value){
  object@zfreq <- value
  object
})

setMethod("mcmcParams", "MixtureModel", function(object) object@mcmc.params )

setMethod("iter", "MixtureModel", function(object) iter(mcmcParams(object)))
setMethod("nStarts", "MixtureModel", function(object) nStarts(mcmcParams(object)))

setReplaceMethod("nStarts", "MixtureModel", function(object, value){
  mcmcParams(object)@nstarts <- as.integer(value)
  object
})

setMethod("thin", "MixtureModel", function(object) thin(mcmcParams(object)))
setMethod("burnin", "MixtureModel", function(object) burnin(mcmcParams(object)))

setReplaceMethod("burnin", "MixtureModel", function(object, value){
  burnin(mcmcParams(object)) <- value
  object
})

setReplaceMethod("iter", "MixtureModel", function(object, force=FALSE, value){
  mp <- mcmcParams(object)
  iter(mp) <- value
  mcmcParams(object, force=force) <- mp
  object
})

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

##
## Any update of the mcmc parameters triggers an update of the chains.
## This is important because it initializes the chains to be of the
## right size. **An important exception: setting iter to zero or a
## negative value -does not- trigger an update of the chains -- chains
## are kept as is
##
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
        mcmc_chains <- mcmcChains(object)[index, ]
      }
      mcmcChains(object) <- mcmc_chains
      return(object)
    }
  }
  ## if we've got to this point, it must be safe to update mcmc.params
  ## (i.e., size of chains is not effected)
  object@mcmc.params <- value
  object
})



setMethod("zChain", "MixtureModel", function(object) chains(object)@z)
setReplaceMethod("zChain", "MixtureModel", function(object, value){
  mcmcChains(object)@z <- value
  object
})

useModes <- function(object){
  m2 <- object
  theta(m2) <- modes(object)[["theta"]]
  sigma2(m2) <- modes(object)[["sigma2"]]
  tau2(m2) <- modes(object)[["tau2"]]
  nu.0(m2) <- modes(object)[["nu0"]]
  sigma2.0(m2) <- modes(object)[["sigma2.0"]]
  p(m2) <- modes(object)[["mixprob"]]
  zFreq(m2) <- as.integer(modes(object)[["zfreq"]])
  logLik(m2) <- modes(object)[["loglik"]]
  logPrior(m2) <- modes(object)[["logprior"]]
  ##
  ## update z using the modal values from above
  ##
  if(is(object, "MarginalModel")){
    z(m2) <- .Call("update_z", m2)
  } else {
    z(m2) <- .Call("update_z_batch", m2)
  }
  m2
}


mapModel <- function(model){
  model2 <- restartAtChainIndex(model, argMax(model))
  model2
}



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
