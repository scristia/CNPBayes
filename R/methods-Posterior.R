Posterior <- function(hyperparams=Hyperparameters(),
                      theta=numeric(),
                      sigma2=numeric(),
                      mu=numeric(),
                      tau2=numeric(),
                      nu.0=numeric(),
                      sigma2.0=numeric(),
                      pi=numeric(),
                      data=numeric(),
                      data.mean=numeric(),
                      data.prec=numeric(),
                      z=factor(),
                      logpotential=numeric(),
                      mcmc.chains=McmcChains(),
                      batch=vector()){
  new("Posterior",
      hyperparams=hyperparams,
      theta=theta,
      sigma2=sigma2,
      mu=mu,
      tau2=tau2,
      nu.0=nu.0,
      sigma2.0=sigma2.0,
      pi=pi,
      data=data,
      data.mean=data.mean,
      data.prec=data.prec,
      z=z,
      logpotential=logpotential,
      mcmc.chains=mcmc.chains,
      batch=batch)
}

##setMethod("BatchModel", "SummarizedExperiment", function(object, ncomp, batch){
##  initializePosterior(means(object)[1, ], ncomp, batch)
##})
##
##setMethod("BatchModel", "numeric", function(object, ncomp, batch){
##  initializePosterior(object, ncomp, batch)
##})

setValidity("Posterior", function(object){
  msg <- TRUE
  if(length(theta(object)) > 0){
    if(!identical(sort(theta(object)), theta(object))){
      msg <- "theta must be strictly increasing..."
      return(msg)
    }
  }
})




splitDataByComponent <- function(object){
  k <- nComp(object)
  ylist <- split(y(object), factor(z(object), levels=seq_len(k)))
  ylist
}

calculateComponentMean <- function(object){
  ylist <- splitDataByComonent(object)
  sapply(ylist, mean, na.rm=TRUE)
}

calculateComponentPrec <- function(object){
  ylist <- splitDataByComonent(object)
  1/sapply(ylist, var, na.rm=TRUE)
}


posteriorModels <- function(Y, K, mcmcp){
  postlist <- foreach(k = K, .packages="CNPBayes") %dopar%{
    post <- initializePosterior(Y, ncomp=k)
    mcmcChains(post) <- McmcChains(post, mcmcp)
    post <- posteriorSimulation(post, mcmcp)
  }
  bic <- sapply(postlist, BIC)
  post <- postlist[[which.min(bic)]]
  post
}
