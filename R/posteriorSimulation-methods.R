#' @include methods-MixtureModel.R
NULL

setMethod("runBurnin", "MultiBatchModel", function(object){
  cpp_burnin(object)
})

setMethod("runBurnin", "TrioBatchModel", function(object){
  trios_burnin(object, mcmcParams(object))
})

setMethod("runMcmc", "MultiBatchModel", function(object){
  cpp_mcmc(object)
})

setMethod("runMcmc", "TrioBatchModel", function(object){
  trios_mcmc(object, mcmcParams(object))
})


setMethod("runBurnin", "MultiBatchPooled", function(object){
  burnin_multibatch_pvar(object, mcmcParams(object))
})

setMethod("runMcmc", "MultiBatchPooled", function(object){
  tmp <- mcmc_multibatch_pvar(object, mcmcParams(object))
  tmp
})


posteriorSimulationNoMendel <- function(object){
  chains(object)@is_mendelian[] <- 0L
  post <- burnin_nomendelian_update(object)
  if(!isOrdered(post)) label_switch(post) <- TRUE
  post <- sortComponentLabels(post)
  if( iter(post) < 1 ) return(post)
  post <- mcmc_nomendelian_update(post)
  modes(post) <- computeModes(post)
  if(isOrdered(post)){
    label_switch(post) <- FALSE
    return(post)
  }
  ## not ordered: try additional MCMC simulations
  label_switch(post) <- TRUE
  post <- sortComponentLabels(post)
  ## reset counter for posterior probabilities
  post@probz[] <- 0
  chains(post)@is_mendelian[] <- 0L
  post <- mcmc_nomendelian_update(post)
  modes(post) <- computeModes(post)
  ##mcmcParams(post) <- mp.orig
  if(isOrdered(post)){
    label_switch(post) <- FALSE
    return(post)
  }
  label_switch(post) <- TRUE
  if(params[["warnings"]]) {
    ##
    ## at this point, we've tried to run the twice after burnin and we still
    ## have mixing. Most likely, we are fitting a model with k too big
    ##warning("label switching: model k=", k(post))
  }
  post <- sortComponentLabels(post)
  post
}

.posteriorSimulation2 <- function(post, params=psParams()){
  post <- runBurnin(post)
  if(!isOrdered(post)) label_switch(post) <- TRUE
  post <- sortComponentLabels(post)
  if( iter(post) < 1 ) return(post)
  post <- runMcmc(post)
  modes(post) <- computeModes(post)
  if(isOrdered(post)){
    label_switch(post) <- FALSE
    return(post)
  }
  ## not ordered: try additional MCMC simulations
  label_switch(post) <- TRUE
  post <- sortComponentLabels(post)
  ## reset counter for posterior probabilities
  post@probz[] <- 0
  post <- runMcmc(post)
  modes(post) <- computeModes(post)
  ##mcmcParams(post) <- mp.orig
  if(isOrdered(post)){
    label_switch(post) <- FALSE
    return(post)
  }
  label_switch(post) <- TRUE
  if(params[["warnings"]]) {
    ##
    ## at this point, we've tried to run the twice after burnin and we still
    ## have mixing. Most likely, we are fitting a model with k too big
    ##warning("label switching: model k=", k(post))
  }
  post <- sortComponentLabels(post)
  post
}

posteriorSimulationPooled <- function(object, iter=1000,
                                      burnin=1000,
                                      thin=10,
                                      param_updates){
  if(missing(param_updates)){
    param_updates <- paramUpdates(object)
  }
  mp <- McmcParams(iter=iter, burnin=burnin, thin=thin,
                   param_updates=param_updates)
  mcmcParams(object, force=TRUE) <- mp
  object <- runBurnin(object)
  object <- sortComponentLabels(object)
  if(!iter(object) > 0) return(object)
  object <- runMcmc(object)
  modes(object) <- computeModes(object)
  object <- sortComponentLabels(object)
  object
}



#' @rdname posteriorSimulation-method
#' @aliases posteriorSimulation,MixtureModel-method
setMethod("posteriorSimulation", "MixtureModel", function(object){
  .posteriorSimulation2(object)
})

setMethod("runBurnin", "MultiBatchModel", function(object){
  cpp_burnin(object)
})

setMethod("runBurnin", "MultiBatch", function(object){
  mbm <- as(object, "MultiBatchModel")
  mbm <- runBurnin(mbm, mcmcParams(mbm))
  ## The MBM model does not have the original data
  mb <- as(mbm, "MultiBatch")
  mb <- revertBack(object, mbm)
  mb
})

setMethod("runMcmc", "MultiBatch", function(object){
  mbm <- as(object, "MultiBatchModel")
  mbm <- runMcmc(mbm)
  mb <- revertBack(object, mbm)
  mb
})

setMethod("revertBack", "MultiBatch", function(object, mbm){
  mb <- as(mbm, "MultiBatch")
  ## do not use replacement method
  mb@data <- assays(object)
  ## using specs replacement method has side effects
  mb@specs$number_obs <- nrow(assays(object))
  mb
})

setMethod("revertBack", "MultiBatchP", function(object, mbm){
  mb <- as(mbm, "MultiBatchP")
  ## do not use replacement method
  mb@data <- assays(object)
  ## using specs replacement method has side effects
  mb@specs$number_obs <- nrow(assays(object))
  mb
})

setMethod("posteriorSimulation", "MultiBatch", function(object){
  mbm <- as(object, "MultiBatchModel")
  mbm <- runBurnin(mbm)
  if(!isOrdered(mbm)) label_switch(mbm) <- TRUE
  mbm <- sortComponentLabels(mbm)
  if( iter(mbm) < 1 ) return(mbm)
  mbm <- runMcmc(mbm)
  modes(mbm) <- computeModes(mbm)
  if(isOrdered(mbm)){
    label_switch(mbm) <- FALSE
    mb <- revertBack(object, mbm)
    return(mb)
  }
  ## not ordered: try additional MCMC simulations
  label_switch(mbm) <- TRUE
  mbm <- sortComponentLabels(mbm)
  ## reset counter for posterior probabilities
  mbm@probz[] <- 0
  mbm <- runMcmc(mbm)
  modes(mbm) <- computeModes(mbm)
  if(isOrdered(mbm)){
    label_switch(mbm) <- FALSE
    mb <- revertBack(object, mbm)
    return(mb)
  }
  label_switch(mbm) <- TRUE
  mbm <- sortComponentLabels(mbm)
  mb <- revertBack(object, mbm)
  mb
})

setMethod("posteriorSimulation", "MultiBatchP", function(object){
  mbm <- as(object, "MultiBatchPooled")
  mbm <- runBurnin(mbm)
  if(!isOrdered(mbm)) label_switch(mbm) <- TRUE
  mbm <- sortComponentLabels(mbm)
  if( iter(mbm) < 1 ) return(mbm)
  mbm <- runMcmc(mbm)
  modes(mbm) <- computeModes(mbm)
  if(isOrdered(mbm)){
    label_switch(mbm) <- FALSE
    mb <- revertBack(object, mbm)
    return(mb)
  }
  ## not ordered: try additional MCMC simulations
  label_switch(mbm) <- TRUE
  mbm <- sortComponentLabels(mbm)
  ## reset counter for posterior probabilities
  mbm@probz[] <- 0
  mbm <- runMcmc(mbm)
  modes(mbm) <- computeModes(mbm)
  if(isOrdered(mbm)){
    label_switch(mbm) <- FALSE
    mb <- revertBack(object, mbm)
    return(mb)
  }
  label_switch(mbm) <- TRUE
  mbm <- sortComponentLabels(mbm)
  mb <- revertBack(object, mbm)
  mb
})

setMethod("posteriorSimulation", "MultiBatchList", function(object){
  for(i in seq_along(object)){
    object[[i]] <- posteriorSimulation(object[[i]])
  }
  object
})

setMethod("posteriorSimulation", "list", function(object){
  for(i in seq_along(object)){
    object[[i]] <- posteriorSimulation(object[[i]])
  }
  object
})

#' @rdname posteriorSimulation-method
#' @aliases posteriorSimulation,TrioBatchModel-method
setMethod("posteriorSimulation", "TrioBatchModel", function(object){
  chains(object)@is_mendelian[] <- 0L
  m <- trios_burnin(object, mcmcParams(object))
  m <- trios_mcmc(object, mcmcParams(object))
  if(!isOrdered(m)) label_switch(m) <- TRUE
  m <- sortComponentLabels(m)
  m
})



