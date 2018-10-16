#' @include methods-MixtureModel.R
NULL

setMethod("runBurnin", "MultiBatchModel", function(object){
  cpp_burnin(object, mcmcParams(object))
})

setMethod("runMcmc", "MultiBatchModel", function(object){
  cpp_mcmc(object, mcmcParams(object))
})

setMethod("runBurnin", "MultiBatchPooled", function(object){
  burnin_multibatch_pvar(object, mcmcParams(object))
})

setMethod("runMcmc", "MultiBatchPooled", function(object){
  mcmc_multibatch_pvar(object, mcmcParams(object))
})

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
    warning("label switching: model k=", k(post))
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

#' @rdname posteriorSimulation-method
#' @aliases posteriorSimulation,SingleBatchModel-method
setMethod("posteriorSimulation", "SingleBatchModel", function(object){
  object <- as(object, "MultiBatchModel")
  .posteriorSimulation2(object)
})

#' @rdname posteriorSimulation-method
#' @aliases posteriorSimulation,SingleBatchPooled-method
setMethod("posteriorSimulation", "SingleBatchPooled", function(object){
  object <- as(object, "MultiBatchPooled")
  .posteriorSimulation2(object)
})



setMethod("runBurnin", "MultiBatchModel", function(object){
  cpp_burnin(object, mcmcParams(object))
})

setMethod("runBurnin", "MultiBatch", function(object){
  mbm <- as(object, "MultiBatchModel")
  mbm <- runBurnin(mbm, mcmcParams(mbm))
  ## The MBM model does not have the original data
  mb <- as(mbm, "MultiBatch")
  ## do not use replacement method
  mb@data <- assays(object)
  down_sample(mb) <- down_sample(object)
  mb
})

setMethod("runMcmc", "MultiBatch", function(object){
  mbm <- as(object, "MultiBatchModel")
##  mp <- mcmcParams(mbm)
##  up <- .param_updates(mp)
##  up[ ] <- 0L
##  up[ c("theta", "sigma2") ] <- 1L
##  mp@param_updates <- up
##  mcmcParams(mbm) <- mp
  mbm <- runMcmc(mbm)
  mb <- as(mbm, "MultiBatch")
  ## do not use replacement method
  mb@data <- assays(object)
  down_sample(mb) <- down_sample(object)
  mb
})

setMethod("posteriorSimulation", "MultiBatch", function(object){
  object <- runBurnin(object)
  if(!isOrdered(object)) label_switch(object) <- TRUE
  object <- sortComponentLabels(object)
  if( iter(object) < 1 ) return(object)
  object <- runMcmc(object)
  modes(object) <- computeModes(object)
  if(isOrdered(object)){
    label_switch(object) <- FALSE
    return(object)
  }
  ## not ordered: try additional MCMC simulations
  label_switch(object) <- TRUE
  object <- sortComponentLabels(object)
  ## reset counter for posterior probabilities
  object@probz[] <- 0
  object <- runMcmc(object)
  modes(object) <- computeModes(object)
  if(isOrdered(object)){
    label_switch(object) <- FALSE
    return(object)
  }
  label_switch(object) <- TRUE
  object <- sortComponentLabels(object)
  object
})
