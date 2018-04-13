#' @include methods-MixtureModel.R
NULL

setMethod("runBurnin", "MultiBatchModel", function(object){
  cpp_burnin(object, mcmcParams(object))
})

setMethod("runBurnin", "TrioBatchModel", function(object){
  trios_burnin(object, mcmcParams(object))
})

setMethod("runMcmc", "MultiBatchModel", function(object){
  cpp_mcmc(object, mcmcParams(object))
})

setMethod("runMcmc", "TrioBatchModel", function(object){
  trios_mcmc(object, mcmcParams(object))
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
