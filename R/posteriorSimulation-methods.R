#' @include methods-MixtureModel.R
NULL

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

mcmcMultiBatchPooled <- function(object){
  .posteriorSimulation(object)
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


.posteriorSimulation <- function(post, params=psParams()){
  .Deprecated("see .posteriorSimulation2")
  if(nStarts(post) > 1){
    post <- multipleStarts2(post)
  }
  if(burnin(post) > 0 ){
    post <- runBurnin(post)
  }
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

#' @rdname posteriorSimulation-method
#' @aliases posteriorSimulation,MixtureModel-method
setMethod("posteriorSimulation", "MixtureModel", function(object){
  .posteriorSimulation2(object)
})
