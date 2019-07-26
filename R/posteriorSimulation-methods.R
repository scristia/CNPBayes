#' @include methods-MixtureModel.R
NULL

setMethod("runBurnin", "MultiBatchModel", function(object){
  cpp_burnin(object)
})

setMethod("runBurnin", "TrioBatchModel", function(object){
  trios_burnin(object)
})

setMethod("runMcmc", "MultiBatchModel", function(object){
  cpp_mcmc(object)
})

setMethod("runMcmc", "TrioBatchModel", function(object){
  trios_mcmc(object)
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

anyWarnings <- function(object){
  fl <- flags(object)
  fl2 <- sapply(fl, as.logical)
  any(fl2)
}

.posteriorSimulation2 <- function(post){
  post <- runBurnin(post)
  post <- sortComponentLabels(post)
  post <- runMcmc(post)
  if(!isOrdered(post)) label_switch(post) <- TRUE
  modes(post) <- computeModes(post)
  return(post)
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
  ##browser()
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
  mbm <- sortComponentLabels(mbm)
  mbm <- runMcmc(mbm)
  if(!isOrdered(mbm)) label_switch(mbm) <- TRUE
  mbm <- sortComponentLabels(mbm)
  mb <- revertBack(object, mbm)
  summaries(mb) <- summarizeModel(mb)
  return(mb)
})

setGeneric("mcmc_homozygous", function(object) standardGeneric("mcmc_homozygous"))

setGeneric("burnin_homozygous", function(object) standardGeneric("burnin_homozygous"))

setMethod("burnin_homozygous", "MultiBatchModel", function(object) {
  mb_homozygous_burnin(object)
})

setMethod("burnin_homozygous", "MultiBatchPooled", function(object) {
  mbp_homozygous_burnin(object)
})

setMethod("burnin_homozygous", "MultiBatch", function(object) {
  mbm <- as(object, "MultiBatchModel")
  mbm <- mb_homozygous_burnin(mbm)
  mbm <- sortComponentLabels(mbm)
  mb <- revertBack(object, mbm)
  mb
})

setMethod("burnin_homozygous", "MultiBatchP", function(object) {
  mbm <- as(object, "MultiBatchPooled")
  mbm <- mbp_homozygous_burnin(mbm)
  mbm <- sortComponentLabels(mbm)
  mb <- revertBack(object, mbm)
  mb
})

setMethod("mcmc_homozygous", "MultiBatchModel", function(object) {
  mb_homozygous_mcmc(object)
})

setMethod("mcmc_homozygous", "MultiBatchPooled", function(object) {
  mbp_homozygous_mcmc(object)
})

setMethod("mcmc_homozygous", "MultiBatch", function(object){
  mbm <- as(object, "MultiBatchModel")
  mbm <- burnin_homozygous(mbm)
  ##  frac <- min(mbm@.internal.counter/iter(mbm), 1)
  ##  if(frac > 0.9){
  ##    stop("Missing component")
  ##    mb <- revertBack(object, mbm)
  ##    summaries(mb) <- summarizeModel(mb)
  ##  }
  ##mbm <- sortComponentLabels(mbm)
  mbm <- mcmc_homozygous(mbm)
  if(!isOrdered(mbm)) label_switch(mbm) <- TRUE
  ##mbm <- sortComponentLabels(mbm)
  mb <- revertBack(object, mbm)
  summaries(mb) <- summarizeModel(mb)
  return(mb)
})

setMethod("mcmc_homozygous", "MultiBatchP", function(object){
  mbm <- as(object, "MultiBatchPooled")
  mbm <- burnin_homozygous(mbm)
  ##mbm <- sortComponentLabels(mbm)
  mbm <- mcmc_homozygous(mbm)
  if(!isOrdered(mbm)) label_switch(mbm) <- TRUE
  ##mbm <- sortComponentLabels(mbm)
  mb <- revertBack(object, mbm)
  summaries(mb) <- summarizeModel(mb)
  return(mb)
})

setMethod("posteriorSimulation", "MultiBatchP", function(object){
  if(any(is.na(dataMean(object)))){
    summaries(object)[["data.mean"]] <- computeMeans(object)
    summaries(object)[["data.prec"]] <- computePrec(object)
  }
  mbm <- as(object, "MultiBatchPooled")
  mbm@.internal.counter <- 0L
  mbm <- runBurnin(mbm)
  mbm <- sortComponentLabels(mbm)
  if(mbm@.internal.counter / burnin(mbm) > 0.9){
    ## there were no valid updates for the z chain
    ## no point in continuing
    warning("Burnin failed to find valid updates for z. Stopping MCMC after burnin.")
    mb <- revertBack(object, mbm)
    summaries(mb) <- summarizeModel(mb)
    return(mb)
  }
  mbm <- runMcmc(mbm)
  label_switch(mbm) <- !isOrdered(mbm)
  mbm <- sortComponentLabels(mbm)
  mb <- revertBack(object, mbm)
  summaries(mb) <- summarizeModel(mb)
  return(mb)
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
  object <- runBurnin(object)
  if(!isOrdered(object)) label_switch(object) <- TRUE
  object <- sortComponentLabels(object)
  object <- runMcmc(object)
  object
})
