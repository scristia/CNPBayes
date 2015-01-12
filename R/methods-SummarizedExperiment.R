#' @export
setMethod("copyNumber", "SummarizedExperiment", function(object, ...){
  assays(object)[["medr"]]/1000
})

setMethod("fitMixtureModels", "SummarizedExperiment", function(object, mcmcp, K=1:5){
  cn <- copyNumber(object)[1, ]
  object$plate <- collapseBatch(object)
  freq <- table(object$plate)
##  if(any(freq < 25)){
##    message("Removing ", sum(freq < 25), " plate(s) with fewer than 25 observations")
##    batchid <- names(freq)[freq < 25]
##    object <- object[, !object$plate %in% batchid]
##    cn <- copyNumber(object)[1, ]
##  }
  message("Fitting ", length(K), " mixture models")
  fit <- foreach(j = seq_along(K)) %do% {
    cat(".")
    kk <- K[j]
    params <- ModelParams("batch", y=cn, k=kk,
                          batch=object$plate,
                          mcmc.params=mcmcp)
    modelk <- initializeModel(params)
    modelk <- posteriorSimulation(modelk, mcmcp)
    modelk
  }
  fit
})



## Allows restarting a chain
##setMethod("fitMixtureModels", "BatchModel", function(object, mcmcp, K=1:5){
##  ##
##  ## K is ignored
##  ##
##  params <- ModelParams("batch", y=y(object),
##                        k=2,
##                        batch=batch(object),
##                        mcmc.params=mcmcp)
##  model <- initializeModel(params)
##  message("Defining batch variable")
##  model <- collapseBatch(model, mcmcp)
##  message("Fitting ", length(K), " mixture models")
##  fit <- foreach(j = seq_along(K)) %do% {
##    cat(".")
##    kk <- K[j]
##    params <- ModelParams("batch", y=cn, k=kk,
##                          batch=batch(model),
##                          mcmc.params=mcmcp)
##    modelk <- initializeModel(params)
##    modelk <- posteriorSimulation(modelk, mcmcp)
##    modelk
##  }
##  fit
##})
