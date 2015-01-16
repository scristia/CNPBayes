#' @export
setMethod("copyNumber", "SummarizedExperiment", function(object, ...){
  assays(object)[["medr"]]/1000
})

setMethod("fitMixtureModels", "SummarizedExperiment", function(object, mcmcp, K=1:5){
  cn <- copyNumber(object)[1, ]
  object$plate <- collapseBatch(object)
  freq <- table(object$plate)
  message("Fitting ", length(K), " mixture models")
  fit <- foreach(j = seq_along(K)) %do% {
    cat(".")
    kk <- K[j]
    mp <- mcmcp[j]
    params <- ModelParams("batch", y=cn, k=kk,
                          batch=object$plate,
                          mcmc.params=mp)
    modelk <- initializeBatchModel(params)
    modelk <- suppressMessages(posteriorSimulation(modelk, mp))
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
