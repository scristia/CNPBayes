#' @export
ModelParams <- function(type=c("marginal", "batch"),
                        y=numeric(),
                        k,
                        batch,
                        mcmc.params){
  if(missing(mcmc.params)){
    mcmc.params <- McmcParams(iter=1000, burnin=0)
  }
  y <- y[!is.na(y)]
  if(missing(batch)) batch <- rep("A", length(y))
  new("ModelParams", type=match.arg(type), data=y, k=k, batch=batch, mcmc.params=mcmc.params)
}

type <- function(object) object@type

setMethod("k", "ModelParams", function(object) object@k)

setMethod("y", "ModelParams", function(object) object@data)

#' @export
setMethod("batch", "ModelParams", function(object) object@batch)

mcmcParams <- function(object) object@mcmc.params

N <- function(object) length(y(object))

setMethod("show", "ModelParams", function(object){
  cat("An object of class 'ModelParams'\n")
  cat("      type: ", type(object), "\n")
  cat("         k: ", k(object), "\n")
  cat("n. batches: ", nBatch(object), "\n")
  cat("See mcmcParams() for thinning, burnin, and number of iterations after burnin.\n")
})
