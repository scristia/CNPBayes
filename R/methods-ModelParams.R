#' Create an object of class 'ModelParams' to specify data, parameters, etc.
#'
#' @param type marginal or batch
#' @param y the data
#' @param k the number of components
#' @param batch a vector specifying batch labels
#' @param mcmc.params an object of class 'McmcParams'
#' @export
ModelParams <- function(type=c("marginal", "batch"),
                        y=numeric(),
                        k=2L,
                        batch,
                        mcmc.params=McmcParams(iter=1000, burnin=0)){
  if(missing(batch)){
    batch <- factor(rep("A", length(y)))
  } else {
    batch <- factor(batch)
  }
  new("ModelParams", type=match.arg(type), data=y, k=k, 
      batch=batch, mcmc.params=mcmc.params)
}

type <- function(object) object@type

#' @rdname k-method
#' @aliases k,ModelParams-method
setMethod("k", "ModelParams", function(object) object@k)

#' @rdname y-method
#' @aliases y,ModelParams-method
setMethod("y", "ModelParams", function(object) object@data)

#' @rdname batch-method
#' @aliases batch,ModelParams-method
setMethod("batch", "ModelParams", function(object) object@batch)

N <- function(object) length(y(object))

setMethod("show", "ModelParams", function(object){
  cat("An object of class 'ModelParams'\n")
  cat("      type: ", type(object), "\n")
  cat("         k: ", k(object), "\n")
  cat("n. batches: ", nBatch(object), "\n")
  cat("See mcmcParams() for thinning, burnin, and number of iterations after burnin.\n")
})
