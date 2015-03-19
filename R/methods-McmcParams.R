#' @export
McmcParams <- function(iter=1000L, burnin=0L, thin, nStarts=1, nStartIter=200, checkLabels=FALSE){
  if(missing(thin)) thin <- rep(1L, length(iter))
  new("McmcParams", iter=iter, burnin=burnin, thin=thin, nstarts=nStarts, nstart_iter=nStartIter,
      check_labels=checkLabels)
}


setMethod("burnin", "McmcParams", function(object)  object@burnin)

setReplaceMethod("burnin", "McmcParams", function(object,value){
  object@burnin <- value
  object
})


thin <- function(object) object@thin
iter <- function(object) object@iter
savedIterations <- function(object)iter(object)/(max(thin(object), 1))


setMethod("show", "McmcParams", function(object){
  cat("An object of class 'McmcParams'\n")
  cat("   iterations:", paste(iter(object), collapse=","), "\n")
  cat("   burnin    :", paste(burnin(object), collapse=","),  "\n")
  cat("   thin      :", paste(thin(object), collapse=","), "\n")
  cat("   n starts  :", nStarts(object), "\n")
})


setMethod("[", "McmcParams", function(x, i, j, ..., drop=FALSE){
  ## allow one to pass a vector for burnin, thin, and iter slots
  if(length(iter(x))==1) return(x)
  if(!missing(i)){
    x@iter <- iter(x)[i]
    x@burnin <- burnin(x)[i]
    x@thin <- thin(x)[i]
  }
  x
})

setValidity("McmcParams", function(object){
  msg <- TRUE
  if(length(thin(object)) != length(burnin(object)) || length(thin(object)) != length(iter(object))){
    msg <- "thin, burnin, and iter vectors must be the same length"
    return(msg)
  }
  msg
})

setMethod("nStarts", "McmcParams", function(object) object@nstarts)
setMethod("checkLabels", "McmcParams", function(object) object@check_labels)
setMethod("nStartIter", "McmcParams", function(object) object@nstart_iter)

#' @export
mcmcpList <- function(test=FALSE, iter=c(500, 3000, 3000, 3000),
                      nStarts=20, nStartIter=100){
  if(test){
    mcmcp.list <- list(McmcParams(iter=rep(2, 4),
                                  burnin=rep(1,4),
                                  thin=rep(1,4),
                                  nStarts=2,
                                  nStartIter=2),
                       McmcParams(iter=2, burnin=1, thin=1))
  } else {
    mcmcp.list <- list(McmcParams(iter=iter,
                                  burnin=as.integer(pmax(1, iter/10)),
                                  thin=as.integer(pmax(1, iter/1000)),
                                  nStarts=nStarts,
                                  nStartIter=nStartIter),
                       McmcParams(iter=1000, burnin=100, thin=1))
  }
  mcmcp.list
}
