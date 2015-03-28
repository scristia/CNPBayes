#' @export
McmcParams <- function(iter=1000L, burnin=0L, thin, nStarts=1,
                       nStartIter=200, checkLabels=FALSE,
                       param_updates=.param_updates()){
  if(missing(thin)) thin <- rep(1L, length(iter))
  new("McmcParams", iter=as.integer(iter),
      burnin=as.integer(burnin),
      thin=as.integer(thin),
      nstarts=as.integer(nStarts),
      nstart_iter=as.integer(nStartIter),
      check_labels=checkLabels,
      param_updates=param_updates)
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
  up <- paramUpdates(object)
  if(!identical(names(up), names(.param_updates()))){
    msg <- "vector for slot param_updates should have same names as .param_updates()"
    return(msg)
  }
})

setMethod("nStarts", "McmcParams", function(object) object@nstarts)

setReplaceMethod("nStarts", "McmcParams", function(object, value){
  object@nstarts <- as.integer(value)
  object
})

setMethod("checkLabels", "McmcParams", function(object) object@check_labels)
setMethod("nStartIter", "McmcParams", function(object) object@nstart_iter)

#' @export
mcmcpList <- function(test=FALSE, iter=c(500, 3000, 3000, 3000),
                      nStarts=20, nStartIter=100){
  if(test){
    mcmcp.list <- list(McmcParams(iter=rep(2L, 4),
                                  burnin=rep(1L,4),
                                  thin=rep(1L,4),
                                  nStarts=2L,
                                  nStartIter=2L),
                       McmcParams(iter=2L, burnin=1L, thin=1L))
  } else {
    mcmcp.list <- list(McmcParams(iter=as.integer(iter),
                                  burnin=as.integer(pmax(1, iter/10)),
                                  thin=as.integer(pmax(1, iter/1000)),
                                  nStarts=as.integer(nStarts),
                                  nStartIter=as.integer(nStartIter)),
                       McmcParams(iter=1000L, burnin=100L, thin=1L))
  }
  mcmcp.list
}


setMethod("paramUpdates", "McmcParams", function(x){
  x@param_updates
})


.param_updates <- function(x){
  x <- setNames(rep(1L, 8),
                c("theta", "sigma2", "p", "mu", "tau2", "nu.0", "sigma2.0", "z"))
  x
}

setReplaceMethod("paramUpdates", "McmcParams", function(x, value){
  x@param_updates <- value
  x
})

setReplaceMethod("iter", "McmcParams", function(object, value){
  object@iter <- value
  object
})
