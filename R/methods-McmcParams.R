#' Create an object of class 'McmcParams' to specify iterations, burnin, etc.
#'
#' @examples
#'      mp <- McmcParams(iter=100, burnin=10)
#' @param iter number of iterations
#' @param burnin number of burnin iterations
#' @param thin thinning interval
#' @param nStarts number of chains to run
#' @param param_updates labeled vector specifying whether each parameter is to be updated (1) or not (0).
#' @return An object of class 'McmcParams'
#' @export
McmcParams <- function(iter=1000L, burnin=0L, thin, nStarts=1,
                       param_updates=.param_updates()){
  if(missing(thin)) thin <- rep(1L, length(iter))
  new("McmcParams", iter=as.integer(iter),
      burnin=as.integer(burnin),
      thin=as.integer(thin),
      nstarts=as.integer(nStarts),
      param_updates=param_updates)
}


#' @rdname burnin-method
#' @aliases burnin,McmcParams-method
setMethod("burnin", "McmcParams", function(object)  object@burnin)

#' @rdname burnin-method
#' @aliases burnin<-,McmcParams-method
setReplaceMethod("burnin", "McmcParams", function(object,value){
  object@burnin <- value
  object
})

#' @rdname thin-method
#' @aliases thin,McmcParams-method
setMethod("thin", "McmcParams", function(object) object@thin)

#' @rdname iter-method
#' @aliases iter,McmcParams-method
setMethod("iter", "McmcParams", function(object) object@iter)

setMethod("show", "McmcParams", function(object){
  cat("An object of class 'McmcParams'\n")
  cat("   iterations:", paste(iter(object), collapse=","), "\n")
  cat("   burnin    :", paste(burnin(object), collapse=","),  "\n")
  cat("   thin      :", paste(thin(object), collapse=","), "\n")
  cat("   n starts  :", nStarts(object), "\n")
})


#' extract iter, thin, and burnin for a specified instance of McmcParams
#'
#' Allows a user to pass a vector for burnin, thin, and iter.
#' @aliases [,McmcParams-method [,McmcParams,ANY-method
#' @return An object of class 'McmcParams'
#' @docType methods
#' @rdname extract-methods
setMethod("[", "McmcParams", function(x, i, j, ..., drop=FALSE){
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

#' @rdname nStarts-method
#' @aliases nStarts,McmcParams-method
setMethod("nStarts", "McmcParams", function(object) object@nstarts)

#' @rdname nStarts-method
#' @aliases nStarts<-,McmcParams-method
setReplaceMethod("nStarts", "McmcParams", function(object, value){
  object@nstarts <- as.integer(value)
  object
})

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

#' @rdname iter-method
#' @aliases iter<-,McmcParams-method
setReplaceMethod("iter", "McmcParams", function(object, force=FALSE, value){
  object@iter <- value
  object
})
