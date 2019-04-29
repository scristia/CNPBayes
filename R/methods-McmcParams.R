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
McmcParams <- function(iter=1000L,
                       burnin=0L,
                       thin=1L,
                       nStarts=1L,
                       param_updates=.param_updates(),
                       min_GR=1.2,
                       min_effsize=round(1/3*iter, 0),
                       max_burnin=32000,
                       min_chains=1){
  if(missing(thin)) thin <- rep(1L, length(iter))
  new("McmcParams", iter=as.integer(iter),
      burnin=as.integer(burnin),
      thin=as.integer(thin),
      nstarts=as.integer(nStarts),
      param_updates=param_updates,
      min_GR=min_GR,
      min_effsize=min_effsize,
      max_burnin=max_burnin,
      min_chains=min_chains)
}


#' @rdname burnin-method
#' @aliases burnin,McmcParams-method
setMethod("burnin", "McmcParams", function(object)  object@burnin)


setMethod("min_GR", "McmcParams", function(object)  object@min_GR)


setMethod("min_chains", "McmcParams", function(object)  object@min_chains)

setMethod("max_burnin", "McmcParams", function(object)  object@max_burnin)

setReplaceMethod("max_burnin", "McmcParams", function(object, value){
  object@max_burnin <- as.integer(value)
  object
})

setMethod("min_effsize", "McmcParams", function(object)  object@min_effsize)

#' @rdname burnin-method
#' @aliases burnin<-,McmcParams-method
setReplaceMethod("burnin", "McmcParams", function(object,value){
  object@burnin <- value
  object
})

#' @rdname thin-method
#' @aliases thin,McmcParams-method
setMethod("thin", "McmcParams", function(object) object@thin)

setReplaceMethod("thin", c("McmcParams", "numeric"), function(object, value){
  object@thin <- value
  object
})

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
##  if(nStarts(object) < min_chains(object)){
##    msg <- "number of independent starts is less than the mininum number required for assessing convergence"
##  }
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
setReplaceMethod("iter", "McmcParams", function(object, value){
  object@iter <- value
  object
})

setMethod("updateObject", "McmcParams", function(object){
  obj <- callNextMethod(object)
  obj@max_burnin <- 32000L
  obj@min_chains <- 1L
  obj@min_effsize <- round(1/3 * iter(object), 0)
  obj@min_GR <- 1.2
  obj
})
