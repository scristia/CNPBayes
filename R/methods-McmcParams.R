#' @export
McmcParams <- function(iter=1000L, burnin=0L, thin=1L, constrainTheta=TRUE){
  new("McmcParams", iter=iter, burnin=burnin, thin=thin, constrainTheta=constrainTheta)
}

burnin <- function(object) object@burnin
thin <- function(object) object@thin
iter <- function(object) object@iter
constrainTheta <- function(object) object@constrainTheta
savedIterations <- function(object)iter(object)/thin(object)


setMethod("show", "McmcParams", function(object){
  cat("An object of class 'McmcParams'\n")
  cat("   iterations:", iter(object), "\n")
  cat("   burnin    :", burnin(object), "\n")
  cat("   thin      :", thin(object), "\n")
})
