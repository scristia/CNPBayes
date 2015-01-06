#' @include AllClasses.R
HyperparametersBatch <- function(k=0L,
                                 mu.0=0,
                                 tau2.0=1000,
                                 eta.0,
                                 m2.0,
                                 alpha,
                                 beta=0.1, ## mean is 1/10
                                 a=1.8,
                                 b=6){
  if(missing(eta.0) || m2.0){
    ab <- gammaShapeRate2(1, sqrt(1000))
    eta.0 <- ab["a"]
    m2.0 <- ab["b"]
  }
  if(missing(alpha)) alpha <- rep(1, k)
  new("HyperparametersBatch",
      k=as.integer(k),
      mu.0=mu.0,
      tau2.0=tau2.0,
      eta.0=eta.0,
      m2.0=m2.0,
      alpha=alpha,
      beta=beta,
      a=a,
      b=b)
}

HyperparametersMarginal <- function(k=0L,
                                    mu,
                                    tau2,
                                    eta.0,
                                    m2.0,
                                    alpha,
                                    beta=0.1, ## mean is 1/10
                                    a=1.8,
                                    b=6){
  if(missing(eta.0) || m2.0){
    ab <- gammaShapeRate2(1, sqrt(1000))
    eta.0 <- ab["a"]
    m2.0 <- ab["b"]
  }
  if(missing(alpha)) alpha <- rep(1, k)
  if(missing(mu)) mu <- initializeMu(k)
  if(missing(tau2)) tau2 <- rep(1, k)
  new("HyperparametersMarginal",
      k=as.integer(k),
      mu=mu,
      tau2=tau2,
      eta.0=eta.0,
      m2.0=m2.0,
      alpha=alpha,
      beta=beta,
      a=a,
      b=b)
}

Hyperparameters <- function(type, k){
  if(type=="marginal") return(HyperparametersMarginal(k))
  if(type=="batch") return(HyperparametersBatch(k))
}

setValidity("Hyperparameters", function(object){
  msg <- NULL
  if(length(alpha(object)) != k(object)){
    msg <- "alpha should be a numeric vector with length equal to the number of components"
    return(msg)
  }
})


##setMethod("params", "Hyperparameters", function(object) object@params)
setMethod("k", "Hyperparameters", function(object) object@k)

#' @export
mu.0 <- function(object) object@mu.0

#' @export
setMethod("mu", "HyperparametersMarginal", function(object) object@mu)

#' @export
setMethod("tau2", "HyperparametersMarginal", function(object) object@tau2)

#' @export
tau2.0 <- function(object) object@tau2.0

#' @export
eta.0 <- function(object) object@eta.0

#' @export
m2.0 <- function(object) object@m2.0

#' @export
alpha <- function(object) object@alpha

## beta is a base function
#' @export
betas <- function(object) object@beta

#' @export
a <- function(object) object@a

#' @export
b <- function(object) object@b

#' @export
setReplaceMethod("alpha", "Hyperparameters", function(object, value){
  object@alpha <- value
  object
})

setMethod("show", "HyperparametersMarginal", function(object){
  cat("An object of class 'Hyperparameters'\n")
  cat("   k      :", k(object), "\n")
  cat("   mu   :", mu(object), "\n")
  cat("   tau2 :", tau2(object), "\n")
  cat("   eta.0  :", eta.0(object), "\n")
  cat("   m2.0   :", round(m2.0(object), 3), "\n")
  cat("   alpha  :", alpha(object), "\n")
  cat("   beta   :", betas(object), "\n")
  cat("   a      :", a(object), "\n")
  cat("   b      :", b(object), "\n")
})

setMethod("show", "Hyperparameters", function(object){
  cat("An object of class 'Hyperparameters'\n")
  cat("   k      :", k(object), "\n")
  cat("   mu.0   :", mu.0(object), "\n")
  cat("   tau2.0 :", tau2.0(object), "\n")
  cat("   eta.0  :", eta.0(object), "\n")
  cat("   m2.0   :", round(m2.0(object), 3), "\n")
  cat("   alpha  :", alpha(object), "\n")
  cat("   beta   :", betas(object), "\n")
  cat("   a      :", a(object), "\n")
  cat("   b      :", b(object), "\n")
})

setMethod("initializeMu", "numeric", function(object){
  means <- switch(paste0("k", object),
                  k1=0,
                  k2=c(-0.5, 0),
                  k3=c(-2, -0.5, 0),
                  k4=c(-2, -0.5, 0, 0.5),
                  k5=c(-2, -0.5, 0, 0.5, 1),
                  k6=c(-2, -0.5, -0.2, 0.2, 0.5, 1),
                  k7=c(-2, -0.5, -0.2, 0, 0.2, 0.5, 1),
                  NULL)
  if(is.null(means)) stop("k needs to be 1-7")
  means
})

##setMethod("initializeTau2", "numeric", function(object){
##  means <- switch(paste0("k", object),
##                  k1=1,
##                  k2=c(-0.5, 0),
##                  k3=c(-2, -0.5, 0),
##                  k4=c(-2, -0.5, 0, 0.5),
##                  k5=c(-2, -0.5, 0, 0.5, 1),
##                  k6=c(-2, -0.5, -0.2, 0.2, 0.5, 1),
##                  k7=c(-2, -0.5, -0.2, 0, 0.2, 0.5, 1),
##                  NULL)
##  if(is.null(means)) stop("k needs to be 1-7")
##  means
##})
