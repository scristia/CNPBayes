#' @include AllClasses.R

.parameterizeGammaByMeanSd <- function(mn, sd){
  rate <- mn/sd^2
  shape <- mn^2/sd^2
  setNames(c(shape, rate), c("shape", "rate"))
}


#' Quantiles, shape, and rate of the prior for the inverse of tau2 (the precision)
#'
#' The precision prior for tau2 in the hiearchical model is given by
#' gamma(shape, rate). The shape and rate are a function of the
#' hyperparameters eta.0 and m2.0.  Specifically, shape=1/2*eta.0 and
#' the rate=1/2*eta.0*m2.0.  Quantiles for this distribution and the
#' shape and rate can be obtained by specifying the hyperparameters
#' eta.0 and m2.0, or alternatively by specifying the desired mean and
#' standard deviation of the precisions.
#' @param eta.0 hyperparameter for precision
#' @param m2.0 hyperparameter for precision
#' @param mn mean of precision
#' @param sd standard deviation of precision
#' @return a list with elements 'quantiles', 'eta.0', 'm2.0', 'mean', and 'sd'
#'
#'
#' @examples
#' results <- qInverseTau2(mn=100, sd=1)
#' precision.quantiles <- results$quantiles
#' sd.quantiles <- sqrt(1/precision.quantiles)
#' results$mean
#' results$sd
#' results$eta.0
#' results$m2.0
#'
#' results2 <- qInverseTau2(eta.0=1800, m2.0=100)
#'
#' ## Find quantiles from the default set of hyperparameters
#' hypp <- Hyperparameters(type="batch")
#' results3 <- qInverseTau2(eta.0(hypp), m2.0(hypp))
#' default.precision.quantiles <- results3$quantiles
#' @export
qInverseTau2 <- function(eta.0=1800, m2.0=100, mn, sd){
  if(!missing(mn) && !missing(sd)){
    params <- .parameterizeGammaByMeanSd(mn, sd)
    a <- params[["shape"]]
    b <- params[["rate"]]
    eta.0 <- 2*a
    m2.0 <- b/a
  }
##  shape <- 0.5*eta.0
##  rate <- 0.5*eta.0*m2.0
##  mn <- shape/rate
##  sd <- sqrt(shape/rate^2)
  x <- qgamma(seq(0, 1-0.001, 0.001), shape=0.5*eta.0, rate=0.5*eta.0*m2.0)
  list(quantiles=x, eta.0=eta.0, m2.0=m2.0, mean=mn, sd=sd)
}

HyperparametersBatch <- function(k=0L,
                                 mu.0=0,
                                 tau2.0=100,
                                 eta.0=1800,
                                 m2.0=1/60,
                                 alpha,
                                 beta=0.1, ## mean is 1/10
                                 a=1.8,
                                 b=6){
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
                                    mu.0=0,
                                    tau2.0=100,
                                    eta.0=1,
                                    m2.0=0.1,
                                    alpha,
                                    beta=0.1, ## mean is 1/10
                                    a=1.8,
                                    b=6){
  if(missing(alpha)) alpha <- rep(1, k)
  ##if(missing(mu)) mu <- initializeMu(k)
  ##if(missing(tau2)) tau2 <- rep(1, k)
  new("HyperparametersMarginal",
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

setValidity("Hyperparameters", function(object){
  msg <- TRUE
  if(k(object) != alpha(object)){
    msg <- "alpha vector must be the same length as k"
    return(msg)
  }
  msg
})

#' Create an object of class 'Hyperparameters'
#'
#' @examples
#'      hypp <- Hyperparameters("marginal", k=2)
#' @param type specifies 'marginal' or 'batch'
#' @param k number of components
#' @param ... optional parameters, i.e. mu.0, tau2.0, eta.0, etc.
#' @export
Hyperparameters <- function(type="batch", k=2L, ...){
  if(type=="marginal") return(HyperparametersMarginal(k, ...))
  if(type=="batch") return(HyperparametersBatch(k, ...))
}

#' @rdname k-method
#' @aliases k<-,Hyperparemeters-method
setReplaceMethod("k", "Hyperparameters", function(object, value){
  object@k <- as.integer(value)
  object@alpha <- rep(1L, value)
  object
})

setValidity("Hyperparameters", function(object){
  msg <- NULL
  if(length(alpha(object)) != k(object)){
    msg <- "alpha should be a numeric vector with length equal to the number of components"
    return(msg)
  }
})

setMethod("k", "Hyperparameters", function(object) object@k)
setMethod("alpha", "Hyperparameters", function(object) object@alpha)

## beta is a base function
betas <- function(object) object@beta
a <- function(object) object@a
b <- function(object) object@b

setReplaceMethod("alpha", "Hyperparameters", function(object, value){
  object@alpha <- value
  object
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
  rnorm(k(object), mu.0(object), tau.0(object))
})
