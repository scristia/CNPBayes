#' @include AllClasses.R

.parameterizeGammaByMeanSd <- function(mn, sd){
  rate <- mn/sd^2
  shape <- mn^2/sd^2
  setNames(c(shape, rate), c("shape", "rate"))
}


qInverseTau2 <- function(eta.0=1800, m2.0=100, mn, sd){
  if(!missing(mn) && !missing(sd)){
    params <- .parameterizeGammaByMeanSd(mn, sd)
    a <- params[["shape"]]
    b <- params[["rate"]]
    eta.0 <- 2*a
    m2.0 <- b/a
  }
  shape <- 0.5*eta.0
  rate <- 0.5*eta.0*m2.0
  mn <- shape/rate
  sd <- sqrt(shape/rate^2)
  x <- qgamma(seq(0, 1-0.001, 0.001), shape=0.5*eta.0, rate=0.5*eta.0*m2.0)
  x <- x[is.finite(x) & x > 0]
  list(quantiles=x, eta.0=eta.0, m2.0=m2.0, mean=mn, sd=sd)
}


#' Create an object of class 'HyperparametersMultiBatch' for the
#' batch mixture model
#'
#' @param k  length-one integer vector specifying number of components
#' (typically 1 <= k <= 4)
#' @param mu.0 length-one numeric vector of the of the normal prior
#' for the component means.
#' @param tau2.0 length-one numeric vector of the variance for the normal
#' prior of the component means
#' @param eta.0 length-one numeric vector of the shape parameter for
#' the Inverse Gamma prior of the component variances, tau2_h.  The
#' shape parameter is parameterized as 1/2 * eta.0.  In the batch
#' model, tau2_h describes the inter-batch heterogeneity of means for
#' component h.
#' @param m2.0 length-one numeric vector of the rate parameter for the
#' Inverse Gamma prior of the component variances, tau2_h.  The rate
#' parameter is parameterized as 1/2 * eta.0 * m2.0.  In the batch
#' model, tau2_h describes the inter-batch heterogeneity of means for
#' component h.
#' @param alpha length-k numeric vector of the shape parameters for
#' the dirichlet prior on the mixture probabilities
#' @param beta length-one numeric vector for the parameter of the
#' geometric prior for nu.0 (nu.0 is the shape parameter of the
#' Inverse Gamma sampling distribution for the component-specific
#' variances. Together, nu.0 and sigma2.0 model inter-component
#' heterogeneity in variances.).  beta is a probability and must be
#' in the interval [0,1].
#' @param a length-one numeric vector of the shape parameter for the
#' Gamma prior used for sigma2.0 (sigma2.0 is the shape parameter of
#' the Inverse Gamma sampling distribution for the component-specific
#' variances).
#' @param b a length-one numeric vector of the rate parameter for the
#' Gamma prior used for sigma2.0 (sigma2.0 is the rate parameter of
#' the Inverse Gamma sampling distribution for the component-specific
#' variances)
#' @param dfr length-one numeric vector for t-distribution degrees of freedom
#' @return An object of class HyperparametersBatch
HyperparametersMultiBatch <- function(k=3L,
                                      mu.0=0,
                                      tau2.0=0.4,
                                      eta.0=32,
                                      m2.0=0.5,
                                      alpha,
                                      beta=0.1, ## mean is 1/10
                                      a=1.8,
                                      b=6,
                                      dfr=100){
  if(missing(alpha)) alpha <- rep(1, k)
  new("HyperparametersMultiBatch",
      k=as.integer(k),
      mu.0=mu.0,
      tau2.0=tau2.0,
      eta.0=eta.0,
      m2.0=m2.0,
      alpha=alpha,
      beta=beta,
      a=a,
      b=b,
      dfr=dfr)
}


#' Create an object of class 'HyperparametersSingleBatch' for the
#' single batch mixture model
#'
#' @param k  length-one integer vector specifying number of components
#' (typically 1 <= k <= 4)
#' @param mu.0  length-one numeric vector of the mean for the normal
#' prior of the component means
#' @param tau2.0 length-one numeric vector of the variance for the normal
#' prior of the component means
#' @param eta.0 length-one numeric vector of the shape parameter for
#' the Inverse Gamma prior of the component variances.  The shape
#' parameter is parameterized as 1/2 * eta.0.
#' @param m2.0 length-one numeric vector of the rate parameter for
#' the Inverse Gamma prior of the component variances.  The rate
#' parameter is parameterized as 1/2 * eta.0 * m2.0.
#' @param alpha length-k numeric vector of the shape parameters for
#' the dirichlet prior on the mixture probabilities
#' @param beta length-one numeric vector for the parameter of the
#' geometric prior for nu.0 (nu.0 is the shape parameter of the
#' Inverse Gamma sampling distribution for the component-specific
#' variances).  beta is a probability and must be in the interval
#' [0,1].
#' @param a length-one numeric vector of the shape parameter for the
#' Gamma prior used for sigma2.0 (sigma2.0 is the shape parameter of
#' the Inverse Gamma sampling distribution for the component-specific
#' variances)
#' @param b a length-one numeric vector of the rate parameter for the
#' Gamma prior used for sigma2.0 (sigma2.0 is the rate parameter of
#' the Inverse Gamma sampling distribution for the component-specific
#' variances)
#' @param dfr length-one numeric vector for t-distribution degrees of freedom
#'
#' @return An object of class HyperparametersSingleBatch
HyperparametersSingleBatch <- function(k=0L,
                                    mu.0=0,
                                    tau2.0=0.4,
                                    eta.0=32,
                                    m2.0=0.5,
                                    alpha,
                                    beta=0.1, ## mean is 1/10
                                    a=1.8,
                                    b=6,
                                    dfr=100){
  if(missing(alpha)) alpha <- rep(1, k)
  ##if(missing(tau2)) tau2 <- rep(1, k)
  new("HyperparametersSingleBatch",
      k=as.integer(k),
      mu.0=mu.0,
      tau2.0=tau2.0,
      eta.0=eta.0,
      m2.0=m2.0,
      alpha=alpha,
      beta=beta,
      a=a,
      b=b,
      dfr=dfr)
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
#' @param type specifies 'marginal' or 'batch'
#' @param k number of components
#' @param ... additional arguments are names of the slots of the Hyperparameters class, including `mu.0`, `tau2.0`, `eta.0`, `m2.0`, `alpha`, `beta`, `a`, `b`, and `dfr`).
#' @return An object of class HyperparametersMarginal or HyperparametersBatch
#'
#' @details
#' Additional hyperparameters can be passed to the HyperparametersMarginal and HyperparametersBatch models.
#'
#' @rdname Hyperparameters
#' @export
Hyperparameters <- function(type="batch", k=2L, ...){
  if(type=="marginal") return(HyperparametersSingleBatch(k, ...))
  if(type=="batch") return(HyperparametersMultiBatch(k, ...))
  ##if(type=="trios") return(HyperparametersTrios(k, ...))
}

#' Accessor for number of mixture components
#' 
#' @rdname k-method
#' @aliases k<-,Hyperparameters-method
#' @param object A Hyperparameters instance
#' @param value An integer specifying the number of mixture components for the mixture model
setReplaceMethod("k", "Hyperparameters", function(object, value){
  object@k <- as.integer(value)
  object@alpha <- rep(1L, value)
  object
})

## #' @rdname k-method
## #' @aliases k<-,HyperparametersTrios-method
## setReplaceMethod("k", "HyperparametersTrios", function(object, value){
##   object@k <- as.integer(value)
##   object@alpha <- rep(1L, value)
##   object
## })

setValidity("Hyperparameters", function(object){
  msg <- NULL
  if(length(alpha(object)) != k(object)){
    msg <- "alpha should be a numeric vector with length equal to the number of components"
    return(msg)
  }
})

#' @rdname k-method
#' @export
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

# Create a list of hyperparameter objects for each of the four mixture model implmentations
hpList <- function(...){
  sbp <- sb <- Hyperparameters(...)
  mb <- mbp <- HyperparametersMultiBatch(...)
  ##tbm <- HyperparametersTrios(...)
  list(SB=sb,
       MB=mb,
       SBP=sbp,
       MBP=mbp)
}

setMethod("dfr", "Hyperparameters", function(object) object@dfr)

setReplaceMethod("dfr", c("Hyperparameters", "numeric"),
                 function(object, value){
                   object@dfr <- value
                   object
                 })
