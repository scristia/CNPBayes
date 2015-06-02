#' @include AllClasses.R

gammaShapeRate2 <- function(mn, sd){
  ##1/2*a   = shape
  ##1/2*a*b = rate
  ##shape/rate   = mn (1)
  ##shape/rate^2 = sd (2)
  ##shape = mn * rate
  ##mn * rate / rate^2 = sd
  ##mn/rate = sd
  ##mn/sd = rate
  rate <- mn/sd
  shape <- mn*rate
  ##a=shape*2
  ##1/2*(shape*2)*b = rate
  ## shape*b=rate
  ## b=rate/shape
  a <- shape*2
  b <- rate/shape
  setNames(c(a, b), c("a", "b"))
}

gammaShapeRate3 <- function(eta.0, m2.0){
  eta.0 <- 1800/100
  ##eta.0 <- 180
  m2.0 <- 1/(60/100)
  a <- 0.5*eta.0
  b <- 0.5*eta.0*m2.0
  x <- rgamma(1000, a, rate=1/b)
  ##  ##shape/rate   = mn (1)
  ##  ##shape/rate^2 = sd (2)
  ##  ##shape = mn * rate
  ##  ##mn * rate / rate^2 = sd
  ##  ##mn/rate = sd
  ##  ##mn/sd = rate
  ##  rate <- mn/sd
  ##  shape <- mn*rate
  ##  setNames(c(shape, rate), c("shape", "rate"))
  x
}

#' @export
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

#' @export
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

##setGeneric("Hyperparameters", function(object, k, ...) standardGeneric("Hyperparameters"))
##
##setMethod("Hyperparameters", "character", function(object, k, ...){
##  if(object=="marginal") return(HyperparametersMarginal(k, ...))
##  if(object=="batch") return(HyperparametersBatch(k, ...))
##})
##
##setMethod("Hyperparameters", "Hyperparameters", function(object, k, ...){
##  if(object=="marginal") return(HyperparametersMarginal(k, ...))
##  if(object=="batch") return(HyperparametersBatch(k, ...))
##})

#' @export
Hyperparameters <- function(type="batch", k=2L, ...){
  if(type=="marginal") return(HyperparametersMarginal(k, ...))
  if(type=="batch") return(HyperparametersBatch(k, ...))
}

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


##setMethod("params", "Hyperparameters", function(object) object@params)
setMethod("k", "Hyperparameters", function(object) object@k)


# #
# setMethod("mu", "HyperparametersMarginal", function(object) object@mu)

# #
# setMethod("tau2", "HyperparametersMarginal", function(object) object@tau2)


setMethod("alpha", "Hyperparameters", function(object) object@alpha)

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

##setMethod("show", "HyperparametersMarginal", function(object){
##  cat("An object of class 'Hyperparameters'\n")
##  cat("   k      :", k(object), "\n")
##  ##cat("   mu   :", mu(object), "\n")
##  ##cat("   tau2 :", tau2(object), "\n")
##  cat("   eta.0  :", eta.0(object), "\n")
##  cat("   m2.0   :", round(m2.0(object), 3), "\n")
##  cat("   alpha  :", alpha(object), "\n")
##  cat("   beta   :", betas(object), "\n")
##  cat("   a      :", a(object), "\n")
##  cat("   b      :", b(object), "\n")
##})

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
##  means <- switch(paste0("k", object),
##                  k1=0,
##                  k2=c(-0.5, 0),
##                  k3=c(-2, -0.5, 0),
##                  k4=c(-2, -0.5, 0, 0.5),
##                  k5=c(-2, -0.5, 0, 0.5, 1),
##                  k6=c(-2, -0.5, -0.2, 0.2, 0.5, 1),
##                  k7=c(-2, -0.5, -0.2, 0, 0.2, 0.5, 1),
##                  NULL)
##  if(is.null(means)) stop("k needs to be 1-7")
##  means
})

setMethod("HyperParameterList", "HyperparametersMarginal", function(hypp, K){
  hyplist <- foreach(k = K) %do% {
    HyperparametersMarginal(k=k, mu.0=mu.0(hypp),
                            tau2.0=tau2.0(hypp),
                            eta.0=eta.0(hypp),
                            m2.0=m2.0(hypp),
                            beta=betas(hypp),
                            a=a(hypp),
                            b=b(hypp))
  }
  hyplist
})

setMethod("HyperParameterList", "HyperparametersBatch", function(hypp, K){
  hyplist <- foreach(k = K) %do% {
    HyperparametersBatch(k=k, mu.0=mu.0(hypp),
                         tau2.0=tau2.0(hypp),
                         eta.0=eta.0(hypp),
                         m2.0=m2.0(hypp),
                         beta=betas(hypp),
                         a=a(hypp),
                         b=b(hypp))
  }
  hyplist
})

setMethod("ModelParamList", "HyperparametersMarginal", function(hypp, K, data, mcmcp, batch){
  plist <- foreach(k=K) %do% {
    ModelParams("marginal", y=data, k=k,
                mcmc.params=mcmcp[k])
  }
  plist
})

setMethod("ModelParamList", "HyperparametersBatch", function(hypp, K, data, mcmcp, batch){
  plist <- foreach(k=K) %do% {
    ModelParams("batch", y=data, k=k,
                batch=batch,
                mcmc.params=mcmcp[k])
  }
  plist
})
