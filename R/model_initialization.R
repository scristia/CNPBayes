constructModel <- function(type, data, k, batch, hypp){
  nbatch <- length(unique(batch))
  if(k > 1){
    if(type=="batch"){
      object <- BatchModel(data, k, batch, hypp)
    } else{
      object <- MarginalModel(data, k, batch, hypp)
    }
  } else {
    if(type=="batch"){
      object <- UnivariateBatchModel(data, 1, batch, hypp)
    } else  {
      object <- UnivariateMarginalModel(data, 1 , batch, hypp)
    }
  }
  object
}

setMethod("initializeTheta", "BatchModel", function(object){
  th <- matrix(NA, nBatch(object), k(object))
  for(j in seq_len(k(object))){
    th[, j] <- rnorm(nBatch(object), mu(object)[j], tau(object)[j])
  }
  th
})

setMethod("initializeSigma2", "BatchModel", function(object){
  s2 <- 1/rgamma(nBatch(object)*k(object), shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  s2 <- matrix(s2, nBatch(object), k(object))
  s2
})

setMethod("initializeModel", "ModelParams", function(params, hypp){
  object <- constructModel(type(params),
                           data=y(params),
                           k=k(params),
                           batch=batch(params),
                           hypp=hypp)
  object <- startingValues(object, params)
  object
})

#' @export
initializeBatchModel <- function(params, zz, hypp){
  object <- constructModel(type(params), data=y(params), k=k(params),
                           batch=batch(params),
                           hypp=hypp)
  startingValues(object, params, zz)
}


setMethod("startingValues", "MarginalModel", function(object, params){
  ##hypp <- hyperParams(object) <- Hyperparameters(type(params),
  ##k=k(params))
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  mu(object) <- rnorm(1, mu.0(object), sigma.0(object))
  tau2(object) <- 1/rgamma(1, shape=1/2*eta.0(object), rate=1/2*eta.0(object)*m2.0(object))
  theta(object) <- rnorm(k(object), mu(object), tau(object))
  sigma2(object) <- 1/rgamma(k(object), shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  p(object) <- as.numeric(rdirichlet(1, alpha(hypp))) ## rows are platform, columns are components
  zz <- simulateZ(length(y(object)), p(object))
  z(object) <- factor(zz, levels=unique(sort(zz)))
  dataPrec(object) <- 1/computeVars(object)
  dataMean(object) <- computeMeans(object)
  logpotential(object) <- computePotential(object)
  logLik(object) <- computeLoglik(object)
  mcmcChains(object) <- McmcChains(object, mcmcParams(params))
  object
})

setMethod("startingValues", "BatchModel", function(object, params, zz){
  ##  object <- constructModel(type(params), data=y(params), k=k(params),
  ##                           batch=batch(params))
  ##  hypp <- hyperParams(object) <- Hyperparameters(type(params), k=k(params))
  ##  hyperParams(object) <- hypp
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  mu(object) <- rnorm(k(object), mu.0(object), sigma.0(object))
  tau2(object) <- 1/rgamma(k(object), shape=1/2*eta.0(object), rate=1/2*eta.0(object)*m2.0(object))
  nB <- nBatch(object)
  theta(object) <- initializeTheta(object)
  sigma2(object) <- initializeSigma2(object)
  p(object) <- as.numeric(rdirichlet(1, alpha(hypp))) ## rows are platform, columns are components
  ##p(object) <- rdirichlet(nBatch(object), alpha(hypp))
  if(missing(zz)){
    zz <- simulateZ(length(y(params)), p(object))
    z(object) <- factor(zz, levels=unique(sort(zz)))
  } else z(object) <- zz
  dataPrec(object) <- 1/computeVars(object)
  dataMean(object) <- computeMeans(object)
  logpotential(object) <- computePotential(object)
  mcmcChains(object) <- McmcChains(object, mcmcParams(params))
  object
})

setMethod("initializeTheta", "numeric", function(object){
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
