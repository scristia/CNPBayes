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

setMethod("startingValues", "MarginalModel", function(object){
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  ##mu(object) <- rnorm(1, mu.0(object), tau.0(object))
  mu(object) <- mean(y(object))
  ##tau2(object) <- 1/rgamma(1, shape=1/2*eta.0(object),
  ##rate=1/2*eta.0(object)*m2.0(object))
  tau2(object) <- var(y(object))*10
  theta(object) <- rnorm(k(object), mu(object), tau(object))
  sigma2(object) <- 1/rgamma(k(object), shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  p(object) <- as.numeric(rdirichlet(1, alpha(hypp))) ## rows are
  ## platform, columns are components
  if(length(y(object)) > 0){
    zz <- as.integer(simulateZ(length(y(object)), p(object)))
    cnt <- 1
    while (length(table(zz)) < k(object)) {
      if (cnt > 10) {
        stop("Too few observations or too many components.")
      }
      p(object) <- as.numeric(rdirichlet(1, alpha(hypp))) ## rows are
      zz <- as.integer(simulateZ(length(y(object)), p(object)))
      cnt <- cnt + 1
    }
  } else zz <- integer()
  z(object) <- zz
  zFreq(object) <- as.integer(table(z(object)))
  dataPrec(object) <- 1/computeVars(object)
  dataMean(object) <- computeMeans(object)
  log_lik(object) <- computeLoglik(object)
  logPrior(object) <- computePrior(object)
  chains(object) <- McmcChains(object)
  object
})

setMethod("startingValues", "BatchModel", function(object){
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  mu(object) <- rnorm(k(object), mu.0(object), tau.0(object))
  tau2(object) <- 1/rgamma(k(object), shape=1/2*eta.0(object), rate=1/2*eta.0(object)*m2.0(object))
  nB <- nBatch(object)
  th <- initializeTheta(object)
  if(is(object, "UnivariateBatchModel") || ncol(th) == 1){
    theta(object) <- th
  } else {
    theta(object) <- t(apply(th, 1, sort))
  }
  sigma2(object) <- initializeSigma2(object)
  p(object) <- as.numeric(rdirichlet(1, alpha(hypp))) ## rows are platform, columns are components
  ##p(object) <- rdirichlet(nBatch(object), alpha(hypp))
  ##if(missing(zz)){
  if(length(y(object)) > 0){
    zz <- simulateZ(length(y(object)), p(object))
    z(object) <- as.integer(factor(zz))

  } else zz <- as.integer(factor(levels=seq_len(k(hypp))))
  zFreq(object) <- as.integer(table(zz))
  if(length(y(object)) > 0){
    dataMean(object) <- computeMeans(object)
    dataPrec(object) <- computePrec(object)
    log_lik(object) <- computeLoglik(object)
    logPrior(object) <- computePrior(object)
  }
  probz(object) <- .computeProbZ(object)
  chains(object) <- McmcChains(object)
  object
})
