BatchModel <- function(data, k, batch){
  mcmc.chains <- McmcChains()
  B <- length(unique(batch))
  new("BatchModel",
      hyperparams=HyperparametersBatch(k=k),
      theta=matrix(NA, B, k),
      sigma2=matrix(NA, B, k),
      mu=numeric(k),
      tau2=numeric(k),
      nu.0=numeric(1),
      sigma2.0=numeric(1),
      pi=numeric(k),
      data=data,
      data.mean=matrix(NA, B, k),
      data.prec=matrix(NA, B, k),
      z=factor(numeric(length(data))),
      probz=matrix(0, length(data), k),
      logpotential=numeric(1),
      mcmc.chains=mcmc.chains,
      batch=batch,
      hwe=numeric(),
      theta_order=numeric(B*k))
}

##BatchModelNoHom <- function(data, k, batch){
##  mcmc.chains <- McmcChains()
##  B <- length(unique(batch))
##  new("BatchModelNoHom",
##      hyperparams=HyperparametersBatch(k=k),
##      theta=matrix(NA, B, k),
##      sigma2=matrix(NA, B, k),
##      mu=numeric(k),
##      tau2=numeric(k),
##      nu.0=numeric(1),
##      sigma2.0=numeric(1),
##      pi=numeric(k),
##      data=data,
##      data.mean=matrix(NA, B, k),
##      data.prec=matrix(NA, B, k),
##      z=factor(numeric(length(data))),
##      probz=matrix(0, length(data), k),
##      logpotential=numeric(1),
##      mcmc.chains=mcmc.chains,
##      batch=batch,
##      hwe=numeric())
##}
##
##BatchModelPlusHom <- function(data, k, batch){
##  mcmc.chains <- McmcChains()
##  B <- length(unique(batch))
##  new("BatchModelPlusHom",
##      hyperparams=HyperparametersBatch(k=k),
##      theta=matrix(NA, B, k),
##      sigma2=matrix(NA, B, k),
##      mu=numeric(k),
##      tau2=numeric(k),
##      nu.0=numeric(1),
##      sigma2.0=numeric(1),
##      pi=numeric(k),
##      data=data,
##      data.mean=matrix(NA, B, k),
##      data.prec=matrix(NA, B, k),
##      z=factor(numeric(length(data))),
##      probz=matrix(0, length(data), k),
##      logpotential=numeric(1),
##      mcmc.chains=mcmc.chains,
##      batch=batch,
##      hwe=numeric())
##}

##
## Multiple batches, but only 1 component
##
UnivariateBatchModel <- function(data, k=1, batch){
  mcmc.chains <- McmcChains()
  B <- length(unique(batch))
  new("UnivariateBatchModel",
      hyperparams=HyperparametersBatch(k=1),
      theta=matrix(NA, B, 1),
      sigma2=matrix(NA, B, 1),
      mu=numeric(k),
      tau2=numeric(k),
      nu.0=numeric(1),
      sigma2.0=numeric(1),
      ##pi=numeric(k),
      pi=matrix(NA, B, 1),
      data=data,
      data.mean=matrix(NA, B, 1),
      data.prec=matrix(NA, B, 1),
      z=factor(numeric(length(data))),
      probz=matrix(1, length(data), 1),
      logpotential=numeric(1),
      mcmc.chains=mcmc.chains,
      batch=batch,
      hwe=numeric())
}
##
##setMethod("startingValues", "BatchModel", function(object){
##  ##
##  ## Initialize all batches to the same starting values by ignoring
##  ## batch
##  ##
##  ##browser()
##  hypp <- hyperParams(object)
##  if(FALSE){
##    tmp.file <- tempfile()
##    sink(tmp.file)
##    mmfit <- normalmixEM(y(object), arbvar = FALSE, epsilon = 1e-03, k=k(hypp), maxit=2000)
##    sink()
##    unlink(tmp.file)
##    mus <- mmfit$mu
##    vars <- (mmfit$sigma[order(mmfit$mu)])^2
##    B <- nBatch(object)
##    mus <- matrix(mus, B, k(object), byrow=TRUE)
##    vars <- matrix(vars, B, k(object), byrow=TRUE)
##    rownames(vars) <- rownames(mus) <- uniqueBatch(object)
##    theta(object) <- mus
##    sigma2(object) <- vars
##  }
##  if(TRUE){
##    if(all(is.na(theta(object)))){
##      theta(object) <- initializeTheta(object)
##    }
##    hypp <- hyperParams(object)
##    B <- nBatch(object)
##    sigma2(object) <- matrix(mad(as.numeric(y(object)), na.rm=TRUE)^2, B, k(object))
##    rownames(theta(object)) <- rownames(sigma2(object)) <- uniqueBatch(object)
##  }
##  object
##})

setMethod("initializeSigma2.0", "BatchModel", function(object){
  hypp <- hyperParams(object)
  sum(alpha(hypp)*colMeans(sigma2(object)))/sum(alpha(hypp))
})


setMethod("posteriorMultinomial", "BatchModel", function(object){
  .multinomial_probs <- .multBatch(y(object),
                                   theta(object),
                                   sqrt(sigma2(object)),
                                   p(object))
})


##
## z has length y.  Each observation is a sample.
##
setMethod("updateZ", "BatchModel", function(object){
  ##plist <- posteriorMultinomial(object)
  P <- posteriorMultinomial(object)
  zz <- .updateZ(P)
  ##zz <- simulateZ(length(y(object)), p)
  ##  zz <- rep(NA, length(y(object)))
  ##  ub <- uniqueBatch(object)
  ##  for(b in seq_along(plist)){
  ##    zz[batch(object)==ub[b]] <- .updateZ(plist[[b]])
  ##  }
  factor(zz, levels=seq_len(k(object)))
})


##setMethod("posteriorMultinomial", "BatchModel", function(object){
##  B <- nBatch(object)
##  plist <- vector("list", B)
##  ub <- uniqueBatch(object)
##  for(b in seq_len(B)){
##    plist[[b]] <- .multBatchSpecific(y(object)[batch(object)==ub[b]],
##                                     theta(object)[b, ],
##                                     sigma(object)[b, ],
##                                     p(object)[b, ])
##  }
##  return(plist)
##})


setMethod("posteriorMultinomial", "UnivariateBatchModel", function(object) return(1))

#' @export
uniqueBatch <- function(object) unique(batch(object))

nBatch <- function(object) length(uniqueBatch(object))

.multBatch <- function(y, theta, sd, pi){
  K <- seq_len(length(pi))
  B <- nrow(theta)
  tmp <- matrix(NA, length(y), B)
  numerator <- list()
  for(j in K){
    for(b in seq_len(B)){
      tmp[, b] <- pi[j]*dnorm(y, theta[b, j], sd[b, j])
    }
    numerator[[j]] <- rowSums(tmp)
  }
  numerator <- do.call(cbind, numerator)
  denominator <- rowSums(numerator)
  mix.probs <- numerator/denominator
  ## mix.probs has dimension N x K
  mix.probs
}

## y: length n_b vector  (number of samples in batch b)
## theta: length K vector for batch b
## sd:  length K vector for batch b
## pi: length K vector for batch b.
## Returns:  n_b x K matrix for batch b
.multBatchSpecific <- function(y, theta, sd, pi){
  K <- seq_len(length(pi))
  result <- matrix(NA, length(y), length(theta))
  for(j in K){
    result[, j] <- pi[j]*dnorm(y, theta[j], sd[j])
  }
  mix.probs <- result/rowSums(result)
  mix.probs
}

setMethod("initializeSigma2.0", "BatchModel", function(object){
  hypp <- hyperParams(object)
  ##
  ## In the marginal model, sigma2.0 is initialized as a weighted
  ## average (weights are number of observations in each component)
  ##
  ## In the batch model, average the variances across batch and then
  ## average across components.
  sigma2.bar <- colMeans(sigma2(object))
  sum(alpha(hypp)*sigma2.bar)/sum(alpha(hypp))
})

setMethod("computeMeans", "BatchModel", function(object){
  ubatch <- uniqueBatch(object)
  B <- factor(batch(object), levels=ubatch)
  ybatch <- split(y(object), B)
  zbatch <- split(z(object), B)
  ymeans <- foreach(y=ybatch, z=zbatch, .combine='rbind') %do%{
    mns <- sapply(split(y, z), mean, na.rm=TRUE)
    mns
  }
  rownames(ymeans) <- ubatch
  ymeans
})

setMethod("computeVars", "BatchModel", function(object){
  ubatch <- uniqueBatch(object)
  B <- factor(batch(object), levels=ubatch)
  ybatch <- split(y(object), B)
  zbatch <- split(z(object), B)
  s2.0 <- sigma2.0(object)
  yvars <- foreach(y=ybatch, z=zbatch, .combine='rbind') %do%{
    v <- sapply(split(y, z), var)
    v
  }
  rownames(yvars) <- ubatch
  yvars
})


setMethod("computePotential", "BatchModel", function(object){
  hypp <- hyperParams(object)
  K <- k(hypp)
  yy <- y(object)
  zz <- z(object)
  thetas <- theta(object)
  sigma2s <- sigma2(object)
  mus <- mu(object)
  tau2s <- tau2(object)
  pp <- p(object)
  p.mu <- dnorm(mus, mu.0(hypp), sqrt(tau2.0(hypp)))
  p.sigma2.0 <- dgamma(sigma2.0(object), a(hypp), b(hypp))
  p.nu.0 <- dgeom(as.integer(nu.0(object)), betas(hypp))
  p.theta <- matrix(NA, nrow(thetas), ncol(thetas))
  for(j in seq_len(ncol(p.theta))){
    p.theta[, j] <- dnorm(thetas[, j], mus[j], sqrt(tau2s[j]))
  }
  p.sigma2 <- dgamma(1/sigma2s, 1/2*nu.0(object), 1/2*nu.0(object)*sigma2.0(object))
  pot <- list()
  ##
  ##
  ##
  sigmas <- sqrt(sigma2s)
  p.y <- rep(NA, length(yy))
  ylist <- split(yy, batch(object))
  ##rownames(thetas) <- rownames(sigmas) <- uniqueBatch(object)
  for(b in uniqueBatch(object)){
    dat <- ylist[[b]]
    pot <- matrix(NA, length(dat), K)
    for(i in seq_len(K)){
      pot[, i] <- pp[i]*dnorm(dat, thetas[b, i], sigmas[b, i])
    }
    p.y[batch(object)==b] <- rowSums(pot)
  }
  total_pot <- sum(log(p.y)) + sum(log(p.theta)) + sum(log(p.sigma2)) + sum(log(p.mu)) + log(p.nu.0) + log(p.sigma2.0)
  total_pot
})


setMethod("show", "BatchModel", function(object){
  ##callNextMethod()
  cat("An object of class 'BatchModel'\n")
  cat("     n. obs      :", length(y(object)), "\n")
  cat("     n. batches  :", nBatch(object), "\n")
  cat("     k           :", k(object), "\n")
  cat("     nobs/batch  :", table(batch(object)), "\n")
})


setMethod("batchCorrect", "BatchModel", function(object){
  B <- factor(batch(object), levels=unique(batch(object)))
  yy <- dat(object)
  names(yy) <- seq_along(yy)
  ybatch <- split(yy, B)
  zbatch <- split(z(object), B)
  thetas <- theta(object)
  sigmas <- sqrt(sigma2(object))
  i <- NULL
  ystar <- foreach(i=seq_along(ybatch), .combine='c') %do% {
    yy <- ybatch[[i]]
    zz <- zbatch[[i]]
    s <- sigmas[i, ]
    m <- thetas[i, ]
    (yy - m[zz])/s[zz]
  }
  ystar[names(yy)]
})

##.plotBatch <- function(object, use.current=FALSE, show.batch=TRUE, ...){
##  ##  browser()
##  pi <- p(object)
##  L <- length(y(object))
##  hist(y(object), breaks=L/50, col="gray", border="gray", freq=FALSE, ...)
##  if(!use.current){
##    mc <- mcmcChains(object)
##    thetas <- matrix(colMeans(theta(mc)), nBatch(object), k(object))
##    sds <- matrix(colMeans(sigma(mc)), nBatch(object), k(object))
##    rownames(thetas) <- rownames(sds) <- uniqueBatch(object)
##    pi <- colMeans(pic(object))
##  } else {
##    ## use current value
##    thetas <- theta(object)
##    sds <- sigma(object)
##  }
##  cols <- brewer.pal(max(k(object), 3),  "Set1")
##  xx <- seq(min(observed(object)), max(observed(object)),  length.out=10e3)
##  P <- matrix(pi, length(xx), k(object), byrow=TRUE)
##  B <- batch(object)
##  marginal.prob <- matrix(NA, length(xx), k(object))
##  for(j in seq_len(k(object))){
##    p.cummulative <- matrix(NA, length(xx), nBatch(object))
##    k <- 1
##    for(b in uniqueBatch(object)){
##      ##p.x <- dnorm(xx[B==b], mean=thetas[b, j], sd=sds[b, j])
##      p.x <- dnorm(xx, mean=thetas[b, j], sd=sds[b, j])
##      if(show.batch) lines(xx, mean(B==b)*pi[j]*p.x, col=cols[j], lwd=2)
##      ##p.cummulative[, k] <- dnorm(xx, mean=thetas[b, j], sd=sds[b, j])
##      p.cummulative[, k] <- p.x
##      k <- k+1
##    }
##    pbatch <- table(batch(object))/L
##    pbatch <- matrix(pbatch, length(xx), nBatch(object), byrow=TRUE)
##    pcum <- rowSums(pbatch * p.cummulative)
##    ##lines(xx, pcum, col="gray", lwd=2)
##    marginal.prob[, j] <- pcum
##  }
##  marginal.cum.prob <- rowSums(P*marginal.prob)
##  limits <- list(range(y(object), na.rm=TRUE), range(marginal.cum.prob, na.rm=TRUE))
##  lines(xx, marginal.cum.prob, col="black", lwd=2)
##}

setGeneric("getThetaOrder", function(object) standardGeneric("getThetaOrder"))

setMethod("getThetaOrder", "MarginalModel", function(object) order(thetaMean(object)))
setMethod("getThetaOrder", "BatchModel", function(object){
  .getThetaOrdering(object)
})

.getThetaOrdering <- function(object){
  ##th <- thetaMean(object)
  th <- matrix(colMeans(thetac(object)), nBatch(object), k(object))
  ##ix <- object@theta_order
  ix <- matrix(NA, nBatch(object), k(object))
  index <- matrix(seq_len(nBatch(object)*k(object)), nBatch(object), k(object))
  for(j in 1:nrow(th)){
    ix[j, ] <- order(th[j,])
    index[j, ] <- index[j, ix[j, ]]
  }
  index <- as.numeric(index)
  index
}

orderCols <- function(th, x){
  ix <- t(apply(th, 1, order))
  for(j in 1:nrow(th)){
    k <- ix[j, ]
    x[j, ] <- x[j, k]
  }
  x
}

setMethod("orderTheta", "MixtureModel", function(object) object@theta_order)


setReplaceMethod("orderTheta", "MixtureModel", function(object, value) {
  object@theta_order <- value
  object
})

##
##
## use empirical, batch=specific mixing probabilities
##
##
.plotBatch <- function(object, use.current=FALSE, show.batch=TRUE, ...){
  L <- length(y(object))
  hist(object)
  xx <- seq(min(observed(object)), max(observed(object)),  length.out=10e3)
  if(!use.current){
    thetas <- thetaMean(object)
    sds <- sigmaMean(object)
    sds <- orderCols(thetas, sds)
    P <- pMean(object)
    P <- orderCols(thetas, P)
    thetas <- orderCols(thetas, thetas)
  } else {
    thetas <- theta(object)
    sds <- sigma(object)
    P <- p(object)
  }
  marginal <- matrix(NA, length(xx), nBatch(object)*k(object))
  cols <- brewer.pal(max(k(object), 3),  "Set1")
  B <- batch(object)
  marginal.prob <- matrix(NA, length(xx), k(object))
  batchPr <- table(batch(object))/length(y(object))
  m <- 1
  for(j in seq_len(k(object))){
    for(b in uniqueBatch(object)){
      p.x <- dnorm(xx, mean=thetas[b, j], sd=sds[b, j])
      p.x <- batchPr[b] * P[b, j] * p.x
      if(show.batch) lines(xx, p.x, col=cols[j], lwd=2)
      marginal[, m] <- p.x
      m <- m+1
    }
  }
  marginal.cum.prob <- rowSums(marginal)
  limits <- list(range(y(object), na.rm=TRUE), range(marginal.cum.prob, na.rm=TRUE))
  lines(xx, marginal.cum.prob, col="black", lwd=2)
}

setMethod("plot", "BatchModel", function(x, y, use.current=FALSE, show.batch=TRUE, ...){
  .plotBatch(x, use.current, show.batch, ...)
})

setMethod("simulateY", "BatchModel", function(object){
##  browser()
  zz <- setNames(z(object), seq_along(z(object)))
  B <- batch(object)
  zbatch <- split(zz, B)
  sigmas <- sqrt(sigma2(object))
  thetas <- theta(object)
  ysim <- foreach(b = uniqueBatch(object), .combine="c") %do% {
    m <- thetas[b, ]
    s <- sigmas[b, ]
    Z <- zbatch[[b]]
    setNames(rnorm(length(Z), mean=m[Z], s=s[Z]), names(Z))
  }
  ysim[names(zz)]
})

setMethod("moveChain", "BatchModel", function(object, s){
  mcmc <- mcmcChains(object)
  theta(mcmc)[s, ] <- as.numeric(theta(object))
  sigma2(mcmc)[s, ] <- as.numeric(sigma2(object))
  p(mcmc)[s, ] <- p(object)
  mu(mcmc)[s, ] <- mu(object)
  tau2(mcmc)[s, ] <- tau2(object)
  nu.0(mcmc)[s] <- nu.0(object)
  sigma2.0(mcmc)[s] <- sigma2.0(object)
  logpotential(mcmc)[s] <- logpotential(object)
  mcmcChains(object) <- mcmc
  object
})

##
## TODO: pass arguments to .updateThetaBatch to make it clear what
## parameters the theta update depends on
##
setMethod("updateTheta", "BatchModel", function(object) {
  .updateThetaBatch(object)
})

.updateThetaBatch <- function(object){
  ##  if(constrain==2) browser()
  ##  if(constrain!=2) constrain <- TRUE
  tau2.tilde <- 1/tau2(object)
  sigma2.tilde <- 1/sigma2(object)
  K <- k(object)
  ## should be a vector of length K
  tau2.n.tilde <- rep(NA, K)
  n.hp <- tablez(object)
  ##
  ## Guard against zero-components
  ##
  n.hp <- pmax(n.hp, 1)
  ##
  ## mu and tau2 are not batch-specific
  tau2.tilde <- matrix(tau2.tilde, nBatch(object), k(object), byrow=TRUE)
  ##
  ## Make mu the same dimension to make the arithmetic obvious
  ##
  mus <- matrix(mu(object), nBatch(object), k(object), byrow=TRUE)
  ##
  tau2.n.tilde <- tau2.tilde + n.hp * sigma2.tilde
  tau2.n <- 1/tau2.n.tilde
  tau.n <- sqrt(tau2.n)
  ##
  denom <- tau2.tilde + n.hp*sigma2.tilde
  w1 <- tau2.tilde/denom
  w2 <- n.hp*sigma2.tilde/denom
  ##
  ## when a component has 0 observations, mu.n should just be w1*mu
  ##
  ybar <- dataMean(object)
  if(any(is.nan(ybar))){
    ybar[is.nan(ybar)] <- 0
  }
  mu.n <- w1*mus + w2*ybar
  rownames(tau.n) <- rownames(mu.n) <- uniqueBatch(object)
  ##
  ##  If there are very few observations, we will be sampling from a
  ##  normal distribution with mean equal to the marginal mean. This
  ##  will result in thetas that do not satisfy the order constraints
  ##
  ##
  thetas <- matrix(NA, nBatch(object), k(object))
  rownames(thetas) <- uniqueBatch(object)
  for(b in uniqueBatch(object)){
    thetas[b, ] <- rnorm(K, mu.n[b, ], tau.n[b, ])
  }
  if(any(is.na(thetas))) stop("NAs in thetas")
  return(thetas)
}

truncNormLower <- function(theta.last, epsilon){
  a <- theta.last[1] + epsilon
}
truncNormUpper <- function(theta.last, epsilon, iequalsk){
  if(!iequalsk){
    b <- theta.last - epsilon
  } else b <- 5
  b
}

.update_equal_s2 <- function(object){
  sigma2.current <- sigma2(object)
  n.hp <- table(batch(object))[uniqueBatch(object)]
  ##
  ## guard against zeros
  ##
  nu.n <- nu.0(object) + n.hp
  thetas <- theta(object)
  ##
  ## assume variance for all components is the same
  ##
  ss <- sumSquares2(object)
  ##
  ## weighted average of sums of squares
  ##
  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- setNames(rep(NA, nBatch(object)), uniqueBatch(object))
  for(b in uniqueBatch(object)){
    sigma2.h.tilde[b] <- rgamma(1, shape=shape[b], rate=rate[b])
  }
  sigma2.h <- 1/sigma2.h.tilde
  v <- var(y(object), na.rm=TRUE) + 0.05
  if(any(sigma2.h > v)) {
    sigma2.current <- sigma2.current[, 1]
    tmp <- tryCatch(sigma2.h[sigma2.h > v] <- sigma2.current[sigma2.h > v], error=function(e) NULL)
    if(is.null(tmp)) browser()
  }
  ##
  ## return matrix of original dimension
  ##
  s2 <- matrix(sigma2.h, nBatch(object), k(object))
  rownames(s2) <- uniqueBatch(object)
  s2
}

##
## assumes all variance components are the same
## - (the UnivariateBatchModel should also be able to use this update)
##setMethod("updateSigma2", "BatchModelNoHom", function(object){
##  .update_equal_s2(object)
##})
##
####
#### Allows the first component to have a different variance, but
#### currently no constraint on the relationship between the first
#### component variance and the other components
####
##setMethod("updateSigma2", "BatchModelPlusHom", function(object){
##  .update_nonzero_s2(object)
##})
##
##.update_nonzero_s2 <- function(object){
##  sigma2.current <- sigma2(object)
##  nz <- nonZeroCopynumber(object)
##  if(length(unique(nz)) ==1){
##    ## guard against zeros
##    if(all(nz > 0)){
##      nz[which.min(y(object))] <- 0
##    } else nz[which.max(y(object))] <- 1
##  }
##  n.hp <- table(batch(object), nz)
##  n.hp <- n.hp[uniqueBatch(object), ]
##  ##
##  ## guard against zeros
##  ##
##  n.hp <- pmax(n.hp, 1)
##  nu.n <- nu.0(object) + n.hp
##  thetas <- theta(object)
##  ##
##  ## assume variance for copy number 1-k is the same
##  ##
##  ss <- sumSquares(object)
##  ##
##  ## weighted average of sums of squares
##  ##
##  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
##  shape <- 1/2*nu.n
##  rate <- shape*sigma2.nh
##  sigma2.h.tilde <- matrix(NA, nBatch(object), 2)
##  rownames(sigma2.h.tilde) <- rownames(thetas)
##  for(b in uniqueBatch(object)){
##    sigma2.h.tilde[b, ] <- rgamma(2, shape=shape[b, ], rate=rate[b, ])
##  }
##  sigma2.h <- 1/sigma2.h.tilde
##  v <- var(y(object), na.rm=TRUE) + 0.05
##  if(any(sigma2.h > v)) {
##    sigma2.current <- sigma2.current[, 1:2]
##    tmp <- tryCatch(sigma2.h[sigma2.h > v] <- sigma2.current[sigma2.h > v], error=function(e) NULL)
##    if(is.null(tmp)) browser()
##  }
##  ##
##  ## return matrix of original dimension
##  ##
##  s2 <- cbind(sigma2.h[, 1], matrix(sigma2.h[, 2], nBatch(object), k(object)-1))
##  s2
##}

sumSquares <- function(object){
  ss <- matrix(NA, nBatch(object), k(object))
  rownames(ss) <- uniqueBatch(object)
  B <- batch(object)
  thetas <- theta(object)
  yy <- y(object)
  zz <- z(object)
  for(b in uniqueBatch(object)){
    y <- yy[B==b]
    cn <- zz[B==b]
    ##nonzero <- nz[B==b]
    m <- thetas[b, ]
    ## This could be tricky in C.  It works in R because of the factor to an integer:
    ##  as.integer(factor(c(1, 3), levels=c("1", "2", "3"))) evaluates to 1,3
    m <- m[as.integer(cn)]
    squares <- (y - m)^2
    ss[b, ] <- sapply(split(squares, cn), sum)
  }
  ss
}

.update_sigma2 <- function(object){
  sigma2.current <- sigma2(object)
  n.hp <- tablez(object)
  ##
  ## guard against zeros
  ##
  n.hp <- pmax(n.hp, 1)
  nu.n <- nu.0(object) + n.hp
  ss <- sumSquares(object)
  ##
  ## Zeros in sums of squares occurs for batches with no observations
  ##
  ## should handle this by polymorphism
  if(k(object) == 1) ss <- ss[, 1, drop=FALSE]
  ##
  ## weighted average of sums of squares
  ##
  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- matrix(NA, nBatch(object), k(object))
  rownames(sigma2.h.tilde) <- uniqueBatch(object)
  for(b in uniqueBatch(object)){
    sigma2.h.tilde[b, ] <- rgamma(k(object), shape=shape[b, ], rate=rate[b, ])
  }
  sigma2.h <- 1/sigma2.h.tilde
  stopif(any(is.nan(sigma2.h)))
  sigma2.h
}

## special case when there is only one component
setMethod("updateSigma2", "UnivariateBatchModel", function(object){
  .update_sigma2(object)
})

setMethod("updateSigma2", "BatchModel", function(object){
  .update_sigma2(object)
})

nonZeroCopynumber <- function(object) as.integer(as.integer(z(object)) > 1)



## This is a more parsimonious model.  There are only 2 variance
## estimates for each batch: the variance of the first component and
## the variance of components k>1
##.updateSigma2Batch_2 <- function(object){
##
##}

sumSquares2 <- function(object){
  ss <- setNames(rep(NA, nBatch(object)), uniqueBatch(object))
  B <- batch(object)
  thetas <- theta(object)
  ##nz <- nonZeroCopynumber(object)
  for(b in uniqueBatch(object)){
    yy <- y(object)[B==b]
    zz <- z(object)[B==b]
    ##nonzero <- nz[B==b]
    m <- thetas[b, ]
    m <- m[as.integer(zz)]
    ss[b] <- sum((yy - m)^2)
  }
  ss
}

##
## If K is 2, assume there is no homozygous deletion component
##  Constrain the variances to be the same so that one component will
##  not have heavier tails and capture the outliers
##.updateSigma2Batch_samevar <- function(object){
##
##}

setMethod("updateMu", "BatchModel", function(object){
  .updateMuBatch(object)
})

##.updateMu <- function(tau2.0, tau2, k, z, theta, mu.0){
.updateMuBatch <- function(object){
  hypp <- hyperParams(object)
  tau2.0.tilde <- 1/tau2.0(hypp)
  tau2.tilde <- 1/tau2(object)
  P <- nBatch(object)
  tau2.P.tilde <- tau2.0.tilde + P*tau2.tilde
  n.h <- tablez(object)
  ## guard against components with zero observations
  n.h <- pmax(n.h, 1)
  thetas <- theta(object)
  ##
  ## weights for within-component average of thetas
  w1 <- tau2.0.tilde/(tau2.0.tilde + P*tau2.tilde)
  w2 <- P*tau2.tilde/(tau2.0.tilde + P*tau2.tilde)
  ##
  ## average thetas, giving more weight to batches with more
  ## observations
  ##
  theta.bar <- colSums(n.h*thetas)/colSums(n.h)
  ## when the prior is zero, mu.P is shrunk towards zero
  mu.P <- w1*mu.0(hypp) + w2*theta.bar
  mu.P
}

setMethod("updateTau2", "BatchModel", function(object){
  .updateTau2Batch(object)
})

##.updateTau2Batch <- function(eta.0, m2.0, theta, mu, k){
.updateTau2Batch <- function(object){
  hypp <- hyperParams(object)
  P <- nBatch(object)
  eta.P <- eta.0(hypp)+P
  mus <- mu(object)
  mus <- matrix(mus, P, k(object), byrow=TRUE)
  thetas <- theta(object)
  s2.P <- colSums((thetas-mus)^2)
  m2.P <- 1/eta.P * (eta.0(hypp) * m2.0(hypp) + s2.P)
  tau2 <- 1/rgamma(k(object), shape=1/2 * eta.P, rate=1/2 * eta.P * m2.P)
  tau2
}

setMethod("updateSigma2.0", "BatchModel", function(object){
  .updateSigma2.0Batch(object)
})

##.updateSigma2.0 <- function(a, b, nu.0, sigma2.h, k){
.updateSigma2.0Batch <- function(object){
  hypp <- hyperParams(object)
  P <- nBatch(object)
  sigma2s <- as.numeric(sigma2(object))
  a.k <- a(hypp) + 1/2 * (k(object) * P) * nu.0(object)
  b.k <- b(hypp) + 1/2 * sum(1/sigma2s)
  sigma2.0 <- rgamma(1, shape=a.k, rate=b.k)
  stopifnot(sigma2.0 > 0)
  stopif(is.nan(sigma2.0))
  sigma2.0
}

setMethod("updateNu.0", "BatchModel", function(object){
  .updateNu.0Batch(object)
})

##.updateNu.0Batch <- function(NUMAX=100, beta, sigma2.0, sigma2.h, nu.0, k){
.updateNu.0Batch <- function(object){
  NUMAX <- 100
  hypp <- hyperParams(object)
  x <- seq_len(NUMAX)
  k <- k(object)
  P <- nBatch(object)
  sigma2s <- as.numeric(sigma2(object))
  lpnu0 <- (k*P) * (0.5 * x * log(sigma2.0(object) * x/2)-lgamma(x/2)) +
    (x/2 - 1) * sum(log(1/sigma2s)) +
      -x * (betas(hypp) + 0.5 * sigma2.0(object) * sum(1/sigma2s))
  prob <- exp(lpnu0 - max(lpnu0))
  nu0 <- sample(x, 1, prob=prob)
  nu0
}

#' @export
setMethod("mu", "BatchModel", function(object) object@mu)

#' @export
setMethod("tau2", "BatchModel", function(object) object@tau2)

setReplaceMethod("tau2", "BatchModel", function(object, value){
  object@tau2 <- value
  object
})

setReplaceMethod("mu", "BatchModel", function(object, value){
  object@mu <- value
  object
})



##setMethod("initializeSigma2", "UnivariateBatchModel", function(object){
##  tmp <- t(replicate(nBatch(object), 1/rgamma(k(object), shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))))
##  matrix(tmp, ncol=1)
##})

##setMethod("initializeTheta", "UnivariateBatchModel", function(object){
##  tmp <- t(replicate(nBatch(object), sort(rnorm(k(object), mu(object), tau2(object)))))
##  tmp <- matrix(tmp, ncol=1)
##  tmp
##})

##componentCapturesTails <- function(object){
##  ix <- order(y(object))
##  if(any(head(z(object)[ix]) %in% tail(z(object)[ix]))){
##    is_outlier <- TRUE
##  } else is_outlier <- FALSE
##  is_outlier
##}

setMethod("bic", "BatchModel", function(object, ...){
  if(k(object) > 1){
    object <- updateWithPosteriorMeans(object)
  }
  ## K: number of free parameters to be estimated
  ##   - component and batch-specific parameters:  theta, sigma2  ( k(model) * nBatch(model))
  ##   - mixing probabilities: (k-1)*nBatch
  ##   - component-specific parameters: mu, tau2                 2 x k(model)
  ##   - length-one parameters: sigma2.0, nu.0                   +2
  K <- 2*k(object)*nBatch(object) + nBatch(object)*(k(object)-1) + 2*k(object) + 2
  n <- length(y(object))
  bicstat <- -2*logpotential(object) + K*(log(n) - log(2*pi))
  bicstat
})


##setMethod("bic", "BatchModelNoHom", function(object, ...){
##  if(k(object) > 1){
##    object <- updateWithPosteriorMeans(object)
##  }
##  ## K: number of free parameters to be estimated
##  ##   - component and batch-specific parameters:  theta  ( k(model) * nBatch(model))
##  ##   - 1 variance estimate for each batch:  nBatch(model)
##  ##   - component-specific parameters: mu, tau2                 2 x k(model)
##  ##   - length-one parameters: sigma2.0, nu.0                   +2
##  ##   - mixture probs:  +3
##  K <- k(object)*nBatch(object) + nBatch(object) +  2*k(object) + 2 + 3
##  ## Experimental: extra penalty
##  K <- K + 2
##  n <- length(y(object))
##  bicstat <- -2*logpotential(object) + K*(log(n) - log(2*pi))
##  bicstat
##})

setMethod("theta", "BatchModel", function(object) {
  b <- object@theta
  b <- matrix(b, nBatch(object), k(object))
  rownames(b) <- uniqueBatch(object)
  b
})

setReplaceMethod("theta", "BatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@theta <- value
  object
})

setReplaceMethod("p", "BatchModel", function(object, value){
  ##rownames(value) <- uniqueBatch(object)
  object@pi <- value
  object
})

setReplaceMethod("sigma2", "BatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})

setMethod("sigma2", "BatchModel", function(object) {
  s2 <- object@sigma2
  s2 <- matrix(s2, nBatch(object), k(object))
  rownames(s2) <- uniqueBatch(object)
  s2
})

setMethod("reorderComponents", "BatchModel", function(object){
  thetas <- theta(object)
  sigma2s <- sigma2(object)
  for(i in seq_len(nrow(thetas))){
    x <- thetas[i, ]
    j <- order(x)
    thetas[i, ] <- x[j]
    sigma2s[i, ] <- sigma2s[i, j]
    if(i == 1){
      p(object) <- p(object)[j]
    }
    dataMean(object)[i, ] <- dataMean(object)[i, j]
    dataPrec(object)[i, ] <- dataPrec(object)[i, j]
  }
  theta(object) <- thetas
  sigma2(object) <- sigma2s
  ##ix <- order(theta(object)[1, ])
  ##theta(object) <- theta(object)[, ix]
  ##sigma2(object) <- sigma2(object)[, ix]
  ##p(object) <- p(object)[ix]
  zz <- z(object)
  zz <- factor(as.integer(factor(zz, levels=j)), levels=seq_len(k(object)))
  z(object) <- zz
  ##  dataMean(object) <- dataMean(object)[, ix]
  ##  dataPrec(object) <- dataPrec(object)[, ix]
  object
})


setMethod("updateWithPosteriorMeans", "BatchModel", function(object){
  mc <- mcmcChains(object)
  theta(object) <- matrix(colMeans(theta(mc)), nBatch(object), k(object))
  sigma2(object) <- matrix(colMeans(sigma2(mc)), nBatch(object), k(object))
  p(object) <- matrix(colMeans(p(mc)), nBatch(object), k(object))
  nu.0(object) <- median(nu.0(mc))
  mu(object) <- colMeans(mu(mc))
  tau2(object) <- colMeans(tau2(mc))
  sigma2.0(object) <- mean(sigma2.0(object))
  logpotential(object) <- computePotential(object)
  z(object) <- factor(map(object), levels=seq_len(k(object)))
  object
})

setAs("BatchModel", "SummarizedExperiment", function(from, to){
  cnmat <- matrix(y(from), 1, length(y(from)))
  cnmat <- oligoClasses::integerMatrix(cnmat, 1000)
  se <- SummarizedExperiment(assays=SimpleList(medr=cnmat),
                             colData=DataFrame(plate=batch(from)))
  se
})

#' @export
setMethod("sort", "BatchModel", function(x, decreasing=FALSE, ...){
  mc <- mcmcChains(x)
  pot <- logpotential(mc)
  index <- which.max(pot)
  if(FALSE){
    modal.params <- list(p=pic(x)[index, ],
                         theta=thetac(x)[index, ],
                         sigma=sigmac(x)[index, ])
    ##
    ## TODO: foreach iteration, order so that the distance to the modal
    ## parameters is minimized
    ##
  }
  thetas <- matrix(theta(mc)[index, ], nBatch(x), k(x))
  thetas <- thetas[1, ]
  if(identical(thetas, sort(thetas))){
    ## nothing to do
    return(x)
  }
  B <- nBatch(x); K <- k(x)
  cn <- order(thetas)

  ## figure out the appropriate indices to sort
  tmp <- matrix(seq_len(B*K), B, K)
  tmp <- tmp[, cn]
  ix <- as.numeric(tmp)

  theta(mc) <- theta(mc)[, ix]
  theta(x) <- theta(x)[, cn]

  sigma2(mc) <- sigma2(mc)[, ix]
  sigma2(x) <- sigma2(x)[, cn]

  p(mc) <- p(mc)[, cn]
  p(x) <- p(x)[cn]

  mu(mc) <- mu(mc)[, cn]
  tau2(mc) <- tau2(mc)[, cn]
  mu(x) <- mu(x)[cn]
  tau2(x) <- tau2(x)[cn]

  probz(x) <- probz(x)[, cn]

  zz <- as.integer(z(x))
  z(x) <- factor(as.integer(factor(zz, levels=cn)), levels=sort(unique(zz)))
  dataMean(x) <- dataMean(x)[, cn]
  dataPrec(x) <- dataPrec(x)[, cn]
  mcmcChains(x) <- mc
  x
})

.computeDistOneBatch <- function(th, s2, P, mus, tau2s, modes, param.sds){
  thetad <- .absoluteDistance(th, modes[["theta"]])
  sigma2d <- .absoluteDistance(s2, modes[["sigma2"]])
  pid <- .absoluteDistance(P, modes[["mixprob"]])
  mud <- .absoluteDistance(mus, modes[["mu"]])
  tau2d <- .absoluteDistance(tau2s, modes[["tau2"]])
  ##
  ## sum the distances for each parameter matrix and standardize the
  ## total distance by the standard deviation of the modal parameter
  ## estimates
  ##
  thd <- rowSums(thetad)/param.sds[["theta"]]
  s2d <- rowSums(sigma2d)/param.sds[["sigma2"]]
  pid <- rowSums(pid)/param.sds[["mixprob"]]
  mud <- rowSums(mud)/param.sds[["mu"]]
  tau2d <- rowSums(tau2d)/param.sds[["tau2"]]
  if(param.sds[["sigma2"]] > 0){
    tot <- thd+s2d+pid+mud+tau2d
  } else tot <- thd+pid+mud+tau2d
  tot
}

## compute distance for a given permutation of columns
.computeDistanceOnePerm <- function(mc, column.permutation, modes){
  param.sds <- sapply(modes, sd)
  .computeDistOneBatch(th=theta(mc)[, column.permutation],
                       s2=sigma2(mc)[, column.permutation],
                       mus=mu(mc)[, column.permutation],
                       tau2s=tau2(mc)[, column.permutation],
                       P=p(mc)[, column.permutation],
                       modes=modes,
                       param.sds=param.sds)
}

setMethod("computeDistance", "BatchModel", function(object){
  modal.params <- modes(object)
  ##param.sds <- sapply(modal.params, function(x) sd(x[1,]))
  mc <- mcmcChains(object)
  th <- theta(mc)
  s2 <- sigma2(mc)
  ## ix is the number of possible orderings for each batch
  ix <- permutations(k(object), k(object))##
  nr <- nrow(th)
  nc <- nrow(ix)
  Dlist <- vector("list", nBatch(bmodel))
  ## iterate over batches
  for(b in seq_len(nBatch(object))){
    mc2 <- mc
    batch.index <- seq(b, nBatch(object)*k(object), by=nBatch(object))
    theta(mc2) <- th[, batch.index]
    sigma2(mc2) <- s2[, batch.index]
    D <- matrix(NA, nr, nc)
    m.params <- list(theta=modal.params[["theta"]][b,],
                     sigma2=modal.params[["sigma2"]][b,],
                     mu=modal.params[["mu"]],
                     tau2=modal.params[["tau2"]],
                     mixprob=modal.params[["mixprob"]])
    ## iterate over all possible permutations
    for(j in 1:nrow(ix)){
      J <- ix[j, ]
      D[, j] <- .computeDistanceOnePerm(mc=mc2, column.permutation=J,
                                        modes=m.params)
    }
    Dlist[[b]] <- D
  }
  Dlist
})

setMethod("switchLabels", "BatchModel", function(object){
  Dlist <- computeDistance(object)
  mc <- mcmcChains(object)
  warn <- FALSE
  for(b in seq_along(Dlist)){
    D <- Dlist[[b]]
    ordering_index <- apply(D, 1, which.min)
    if(all(ordering_index == 1)) next()
    warn <- TRUE
    batch.index <- seq(b, nBatch(object)*k(object), by=nBatch(object))
    mc2 <- mc
    theta(mc2) <- theta(mc)[, batch.index]
    sigma2(mc2) <- sigma2(mc)[, batch.index]
    perms <- permutations(k(object), k(object))
    tab <- as.integer(names(table(ordering_index)))
    tab <- tab[tab!=1]
    for(i in seq_along(tab)){
      mcmc.index <- which(ordering_index == tab[i])
      j <- perms[tab[i], ]
      ## rewrite batch.index in mc from the permuted index in mc2
      theta(mc)[mcmc.index, batch.index] <- theta(mc2)[mcmc.index, j]
      sigma2(mc)[mcmc.index, batch.index] <- sigma2(mc2)[mcmc.index, j]
      p(mc)[mcmc.index, ] <- p(mc)[mcmc.index, j]
      mu(mc)[mcmc.index,] <- mu(mc)[mcmc.index, j]
      tau2(mc)[mcmc.index,] <- tau2(mc)[mcmc.index, j]
    }
    mcmcChains(object) <- mc
  }
  if(warn) warning("Label switching occurred. Posterior probabilities for z may be incorrect")
  object
})


.computeModesBatch <- function(object){
  mc <- mcmcChains(object)
  th <- theta(mc)
  nr <- nrow(th)
  nc <- ncol(th)
  pot <- logpotential(mc)
  i <- which.max(pot)
  nb <- nBatch(object)
  kk <- k(object)
  thetamax <- matrix(theta(mc)[i, ], nb, kk)
  sigma2max <- matrix(sigma2(mc)[i, ], nb, kk)
  pmax <- p(mc)[i, ]
  mumax <- mu(mc)[i, ]
  tau2max <- tau2(mc)[i,]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mumax,
                tau2=tau2max)
  modes
}

setMethod("computeModes", "BatchModel", function(object){
  .computeModesBatch(object)
})

setMethod("tracePlot", "BatchModel", function(object, name, ...){
  ilist <- foreach(j=1:nBatch(object)) %do% seq(j, nBatch(object)*k(object), nBatch(object))
  uB <- uniqueBatch(object)
  if(name=="theta"){
    ##op <- par(mfrow=c(3, 3), las=1)
    foreach(k=1:nBatch(object)) %do% {
      plot.ts(thetac(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k], ...)
    }
    ##par(op)
  }
  if(name=="sigma"){
    ##op <- par(mfrow=c(nBatch(object)/3, 3), las=1)
    foreach(k=1:nBatch(object)) %do% {
      plot.ts(sigmac(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k],...)
    }
    ##par(op)
  }
  if(name=="p"){
    ##op <- par(mfrow=c(1, k(object)), las=1)
    foreach(k=1:nBatch(object)) %do% {
      plot.ts(pic(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k], ...)
    }
    ##plot.ts(pic(object), col="gray", ...)
    ##par(op)
  }
  if(name=="mu"){
    ##op <- par(mfrow=c(1, k(object)), las=1)
    plot.ts(muc(object),  ...)
    ##par(op)
  }
  if(name=="tau"){
    ##op <- par(mfrow=c(1, k(object)), las=1)
    plot.ts(tauc(object),  ...)
    ##par(op)
  }
})

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})

setMethod("updateZ", "UnivariateBatchModel", function(object){
  factor(rep(1, length(y(object))))
})

setMethod("sigmaMean", "BatchModel", function(object) {
  mns <- colMeans(sigmac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

setMethod("thetaMean", "BatchModel", function(object) {
  mns <- colMeans(thetac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

setMethod("pMean", "BatchModel", function(object) {
  mns <- colMeans(pic(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})
