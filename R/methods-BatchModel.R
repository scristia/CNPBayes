
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
      hwe=numeric())
}

##
## Multiple batches, but only 1 component
##
UnivariateBatchModel <- function(data, k=1, batch){
  mcmc.chains <- McmcChains()
  B <- length(unique(batch))
  new("UnivariateBatchModel",
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
      probz=matrix(1, length(data), 1),
      logpotential=numeric(1),
      mcmc.chains=mcmc.chains,
      batch=batch,
      hwe=numeric())
}

setMethod("startingValues", "BatchModel", function(object){
  ##
  ## Initialize all batches to the same starting values by ignoring
  ## batch
  ##
  ##browser()
  hypp <- hyperParams(object)
  if(FALSE){
    tmp.file <- tempfile()
    sink(tmp.file)
    mmfit <- normalmixEM(y(object), arbvar = FALSE, epsilon = 1e-03, k=k(hypp), maxit=2000)
    sink()
    unlink(tmp.file)
    mus <- mmfit$mu
    vars <- (mmfit$sigma[order(mmfit$mu)])^2
    B <- nBatch(object)
    mus <- matrix(mus, B, k(object), byrow=TRUE)
    vars <- matrix(vars, B, k(object), byrow=TRUE)
    rownames(vars) <- rownames(mus) <- uniqueBatch(object)
    theta(object) <- mus
    sigma2(object) <- vars
  }
  if(TRUE){
    if(all(is.na(theta(object)))){
      theta(object) <- initializeTheta(object)
    }
    hypp <- hyperParams(object)
    B <- nBatch(object)
    sigma2(object) <- matrix(mad(as.numeric(y(object)), na.rm=TRUE)^2, B, k(object))
    rownames(theta(object)) <- rownames(sigma2(object)) <- uniqueBatch(object)
  }
  object
})

setMethod("initializeSigma2.0", "BatchModel", function(object){
  hypp <- hyperParams(object)
  sum(alpha(hypp)*colMeans(sigma2(object)))/sum(alpha(hypp))
})


setMethod("posteriorMultinomial", "BatchModel", function(object){
  .multinomial_probs <- .posteriorMultinomialBatch(y(object),
                                                   theta(object),
                                                   sqrt(sigma2(object)),
                                                   p(object))

})


setMethod("posteriorMultinomial", "UnivariateBatchModel", function(object) return(1))

#' @export
uniqueBatch <- function(object) unique(batch(object))

nBatch <- function(object) length(uniqueBatch(object))


.posteriorMultinomialBatch <- function(y, theta, sd, pi){
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
  ## all(rowSums(mix.probs) == 1)
  mix.probs
}

##setMethod("updateMu", "BatchModel", function(object){
##  hypp <- hyperParams(object)
##  mu.k <- .updateMuBatch(tau2.0(hypp), tau2(object), k(object), z(object),
##                         theta(object), mu.0(hypp))
##  mu.k
##})
##
##
##
##.updateMuBatch <- function(tau2.0, tau2, k, z, theta, mu.0){
##  P <- nrow(theta)
##  tau2.0.tilde <- 1/tau2.0
##  tau2.tilde <- 1/tau2
##  tau2.P.tilde <- tau2.0.tilde + P*tau2.tilde
##  nn <- table(z)
##  ## average across batch for each component
##  theta.bar.h <- colMeans(theta)
##  mu.h <- tau2.0.tilde/(tau2.P.tilde)*mu.0 +
##    P*tau2.tilde/(tau2.P.tilde)*theta.bar.h
##  mu.h
##}

##setMethod("initializeTau2", "BatchModel", function(object){
##  hypp <- hyperParams(object)
##  tau2s <- 1/rgamma(1, shape=1/2*eta.0(hypp), rate=1/2*eta.0(hypp)*m2.0(hypp))
##  ## assume initially that the taus are the same
##  rep(tau2s, nBatch(object))
##})

##
## Within-component heterogeneity of the thetas
##
setMethod("initializeTau2", "BatchModel", function(object){
  s <- mad(y(object), na.rm=TRUE)/4
  rep(s^2,  k(object))
})

setMethod("initializeMu", "BatchModel", function(object) colMeans(theta(object)))

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
  mus <- mu(object)
  ymeans <- foreach(y=ybatch, z=zbatch, .combine='rbind') %do%{
    mns <- sapply(split(y, z), mean)
    mns[is.na(mns)] <- mus[is.na(mns)]
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
    v[is.na(v)] <- s2.0
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
  callNextMethod()
  cat("     n. batches  :", nBatch(object), "\n")
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

## use empirical, batch=specific mixing probabilities
.plotBatch <- function(object, use.current=FALSE, show.batch=TRUE, ...){
  L <- length(y(object))
  hist(y(object), breaks=L/50, col="gray", border="gray", freq=FALSE, ...)
  xx <- seq(min(observed(object)), max(observed(object)),  length.out=10e3)
  if(!use.current){
    zz <- map(object)
    tabz <- table(batch(object), zz)
    ttl <- rowSums(tabz)
    pi <- tabz/ttl
    pi <- pi[uniqueBatch(object), , drop=FALSE]

    thetas <- matrix(colMeans(thetac(object)), nBatch(object), k(object))
    sds <- matrix(colMeans(sigmac(object)), nBatch(object), k(object))
    rownames(thetas) <- rownames(sds) <- uniqueBatch(object)
    P <- matrix(colMeans(pic(object)), length(xx), k(object), byrow=TRUE)
    ##pi <- colMeans(pic(object))
  } else {
    tabz <- table(batch(object), z(object))
    ttl <- rowSums(tabz)
    pi <- tabz/ttl
    pi <- pi[uniqueBatch(object), ]
    ## use current value
    thetas <- theta(object)
    sds <- sigma(object)
    P <- matrix(p(object), length(xx), k(object), byrow=TRUE)
    ##P <- pi[batch(object), ]
  }
  cols <- brewer.pal(max(k(object), 3),  "Set1")
  B <- batch(object)
  marginal.prob <- matrix(NA, length(xx), k(object))
  for(j in seq_len(k(object))){
    p.cummulative <- matrix(NA, length(xx), nBatch(object))
    k <- 1
    for(b in uniqueBatch(object)){
      ##p.x <- dnorm(xx[B==b], mean=thetas[b, j], sd=sds[b, j])
      p.x <- dnorm(xx, mean=thetas[b, j], sd=sds[b, j])
      if(show.batch) lines(xx, mean(B==b)*pi[b, j]*p.x, col=cols[j], lwd=2)
      ##p.cummulative[, k] <- dnorm(xx, mean=thetas[b, j], sd=sds[b, j])
      p.cummulative[, k] <- p.x
      k <- k+1
    }
    pbatch <- table(batch(object))/L
    pbatch <- pbatch[uniqueBatch(object)]
    pbatch <- matrix(pbatch, length(xx), nBatch(object), byrow=TRUE)
    pcum <- rowSums(pbatch * p.cummulative)
    ##lines(xx, pcum, col="gray", lwd=2)
    marginal.prob[, j] <- pcum
  }
  ##browser()
  marginal.cum.prob <- rowSums(P*marginal.prob)
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
setMethod("updateTheta", "BatchModel", function(object, constrain=TRUE) {
  .updateThetaBatch(object, constrain=constrain)
})

.updateThetaBatch <- function(object, constrain=TRUE){
  ##  if(constrain==2) browser()
  ##  if(constrain!=2) constrain <- TRUE
  tau2.tilde <- 1/tau2(object)
  sigma2.tilde <- 1/sigma2(object)
  K <- k(object)
  ## should be a vector of length K
  tau2.n.tilde <- rep(NA, K)
  n.hp <- table(batch(object), z(object))
  ##
  ## Guard against zero-components
  ##
  n.hp <- pmax(n.hp, 1)
  ##
  ## mu and tau2 are not batch-specific
  tau2.tilde <- matrix(tau2.tilde, nBatch(object), k(object), byrow=TRUE)
  mus <- matrix(mu(object), nBatch(object), k(object), byrow=TRUE)
  ##
  tau2.n.tilde <- tau2.tilde + n.hp * sigma2.tilde
  tau2.n <- 1/tau2.n.tilde
  tau.n <- sqrt(tau2.n)
  ##
  denom <- tau2.tilde + n.hp*sigma2.tilde
  w1 <- tau2.tilde/denom
  w2 <- n.hp*sigma2.tilde/denom
  mu.n <- w1*mus + w2*dataMean(object)
  rownames(tau.n) <- rownames(mu.n) <- uniqueBatch(object)
  ##
  ##  If there are very few observations, we will be sampling from a
  ##  normal distribution with mean equal to the marginal mean. This
  ##  will result in thetas that do not satisfy the order constraints
  ##
  ##
  thetas.last <- theta(object)
  ##thetas.last <- matrix(theta(object), nBatch(object), k(object))
  ##rownames(thetas.last) <- rownames(mu.n)
  epsilon <- 1/1000
  thetas <- matrix(NA, nrow(thetas.last), ncol(thetas.last))
  rownames(thetas) <- rownames(thetas.last)
  ##
  if(!constrain){
    for(B in uniqueBatch(object)){
      thetas[B, ] <- rnorm(K, mu.n[B, ], tau.n[B, ])
    }
    return(thetas)
  }
  ## constrain thetas:
  for(B in uniqueBatch(object)){
    tmp <- rnorm(K, mu.n[B, ], tau.n[B, ])
    if(identical(sort(tmp), tmp)) {
      thetas[B, ] <- tmp
      next()
    }
    if(!constrain){
      thetas[B, ] <- tmp
      next()
    }
    thetas[B, 1] <- rtruncnorm(1, a=-Inf,
                               b=thetas.last[B, 2]-epsilon,
                               mean=thetas.last[B, 1], sd=tau.n[B, 1])
    for(i in 2:K){
      a <- truncNormLower(thetas.last[B, i-1], epsilon)
      b <- truncNormUpper(thetas.last[B, i+1], epsilon, i==K)
      if(a > b){
        epsilon <- epsilon/100
        a <- truncNormLower(thetas.last[B, i-1], epsilon)
        b <- truncNormUpper(thetas.last[B, i+1], epsilon, i==K)
      }
      thetas[B, i] <- rtruncnorm(1, a=a, b=b,
                                 mean=mu.n[B, i],
                                 sd=tau.n[B, i])
    }
  }
  ##stopif(any(is.na(thetas)))
  if(any(is.na(thetas))) browser()
  thetas
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


setMethod("updateSigma2", "BatchModel", function(object) {
##  .updateSigma2(split(y(object),
##                      z(object)),
##                theta(object),
##                nu.0(object),
##                sigma2.0(object),
  ##                tablez(object))
  .updateSigma2Batch_2(object)
  ##.updateSigma2Batch(object)
})

setMethod("updateSigma2", "UnivariateBatchModel", function(object) {
  .updateSigma2Batch(object)
})

##.updateSigma2 <- function(data.list, thetas, nu.0, sigma2.0, n.h){
.updateSigma2Batch <- function(object){
  sigma2.current <- sigma2(object)
  n.hp <- table(batch(object), z(object))
  ##
  ## guard against zeros
  ##
  n.hp <- pmax(n.hp, 1)
  n.hp <- n.hp[uniqueBatch(object), , drop=FALSE]
  ##
  ##
  nu.n <- nu.0(object) + n.hp
  ##k <- length(nu.n)
  thetas <- theta(object)
  rownames(thetas) <- uniqueBatch(object)
  B <- batch(object)
  ss <- sumSquares(object)
  if(k(object) == 1) ss <- ss[, 1, drop=FALSE]
  ##
  ## weighted average of sums of squares
  ##
  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- matrix(NA, nBatch(object), k(object))
  rownames(sigma2.h.tilde) <- rownames(thetas)
  for(b in uniqueBatch(object)){
    sigma2.h.tilde[b, ] <- rgamma(k(object), shape=shape[b, ], rate=rate[b, ])
  }
  ##tmp <- rgamma(1000, shape=1/2*nu.n[1], rate=1/2*nu.n[1]*sigma2.nh[1])
  sigma2.h <- 1/sigma2.h.tilde
  stopif(any(is.nan(sigma2.h)))
  v <- var(y(object), na.rm=TRUE) + 0.01
  if(any(sigma2.h > v)) {
    sigma2.h[sigma2.h > v] <- sigma2.current[sigma2.h > v]
  }
  sigma2.h
}

nonZeroCopynumber <- function(object) as.integer(as.integer(z(object)) > 1)

sumSquares <- function(object){
  ss <- matrix(NA, nBatch(object), 2)
  rownames(ss) <- uniqueBatch(object)
  B <- batch(object)
  thetas <- theta(object)
  nz <- nonZeroCopynumber(object)
  for(b in uniqueBatch(object)){
    yy <- y(object)[B==b]
    zz <- z(object)[B==b]
    nonzero <- nz[B==b]
    m <- thetas[b, ]
    m <- m[as.integer(zz)]
    squares <- (yy - m)^2
    ss[b, ] <- sapply(split(squares, nonzero), sum)
  }
  ss
}
## This is a more parsimonious model.  There are only 2 variance
## estimates for each batch: the variance of the first component and
## the variance of components k>1
.updateSigma2Batch_2 <- function(object){
  sigma2.current <- sigma2(object)
  nz <- nonZeroCopynumber(object)
  if(length(unique(nz)) ==1){
    ## guard against zeros
    if(all(nz > 0)){
      nz[which.min(y(object))] <- 0
    } else nz[which.max(y(object))] <- 1
  }
  ##n.hp <- table(batch(object), z(object))
  n.hp <- table(batch(object), nz)
  n.hp <- n.hp[uniqueBatch(object), ]
  ##
  ## guard against zeros
  ##
  n.hp <- pmax(n.hp, 1)
  ##
  ##
  nu.n <- nu.0(object) + n.hp
  ##k <- length(nu.n)
  thetas <- theta(object)
  ##
  ## assume variance for copy number 1-k is the same
  ##
  ss <- sumSquares(object)
  ##
  ## weighted average of sums of squares
  ##
  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- matrix(NA, nBatch(object), 2)
  rownames(sigma2.h.tilde) <- rownames(thetas)
  for(b in uniqueBatch(object)){
    sigma2.h.tilde[b, ] <- rgamma(2, shape=shape[b, ], rate=rate[b, ])
  }
  ##tmp <- rgamma(1000, shape=1/2*nu.n[1], rate=1/2*nu.n[1]*sigma2.nh[1])
  sigma2.h <- 1/sigma2.h.tilde
  ##stopif(any(is.nan(sigma2.h)))
  v <- var(y(object), na.rm=TRUE) + 0.05
  if(any(sigma2.h > v)) {
    sigma2.current <- sigma2.current[, 1:2]
    tmp <- tryCatch(sigma2.h[sigma2.h > v] <- sigma2.current[sigma2.h > v], error=function(e) NULL)
    if(is.null(tmp)) browser()
  }
  ##
  ## return matrix of original dimension
  ##
  s2 <- cbind(sigma2.h[, 1], matrix(sigma2.h[, 2], nBatch(object), k(object)-1))
  s2
}


setMethod("updateMu", "BatchModel", function(object){
  .updateMuBatch(object)
})

##.updateMu <- function(tau2.0, tau2, k, z, theta, mu.0){
.updateMuBatch <- function(object){
##  browser()
  hypp <- hyperParams(object)
  tau2.0.tilde <- 1/tau2.0(hypp)
  tau2.tilde <- 1/tau2(object)
  P <- nBatch(object)
  tau2.P.tilde <- tau2.0.tilde + P*tau2.tilde
  n.h <- table(batch(object), z(object))
  ## guard against components with zero observations
  n.h <- pmax(n.h, 1)
  thetas <- theta(object)
  ##
  ## between-batch average of thetas
  theta.bar <- colSums(n.h*thetas)/colSums(n.h)
  mu.P <- tau2.0.tilde/(tau2.0.tilde + P*tau2.tilde)*mu.0(hypp) +
    P*tau2.tilde/(tau2.0.tilde + P*tau2.tilde)*theta.bar
  stopif(any(is.nan(mu.P)))
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
  ##s2.P <- sum((theta-mu)^2)
  m2.P <- 1/eta.P * (eta.0(hypp) * m2.0(hypp) + s2.P)
  tau2 <- 1/rgamma(k(object), shape=1/2 * eta.P, rate=1/2 * eta.P * m2.P)
  stopif(is.nan(tau2))
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


setMethod("initializeTheta", "BatchModel", function(object){
  th <- matrix(initializeTheta(k(object)), nBatch(object), k(object), byrow=TRUE)
  rownames(th) <- uniqueBatch(object)
  th
})

setMethod("bic", "BatchModel", function(object, ...){
  if(k(object) > 1){
    object <- updateWithPosteriorMeans(object)
  }
  ## K: number of free parameters to be estimated
  ##   - component and batch-specific parameters:  theta  ( k(model) * nBatch(model))
  ##   - 2 variance estimates for each batch:  2*nBatch(model)
  ##   - component-specific parameters: mu, tau2                 2 x k(model)
  ##   - length-one parameters: sigma2.0, nu.0                   +2
  ##   - mixture probs:  +3
  K <- k(object)*nBatch(object) + 2*nBatch(object) +  2*k(object) + 2 + 3
  n <- length(y(object))
  -2*logpotential(object) + K*(log(n) - log(2*pi))
})

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
  p(object) <- colMeans(p(mc))
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
