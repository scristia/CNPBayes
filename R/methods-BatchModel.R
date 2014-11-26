BatchModel <- function(data, k, batch){
  mcmc.chains=McmcChains()
  B <- length(unique(batch))
  new("BatchModel",
      hyperparams=Hyperparameters(k=k),
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
      logpotential=numeric(1),
      mcmc.chains=mcmc.chains,
      batch=batch)
}

##setMethod("startingValues", "BatchModel", function(object){
##  hypp <- hyperParams(object)
##  B <- batch(object)
##  ylist <- split(y(object), batch(object))
##  foreach(i=seq_along(ylist)) %do% {
##    dat <- ylist[[i]]
##    tmp.file <- tempfile()
##    sink(tmp.file)
##    mmfit <- normalmixEM(dat, arbvar = FALSE,
##                         epsilon = 1e-03,
##                         k=k(hypp), maxit=20)
##    sink()
##    unlink(tmp.file)
##    theta(object)[i, ] <- mmfit$mu[order(mmfit$mu)]
##    sigma2(object)[i, ] <- (mmfit$sigma[order(mmfit$mu)])^2
##  }
##  object
##})

setMethod("startingValues", "BatchModel", function(object){
  ##
  ## Initialize all batches to the same starting values by ignoring
  ## batch
  ##
  hypp <- hyperParams(object)
  tmp.file <- tempfile()
  sink(tmp.file)
  mmfit <- normalmixEM(y(object), arbvar = FALSE,
                       epsilon = 1e-03, k=k(hypp), maxit=20)
  sink()
  unlink(tmp.file)
  B <- length(unique(batch(object)))
  theta(object) <- matrix(mmfit$mu[order(mmfit$mu)], B, k(hypp), byrow=TRUE)
  sigma2(object) <- matrix((mmfit$sigma[order(mmfit$mu)])^2, B, k(hypp), byrow=TRUE)
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

setMethod("initializeTau2", "BatchModel", function(object){
  tau2 <- colSds(theta(object))^2
  tau2[tau2 < 0.1] <- 0.1
  tau2
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
  B <- factor(batch(object), levels=unique(batch(object)))
  ybatch <- split(y(object), B)
  zbatch <- split(z(object), B)
  mus <- mu(object)
  ymeans <- foreach(y=ybatch, z=zbatch, .combine='rbind') %do%{
    mns <- sapply(split(y, z), mean)
    mns[is.na(mns)] <- mus[is.na(mns)]
    mns
  }
  ymeans
})

setMethod("computeVars", "BatchModel", function(object){
  B <- factor(batch(object), levels=unique(batch(object)))
  ybatch <- split(y(object), B)
  zbatch <- split(z(object), B)
  s2.0 <- sigma2.0(object)
  yvars <- foreach(y=ybatch, z=zbatch, .combine='rbind') %do%{
    v <- sapply(split(y, z), var)
    v[is.na(v)] <- s2.0
    v
  }
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
  p.nu.0 <- dgeom(nu.0(object), beta(hypp))
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
  ystar <- foreach(i=seq_along(ybatch), .combine='c') %do% {
    yy <- ybatch[[i]]
    zz <- zbatch[[i]]
    s <- sigmas[i, ]
    m <- thetas[i, ]
    (yy - m[zz])/s[zz]
  }
  ystar[names(yy)]
})

setMethod("plot", "BatchModel", function(x, y, use.current=FALSE, ...){
  object <- x
  hist(object)
  ##yy <- batchCorrect(x)
  pi <- p(object)
  mc <- mcmcChains(object)
  if(!use.current){
    thetas <- matrix(colMeans(theta(mc)), nBatch(object), k(object))
    sds <- matrix(colMeans(sqrt(sigma2(mc))), nBatch(object), k(object))
  } else {
    ## use current value
    thetas <- theta(object)
    sds <- sqrt(sigma2(object))
  }
  cols <- brewer.pal(max(k(object), 3),  "Set1")
  xx <- seq(min(observed(object)), max(observed(object)),  length.out=10e3)
  B <- batch(object)
  for(j in seq_along(pi)){
    for(b in seq_along(uniqueBatch(object))){
      p.x <- dnorm(xx[B==b], mean=thetas[b, j], sd=sds[b, j])
      lines(xx[B==b], mean(B==b)*pi[j]*p.x, col=cols[j], lwd=2)
    }
  }
})

setMethod("simulateY", "BatchModel", function(object){
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
  tau2.tilde <- 1/tau2(object)
  sigma2.tilde <- 1/sigma2(object)
  K <- k(object)
  ## should be a vector of length K
  tau2.n.tilde <- rep(NA, K)
  n.hp <- table(batch(object), z(object))
  ##
  ## mu and tau2 are not batch-specific
  tau2.tilde <- matrix(tau2.tilde, nBatch(object), k(object), byrow=TRUE)
  mus <- matrix(mu(object), nBatch(object), k(object), byrow=TRUE)
  ##
  tau2.n.tilde <- tau2.tilde + n.hp * sigma2.tilde
  tau2.n <- 1/tau2.n.tilde
  ##
  denom <- tau2.tilde + n.hp*sigma2.tilde
  w1 <- tau2.tilde/denom
  w2 <- n.hp*sigma2.tilde/denom
  mu.n <- w1*mus + w2*dataMean(object)
  ##
  ##  If there are very few observations, we will be sampling from a
  ##  normal distribution with mean equal to the marginal mean. This
  ##  will result in thetas that do not satisfy the order constraints
  ##
  ##
  tau.n <- sqrt(tau2.n)
  thetas.last <- theta(object)
  epsilon <- 0.01
  thetas <- foreach(b = uniqueBatch(object), .combine="rbind") %do% {
    thetas <- rnorm(K, mu.n[b, ], tau.n[b, ])
    stopif(any(is.na(thetas)))
    if(identical(sort(thetas), thetas))  return(thetas)
    ##
    ## sorted thetas are not the same
    ##
    ## - constrain updates according to theta.last
    thetas[1] <- rtruncnorm(1, a=-Inf,
                            b=thetas.last[b, 2]-epsilon,
                            mean=thetas.last[b, 1], sd=tau.n[b, 1])
    for(i in 2:K){
      a <- thetas.last[b, i-1] + epsilon
      if(i < K){
        b <- thetas.last[b, i+1] - epsilon
      } else b <- 5
      ## simulate from a truncated normal
      thetas[i] <- rtruncnorm(1, a=a, b=b,
                              mean=mu.n[b, i],
                              sd=tau.n[b, i])
    }
    thetas
  }
  stopif(any(is.na(thetas)))
  thetas
}


setMethod("updateSigma2", "BatchModel", function(object) {
##  .updateSigma2(split(y(object),
##                      z(object)),
##                theta(object),
##                nu.0(object),
##                sigma2.0(object),
  ##                tablez(object))
  .updateSigma2Batch(object)
})

##.updateSigma2 <- function(data.list, thetas, nu.0, sigma2.0, n.h){
.updateSigma2Batch <- function(object){
  n.hp <- table(batch(object), z(object))
  nu.n <- nu.0(object) + n.hp
  ##k <- length(nu.n)
  thetas <- theta(object)
  B <- batch(object)
  ss <- matrix(NA, nBatch(object), k(object))
  for(b in uniqueBatch(object)){
    yy <- y(object)[B==b]
    zz <- z(object)[B==b]
    m <- thetas[b, ]
    m <- m[zz]
    squares <- (yy - m)^2
    ss[b, ] <- sapply(split(squares, zz), sum)
  }
  ##
  ## weighted average of sums of squares
  ##
  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- matrix(NA, nBatch(object), k(object))
  for(b in uniqueBatch(object)){
    sigma2.h.tilde[b, ] <- rgamma(k(object), shape=shape[b, ], rate=rate[b, ])
  }
  ##tmp <- rgamma(1000, shape=1/2*nu.n[1], rate=1/2*nu.n[1]*sigma2.nh[1])
  sigma2.h <- 1/sigma2.h.tilde
  stopif(any(is.nan(sigma2.h)))
  sigma2.h
}


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
  n.h <- table(batch(object), z(object))
  thetas <- theta(object)
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
      -x * (beta(hypp) + 0.5 * sigma2.0(object) * sum(1/sigma2s))
  prob <- exp(lpnu0 - max(lpnu0))
  nu0 <- sample(x, 1, prob=prob)
  nu0
}
