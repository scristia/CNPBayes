
setMethod("posteriorMultinomial", "UnivariateBatchModel",
          function(object) return(1))

setMethod("posteriorMultinomial", "BatchModel", function(object){
  .multinomial_probs <- .multBatch(y(object),
                                   theta(object),
                                   sqrt(sigma2(object)),
                                   p(object))
})


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



setMethod("updateSigma2", "BatchModel", function(object){
  .update_sigma2(object)
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

## special case when there is only one component
setMethod("updateSigma2", "UnivariateBatchModel", function(object){
  .update_sigma2(object)
})

setMethod("updateSigma2", "BatchModel", function(object){
  .update_sigma2(object)
})

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


##
## z has length y.  Each observation is a sample.
##
setMethod("updateZ", "BatchModel", function(object){
  P <- posteriorMultinomial(object)
  zz <- .updateZ(P)
  factor(zz, levels=seq_len(k(object)))
})

setMethod("updateZ", "UnivariateBatchModel", function(object){
  factor(rep(1, length(y(object))))
})


##
##.update_equal_s2 <- function(object){
##  sigma2.current <- sigma2(object)
##  n.hp <- table(batch(object))[uniqueBatch(object)]
##  ##
##  ## guard against zeros
##  ##
##  nu.n <- nu.0(object) + n.hp
##  thetas <- theta(object)
##  ##
##  ## assume variance for all components is the same
##  ##
##  ss <- sumSquares2(object)
##  ##
##  ## weighted average of sums of squares
##  ##
##  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
##  shape <- 1/2*nu.n
##  rate <- shape*sigma2.nh
##  sigma2.h.tilde <- setNames(rep(NA, nBatch(object)), uniqueBatch(object))
##  for(b in uniqueBatch(object)){
##    sigma2.h.tilde[b] <- rgamma(1, shape=shape[b], rate=rate[b])
##  }
##  sigma2.h <- 1/sigma2.h.tilde
##  v <- var(y(object), na.rm=TRUE) + 0.05
##  if(any(sigma2.h > v)) {
##    sigma2.current <- sigma2.current[, 1]
##    tmp <- tryCatch(sigma2.h[sigma2.h > v] <- sigma2.current[sigma2.h > v], error=function(e) NULL)
##    if(is.null(tmp)) browser()
##  }
##  ##
##  ## return matrix of original dimension
##  ##
##  s2 <- matrix(sigma2.h, nBatch(object), k(object))
##  rownames(s2) <- uniqueBatch(object)
##  s2
##}
##

##sumSquares <- function(object){
##  ss <- matrix(NA, nBatch(object), k(object))
##  rownames(ss) <- uniqueBatch(object)
##  B <- batch(object)
##  thetas <- theta(object)
##  yy <- y(object)
##  zz <- z(object)
##  for(b in uniqueBatch(object)){
##    y <- yy[B==b]
##    cn <- zz[B==b]
##    ##nonzero <- nz[B==b]
##    m <- thetas[b, ]
##    ## This could be tricky in C.  It works in R because of the factor to an integer:
##    ##  as.integer(factor(c(1, 3), levels=c("1", "2", "3"))) evaluates to 1,3
##    m <- m[as.integer(cn)]
##    squares <- (y - m)^2
##    ss[b, ] <- sapply(split(squares, cn), sum)
##  }
##  ss
##}


##setMethod("updateMixProbs", "BatchModel", function(object){
##  ##With batch-specific, there's not any borrowing
##  ## of strength between batches even in the higher levels of the
##  ## model.  Unless we have a prior on the parameters for the
##  ## Dirichlet.
##  ##
##  alpha.n <- updateAlpha(object)
##  as.numeric(rdirichlet(1, alpha.n))
##})
