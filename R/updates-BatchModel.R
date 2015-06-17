
.updateZ <- function(p){
  ## generalize to any k, k >= 1
  ##cumP <- t(apply(p, 1, function(x) cumsum(x)))
  cumP <- rowCumsums(p)
  N <- nrow(p)
  u <- runif(N)
  zz <- rep(NA, N)
  zz[u < cumP[, 1]] <- 1
  k <- 2
  while(k <= ncol(p)){
    zz[u < cumP[, k] & u >= cumP[, k-1]] <- k
    k <- k+1
  }
  ##if(any(is.na(zz))) stop("missing values in zz")
  return(zz)
}

setMethod("updateZ", "MarginalModel", function(object){
  zz <- .Call("update_z", object)
  zz
})

setMethod("updateZ", "BatchModel", function(object){
  zz <- .Call("update_z_batch", object)
  zz
})


setMethod("posteriorMultinomial", "UnivariateBatchModel",
          function(object) return(1))

setMethod("posteriorMultinomial", "BatchModel", function(object){
  ##.multBatch(object)
  .Call("update_multinomialPr_batch", object)
})

.multBatch <- function(object){
  B <- batch(object)
  pi <- p(object)
  sds <- as.numeric(sigma(object))
  thetas <- as.numeric(theta(object))
  x <- y(object)
  K <- k(object)
  nb <- rep(batchElements(object), K)
  xx <- rep(x, K)
  thetas <- rep(thetas, nb)
  sds <- rep(sds, nb)
  nr <- length(x)
  pi <- rep(pi, each=nr)
  lik <- pi*dnorm(x, thetas, sds)
  lik <- matrix(lik, nr, K)
  lik/rowSums(lik)
}

.multBatch2 <- function(object){
  B <- batch(object)
  pi <- p(object)
  sigmas <- sigma(object)
  thetas <- theta(object)
  x <- y(object)
  K <- k(object)
  lik <- matrix(NA, length(x), K)
  for(j in seq_len(K)){
    means <- thetas[B, j]
    sds <- sigmas[B, j]
    lik[, j] <- pi[j]*dnorm(x, means, sds)
  }
  lik/rowSums(lik)
}





setMethod("updateMu", "BatchModel", function(object){
  ##.updateMuBatch(object)
  .Call("update_mu_batch", object)
})

##.updateMu <- function(tau2.0, tau2, k, z, theta, mu.0){
.updateMuBatch <- function(object){
  hypp <- hyperParams(object)
  tau2.0.tilde <- 1/tau2.0(hypp)
  tau2.tilde <- 1/tau2(object)
  P <- nBatch(object)
  post.prec <- tau2.0.tilde + P*tau2.tilde
  n.h <- tablez(object)
  ## guard against components with zero observations
  n.h <- pmax(n.h, 1)
  thetas <- theta(object)
  thetas <- t(apply(thetas, 1, sort))
  ##
  ## weights for within-component average of thetas
  w1 <- tau2.0.tilde/post.prec
  w2 <- P*tau2.tilde/post.prec
  ##
  ## average thetas, giving more weight to batches with more
  ## observations
  ##
  theta.bar <- colSums(n.h*thetas)/colSums(n.h)
  ## when the prior is zero, mu.P is shrunk towards zero
  mu.P <- w1*mu.0(hypp) + w2*theta.bar
  s <- sqrt(1/post.prec)
  rnorm(k(object), mu.P, s)
}

.update_sigma2 <- function(object){
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
  ##browser()
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
  .Call("update_sigma2_batch", object)
  ##.update_sigma2(object)
})


##
## TODO: pass arguments to .updateThetaBatch to make it clear what
## parameters the theta update depends on
##
setMethod("updateTheta", "BatchModel", function(object) {
  ##.updateThetaBatch(object)
  .Call("update_theta_batch", object)
})

.updateThetaBatch <- function(object){
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
  w1 <- tau2.tilde/tau2.n.tilde
  w2 <- n.hp*sigma2.tilde/tau2.n.tilde
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
  ##.updateTau2Batch(object)
  .Call("update_tau2_batch", object)
})

setMethod("updateSigma2.0", "BatchModel", function(object){
  ##.updateSigma2.0Batch(object)
  .Call("update_sigma20_batch", object)
})

##.updateSigma2.0 <- function(a, b, nu.0, sigma2.h, k){
.updateSigma2.0Batch <- function(object){
  hypp <- hyperParams(object)
  P <- nBatch(object)
  sigma2s <- as.numeric(sigma2(object))
  prec <- sum(1/sigma2s)
  a.k <- a(hypp) + 1/2 * (k(object) * P) * nu.0(object)
  b.k <- b(hypp) + 1/2 * sum(1/sigma2s)
  sigma2.0 <- rgamma(1, shape=a.k, rate=b.k)
  stopifnot(sigma2.0 > 0)
  stopif(is.nan(sigma2.0))
  sigma2.0
}

setMethod("updateNu.0", "BatchModel", function(object){
  ##.updateNu.0Batch(object)
  .Call("update_nu0_batch", object)
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
  ##return(lpnu0)
  prob <- exp(lpnu0 - max(lpnu0))
  nu0 <- sample(x, 1, prob=prob)
  nu0
}

##
## z has length y.  Each observation is a sample.
##
setMethod("updateZ", "BatchModel", function(object){
  P <- posteriorMultinomial(object)
  zz <- .updateZ(P)
  factor(zz, levels=seq_len(k(object)))
  ##as.integer(zz)
})

setMethod("updateZ", "UnivariateBatchModel", function(object){
  factor(rep(1, length(y(object))))
})
