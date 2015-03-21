loglik.normmix <-
  function(r, mixture, K, burnin=1) {
    loglike <- dnormmix(r, mixture, K, burnin=1, log=TRUE)
    return(sum(loglike))
  }


.computeLoglikBatch <- function(object){
  ll.data <- logLikData(object)
  ll.data
}

setMethod("computeLoglik", "BatchModel", function(object){
  .computeLoglikBatch(object)
})


setMethod("computeLoglik", "MarginalModel", function(object){
  .computeLoglik(object)
})


setMethod("computePotential", "BatchModel", function(object){
  computeLogLikxPrior(object)
})

setMethod("computeLogLikxPrior", "MixtureModel", function(object){
  ##   log.prior <- computePrior(object)
  ##   loglik <- computeLoglik(object)
  ##   loglik + log.prior
  .Call("compute_llxprior", object)
})


logLikData <- function(object){
  B <- batch(object)
  mn <- theta(object)
  ss <- sigma(object)
  x <- y(object)
  tabz <- table(B, z(object))
  P <- tabz/rowSums(tabz)
  ## Vectorize
  lk <- k(object)
  xx <- rep(x, lk)
  nb <- rep(batchElements(object), lk)
  means <- rep(as.numeric(mn), nb)
  sds <- rep(as.numeric(ss), nb)
  p <- rep(as.numeric(P), nb)
  lik <- p*dnorm(xx, means, sds)
  lik <- matrix(lik, length(x), lk)
  lik <- rowSums(lik)
  sum(log(lik))
}

logLikPhi <- function(object){
  thetas <- theta(object)
  mus <- mu(object)
  tau2s <- tau2(object)
  sigma2s <- sigma2(object)
  ##
  ## Vectorize
  ##
  nr <- nrow(thetas)
  nc <- ncol(thetas)
  mus <- rep(mus, each=nr)
  taus <- rep(sqrt(tau2s), each=nr)
  thetas <- as.numeric(thetas)
  p.theta <- dnorm(thetas, mus, taus)
  p.theta <- matrix(p.theta, nr, nc)
  rownames(p.theta) <- uniqueBatch(object)
  p.sigma2 <- dgamma(1/sigma2s, shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  sum(log(p.theta)) + sum(log(p.sigma2))
}

setMethod("computePrior", "BatchModel", function(object){
  hypp <- hyperParams(object)
  K <- k(hypp)
  tau2s <- tau2(object)
  mus <- mu(object)
  p.mu <- dnorm(mus, mu.0(hypp), sqrt(tau2.0(hypp)))
  p.sigma2.0 <- dgamma(sigma2.0(object), shape=a(hypp), rate=b(hypp))
  p.nu.0 <- dgeom(as.integer(nu.0(object)), betas(hypp))
  sum(log(p.mu)) + log(p.sigma2.0) + log(p.nu.0)
})



.loglikMarginal <- function(object){
  x <- y(object)
  nr <- length(x)
  pp <- p(object)
  K <- k(object)
  x <- rep(x, K)
  p <- rep(p(object), each=nr)
  thetas <- rep(theta(object), each=nr)
  sigmas <- rep(sigma(object), each=nr)
  lik <- matrix(p*dnorm(x, thetas, sigmas),
                nr, K)
  sum(log(rowSums(lik)))
}

.loglikPhiMarginal <- function(object){
  thetas <- theta(object)
  mus <- mu(object)
  tau2s <- tau2(object)
  sigma2s <- sigma2(object)
  p.theta <- dnorm(thetas, mus, sqrt(tau2s))
  p.sigma2 <- dgamma(1/sigma2s, shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  sum(log(p.theta)) + sum(log(p.sigma2))
}

.computeLoglik <- function(object){
  ##ll.data <- .loglikMarginal(object)
  ll.data <- .Call("loglik", object)
  ll.data
}

setMethod("computePrior", "MarginalModel", function(object){
  ##  .compute_prior_marginal(object)
  .Call("compute_priorPr", object)
})

.compute_prior_marginal <- function(object){
  hypp <- hyperParams(object)
  K <- k(hypp)
  mus <- mu(object)
  p.sigma2.0 <- dgamma(sigma2.0(object), shape=a(hypp), rate=b(hypp))
  p.nu.0 <- dgeom(nu.0(object), betas(hypp))
  p.mu <- dnorm(mus, mu.0(hypp), sqrt(tau2.0(hypp)))
  log(p.mu) + log(p.sigma2.0) + log(p.nu.0)
}
