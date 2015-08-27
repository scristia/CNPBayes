
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
  zz <- update_z(object)
  zz
})

setMethod("updateZ", "BatchModel", function(object){
  zz <- update_z_batch(object)
  zz
})


setMethod("posteriorMultinomial", "UnivariateBatchModel",
          function(object) return(1))

setMethod("posteriorMultinomial", "BatchModel", function(object){
  ##.multBatch(object)
  update_multinomialPr_batch(object)
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
