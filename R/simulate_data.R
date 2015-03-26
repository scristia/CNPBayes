#' @export
simulateBatchData <- function(N=2500, p, theta, sds, batch, zz){
  if(missing(batch)) {
    batch <- rep(1L, N)
  } else {
    batch <- as.integer(factor(batch))
  }
  if(missing(zz)) {
    zz <- simulateZ(N, p)
  }
  yy <- rep(NA, N)
  ub <- unique(batch)
  rownames(theta) <- rownames(sds) <- ub
  for(b in ub){
    index <- which(batch==b)
    nn <- length(index)
    cn <- zz[index]
    mu <- theta[b, ]
    s <- sds[b, ]
    yy[index] <- rnorm(nn, mu[cn], s[cn])
  }
  ##browser()
  ix <- order(batch)
  object <- BatchModel(yy, batch=batch, k=ncol(theta))
  z(object) <- as.integer(factor(zz))[ix]
  ##
  ## Must initialize the slots independently (must point to different
  ## locations in memory)
  ##
  theta(object) <- computeMeans(object)
  dataMean(object) <- computeMeans(object)
  sigma2(object) <- computeVars(object)
  dataPrec(object) <- 1/computeVars(object)
  p(object) <- as.numeric(table(z(object))/N)
  logLik(object) <- computeLoglik(object)
  object
}


#' @export
simulateData <- function(N, p, theta, sds){
  zz <- simulateZ(N, p)
  y <- rnorm(N, theta[zz], sds[zz])
  ##params <- ModelParams("marginal", y=y, k=length(theta))
  ##object <- initializeModel(params)
  object <- MarginalModel(data=y, k=length(theta))
  z(object) <- as.integer(factor(zz, levels=unique(sort(zz))))
  p(object) <- p
  theta(object) <- as.numeric(sapply(split(y(object), z(object)), mean))
  sigma2(object) <- as.numeric(sapply(split(y(object), z(object)), var))
  p(object) <- as.numeric(sapply(split(y(object), z(object)), length)/length(z(object)))
  mu(object) <- mean(theta(object))
  tau2(object) <- var(theta(object))
  logLik(object) <- computeLoglik(object)
  logPrior(object) <- computePrior(object)
  ##logpotential(object) <- computePotential(object)
  object
}

simulateZ <- function(N, p){
  P <- rdirichlet(N, p)
  cumP <- t(apply(P, 1, cumsum))
  u <- runif(N)
  zz <- rep(NA, N)
  zz[u < cumP[, 1]] <- 1
  k <- 2
  while(k <= ncol(P)){
    zz[u < cumP[, k] & u >= cumP[, k-1]] <- k
    k <- k+1
  }
  zz
}

cumProbs <- function(p, k){
  pcum <- list()
  cols <- 2:(k-1)
  for(j in seq_along(cols)){
    g <- cols[j]
    pcum[[j]] <- rowSums(p[, 1:g, drop=FALSE])
  }
  pcum2 <- cbind(p[, 1], do.call(cbind, pcum), 1)
}
