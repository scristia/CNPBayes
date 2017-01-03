#' @include AllClasses.R
NULL

setMethod("initializeTheta", "BatchModel", function(object){
  th <- matrix(NA, nBatch(object), k(object))
  for(j in seq_len(k(object))){
    th[, j] <- rnorm(nBatch(object), mu(object)[j], tau(object)[j])
  }
  th
})

.sigma2_batch <- function(object){
  s2 <- 1/rgamma(nBatch(object)*k(object),
                 shape=1/2*nu.0(object),
                 rate=1/2*nu.0(object)*sigma2.0(object))
  s2 <- matrix(s2, nBatch(object), k(object))
  s2
}

setMethod("initializeSigma2", "BatchModel", function(object){
  .sigma2_batch(object)
})

simulateThetas <- function(y, K){
  prob.kmeans <- runif(1, 0, 1)
  if(prob.kmeans < 0.5){
    ## simulate thetas based on a priori expectation
    homo.del <- runif(1, -5, -0.5)
    hemi.del <- rnorm(1, -0.5, 0.5)
    diploid <- rnorm(1, 0, 0.5)
    dupl <- rnorm(1, 0.5, 0.5)
    twocopy.gain <- rnorm(1, 1, 0.5)
    thetas <- c(homo.del, hemi.del, diploid, dupl, twocopy.gain)
    thetas <- sample(thetas, K, replace=FALSE)
  } else {
    ## simulate starting values near kmean centers
    thetas <- kmeans(y, K)$centers[, 1]
    thetas <- thetas + rnorm(K, thetas, 0.1)
  }
  sort(thetas)
}


.init_sb2 <- function(object){
  hypp <- hyperParams(object)
  if(length(y(object))==0) return(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  ys <- y(object)
  ys <- sample(ys, length(ys), replace=TRUE)
  K <- k(object)
  mc <- Mclust(ys, G=K)
  params <- mc$parameters
  thetas <- params$mean
  s2s <- params$variance$sigmasq
  ps <- params$pro
  ##zs <- as.integer(mc$classification)
  zz <- as.integer(simulateZ(length(y(object)), p(object)))
  cnt <- 1
  while (length(table(zz)) < K) {
    if (cnt > 10) {
      stop("Too few observations or too many components.")
    }
    p(object) <- as.numeric(rdirichlet(1, alpha(hypp))) ## rows are
    zz <- as.integer(simulateZ(length(y(object)), p(object)))
    cnt <- cnt + 1
  }
  mu(object) <- mean(thetas)
  if(K > 1){
    tau2(object) <- var(thetas)
  } else tau2(object) <- var(y(object))*10
  p(object) <- ps
  z(object) <- zz
  theta(object) <- thetas
  sigma2(object) <- s2s
  zFreq(object) <- as.integer(table(zz))
  dataPrec(object) <- 1/computeVars(object)
  dataMean(object) <- computeMeans(object)
  log_lik(object) <- computeLoglik(object)
  logPrior(object) <- computePrior(object)
  chains(object) <- McmcChains(object)
  object
}

.init_sb1 <- function(object){
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  ##mu(object) <- rnorm(1, mu.0(object), tau.0(object))
  mu(object) <- mean(y(object))
  ##tau2(object) <- 1/rgamma(1, shape=1/2*eta.0(object),
  ##rate=1/2*eta.0(object)*m2.0(object))
  tau2(object) <- var(y(object))*10
  theta(object) <- rnorm(k(object), mu(object), tau(object))
  sigma2s <- 1/rgamma(k(object),
                      shape=1/2*nu.0(object),
                      rate=1/2*nu.0(object)*sigma2.0(object))
  sigma2(object) <- sigma2s
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
}

setMethod("startingValues", "MarginalModel", function(object){
  ##object <- .init_sb1(object)
  object <- .init_sb2(object)
  object
})

.init_kmeans <- function(object){
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  B <- uniqueBatch(object)
  nB <- nBatch(object)
  K <- k(object)
  ##T <- theta(object)
  ##S <- sigma(object)
  T <- S <- matrix(NA, nB, K)
  P <- S
  zlist <- vector("list", length(B))
  ##browser()
  for(i in seq_along(B)){
    is.batch <- batch(object)==B[i]
    j <- which(is.batch)
    y.batch <- y(object)[j]
    ys <- sample(y.batch, length(y.batch), replace=TRUE)
    ##ys <- ys + rnorm(length(ys), 0, 0.1)
    km <- kmeans(ys, K)
    thetas <- as.numeric(km$centers[, 1])
    ix <- order(thetas)
    T[i, ] <- thetas[ix]
    ##zs <- as.integer(kmeans(y.batch, centers=thetas[ix])$cluster)
    ps <- as.numeric(rdirichlet(1, alpha(hypp)))
    zs <- as.integer(simulateZ(length(y.batch), ps))
    ##zs <- as.integer(factor(km$cluster, levels=ix))
    ##ps <- (table(zs)/length(zs))
    if(length(ps) != K) browser()
    P[i, ] <- ps
    sds <- sapply(split(ys, zs), sd)
    if(length(sds) < K){
      sds <- sqrt(initializeSigma2(object))[1, ]
    }
    S[i, ] <- sds
    zlist[[i]] <- zs
  }
  zz <- unlist(zlist)
  p <- colMeans(P)
  mu(object) <- colMeans(T)
  tau2(object) <- colVars(T)
  p(object) <- p
  z(object) <- zz
  theta(object) <- T
  sigma2(object) <- S^2
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
}

.init_priors <- function(object){
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  K <- k(object)
  B <- nBatch(object)
  homo.del <- runif(B, -5, -0.5)
  hemi.del <- rnorm(B, -0.5, 0.5)
  diploid <- rnorm(B,  0, 0.5)
  dupl <- rnorm(B, 0.5, 0.5)
  twocopy.gain <- rnorm(B, 1, 0.5)
  thetas <- cbind(homo.del, hemi.del, diploid, dupl, twocopy.gain)
  if(K > 5){
    replace <- TRUE
  } else replace <- FALSE
  thetas <- thetas[, sample(1:5, K, replace=replace), drop=FALSE]
  for(i in 1:B){
    thetas[i, ] <- sort(thetas[i, ])
  }
  sigma2s <- initializeSigma2(object)
  sigma2(object) <- sigma2s
  mu(object) <- colMeans(thetas)
  tau2(object) <- colVars(thetas)
  p(object) <- as.numeric(rdirichlet(1, alpha(hypp)))
  if(length(y(object)) > 0){
    zz <- simulateZ(length(y(object)), p(object))
    z(object) <- as.integer(factor(zz))
  } else {
    zz <- as.integer(factor(levels=seq_len(k(hypp))))
    z(object) <- zz
  }
  if(length(y(object)) > 0){
    dataMean(object) <- computeMeans(object)
    dataPrec(object) <- computePrec(object)
    log_lik(object) <- computeLoglik(object)
    logPrior(object) <- computePrior(object)
  }
  probz(object) <- .computeProbZ(object)
  chains(object) <- McmcChains(object)
  object
}

.init_mclust <- function(object){
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  B <- uniqueBatch(object)
  K <- k(object)
  nB <- nBatch(object)
  P <- T <- S <- matrix(NA, nB, K)
  ylist <- split(y(object), batch(object))
  n.b <- elementNROWS(ylist)
  zlist <- vector("list", length(B))
  for(i in seq_along(B)){
    ys <- ylist[[i]]
    ys <- sample(ys, length(ys), replace=TRUE)
    mc <- Mclust(ys, G=K)
    params <- mc$parameters
    thetas <- params$mean
    sds <- sqrt(params$variance$sigmasq)
    ps <- params$pro
    T[i, ] <- thetas
    zs <- as.integer(simulateZ(length(ylist[[i]]), ps))
    ##zz <- simulateZ(length(y(object)), p(object))
    ##zz <- simulateZ(length(ys), p(object))
    ##zs <- as.integer(kmeans(ylist[[i]], centers=thetas)$cluster)
    ##zs <- as.integer(mc$classification)
    P[i, ] <- ps
    S[i, ] <- sds
    zlist[[i]] <- zs
  }
  zz <- unlist(zlist)
  ## *should do a weighted average since batches are unequally sized
  p <- colMeans(P)
  mu(object) <- colMeans(T)
  tau2(object) <- colVars(T)
  p(object) <- p
  z(object) <- zz
  theta(object) <- T
  sigma2(object) <- S^2
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
}

.init_batchmodel2 <- function(object){
  u <- runif(1, 0, 1)
  if(u < 1/3){
    object <- .init_priors(object)
  }
  if(u < 2/3 & u >= 1/3){
    object <- .init_kmeans(object)
  }
  if(u >= 2/3){
    object <- .init_mclust(object)
  }
  object
}

.init_batchmodel3 <- function(object){
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  NB <- nBatch(object)
  K <- k(object)
  ylist <- split(y(object), batch(object))
  thetas.list <- lapply(ylist, simulateThetas, K=K)
  thetas <- do.call(rbind, thetas.list)
  theta(object) <- thetas
  mu(object) <- colMeans(thetas)
  tau2(object) <- colVars(thetas)
  nB <- nBatch(object)
  p(object) <- as.numeric(rdirichlet(1, alpha(hypp)))
  if(length(y(object)) > 0){
    zz <- simulateZ(length(y(object)), p(object))
    z(object) <- as.integer(factor(zz))
  } else {
    zz <- as.integer(factor(levels=seq_len(k(hypp))))
    z(object) <- zz
  }
  zFreq(object) <- as.integer(table(zz))
  zlist <- split(zz, batch(object))
  sigma2s <- initializeSigma2(object)
  sigma2(object) <- sigma2s
  if(length(y(object)) > 0){
    dataMean(object) <- computeMeans(object)
    dataPrec(object) <- computePrec(object)
    log_lik(object) <- computeLoglik(object)
    logPrior(object) <- computePrior(object)
  }
  probz(object) <- .computeProbZ(object)
  chains(object) <- McmcChains(object)
  object
}

.init_batchmodel <- function(object){
  hypp <- hyperParams(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  mu(object) <- rnorm(k(object), mu.0(object), tau.0(object))
  tau2(object) <- 1/rgamma(k(object), shape=1/2*eta.0(object),
                           rate=1/2*eta.0(object)*m2.0(object))
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
}

setMethod("startingValues", "BatchModel", function(object){
  ##.init_batchmodel(object)
  if(length(y(object)) == 0) return(object)
  .init_batchmodel2(object)
})
