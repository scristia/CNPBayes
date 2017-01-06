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

tau2Hyperparams <- function(thetas){
  tau2hat <- var(thetas)
  itau2 <- 1/tau2hat
  qInverseTau2(mn=itau2, sd=0.001)
}

.initialize_z_singlebatch <- function(y, centers){
  kmeans(y, centers=centers)$cluster
}


.init_sb2 <- function(object){
  hypp <- hyperParams(object)
  if(length(y(object))==0) return(object)
  sigma2.0(object) <- rgamma(1, a(hypp), b(hypp))
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  ys <- y(object)
  ys <- sample(ys, length(ys), replace=TRUE)
  K <- k(object)
  mc <- tryCatch(Mclust(ys, G=K), warning=function(w) NULL)
  if(is.null(mc)){
    stop("Trouble selecting starting values with mclust.")
  }
  params <- mc$parameters
  thetas <- params$mean
  if(K > 1){
    tau2.hp <- tau2Hyperparams(thetas)
    eta.0(hypp) <- tau2.hp$eta.0
    m2.0(hypp) <- tau2.hp$m2.0
    hyperParams(object) <- hypp
  }
  s2s <- params$variance$sigmasq
  ps <- params$pro
  ##zs.bootstrap <- as.integer(mc$classification)
  ##splits <- split(ys, zs.bootstrap)
  ##mins <- sort(sapply(splits, min))
  yy <- y(object)
  if(K > 1){
    zz <- .initialize_z_singlebatch(yy, thetas)
  } else {
    zz <- rep(1L, length(yy))
  }
  ##zz <- as.integer(simulateZ(length(y(object)), p(object)))
  if(FALSE){
    cnt <- 1
    while (length(table(zz)) < K) {
      if (cnt > 10) {
        stop("Too few observations or too many components.")
      }
      p(object) <- as.numeric(rdirichlet(1, alpha(hypp))) ## rows are
      zz <- as.integer(simulateZ(length(y(object)), p(object)))
      cnt <- cnt + 1
    }
  }
  mu(object) <- mean(thetas)
  if(K > 1){
    tau2(object) <- var(thetas)
  } else tau2(object) <- var(y(object))*10
  p(object) <- ps
  z(object) <- zz
  theta(object) <- thetas
  if(length(s2s) != length(thetas)){
    s2s <- rep(s2s[1], length(thetas))
  }
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
  counter <- 0; model <- NULL
  while(counter < 3 && is.null(model)){
    ## when k is too large, kmeans for finding cluster centers will not work
    ## -- try a couple of times before quitting
    model <- tryCatch(.init_sb2(object), error=function(e) NULL)
    counter <- counter + 1
  }
  K <- k(object)
  while(is.null(model) & K > 1){
    K <- K - 1
    k(object) <- K
    model <- tryCatch(.init_sb2(object), error=function(e) NULL)
  }
  if(is.null(model)) stop("No good initial values identified")
  model
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
  for(i in seq_along(B)){
    is.batch <- batch(object)==B[i]
    j <- which(is.batch)
    y.batch <- y(object)[j]
    ys <- sample(y.batch, length(y.batch), replace=TRUE)
    ##ys <- ys + rnorm(length(ys), 0, 0.1)
    if(K > 1){
      km <- kmeans(ys, K)
      thetas <- as.numeric(km$centers[, 1])
      ##zs <- as.integer(kmeans(y.batch, centers=thetas[ix])$cluster)
    } else {
      thetas <- mean(ys)
    }
    ix <- order(thetas)
    T[i, ] <- thetas[ix]
    ps <- as.numeric(rdirichlet(1, alpha(hypp)))
    if(K > 1){
      zs <- kmeans(y.batch, centers=thetas[ix])$cluster
    } else zs <- rep(1L, length(y.batch))
    ##zs <- as.integer(simulateZ(length(y.batch), ps))
    ##zs <- as.integer(factor(km$cluster, levels=ix))
    ##ps <- (table(zs)/length(zs))
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
  if(K > 1){
    tau2.hypp <- tau2HyperparamsBatch(T)
    eta.0(hypp) <- tau2.hypp$eta.0
    m2.0(hypp) <- tau2.hypp$m2.0
    hyperParams(object) <- hypp
  }
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
  theta(object) <- thetas
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

tau2HyperparamsBatch <- function(thetas){
  tau2s <- colVars(thetas)
  itau2s <- 1/tau2s
  med <- median(itau2s)
  qInverseTau2(mn=med, sd=med/100)
}

##tau2HyperparamsBatch <- function(tau2s){
##  mn <- 1/var(tau2s)
##  qInverseTau2(mn=mn, sd=0.001)
##}

.initialize_z <- function(y, centers){
  if(length(centers) > 1){
    z <- tryCatch(kmeans(y, centers=centers)$cluster, error=function(e) NULL)
    if(is.null(z)) stop("Empty clusters because k is too large")
  } else {
    z <- rep(1L, length(y))
  }
  z
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
    mc <- tryCatch(Mclust(ys, G=K), warning=function(w) NULL)
    if(is.null(mc)) stop("Error initializing values for batch model")
    params <- mc$parameters
    thetas <- params$mean
    sds <- sqrt(params$variance$sigmasq)
    ps <- params$pro
    T[i, ] <- thetas
    zlist[[i]] <- .initialize_z(ylist[[i]], thetas)
    P[i, ] <- ps
    S[i, ] <- sds
  }
  ##
  ## since batches are not ordered, unlist(z) still matches the original y's
  ##
  zz <- unlist(zlist)
  ## *should do a weighted average since batches are unequally sized
  p <- colMeans(P)
  mu(object) <- colMeans(T)
  tau2(object) <- colVars(T)
  if(K > 1){
    tau2.hypp <- tau2HyperparamsBatch(T)
    eta.0(hypp) <- tau2.hypp$eta.0
    m2.0(hypp) <- tau2.hypp$m2.0
    hyperParams(object) <- hypp
  }
  p(object) <- p
  z(object) <- zz
  ## since batches are required to be in batch order, this is not necessary
  ##y(object) <- unlist(ylist)
  ##batch(object) <- unlist(split(batch(object), batch(object)))
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
  if(u <= 0.25){
    ##    object <- .init_priors(object)
    ##    if(!is.finite(log_lik(object)) || is.na(log_lik(object))){
    ##      object <- .init_mclust(object)
    ##    }
    ##  }
    ##  if(u < 0.5 & u >= 0.25){
    counter <- 0; model <- NULL; not.valid <- FALSE
    while(counter < 3 && is.null(model) || not.valid){
      model <- tryCatch(.init_kmeans(object), error=function(e) NULL)
      if(!is.null(model)){
        not.valid <- !is.finite(log_lik(model)) || is.na(log_lik(model))
      }
      counter <- counter+1
    }
    if(is.null(model)) stop("No good initial values identified")
  }
  if(u > 0.25){
    counter <- 0; model <- NULL; not.valid <- FALSE
    while(counter < 3 && is.null(model) || not.valid){
      model <- tryCatch(.init_mclust(object), error=function(e) NULL)
      if(!is.null(model)){
        not.valid <- !is.finite(log_lik(model)) || is.na(log_lik(model))
      }
      counter <- counter+1
    }
    if(is.null(model)) stop("No good initial values identified")
  }
  model
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
