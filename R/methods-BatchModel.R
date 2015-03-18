setMethod("[", "BatchModel", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    y(object) <- y(object)[i]
    z(object) <- z(object)[i]
    batch(objct) <- batch(object)[i]
  }
  object
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


BatchModel <- function(data, k, batch, hypp){
  mcmc.chains <- McmcChains()
  B <- length(unique(batch))
  if(missing(hypp)) hypp <- HyperparametersBatch(k=k)
  new("BatchModel",
      hyperparams=hypp,
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
      loglik=numeric(1),
      mcmc.chains=mcmc.chains,
      batch=batch,
      hwe=numeric(),
      theta_order=numeric(B*k),
      m.y=numeric(1))
}

setMethod("bic", "BatchModel", function(object, ...){
  if(k(object) > 1){
    object <- updateWithPosteriorMeans(object)
  }
  ## K: number of free parameters to be estimated
  ##   - component and batch-specific parameters:  theta, sigma2  ( k(model) * nBatch(model))
  ##   - mixing probabilities: (k-1)*nBatch
  ##   - component-specific parameters: mu, tau2                 2 x k(model)
  ##   - length-one parameters: sigma2.0, nu.0                   +2
  K <- 2*k(object)*nBatch(object) + (k(object)-1) + 2*k(object) + 2
  n <- length(y(object))
  bicstat <- -2*logpotential(object) + K*(log(n) - log(2*pi))
  bicstat
})

setMethod("collapseBatch", "BatchModel", function(object){
  N <- choose(nBatch(object), 2)
  cond2 <- TRUE
  while(N > 1 && cond2){
    cat('.')
    B <- batch(object)
    batch(object) <- .collapseBatch(y(object), batch(object))
    cond2 <- !identical(B, batch(object))
    N <- nBatch(object)
  }
  makeUnique(batch(object))
})


setMethod("computeLoglik", c("BatchModel", "missing"), function(object, psi){
  .computeLoglikBatch(object)
})

logLikData <- function(object){
  lik <- rep(NA, length(y(object)))
  B <- batch(object)
  ub <- unique(B)
  mn <- theta(object)
  ss <- sigma(object)
  x <- y(object)
  tabz <- table(batch(object), z(object))
  tabz <- tabz/rowSums(tabz)
  tabz <- tabz[uniqueBatch(object), , drop=FALSE]
  P <- tabz
  lik <- matrix(NA, length(x), k(object))
  browser()
  for(j in seq_len(k(object))){
    for(b in ub){
      p.x <- P[b, j] * dnorm(x[B==b], mean=mn[b, j], sd=ss[b, j])
      lik[B==b, j] <- p.x
    }
  }
  lik <- rowSums(lik)
  sum(log(lik))
}

logLikPhi <- function(object){
  thetas <- theta(object)
  mus <- mu(object)
  tau2s <- tau2(object)
  sigma2s <- sigma2(object)
  p.theta <- matrix(NA, nrow(thetas), ncol(thetas))
  for(j in seq_len(ncol(p.theta))){
    p.theta[, j] <- dnorm(thetas[, j], mus[j], sqrt(tau2s[j]))
  }
  p.sigma2 <- dgamma(1/sigma2s, shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  sum(log(p.theta)) + sum(log(p.sigma2))
}

.computeLoglikBatch <- function(object){
  ll.data <- logLikData(object)
  ll.phi <- logLikPhi(object)
  ll.data + ll.phi
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


##setMethod("computeLoglik", c("BatchModel", "McmcChains"), function(object, psi){
##  lik <- rep(NA, length(y(object)))
##  B <- batch(object)
##  ub <- unique(B)
##  ## psi is fixed, but not z
##  mn <- theta(psi)
##  ss <- sigma(psi)
##  pp <- p(psi)
##  yy <- y(object)
##  zz <- as.integer(z(object))
##  for(b in ub){
##    i <- B==b
##    m <- mn[b, ]
##    s <- ss[b, ]
##    z <- zz[i]
##    y <- yy[i]
##    lik[i] <- pp[z]*dnorm(y, m[z], s[z])
##  }
##  log(lik)
##})


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


.computeModesBatch <- function(object){
  i <- argMax(object)
  mc <- mcmcChains(object)
  B <- nBatch(object)
  K <- k(object)
  thetamax <- matrix(theta(mc)[i, ], B, K)
  sigma2max <- matrix(sigma2(mc)[i, ], B, K)
  pmax <- p(mc)[i, ]
  mumax <- mu(mc)[i, ]
  tau2max <- tau2(mc)[i,]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mumax,
                tau2=tau2max,
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i])
  modes
}

setMethod("computeModes", "BatchModel", function(object){
  .computeModesBatch(object)
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
##  hypp <- hyperParams(object)
##  K <- k(hypp)
##  yy <- y(object)
##  zz <- z(object)
##  thetas <- theta(object)
##  sigma2s <- sigma2(object)
##  mus <- mu(object)
##  tau2s <- tau2(object)
##  pp <- p(object)
##
##  p.mu <- dnorm(mus, mu.0(hypp), sqrt(tau2.0(hypp)))
##  p.sigma2.0 <- dgamma(sigma2.0(object), shape=a(hypp), rate=b(hypp))
##  p.nu.0 <- dgeom(as.integer(nu.0(object)), betas(hypp))
##
##  p.theta <- matrix(NA, nrow(thetas), ncol(thetas))
##  for(j in seq_len(ncol(p.theta))){
##    p.theta[, j] <- dnorm(thetas[, j], mus[j], sqrt(tau2s[j]))
##  }
##  p.sigma2 <- dgamma(1/sigma2s, shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
##  pot <- list()
##  ##
##  ##
##  ##
##  sigmas <- sqrt(sigma2s)
##  p.y <- rep(NA, length(yy))
##  ylist <- split(yy, batch(object))
##  ##rownames(thetas) <- rownames(sigmas) <- uniqueBatch(object)
##  for(b in uniqueBatch(object)){
##    dat <- ylist[[b]]
##    pot <- matrix(NA, length(dat), K)
##    for(i in seq_len(K)){
##      pot[, i] <- pp[i]*dnorm(dat, thetas[b, i], sigmas[b, i])
##    }
##    p.y[batch(object)==b] <- rowSums(pot)
##  }
##  total_pot <- sum(log(p.y)) + sum(log(p.theta)) + sum(log(p.sigma2)) + sum(log(p.mu)) + log(p.nu.0) + log(p.sigma2.0)
  ##  total_pot
  computeLogLikxPrior(object)
})

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



setMethod("initializeSigma2.0", "BatchModel", function(object){
  hypp <- hyperParams(object)
  sum(alpha(hypp)*colMeans(sigma2(object)))/sum(alpha(hypp))
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
  logpotential(mcmc)[s] <- computeLogLikxPrior(object)
  logLik(mcmc)[s] <- computeLoglik(object)
  mcmcChains(object) <- mcmc
  object
})


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



#' @export
setMethod("mu", "BatchModel", function(object) object@mu)

setReplaceMethod("mu", "BatchModel", function(object, value){
  object@mu <- value
  object
})

nBatch <- function(object) length(uniqueBatch(object))

setReplaceMethod("p", "BatchModel", function(object, value){
  object@pi <- value
  object
})

setMethod("pMean", "BatchModel", function(object) {
  mns <- colMeans(pic(object))
  mns
})

setMethod("showMeans", "BatchModel", function(object){
  thetas <- round(theta(object), 2)
##  dimnames(thetas) <- list(paste0("batch", uniqueBatch(object)),
##                           paste0("component", seq_len(k(object))))
  ##
  mns <- c("\n", paste0(t(cbind(thetas, "\n")), collapse="\t"))
  mns <- paste0("\t", mns[2])
  mns <- paste0("\n", mns[1])
  mns
})

setMethod("showSigmas", "BatchModel", function(object){
  sigmas <- round(sqrt(sigma2(object)), 2)
##  dimnames(sigmas) <- list(paste0("batch", uniqueBatch(object)),
##                           paste0("component", seq_len(k(object))))
  sigmas <- c("\n", paste0(t(cbind(sigmas, "\n")), collapse="\t"))
  sigmas <- paste0("\t", sigmas[2])
  sigmas <- paste0("\n", sigmas[1])
  sigmas
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

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})


setMethod("sigmaMean", "BatchModel", function(object) {
  mns <- colMeans(sigmac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

#' @export
setMethod("tau2", "BatchModel", function(object) object@tau2)

setReplaceMethod("tau2", "BatchModel", function(object, value){
  object@tau2 <- value
  object
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



setMethod("thetaMean", "BatchModel", function(object) {
  mns <- colMeans(thetac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})


setMethod("show", "BatchModel", function(object){
  ##callNextMethod()
  cat("An object of class 'BatchModel'\n")
  cat("     n. obs      :", length(y(object)), "\n")
  cat("     n. batches  :", nBatch(object), "\n")
  cat("     k           :", k(object), "\n")
  cat("     nobs/batch  :", table(batch(object)), "\n")
  cat("     loglik (s)  :", round(logLik(object), 1), "\n")
  cat("     loglik x prior:", round(logpotential(object), 1), "\n")
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

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})

#' @export
uniqueBatch <- function(object) unique(batch(object))

##
## Multiple batches, but only 1 component
##
UnivariateBatchModel <- function(data, k=1, batch, hypp){
  mcmc.chains <- McmcChains()
  B <- length(unique(batch))
  if(missing(hypp)) hypp <- HyperparametersBatch(k=1)
  new("UnivariateBatchModel",
      hyperparams=hypp,
      theta=matrix(NA, B, 1),
      sigma2=matrix(NA, B, 1),
      mu=numeric(k),
      tau2=numeric(k),
      nu.0=numeric(1),
      sigma2.0=numeric(1),
      pi=numeric(k),
      ##pi=matrix(NA, B, 1),
      data=data,
      data.mean=matrix(NA, B, 1),
      data.prec=matrix(NA, B, 1),
      z=factor(numeric(length(data))),
      probz=matrix(1, length(data), 1),
      logpotential=numeric(1),
      mcmc.chains=mcmc.chains,
      batch=batch,
      hwe=numeric(),
      m.y=numeric(1))
}

##
###' @export
##setMethod("mu", "BatchModel", function(object) object@mu)
##
###' @export
##setMethod("tau2", "BatchModel", function(object) object@tau2)
##
##setReplaceMethod("tau2", "BatchModel", function(object, value){
##  object@tau2 <- value
##  object
##})
##
##setReplaceMethod("mu", "BatchModel", function(object, value){
##  object@mu <- value
##  object
##})
##
##
##
##setReplaceMethod("theta", "BatchModel", function(object, value){
##  rownames(value) <- uniqueBatch(object)
##  object@theta <- value
##  object
##})
##
##
##
##setReplaceMethod("sigma2", "BatchModel", function(object, value){
##  rownames(value) <- uniqueBatch(object)
##  object@sigma2 <- value
##  object
##})
##
##setMethod("sigma2", "BatchModel", function(object) {
##  s2 <- object@sigma2
##  s2 <- matrix(s2, nBatch(object), k(object))
##  rownames(s2) <- uniqueBatch(object)
##  s2
##})
##
##setMethod("reorderComponents", "BatchModel", function(object){
##  thetas <- theta(object)
##  sigma2s <- sigma2(object)
##  for(i in seq_len(nrow(thetas))){
##    x <- thetas[i, ]
##    j <- order(x)
##    thetas[i, ] <- x[j]
##    sigma2s[i, ] <- sigma2s[i, j]
##    if(i == 1){
##      p(object) <- p(object)[j]
##    }
##    dataMean(object)[i, ] <- dataMean(object)[i, j]
##    dataPrec(object)[i, ] <- dataPrec(object)[i, j]
##  }
##  theta(object) <- thetas
##  sigma2(object) <- sigma2s
##  ##ix <- order(theta(object)[1, ])
##  ##theta(object) <- theta(object)[, ix]
##  ##sigma2(object) <- sigma2(object)[, ix]
##  ##p(object) <- p(object)[ix]
##  zz <- z(object)
##  zz <- factor(as.integer(factor(zz, levels=j)), levels=seq_len(k(object)))
##  z(object) <- zz
##  ##  dataMean(object) <- dataMean(object)[, ix]
##  ##  dataPrec(object) <- dataPrec(object)[, ix]
##  object
##})
##
##
##

##
##
##setMethod("sigmaMean", "BatchModel", function(object) {
##  mns <- colMeans(sigmac(object))
##  mns <- matrix(mns, nBatch(object), k(object))
##  rownames(mns) <- uniqueBatch(object)
##  mns
##})
##
##setMethod("thetaMean", "BatchModel", function(object) {
##  mns <- colMeans(thetac(object))
##  mns <- matrix(mns, nBatch(object), k(object))
##  rownames(mns) <- uniqueBatch(object)
##  mns
##})
##
##setMethod("pMean", "BatchModel", function(object) {
##  mns <- colMeans(pic(object))
##  mns
##})

#' @export
batchExperiment <- function(object, outdir, marginaly=TRUE, test=FALSE){
  B <- getFiles(outdir, rownames(object), "batch")
  ## shouldn't batch files be named with respect to the rownames?
  batch.files <- paste0(dirname(model(B)), "/", rownames(object), "_batch.rds")
  mcmcp.list <- mcmcpList()
  hp.list <- HyperParameterList(K=1:4, HyperparametersBatch(tau2.0=1000))
  J <- seq_len(nrow(object)); j <- NULL
  x <- foreach(j = J, .packages=c("CNPBayes", "foreach")) %dopar% {
    cn <- copyNumber(object)[j, ]
    notna <- !is.na(cn)
    bt <- saveBatch(object[j, notna], batch.files[j])
    models <- batchModel(B[j], data=cn[notna],
                         hyp.list=hp.list,
                         mcmcp.list=mcmcp.list,
                         batch=bt, save.it=TRUE, test=test,
                         marginaly=marginaly)
  }
  TRUE
}

#' @export
marginalExperiment <- function(object, outdir, test=FALSE){
  M <- getFiles(outdir, rownames(object), "marginal")
  mcmcp.list <- mcmcpList()
  hp.list <- HyperParameterList(K=1:4, HyperparametersMarginal(tau2.0=1000))
  cn <- copyNumber(object)
  J <- seq_len(nrow(object)); j <- NULL
  x <- foreach(j = J, .packages=c("CNPBayes", "foreach")) %dopar% {
    models <- marginalModel(M[j], data=cn[j, ],
                            hyp.list=hp.list,
                            mcmcp.list=mcmcp.list, save.it=TRUE,
                            test=test)
  }
  TRUE
}
