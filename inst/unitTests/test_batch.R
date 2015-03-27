test_batchEasy <- function(){
  library(oligoClasses)
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             theta=means,
                             sds=sds,
                             p=c(1/5, 1/3, 1-1/3-1/5))
  ##checkTrue(oligoClasses::batch(truth))
  yy <- y(truth)
  checkIdentical(yy, yy[order(batch(truth))])
  mcmcp <- McmcParams(iter=50, burnin=0)
  set.seed(123)
  model <- BatchModel(y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp)
  model <- startAtTrueValues(model, truth)
  checkIdentical(batch(model), batch(truth))
  checkIdentical(y(model), y(truth))
  checkIdentical(theta(model), theta(truth))
  checkIdentical(sigma(model), sigma(truth))
  checkIdentical(p(model), p(truth))
  checkIdentical(z(model), z(truth))
  iter(model) <- 50
  model2 <- CNPBayes:::.runMcmc(model)
  checkEquals(theta(model2), theta(truth), tolerance=0.1)

  if(FALSE){
    iter(model) <- 500
    model2 <- CNPBayes:::.runMcmc(model)
    checkEquals(theta(model2), theta(truth), tolerance=0.1)
  }
  ##paramUpdates(model)[c("nu.0", "sigma2.0")] <- 0L
  ##paramUpdates(model)[c("sigma2.0")] <- 0L
  if(FALSE){
    model <- .Call("mcmc_batch", model, mcmcParams(model))
    set.seed(123)
    .Call("update_sigma20_batch", model)
    set.seed(123)
    .updateSigma2.0Batch(model2)
    ## problem lies with nu.0 or sigma2.0...
  }
  eta.0 <- 100
  m2.0 <- 1/10
  a <- 0.5*eta.0
  b <- 0.5*eta.0*m2.0
  x <- rgamma(1000, a, rate=1/b)
  ix <- 1/x
  hist(sqrt(1/x), breaks=100)
  s <- mean(sqrt(1/x))
  theta <- rnorm(1000, -1, s)
  hist(theta, breaks=100, xlim=c(-2, 1.5))

  eta.0 <- 1800
  m2.0 <- 1/60
  a <- 0.5*eta.0
  b <- 0.5*eta.0*m2.0
  x <- rgamma(1000, a, rate=1/b)
  ix <- 1/x
  hist(sqrt(1/x), breaks=100)
  s <- mean(sqrt(1/x))
  theta <- rnorm(1000, -1, s)
  hist(theta, breaks=100, xlim=c(-2, 1.5))

  mcmcp <- McmcParams(iter=2000, burnin=0)
  model <- BatchModel(y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp,
                      HyperparametersBatch(m2.0=1/60, eta.0=1800, k=3))
  model <- posteriorSimulation(model)
  i <- order(theta(model)[1, ])
  checkEquals(theta(model)[, i], theta(truth), tolerance=0.1)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
    par(mfrow=c(1,3))
    tracePlot(model, "theta", col=1:3)
    plot.ts(muc(model), col=1:3)
  }
}

test_batch_moderate <- function(){
  library(oligoClasses)
  set.seed(100)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95,
                    -0.41, -0.4, -.395,
                    -0.1, 0, 0.05), nbatch, k, byrow=FALSE)
  sds <- matrix(0.15, nbatch, k)
  sds[, 1] <- 0.3
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             p=c(1/10, 1/5, 1-0.1-0.2),
                             theta=means,
                             sds=sds)
  mcmcp <- McmcParams(iter=150, burnin=0)
  hypp <- HyperparametersBatch(m2.0=1/60, eta.0=1800, k=3)
  model <- BatchModel(data=y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp,
                      hypp=hypp)
  model <- startAtTrueValues(model, truth)
  model <- posteriorSimulation(model, mcmcp)
  i <- order(theta(model)[1, ])
  checkEquals(theta(model)[, i], theta(truth), tolerance=0.1)
  ## random starts
  mcmcp <- McmcParams(iter=100, burnin=100, nStarts=20)
  model <- BatchModel(data=y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp,
                      hypp=hypp)
  model <- posteriorSimulation(model)
  i <- order(theta(model)[1, ])
  checkEquals(theta(model)[, i], theta(truth), tolerance=0.1)
  checkEquals(mu(model)[i], mu(truth), tolerance=0.15)
  if(FALSE){
    zz <- as.integer(z(truth))
    ps <- c(mean(zz==1), mean(zz==2), mean(zz==3))
    modelk <- model
    plot.ts(pic(modelk), plot.type="single")
    abline(h=p(truth))
    par(mfrow=c(1,3))
    tracePlot(modelk, "theta", col=1:3)
    abline(h=theta(truth))
    plot.ts(sigmac(modelk), plot.type="single")
    abline(h=sigma(truth))
    op <- par(mfrow=c(1,2),las=1)
    ##trace(.plotBatch, browser)
    ##.plotBatch(truth, use.current=TRUE)
    CNPBayes::plot(truth, use.current=TRUE)
    CNPBayes::plot(modelk, use.current=TRUE)
    par(op)
    mc <- mcmcChains(modelk)
    plot.ts(sigma(mc), col="gray")
    plot.ts(theta(mc), col="gray")
    plot.ts(p(mc), col="gray")
  }
}

hardTruth <- function(prop_comp1=0.005, s=0.3){
  set.seed(1234)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.9, -2, -1.85,
                    -0.45, -0.4, -0.35,
                    -0.1, 0, -0.05), nbatch, k, byrow=FALSE)
  sds <- matrix(s, nbatch, k)
  p1 <- prop_comp1
  p2 <- 20*p1
  p3 <- 1-p1-p2
  N <- 8e3
  truth <- simulateBatchData(N,
                             theta=means,
                             sds=sds,
                             batch=rep(letters[1:3], length.out=N),
                             p=c(p1, p2, p3))
}

test_hard3 <- function(){
  ##
  ## Repeat above, but with smaller simulated variance (sd=0.2)
  ##
  library(oligoClasses)
  library(GenomicRanges)
  truth <- hardTruth(0.005, s=0.1)
  table(z(truth), batch(truth))
  if(FALSE) CNPBayes::plot(truth, use.current=TRUE)
  se <- as(truth, "SummarizedExperiment")
  if(FALSE) hist(oligoClasses::copyNumber(se), breaks=1000, col="gray", border="gray")
  ##
  ## Use defaults
  ##
  mcmcp <- McmcParams(iter=100, burnin=0)
  modelk <- BatchModel(data=y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp,
                       HyperparametersBatch(k=3, m2.0=1/60, eta.0=1800))
  modelk <- startAtTrueValues(modelk, truth)
  mmodel <- posteriorSimulation(modelk, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE)
    ##trace(.plotBatch, browser)
    .plotBatch(mmodel)
    par(op)
    i <- order(theta(mmodel)[1, ])
    theta(mmodel)[, i]
    theta(truth)
    sigma(model)[, i]
    sigma(truth)
    mc <- mcmcChains(mmodel)
    par(mfrow=c(1,3))
    tracePlot(mmodel, "sigma", col=1:3)
    tracePlot(mmodel, "theta", col=1:3)
    plot.ts(sigma(mc), col="gray", )
    plot.ts(theta(mc), col="gray")
    plot.ts(mu(mc), plot.type="single", col=1:3)
    cbind(colMeans(sigma(mc)), as.numeric(sigma(truth)))
  }
  ##mc <- mcmcChains(mmodel)
  pmns <- thetaMean(mmodel)
  j <- order(pmns[1,])
  ps <- sigmaMean(mmodel)[, j]
  pmix <- pMean(mmodel)[j]
  checkEquals(pmns[, j], theta(truth), tolerance=0.04)
  checkEquals(ps, sigma(truth), tolerance=0.15)
  checkEquals(pmix, p(truth), tolerance=0.04)
  ##
  ## With multiple starts we can find the mode
  ##
  mcmcp <- McmcParams(iter=200, burnin=100, nStarts=20)
  modelk <- BatchModel(data=y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp,
                       HyperparametersBatch(k=3, m2.0=1/60, eta.0=1800, tau2.0=1000))
  mmodel <- posteriorSimulation(modelk, mcmcp)
  pmns <- thetaMean(mmodel)
  j <- order(pmns[1,])
  ps <- sigmaMean(mmodel)[, j]
  pmix <- pMean(mmodel)[j]
  checkEquals(pmns[, j], theta(truth), tolerance=0.04)
  checkEquals(ps, sigma(truth), tolerance=0.15)
  checkEquals(pmix, p(truth), tolerance=0.04)
}
