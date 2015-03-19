test_batchEasy <- function(){
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
  checkTrue(is.factor(batch(truth)))
  yy <- y(truth)
  checkIdentical(yy, yy[order(batch(truth))])
  mcmcp <- McmcParams(iter=50, burnin=0)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  checkIdentical(y(params), y(truth))
  checkIdentical(batch(params), batch(truth))

  model <- initializeBatchModel(params)
  checkIdentical(batch(model), batch(truth))
  checkIdentical(y(model), y(truth))

  model <- startAtTrueValues(model, truth)
  checkIdentical(theta(model), theta(truth))
  checkIdentical(sigma(model), sigma(truth))
  checkIdentical(p(model), p(truth))
  checkIdentical(z(model), z(truth))
  model <- posteriorSimulation(model, mcmcp)
  ## Parameters should not stray far from the true values after a
  ## small number of iterations
  checkEquals(theta(model), theta(truth), tolerance=0.1)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    ##abline(v=th, lwd=1)
    par(op)
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
  sds[,1] <- 0.3
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             p=c(1/10, 1/5, 1-0.1-0.2),
                             theta=means,
                             sds=sds)
  mcmcp <- McmcParams(iter=150, burnin=0)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
  model <- startAtTrueValues(model, truth)
  model <- posteriorSimulation(model, mcmcp)
  checkEquals(theta(model), theta(truth), tolerance=0.1) ## fails
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
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  modelk <- initializeBatchModel(params)
  modelk <- startAtTrueValues(modelk, truth)
  mmodel <- posteriorSimulation(modelk, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE)
    CNPBayes::plot(mmodel)
    par(op)
    mc <- mcmcChains(mmodel)
    plot.ts(sigma(mc), col="gray")
    plot.ts(theta(mc), col="gray")
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
}
