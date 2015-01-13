test_batchEasy <- function(){
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             .batch=rep(letters[1:3], length.out=2500),
                             .k=3,
                             .alpha=rep(1, k),
                             means=means,
                             sds=sds)

  ##
  ## Start from true values
  ##
  if(FALSE){
    model <- truth
    if(FALSE) CNPBayes::plot(truth, use.current=TRUE)
    mcmcp <- McmcParams(iter=1500, burnin=100)
    mcmcChains(model) <- McmcChains(model, mcmcp)
    model <- posteriorSimulation(model, mcmcp)
    if(FALSE){
      op <- par(mfrow=c(1,2),las=1)
      plot(truth, use.current=T)
      plot(model)
      par(op)
    }
    mc <- mcmcChains(model)
    pmns <- colMeans(theta(mc))
    checkEquals(pmns, as.numeric(theta(truth)), tolerance=0.03)

    ps <- colMeans(sigma(mc))
    checkEquals(ps, as.numeric(sigma(truth)), tolerance=0.18)

    pmix <- p(truth)
    pm_pmix <- colMeans(p(mc))
    checkEquals(pmix, pm_pmix, tolerance=0.02)
  }
  ##
  ## Start using defaults
  ##
  mcmcp <- McmcParams(iter=1000, burnin=100)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=rep(letters[1:3], length.out=2500),
                        mcmc.params=mcmcp)
  model2 <- initializeModel(params)
  model2 <- posteriorSimulation(model2, mcmcp)
  mc <- mcmcChains(model2)
  pmns <- colMeans(theta(mc))
  checkEquals(pmns, as.numeric(theta(truth)), tolerance=0.02)

  ps <- colMeans(sigma(mc))
  checkEquals(ps, as.numeric(sigma(truth)), tolerance=0.1)

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkEquals(pmix, pm_pmix, tolerance=0.02)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model2)
    par(op)
  }
}

test_selectK_batch_easy <- function(){
  ## Need to replace setParallelization with source code
  ## -investigate BiocParallel
  library(foreach)
  library(GenomicRanges)
  set.seed(1000)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             .batch=rep(letters[1:3], length.out=2500),
                             .k=3,
                             .alpha=rep(1, k),
                             means=means,
                             sds=sds)
  se <- as(truth, "SummarizedExperiment")

  ##
  ## Evaluate at different K
  ##
  mcmcp <- McmcParams(iter=1000, burnin=300)
  ##mcmcp <- McmcParams(iter=5, burnin=3)
  mmodels <- fitMixtureModels(se, mcmcp, K=1:5)
  bicstat <- sapply(mmodels, bic)
  ##
  ## 4-component has higher bic, yet the max a post. estimate has only 3 components
  ## --> inference from map would be the same
  ##
  checkTrue(which.min(bicstat) %in% c(3L, 4L))
  if(FALSE) {
    par(mfrow=c(1,2))
    CNPBayes::plot(truth, use.current=TRUE)
    CNPBayes::plot(mmodels[[4]])
    par(op)
  }
}

test_selectK_batch_moderate <- function(){
  set.seed(100)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95,
                    -0.41, -0.4, -.395,
                    -0.1, 0, 0.05), nbatch, k, byrow=FALSE)
  sds <- matrix(0.15, nbatch, k)
  sds[,1] <- 0.3
  truth <- simulateBatchData(N=2500,
                             means=means,
                             sds=sds,
                             .k=3,
                             .batch=rep(letters[1:3], length.out=2500),
                             .alpha=c(100, 200, 400))
  se <- as(truth, "SummarizedExperiment")
  mcmcp <- McmcParams(iter=1000, burnin=1000)
  params <- ModelParams("batch", y=y(truth), k=3, batch=batch(truth),
                        mcmc.params=mcmcp)
  modelk <- initializeModel(params)
  modelk <- posteriorSimulation(modelk, mcmcp)

  mc <- mcmcChains(modelk)
  pmns <- colMeans(theta(mc))
  checkEquals(pmns, as.numeric(theta(truth)), tolerance=0.03)
  ps <- colMeans(sigma(mc))
  checkEquals(ps, as.numeric(sigma(truth)), tolerance=0.1)
  pmix <- colMeans(p(mc))
  checkEquals(pmix, p(truth), tolerance=0.05)
  if(FALSE){
    zz <- as.integer(z(truth))
    ps <- c(mean(zz==1), mean(zz==2), mean(zz==3))

    plot.ts(pic(modelk), plot.type="single")
    abline(h=p(truth))
    plot.ts(thetac(modelk), plot.type="single")
    abline(h=theta(truth))

    plot.ts(sigmac(modelk), plot.type="single")
    abline(h=sigma(truth))
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE)
    CNPBayes::plot(modelk)
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
  ncomp1 <- ceiling(prop_comp1*1500 )
  ncomp2 <- 500-ncomp1
  ncomp3 <- 1000
  ##trace(simulateBatchData, browser)
  N <- 8e3
  truth <- simulateBatchData(N,
                             means=means,
                             sds=sds,
                             .batch=rep(letters[1:3], length.out=N),
                             .alpha=c(ncomp1, ncomp2, ncomp3))
}


.test_hard1 <- function(){
  truth <- hardTruth(0.005, s=0.3)
  table(z(truth), batch(truth))
  if(FALSE) plot(truth, use.current=TRUE)
  ##
  ## Use defaults
  ##
  mcmcp <- McmcParams(iter=2000, burnin=1000, thin=2)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeModel(params)
  batch(model) <- collapseBatch(model)
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE, xlim=c(-2,1))
    CNPBayes::plot(model, xlim=c(-2,1))
    par(op)
    ## we get the mixture probabilities backwards for comp 2 and 3
  }
  mc <- mcmcChains(model)
  cbind(colMeans(theta(mc)), as.numeric(theta(truth)))
  cbind(colMeans(sigma(mc)), as.numeric(sigma(truth)))
  cbind(colMeans(p(mc)), p(truth))
  ## this fails
  checkEquals(colMeans(p(mc)), p(truth), tolerance=0.1)
}

.test_hard2 <- function(){
  ##
  ## Repeat above, but with smaller simulated variance (sd=0.2)
  ##
  truth <- hardTruth(0.005, s=0.2)
  table(z(truth), batch(truth))
  if(FALSE) CNPBayes::plot(truth, use.current=TRUE)

  ## start at true values
  mcmcp <- McmcParams(iter=1000, burnin=1000, thin=1)
  model <- truth
  model <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model)
  checkEquals(colMeans(theta(mc)), as.numeric(theta(truth)), tolerance=0.06)
  checkEquals(colMeans(sigma(mc)), as.numeric(sigma(truth)), tolerance=0.05)
  checkEquals(colMeans(p(mc)), as.numeric(p(truth)), tolerance=0.17)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE, xlim=c(-2,1))
    CNPBayes::plot(model, xlim=c(-2,1))
    par(op)
    plot.ts(logpotential(mc), col="gray")
  }
  ##
  ## Default starts
  ##
  mcmcp <- McmcParams(iter=1000, burnin=1000, thin=1)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeModel(params)
  batch(model) <- collapseBatch(model)
  model <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE, xlim=c(-2,1))
    CNPBayes::plot(model, xlim=c(-2,1))
    par(op)
    plot.ts(p(mc), col="gray")
    plot.ts(sigma(mc), col="gray")
    plot.ts(mu(mc), col="gray")
    checkEquals(colMeans(theta(mc)), as.numeric(theta(truth)), tolerance=0.06)
    checkEquals(colMeans(sigma(mc)), as.numeric(sigma(truth)), tolerance=0.11)
    checkEquals(colMeans(p(mc)), as.numeric(p(truth)), tolerance=0.1)

    probs <- apply(probz(model), 1, max)
    plot(probs~jitter(as.integer(z(truth)), amount=0.1))
  }
  ## Fails
  checkEquals(colMeans(theta(mc)), as.numeric(theta(truth)), tolerance=0.1)
  ## Fails
  checkEquals(colMeans(sigma(mc)), as.numeric(sigma(truth)), tolerance=0.1)
  ## Fails
  checkEquals(colMeans(p(mc)), as.numeric(p(truth)), tolerance=0.1)
}

test_hard3 <- function(){
  ##
  ## Repeat above, but with smaller simulated variance (sd=0.2)
  ##
  truth <- hardTruth(0.005, s=0.1)
  table(z(truth), batch(truth))
  if(FALSE) CNPBayes::plot(truth, use.current=TRUE)
  se <- as(truth, "SummarizedExperiment")
  hist(oligoClasses::copyNumber(se), breaks=1000, col="gray", border="gray")
  ##
  ## Use defaults
  ##
  mcmcp <- McmcParams(iter=1000, burnin=1000, thin=1)
  ##mcmcp <- McmcParams(iter=1000, burnin=0, thin=1)
  mmodel <- fitMixtureModels(se, mcmcp, K=3)[[1]]
  ##mmodel2 <- posteriorSimulation(mmodel, mcmcp)
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
  mc <- mcmcChains(mmodel)
  checkEquals(colMeans(theta(mc)), as.numeric(theta(truth)), tolerance=0.02)
  checkEquals(colMeans(sigma(mc)), as.numeric(sigma(truth)), tolerance=0.05)
  checkEquals(colMeans(p(mc)), as.numeric(p(truth)), tolerance=0.04)
}

.test_hard4 <- function(){
  truth <- hardTruth(0.01, s=0.15)
  mcmcp <- McmcParams(iter=2000, burnin=1000, thin=2)
  table(z(truth), batch(truth))
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeModel(params)
  batch(model) <- collapseBatch(model)
  model <- posteriorSimulation(model, mcmcp)
  ##
  ## The standard deviation of the second component is nearly twice as
  ## high as the standard deviation of the third component
  ##
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE, xlim=c(-2,1))
    CNPBayes::plot(model, xlim=c(-2,1))
    par(op)
    mc <- mcmcChains(model)
    plot.ts(sigma(mc), col="gray")
    plot.ts(theta(mc), col="gray")
    cbind(colMeans(sigma(mc)), as.numeric(sigma(truth)))

    probs <- probz(model)
    maxprob <- apply(probs, 1, max)
    plot(maxprob~jitter(as.integer(z(truth)), amount=0.1), pch=".")
  }
  mc <- mcmcChains(model)
  checkEquals(colMeans(theta(mc)), as.numeric(theta(truth)), tolerance=0.05)
  checkEquals(colMeans(sigma(mc)), as.numeric(sigma(truth)), tolerance=0.05)
  checkEquals(colMeans(p(mc)), as.numeric(p(truth)), tolerance=0.03)
}




test_unequal_pmix <- function(){
  library(GenomicRanges)
  library(oligoClasses)
  truth <- readRDS("~/Software/CNPBayes/inst/extdata/unequal_mix_model.rds")

  mcmcp <- McmcParams(iter=1000, burnin=1000)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  bmodel <- initializeModel(params)
  bmodel <- posteriorSimulation(bmodel, mcmcp)

  th1 <- as.numeric(theta(truth))
  th2 <- as.numeric(colMeans(thetac(bmodel)))
  checkEquals(th1, th2, tolerance=0.02)

  th1 <- as.numeric(sigma(truth))
  th2 <- as.numeric(colMeans(sigmac(bmodel)))
  checkEquals(th1, th2, tolerance=0.11)

  th1 <- as.numeric(p(truth))
  th2 <- as.numeric(colMeans(pic(bmodel)))
  checkEquals(th1, th2, tolerance=0.01)

  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(bmodel)
    par(op)
  }
}
