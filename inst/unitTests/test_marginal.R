## no reason to test if the other tests are working
test_marginalEasy <- function(){
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3), theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth, use.current=TRUE)

  params <- ModelParams("marginal", y=y(truth), k=3)
  mcmcp <- McmcParams(iter=1000, burnin=500)
  ##trace(initializeModel, browser)
  model <- initializeModel(params)
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  mc <- mcmcChains(model)
  pmns <- colMeans(theta(mc))
  checkEquals(sort(pmns), theta(truth), tolerance=0.03)

  ps <- colMeans(sigma(mc))
  checkEquals(ps[order(pmns)], sigma(truth), tolerance=0.1)

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkEquals(pmix, pm_pmix[order(pmns)], tolerance=0.05)
}




.test_k_too_big <- function(){
  ## when k is too big, chains may have a greater likelihood of crossing
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3), theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth, use.current=TRUE)
  params <- ModelParams("marginal", y=y(truth), k=3)
  mcmcp <- McmcParams(iter=1000, burnin=500)
  params <- ModelParams("marginal", y=y(truth), k=5)
  mcmcp <- McmcParams(iter=1000, burnin=500)
  model <- initializeModel(params)
  model <- posteriorSimulation(model, mcmcp)
}

## no reason to test if the other tests are working
test_normality_assumption <- function(){
  set.seed(1)
  ## does the prior on theta hurt us
  truth <- simulateData(N=2500, p=c(0.01, 1/3, 1-1/3-0.01), theta=c(-5, -0.4, 0),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth, use.current=TRUE)
  params <- ModelParams("marginal", y=y(truth), k=3)
  mcmcp <- McmcParams(iter=1000, burnin=500)
  model <- initializeModel(params)
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  mc <- mcmcChains(model)
  pmns <- colMeans(theta(mc))
  checkEquals(sort(pmns), theta(truth), tolerance=0.03)

  ps <- colMeans(sigma(mc))
  checkEquals(ps[order(pmns)], sigma(truth), tolerance=0.1)

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))[order(pmns)]
  checkEquals(pmix, pm_pmix, tolerance=0.05)
}

test_selectK_easy <- function(){
  ## Need to replace setParallelization with source code
  ## -investigate BiocParallel
  set.seed(1000)
  means <- c(-1, 0, 1)
  sds <- c(0.1, 0.2, 0.2)
  truth <- simulateData(N=2500, p=rep(1/3, 3), theta=means, sds=sds)
  if(FALSE) plot(truth, use.current=TRUE)
  ##
  ## Evaluate at different K
  ##
  mcmcp <- McmcParams(iter=rep(1000, 4), burnin=c(100, 200, 1000, 1000))
  mmodels <- fitMixtureModels(y(truth), mcmcp, K=1:4)
  bicstat <- sapply(mmodels, bic)
  checkIdentical(which.min(bicstat), 3L)
  if(FALSE){
    mc <- mcmcChains(mmodels[[3]])
    plot.ts(sigma(mc), col="gray")
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(mmodels[[3]])
    par(op)
  }
}

test_marginal_Moderate <- function(){
  set.seed(100)
  truth <- simulateData(N=2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.05, 0.1, 0.8))
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=2000, thin=2, burnin=1000)
  params <- ModelParams("marginal", y=y(truth), k=3)
  model <- initializeModel(params)
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    plot.ts(thetac(model))
    plot.ts(sigma(model))
    plot.ts(pic(model))
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  ##
  ## Have to increase the tolerance a bit
  ##
  mc <- mcmcChains(model)
  pmns <- colMeans(theta(mc))
  checkEquals(sort(pmns), theta(truth), tolerance=0.06)

  ps <- colMeans(sigma(mc))
  checkEquals(ps[order(pmns)], sigma(truth), tolerance=0.09)

  pmix <- p(truth)

  pm_pmix <- pic(model)[argMax(model), ]
  checkEquals(pmix, pm_pmix[order(pmns)], tolerance=0.1)
}

test_marginal_hard <- function(){
  ##
  ## The mixture probabilities are slow mixing -- high
  ## autocorrelation.  Much more thinning is probably needed to
  ## adequately explore the space.
  ##
  ##
  ##
  ## Rare components
  ##
  set.seed(2000)
  truth <- simulateData(N=5e3,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.005, 1/10, 1-0.005-1/10))
  if(FALSE) plot(truth, use.current=TRUE)

  mcmcp <- McmcParams(iter=2000, thin=2, burnin=1000)
  params <- ModelParams("marginal", y=y(truth), k=3,
                        mcmc.params=mcmcp)
  modelk <- initializeModel(params)
  model <- posteriorSimulation(modelk, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  ##
  ## Have to increase the tolerance a bit
  ##
  pmns <- thetaMean(model)
  checkEquals(sort(pmns), theta(truth), tolerance=0.08)

  ps <- sigmaMean(model)
  checkEquals(ps[order(pmns)], sigma(truth), tolerance=0.15)

  pmix <- pMean(model)
  checkEquals(p(truth), pmix[order(pmns)], tolerance=0.09)
  ##min(effectiveSize::p(mc)) ## is very low
}

.test_outliers <- function(){
  set.seed(2000)
  truth <- simulateData(N=2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.003, 1/10, 1-0.003-1/10))
  ## heavy-tailed
  ##y(truth)[z==1] <- rt
}



test_selectK_moderate <- function(){
  ##
  ## Could think about having a MarginalPlusHomozygous class and
  ## constraining the variance of one component to be at least the
  ## variance of the other two components.  But then this would make
  ## the components no longer exchangeable...
  ##
  set.seed(1)
  truth <- simulateData(2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(1/10, 1/5, 1-1/10-1/5))
  ##
  ## Evaluate at different K.  models with fewer components
  ## (parameters) require fewer iterations to converge.
  ##
  mcmcp <- McmcParams(iter=c(1000, 1000, 2000, 2000, 2000),
                      burnin=c(100, 200, 1000, 1000, 1000),
                      thin=c(1, 1, 2, 2, 2),
                      nStart=20,
                      nStartIter=150)
  hypp <- HyperparametersMarginal()
  hplist <- HyperParameterList(hypp, K=1:4)
  mplist <- ModelParamList(hypp, K=1:4, data=y(truth),
                           mcmcp=mcmcp)
  modlist <- foreach(hypp=hplist, param=mplist) %do% {
    initializeModel(params=param, hypp=hypp)
  }
  models <- foreach(k=1:4, model=modlist) %do% posteriorSimulation(model, mcmcp[k])
  bicstat <- sapply(models, bic)
  checkTrue(which.min(bicstat) == 3)
  if(FALSE){
    op <- par(mfrow=c(1,3),las=1)
    plot(truth, use.current=T)
    plot(models[[2]])
    plot(models[[3]])
    par(op)
  }
}


test_selectK_hard <- function(){
  set.seed(1234)
  mcmcp <- McmcParams(iter=c(1000, 1000, 2000, 2000, 2000),
                      burnin=c(100, 200, 1000, 1000, 1000),
                      thin=c(1, 1, 2, 2, 2),
                      nStart=20,
                      nStartIter=150)
  truth <- simulateData(2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(1/100, 1/10, 1-0.1-0.001))
  hypp <- HyperparametersMarginal()
  hplist <- HyperParameterList(hypp, K=1:4)
  mplist <- ModelParamList(hypp, K=1:4, data=y(truth),
                           mcmcp=mcmcp)
  modlist <- foreach(hypp=hplist, param=mplist) %do% {
    initializeModel(params=param, hypp=hypp)
  }
  fit <- foreach(k=1:4, model=modlist) %do% posteriorSimulation(model, mcmcp[k])
  ##
  ## Evaluate at different K
  ##
  ##mcmcp <- McmcParams(iter=50000, thin=50, burnin=3000)
  ##params <- ModelParams("marginal", y=y(truth), k=3, mcmc.params=mcmcp)
  ##modelk <- initializeModel(params)
  ##model <- posteriorSimulation(modelk, mcmcp)
  ##fit <- fitMixtureModels(y(truth), mcmcp, K=1:4)
  bicstat <- sapply(fit, bic)
  checkTrue(which.min(bicstat) == 3)
}




.test_cnp360 <- function(){
  library(GenomicRanges)
  se360 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp360.rds")
  r <- assays(se360)[["mean"]][1, ]
  if(FALSE) hist(r, breaks=1000, col="gray", border="gray")

  mcmcp <- McmcParams(iter=2000, burnin=500)
  params <- ModelParams("marginal", y=r, k=3, batch=rep("A", length(r)),
                        mcmc.params=mcmcp)
  model <- initializeModel(params)
  model <- posteriorSimulation(model, mcmcParams(params))
  mc <- mcmcChains(model)
  plot.ts(theta(mc), col="gray")
  plot(model)
  bic(model)
}
