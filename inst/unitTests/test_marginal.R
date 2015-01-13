## no reason to test if the other tests are working
.test_marginalEasy <- function(){
  set.seed(1)
  truth <- simulateData(N=500, .k=3, means=c(-1, 0, 1), sds=c(0.1, 0.1, 0.1))
  model <- truth
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=1000, burnin=100)
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
  checkEquals(pmns, theta(truth), tolerance=0.02)

  ps <- colMeans(sigma(mc))
  checkEquals(ps, sigma(truth), tolerance=0.1)

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkEquals(pmix, pm_pmix, tolerance=0.01)

  pz <- probz(model)
  checkTrue(all(rowSums(pz) == 1))

  ##
  ## HWE
  ##
  zz <- map(model)
  hw.cs <- HardyWeinberg(model)$chisq
  checkIdentical(round(hw.cs, 3), 0.094)
}

test_marginalEasy_default_starts <- function(){
  set.seed(1000)
  truth <- simulateData(N=2500,
                        .k=3,
                        means=c(-1, 0, 1), sds=c(0.1, 0.1, 0.1))
  ##
  ## Start with default starting values
  ##
  mcmcp <- McmcParams(iter=1000, burnin=100)
  params <- ModelParams("marginal", y=y(truth), k=3)
  model2 <- initializeModel(params)
  model2 <- posteriorSimulation(model2, mcmcp)

  mc <- mcmcChains(model2)
  ##plot.ts(theta(mc), col="gray")

  pmns <- colMeans(theta(mc))
  checkEquals(pmns, theta(truth), tolerance=0.02)

  ps <- colMeans(sigma(mc))
  checkEquals(ps, sigma(truth), tolerance=0.06)

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkEquals(pmix, pm_pmix, tolerance=0.025)
}

##.test_wrong_k <- function()## incorrect k
##  mcmcp <- McmcParams(iter=1000, burnin=100, constrainTheta=TRUE)
##  model3 <- initializeModel("marginal", yy=y(truth), k=5)
##  mcmcChains(model3) <- McmcChains(model3, mcmcp)
##  model3 <- posteriorSimulation(model3, mcmcp)
##  mc <- mcmcChains(model3)
##  plot.ts(theta(mc), col="gray")
##  bic5 <- bic(model3)
##  checkTrue(bic3 < bic5)
##}

test_selectK_easy <- function(){
  ## Need to replace setParallelization with source code
  ## -investigate BiocParallel
  library(foreach)
  set.seed(1000)
  means <- c(-1, 0, 1)
  sds <- c(0.1, 0.1, 0.1)
  truth <- simulateData(N=2500, .k=3, means=means, sds=sds)
  if(FALSE) plot(truth, use.current=TRUE)
  ##
  ## Evaluate at different K
  ##
  mcmcp <- McmcParams(iter=1000, burnin=1000)
  mmodels <- fitMixtureModels(y(truth), mcmcp, K=1:5)
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
                        means=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        .alpha=c(100, 200, 400))
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=1000, burnin=1000)

  params <- ModelParams("marginal", y=y(truth), k=3)
  model <- initializeModel(params)
  model <- posteriorSimulation(model, mcmcp)
  model <- sort(model)
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
  ##  setGeneric("sort", "MarginalModel", function(x, decreasing=FALSE, ...){
  mc <- mcmcChains(model)
  pmns <- colMeans(theta(mc))
  checkEquals(pmns, theta(truth), tolerance=0.04)

  ps <- colMeans(sigma(mc))
  checkEquals(ps, sigma(truth), tolerance=0.04)

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkEquals(pmix, pm_pmix, tolerance=0.1)
}

test_bad_starts <- function(){
  set.seed(123)
  ##
  ## Default hyperparameters are not near the true values
  ##
  truth <- simulateData(N=2500,
                        means=c(-1, 0.4, 0.75),
                        sds=c(0.3, 0.1, 0.1),
                        .alpha=c(100, 200, 100))
  mcmcp <- McmcParams(iter=1000, burnin=3000)
  params <- ModelParams("marginal", y=y(truth), k=3, mcmc.params=mcmcp)
  model <- initializeModel(params)
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
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
  checkEquals(pmns, theta(truth), tolerance=0.09)

  ps <- colMeans(sigma(mc))
  checkEquals(ps, sigma(truth), tolerance=0.08)

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkEquals(pmix, pm_pmix, tolerance=0.06)
}




test_marginal_hard <- function(){
  ##
  ## Rare components
  ##
  set.seed(2000)
  truth <- simulateData(N=2500,
                        means=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        .alpha=c(5, 500, 1000))
  model <- truth
  true.sigmas <- sigma(truth)
  if(FALSE) plot(truth, use.current=TRUE)

  mcmcp <- McmcParams(iter=1000, burnin=500)
  params <- ModelParams("marginal", y=y(truth), k=3, mcmc.params=mcmcp)
  modelk <- initializeModel(params)
  model <- posteriorSimulation(modelk, mcmcp)
  model <- sort(model)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  ##
  ## Have to increase the tolerance a bit
  ##
  mc <- mcmcChains(model)
  s <- sigma(mc)
  ps <- colMeans(s)
  checkEquals(ps, sigma(truth), tolerance=0.05)

  pmns <- colMeans(theta(mc))
  checkEquals(pmns, theta(truth), tolerance=0.05)

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkEquals(pmix, pm_pmix, tolerance=0.07)
}



test_selectK_moderate <- function(){
  library(HopkinsHPC)
  setParallelization(7)
  set.seed(1)
  truth <- simulateData(2500,
                        means=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15))
  ##
  ## Evaluate at different K
  ##
  mcmcp <- McmcParams(iter=1000, burnin=500)
  models <- fitMixtureModels(y(truth), mcmcp, K=1:5)
  bicstat <- sapply(models, bic)
  checkTrue(which.min(bicstat) == 3)
}


test_selectK_hard <- function(){
  if(require(doSNOW)){
    library(doSNOW)
    cl <- makeCluster(7, type = "SOCK")
    registerDoSNOW(cl)
  }
  set.seed(1234)
  truth <- simulateData(2500,
                        means=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        .alpha=c(5, 500, 1000))
  ##
  ## Evaluate at different K
  ##
  mcmcp <- McmcParams(iter=2000, burnin=500)
  bicstat <- foreach(k = 1:7, .packages="CNPBayes", .combine="c") %dopar% {
    params <- ModelParams("marginal", k=k, y=y(truth), mcmc.params=mcmcp)
    model <- initializeModel(params)
    model <- posteriorSimulation(model, mcmcp)
    bic(model)
  }
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
