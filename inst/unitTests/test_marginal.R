##test_marginal <- function(){
##  set.seed(1)
##  model <- simulateData(N=2500,
##                        .k=3,
##                        .theta=c(-1, 0, 1),
##                        .sigma=c(0.1, 0.1, 0.1))
##  if(FALSE) plot(model, use.current=TRUE)
##  truth <- model
##  mcmcp <- McmcParams(iter=5, burnin=2)
##  mcmcChains(model) <- McmcChains(model, mcmcp)
##  model <- posteriorSimulation(model, mcmcp)
##  checkTrue(TRUE)
##}

test_marginalEasy <- function(){
  set.seed(1)
  model <- simulateData(N=2500,
                        .k=3,
                        .theta=c(-1, 0, 1),
                        .sigma=c(0.1, 0.1, 0.1))
  if(FALSE) plot(model, use.current=TRUE)
  truth <- model
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
  checkTrue(all.equal(mean(abs(pmns-theta(truth))), 0, 0.01))

  ps <- colMeans(sigma(mc))
  checkTrue(all.equal(mean(abs(ps-sigma(truth))), 0, 0.01))

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkTrue(all.equal(mean(abs(pmix-pm_pmix)), 0, 0.02))
}

test_marginal_Moderate <- function(){
  set.seed(1)
  model <- simulateData(2500,
                        .theta=c(-2, -0.4, 0),
                        .sigma=c(0.3, 0.15, 0.15))
  truth <- model
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=2000, burnin=500)
  mcmcChains(model) <- McmcChains(model, mcmcp)
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
  checkTrue(all.equal(mean(abs(pmns-theta(truth))), 0, 0.01))

  ps <- colMeans(sigma(mc))
  checkTrue(all.equal(mean(abs(ps-sigma(truth))), 0, 0.01))

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkTrue(all.equal(mean(abs(pmix-pm_pmix)), 0, 0.02))
}


test_marginal_hard <- function(){
  ##
  ## Rare components
  ##
  set.seed(1)
  model <- simulateData(2500,
                        .theta=c(-2, -0.4, 0),
                        .sigma=c(0.3, 0.15, 0.15),
                        .alpha=c(5, 500, 1000))
  truth <- model
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=2000, burnin=500)
  mcmcChains(model) <- McmcChains(model, mcmcp)
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
  ##save(pmns, file="pmns.rda")
  checkTrue(all.equal(mean(abs(pmns-theta(truth))), 0, 0.1))

  ##
  ## Variance of component with few observations is higher than the
  ## true variance -- this is reasonable as more weight would go to
  ## the prior
  ps <- colMeans(sigma(mc))
  checkTrue(all.equal(mean(abs(ps-sigma(truth))), 0, 0.02))

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkTrue(all.equal(mean(abs(pmix-pm_pmix)), 0, 0.01))

  logpotential(truth)
  logp <- logpotential(mc)
}

.selectK_easy <- function(){
  library(HopkinsHPC)
  setParallelization(8)
  set.seed(1)
  truth <- simulateData(N=2500,
                        .k=3,
                        .theta=c(-1, 0, 1),
                        .sigma=c(0.1, 0.1, 0.1))
  ##
  ## Evaluate at different K
  ##
  ## if BIC increases for 2 consecutive models, stop
  bic <- rep(NA, 8)
  bic <- foreach(k = 1:8, .packages="CNPBayes", .combine="c") %dopar% {
    model <- initializeModel(y(truth), k=k, batch=batch(truth))
    mcmcp <- McmcParams(iter=2000, burnin=500)
    mcmcChains(model) <- McmcChains(model, mcmcp)
    model <- posteriorSimulation(model, mcmcp)
    BIC(model)
  }
  checkTrue(which.min(bic) == 3)
}

.selectK_moderate <- function(){
  set.seed(1)
  truth <- simulateData(2500,
                        .theta=c(-2, -0.4, 0),
                        .sigma=c(0.3, 0.15, 0.15))
  ##
  ## Evaluate at different K
  ##
  bic <- rep(NA, 8)
  bic <- foreach(k = 1:8, .packages="CNPBayes", .combine="c") %dopar% {
    model <- initializeModel(y(truth), k=k, batch=batch(truth))
    mcmcp <- McmcParams(iter=2000, burnin=500)
    mcmcChains(model) <- McmcChains(model, mcmcp)
    model <- posteriorSimulation(model, mcmcp)
    BIC(model)
  }
  checkTrue(which.min(bic) == 3)
}


.selectK_hard <- function(){
  set.seed(1)
  truth <- simulateData(2500,
                        .theta=c(-2, -0.4, 0),
                        .sigma=c(0.3, 0.15, 0.15),
                        .alpha=c(5, 500, 1000))
  ##
  ## Evaluate at different K
  ##
  bic <- rep(NA, 8)
  bic <- foreach(k = 1:8, .packages="CNPBayes", .combine="c") %dopar% {
    model <- initializeModel(y(truth), k=k, batch=batch(truth))
    mcmcp <- McmcParams(iter=2000, burnin=500)
    mcmcChains(model) <- McmcChains(model, mcmcp)
    model <- posteriorSimulation(model, mcmcp)
    BIC(model)
  }
  save(bic, file="bic.rda")
  checkTrue(which.min(bic) == 3 || which.min(bic) == 2)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)

    ## variances are off
    op <- par(mfrow=c(3,1), las=1)
    for(i in 1:3){
      ## consistent with shrinkage
      plot.ts(sigma2(mc)[,i], col="gray", xaxt="n")
      abline(h=as.numeric(sigma2(truth))[i], col="blue")
      abline(h=mean(sigma2(mc)[, i]), lty=2)
    }
    par(op)

    ## if we started at the true value, would we get the right
    ## variance
    mcmcp <- McmcParams(iter=2000, burnin=500)
    mcmcChains(model) <- McmcChains(model, mcmcp)
    theta(model) <- theta(truth)
    z(model) <- z(truth)
    sigma2(model) <- sigma2(truth)
    model <- posteriorSimulation(model, mcmcp)
    op <- par(mfrow=c(3,1), las=1)
    for(i in 1:3){
      ## consistent with shrinkage
      plot.ts(sigma2(mc)[,i], col="gray", xaxt="n")
      abline(h=as.numeric(sigma2(truth))[i], col="blue")
      abline(h=mean(sigma2(mc)[, i]), lty=2)
    }
    par(op)
    BIC(model)
  }
}
