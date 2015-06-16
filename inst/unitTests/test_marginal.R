checkEquals2 <- function(target, ..., order){
  checkEquals(target[order], ...)
}

  ## no reason to test if the other tests are working
test_marginalEasy <- function(){
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth)
  ##mp <- McmcParams(iter=500, burnin=500)
  mp <- McmcParams(iter=5, burnin=5, nStarts=1)
  model <- MarginalModel(data=y(truth), k=3, mcmc.params=mp)
  model <- CNPBayes:::startAtTrueValues(model, truth)
  model <- posteriorSimulation(model)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth)
    plot(model)
    par(op)
  }
  mc <- chains(model)
  pmns <- colMeans(theta(mc))
  i <- order(pmns)
  checkEquals2(pmns, theta(truth), order=i, tolerance=0.03)

  ps <- colMeans(sigma(mc))
  checkEquals2(ps, sigma(truth), tolerance=0.11, order=i)

  pmix <- p(truth)
  pp <- colMeans(p(mc))
  checkEquals2(pp, pmix, tolerance=0.05, order=i)
  ##
  ## Update slots with the mode of loglik + logprior
  ##
  i <- CNPBayes:::argMax(model)
  checkTrue(i == which.max(logPrior(chains(model)) + log_lik(chains(model))))
  ## Check the modes
  checkIdentical(modes(model)[["theta"]], thetac(model)[i, ])
  checkIdentical(sqrt(modes(model)[["sigma2"]]), sigmac(model)[i, ])
}


test_selectK_easy <- function(){
  library(GenomicRanges)
  set.seed(1000)
  means <- c(-1, 0, 1)
  sds <- c(0.1, 0.2, 0.2)
  truth <- simulateData(N=2500, p=rep(1/3, 3), theta=means, sds=sds)
  ## When the proportion of observations in each state is the same, we
  ## tend to over fit
  x1 <- computeMarginalLik(y(truth), nchains=3, K=1:4)
  m1 <- orderModels(x1)
  checkTrue(k(m1)[1] >= 3)

  set.seed(1000)
  truth <- simulateData(N=2500, p=c(1/4, 1/2, 1-1/2-1/4), theta=means, sds=sds)
  x2 <- computeMarginalLik(y(truth), nchains=3, K=1:4, T=5000, T2=1000, burnin=1000)
  m2 <- orderModels(x2)
  checkTrue(k(m2)[1] >= 3)
   ## T must be >= T2
  checkException(computeMarginalLik(y(truth),
                                    nchains=3,
                                    K=1:4,
                                    T=10, T2=20, burnin=1000))
}

test_marginal_Moderate <- function(){
  set.seed(100)
  truth <- simulateData(N=2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.05, 0.1, 0.8))
  if(FALSE) plot(truth)
  mcmcp <- McmcParams(iter=500, burnin=500, thin=2)
  model <- MarginalModel(y(truth), k=3, mcmc.params=mcmcp)
  model <- CNPBayes:::startAtTrueValues(model, truth)
  model <- posteriorSimulation(model)
  checkEquals(sort(theta(model)), theta(truth), tolerance=0.15)
  if(FALSE){
    plot.ts(thetac(model), plot.type="single")
    plot.ts(sigmac(model), plot.type="single")
    plot.ts(pic(model), plot.type="single")
    op <- par(mfrow=c(1,2),las=1)
    plot(truth)
    plot(model)
    par(op)
  }
  checkEquals(sigma(model)[order(theta(model))], sigma(truth), tolerance=0.15)
  checkEquals(colMeans(pic(model))[order(theta(model))], p(truth), tolerance=0.2)
}

test_marginal_hard <- function(){
  ##
  ## The mixture probabilities are slow mixing -- high
  ## autocorrelation.  Much more thinning is probably needed to
  ## adequately explore the space.
  ##
  ## Rare components
  ##
  set.seed(2000)
  truth <- simulateData(N=5e3,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.005, 1/10, 1-0.005-1/10))
  if(FALSE) plot(truth)
  mcmcp <- McmcParams(iter=1000, burnin=100, thin=10)
  modelk <- MarginalModel(y(truth), k=3, mcmc.params=mcmcp)
  model <- posteriorSimulation(modelk)
  i <- order(theta(model))
  checkEquals(theta(model)[i], theta(truth), tolerance=0.1)
  checkEquals(colMeans(sigmac(model))[i], sigma(truth), tolerance=0.1)
  checkEquals(colMeans(pic(model))[i], p(truth), tolerance=0.15)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth)
    plot(model)
    par(op)
    plot.ts(sigmac(model), col=1:3, plot.type="single")
  }
}
