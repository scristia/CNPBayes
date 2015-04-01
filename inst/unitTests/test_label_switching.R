.test_label_switching <- function(){
  library(devtools)
  library(oligoClasses)
  load_all()
  set.seed(1000)
  truth <- simulateData(N=2500,
                        theta=c(-0.4, 0, 0.3),
                        sds=c(0.1, 0.05, 0.1),
                        p=c(0.01, 0.98, 0.01))
  ll.truth <- computeLoglik(truth)
  lp.truth <- logpotential(truth)
  checkEquals(logpotential(truth), computePrior(truth)+logLik(truth)+.loglikPhiMarginal(truth))
  if(FALSE){
    mcmcp <- McmcParams(iter=1000, burnin=1000, nStarts=20, nStartIter=100)
    params <- ModelParams("marginal", y=y(truth), k=3)
    hypp <- HyperparametersMarginal(k=3, eta.0=1, m2.0=.1, tau2.0=10)
    model <- initializeModel(params, hypp)
    plot(truth, use.current=TRUE)
    ## Prior information that if there are multiple components tau2
    ## should be somewhat big
    model <- posteriorSimulation(model, mcmcp)
    mcmcp <- McmcParams(iter=2000, burnin=1000, thin=2)
    model2 <- posteriorSimulation(model, mcmcp)
    th <- sort(thetaMean(model2))
    checkEquals(th, theta(truth), tolerance=0.3)
    ## mean of 2nd component should be very precise
    checkEquals(th[2], theta(truth)[2], tolerance=0.01)
    plot.ts(tauc(model2), col="gray")
    abline(h=sd(theta(truth)))

    plot.ts(thetac(model2), col="gray", plot.type="single")
    abline(h=theta(truth))

    i <- argMax(model2)
    ii <- which.max(logLikc(model2))
    abline(v=i)
    ## plot will not look good because of  label switching
    plot(model2)
    ## but the plot drawn at the modes should look good
    model.modes <- useModes(model2)
    plot(model.modes, use.current=TRUE)
    ##
    ## How to get a posterior probability for the components at each
    ## observation?  (z's will reflect label switching)
    ##
  }
  mcmcp <- McmcParams(iter=c(500, 5000, 5000, 5000),
                      burnin=c(50, 1000, 1000, 1000),
                      thin=c(1, 5, 5, 5),
                      nStarts=20,
                      nStartIter=200)
  hypp <- HyperparametersMarginal(tau2.0=1000)
  hplist <- HyperParameterList(K=1:4, hypp)
  mplist <- ModelParamList(hypp, K=1:4, data=y(truth), batch=batch(bmodel), mcmcp=mcmcp)
  modlist <- foreach(hypp=hplist, param=mplist) %do% initializeModel(params=param, hypp=hypp)
  mmodels <- foreach(k=1:4, model=modlist) %do% posteriorSimulation(model, mcmcp[k])
  ## Look at other modes, fitting reduced Gibbs
  ## -- use modes of mmodels
  mcmcp <- McmcParams(iter=1000, burnin=100, thin=1)
  x <- lapply(mmodels, computeMarginalPr, mcmcp=mcmcp)
  xx <- do.call(cbind, lapply(x, rowMeans))
  post.range <- unlist(lapply(x, posteriorRange))
  marginals <- computeMarginal(xx)
  marginals[post.range > 5] <- -Inf
  checkIdentical(which.max(marginals), 3L)
  ## Fit batch model
  ##
  ##- the batch model should not do as well as the marginal model
  ##  since there is no batch effect in this simulation.
  ##
  ##- however, each batch should pick up the three modes
  ##
  ## Make up a batch mcmcp <- McmcParams(iter=1000, burnin=100)
  mcmcp <- McmcParams(iter=c(500, 2000, 2000, 2000),
                      burnin=c(50, 100, 100, 100),
                      thin=c(1, 2, 2, 2),
                      nStarts=20,
                      nStartIter=200)
  hypp <- HyperparametersBatch(tau2.0=1000)
  hplist <- HyperParameterList(hypp, K=1:4)
  mplist <- ModelParamList(hypp, K=1:4, data=y(truth), batch=batch(bmodel), mcmcp=mcmcp)
  bmodlist <- foreach(hypp=hplist, param=mplist) %do% initializeBatchModel(params=param, hypp=hypp)
  bmodlist2 <- foreach(k=1:4, model=bmodlist) %do% posteriorSimulation(model, mcmcp[k])
  mcmcp <- McmcParams(iter=1000, burnin=100, thin=1)
  x <- lapply(bmodlist2, computeMarginalPr, mcmcp=mcmcp)
  xx <- do.call(cbind, lapply(x, rowMeans))
  post.range <- unlist(lapply(x, posteriorRange))
  bmarginals <- computeMarginal(xx)
  bmarginals[post.range > 5] <- -Inf
  checkTrue(max(mmarginals) > max(bmarginals))
  table(z(truth), batch(params))
  if(FALSE){
    mcmcp <- McmcParams(iter=1000, burnin=100, thin=1, nStarts=10, nStartIter=100)
    params <- ModelParams("batch", y=y(truth), k=3,
                          batch=rep(letters[1:3], length.out=2500),
                          mcmc.params=mcmcp)
    ##
    ## Does log likelihood approach the true model?
    ##  - no
    hypp <- HyperparametersBatch(k=3, tau2.0=100)
    model <- initializeBatchModel(params, hypp=hypp)
    bmodels <- posteriorSimulation(model, mcmcp)
    plot.ts(c(logLikc(bmodels), logLik(truth)), col="gray")
    abline(h=logLik(truth))

    ll.prior <- logLik(truth) + computePrior(truth)
    plot.ts(c(logpotentialc(bmodels), ll.prior), col="gray")
    abline(h=ll.prior)

    par(mfrow=c(1,3))
    tracePlot(bmodels, "theta")
    plot.ts(tauc(bmodels), plot.type="single", col=1:3)
    plot.ts(muc(bmodels), plot.type="single", col=1:3)

    mp <- McmcParams(iter=5000, burnin=0, thin=5)
    bmodel <- posteriorSimulation(bmodels, mp)
    par(mfrow=c(1,1))
    plot.ts(c(logLikc(bmodel), logLik(truth)), col="gray")
    abline(h=logLik(truth))

    ll.prior <- logLik(truth) + computePrior(truth)
    plot.ts(c(logpotentialc(bmodels), ll.prior), col="gray")
    abline(h=ll.prior)

    par(mfrow=c(1,3))
    tracePlot(bmodels, "theta", col=1:3)
    plot.ts(tauc(bmodels), plot.type="single", col=1:3)
    plot.ts(muc(bmodels), plot.type="single", col=1:3)

    if(FALSE){
      par(mfrow=c(1,3), las=1)
      tracePlot(bmodels, "theta", col=1:3)
    }
  }

}
