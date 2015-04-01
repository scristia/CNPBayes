##
## TODO:
##  - save chain of zfreq
##  - mcmc updates use zfreq instead of tableZ
##  - finish posteriorP , then marginal probabilities
.test_rcpp_model <- function(){
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
                             p=rep(1/3,3))
  ##p=c(1/5, 1/3, 1-1/3-1/5))
  trace(.updateMuBatch, browser)
  .updateMuBatch(truth)
  .Call("update_mu_batch", truth)

  ## Why not order thetabar when updating mu?

  mcmcp <- McmcParams(iter=2, burnin=1)
  ##
  ## selecting eta.0 and m2.0 such that the standard deviation of tau
  ## is small.  tau follows inverse gamma.  1/tau is gamma
  ##  if desired mean of tau is 0.1, desired precision is 10.

  ##(y <- 1/rgamma(1000, 0.5*eta.0, rate=0.5*eta.0*m2.0))
  ##hist(y, breaks=100, xlim=c(0,5))
  model <- BatchModel(y(truth), batch=batch(truth),
                      k=3, mcmc.params=mcmcp)
  ##hypp=HyperparametersBatch(tau2.0=0.1,k=3))

  set.seed(1)
  .updateMuBatch(model)
  set.seed(1)
  .Call("update_mu_batch", model)

  model2 <- startAtTrueValues(model, truth)
  ##trace(.updateMuBatch, browser)
  ##.updateMuBatch(model2)
  ##.Call("update_mu_batch", model)
  ##hypp=HyperparametersBatch(tau2.0=0.001, k=3, eta.0=0.01,
  ##m2.0=0.01))
  ## the mus should be simulated
  iter(model) <- 500
  burnin(model) <- 0
  model <- runBurnin(model)
  model <- runMcmc(model)
  plot.ts(muc(model), col=1:3, plot.type="single")

  .Call("update_mu_batch", model)
  ##paramUpdates(mcmcp)[c("sigma2", "p", "nu.0", "mu", "tau2", "nu.0", "sigma2.0", "z")] <- 0L

  mcmcp <- McmcParams(iter=2, burnin=1)
  ##paramUpdates(mcmcp)["nu.0"] <- 0L
  mcmcParams(model) <- mcmcp
  burnin(model) <- 50
  set.seed(1)
  modelr <- .runBurnin(model)
  set.seed(1)
  modelc <- .Call("mcmc_batch_burnin", model, mcmcParams(model))

  mcmcParams(model) <- mcmcp
  iter(model) <- 1
  modelc <- .Call("mcmc_batch", model, mcmcParams(model))


  .Call("compute_vars_batch", model)
  .Call("compute_prec_batch", model)
  .Call("compute_prec_batch", model)
  1/.computeVarsBatch(model)
  .Call("compute_prec_batch", model)
  .Call("compute_vars_batch", model)

  set.seed(1)
  (zr=.updateMuBatch(model))
  set.seed(1)
  .Call("update_mu_batch", model)

  set.seed(1)
  zr=updateZ(model)
  set.seed(1)
  table(.Call("update_z_batch", model))

  p(modelc)=.Call("update_multinomialPr_batch", model)
  zc=.Call("update_z_batch", modelc)


  set.seed(1)
  .updateTau2Batch(model)
  set.seed(1)
  .Call("update_tau2_batch", model)

  set.seed(1)
  .updateThetaBatch(model)
  set.seed(1)
  .Call("update_theta_batch", model)

  set.seed(1)
  .update_sigma2(model)
  set.seed(1)
  .Call("update_sigma2_batch", model)

  set.seed(1)
  .updateSigma2.0Batch(model)
  set.seed(1)
  .Call("update_sigma20_batch", model)


  ll_r <- system.time(logLikData2(model))
  ll_r <- computeLoglik(model)
  ll_c <- system.time(.Call("compute_loglik_batch", model))

  mu_r <- updateMu(model)
  trace(.updateMuBatch, browser)
  system.time(.updateMuBatch(model))
  system.time(.Call("update_mu_batch", model))
  updateMu(model)

  set.seed(1)
  tmp <- .updateTau2Batch(model)
  set.seed(1)
  .Call("update_tau2_batch", model)

  set.seed(1)
  tmp <- updateSigma2.0(model)
  set.seed(1)
  .Call("update_sigma20_batch", model)

  ##
  ## C does not perfectly duplicate R version.  Isolated the
  ## difference to the sampling scheme used in C
  ##
  set.seed(1)
  updateNu.0(model)
  set.seed(1)
  .Call("update_nu0_batch", model)

  pr <- .Call("update_multinomialPr_batch", model)
  zz <- updateZ(model)

  mn_r <- computeMeans(model)
  mn_c <- .Call("compute_means_batch", model)
  v_r <- computeVars(model)
  v_c <- .Call("compute_vars_batch", model)
  p_c <- .Call("compute_prec_batch", model)

  pr_r <- computePrior(model)
  pr_c <- .Call("compute_logprior_batch", model)

  set.seed(1)
  (ss_r <- updateSigma2(model))
  set.seed(1)
  (ss_c <- .Call("update_sigma2_batch", model))

  (th_r <- updateSigma2(model))

  set.seed(1)
  th_r <- updateTheta(model)
  set.seed(1)
  th_c <- .Call("update_theta_batch", model)


  (th_c <- .Call("update_theta_batch", model))



  ll_r <- logLikData(model)
  ll_r2 <- logLikData2(model)
  ll_c <- .Call("compute_loglik_batch", model)

  load_all()
  x <- as.integer(1:3)
  .Call("sumVector", 1:3)

  tmp <- .Call("tableBatchZ", model) ##k(model), z(model),
  tabz <- table(batch(model), z(model))
  checkEquals(as.integer(tmp), as.integer(tabz))

  r <- y(model)[batch(model)=="c"]


  theta(model)

  model <- .runBurnin(model, mcmcParams(model))
  .runMcmc(model, mcmcParams(model))








  model_c <- .Call("mcmc_marginal_burnin", model1, mcmcp)
  model_c2 <- .Call("mcmc_marginal", model_c, mcmcp)

  computeModes(model_c2)

  coda::effectiveSize(thetac(model_c)[, 1]) ## looks good
  checkEquals(as.numeric(thetac(model_c)), as.numeric(theta(model1)))

  th_c <- .Call("mcmc_marginal", model1)
  th_c <- .Call("update_theta", model1)

  model2 <- startAtTrueValues(model, truth)
  model3 <- posteriorSimulation(model, mcmcp)



  ##
  ## - moving the chain
  ## - burnin
  ## - marginal probabilities
  ##

  theta(model)

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
