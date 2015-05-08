checkEquals2 <- function(target, ..., order){
  checkEquals(target[order], ...)
}

  ## no reason to test if the other tests are working
test_marginalEasy <- function(){
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth, use.current=TRUE)
  ##mp <- McmcParams(iter=500, burnin=500)
  mp <- McmcParams(iter=5, burnin=5, nStarts=1)
  model <- MarginalModel(data=y(truth), k=3, mcmc.params=mp)
  model <- startAtTrueValues(model, truth)
  model <- posteriorSimulation(model)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  mc <- mcmcChains(model)
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
  i <- argMax(model)
  checkTrue(i == which.max(logPrior(chains(model)) + logLik(chains(model))))
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

  truth <- simulateData(N=2500, p=c(1/4, 1/2, 1-1/2-1/4), theta=means, sds=sds)
  x2 <- computeMarginalLik(y(truth), nchains=3, K=1:4)
  m2 <- orderModels(x2)
  zz <- map(m2[[1]])
  tab <- table(zz)
  checkTrue(length(tab)==3)
}

.test_selectK_nobatchEffect <- function(){
  ## Fit batch model when truth is no batch effect
  ## batch found
  library(GenomicRanges)
  set.seed(1000)
  means <- c(-1, 0, 1)
  sds <- c(0.1, 0.2, 0.2)
  truth <- simulateData(N=2500, p=rep(1/3, 3), theta=means, sds=sds)
  se <- as(truth, "SummarizedExperiment")
  library(devtools)
  load_all()
  ##m <- marginal(se, mcmc.params=McmcParams(iter=10, burnin=5, nStarts=1), K=3)

  kmod <- MarginalModel(y(truth), k=3)
  iter(kmod, force=TRUE) <- 500
  kmod <- posteriorSimulation(kmod)
  headtheta <- head(thetac(kmod))

  ##
  ## Marginal likelihood
  ##
  logp_theta <- .Call("marginal_theta", kmod)
  logp <- log(mean(exp(logp_theta)))
  checkTrue(identical(headtheta, head(thetac(kmod))))

  ##  logp_theta <- .Call("p_theta_zpermuted", kmod)
  ##  logp <- log(mean(exp(logp_theta)))
  ##  checkTrue(identical(headtheta, head(thetac(kmod))))

  ## Simulate Z under the reduced model
  kmodZ <- .Call("simulate_z_reduced1", kmod)
  checkTrue(identical(headtheta, head(thetac(kmod))))
  zz <- zChain(kmodZ)

  logp_sigma2 <- .Call("p_sigma2_zpermuted", kmodZ)
  lps2 <- log(mean(exp(logp_sigma2)))
  checkTrue(identical(headtheta, head(thetac(kmod))))

  logp_p <- .Call("p_p_unpermuted", kmod)
  lpp <- log(mean(exp(logp_p)))
  checkTrue(identical(headtheta, head(thetac(kmod))))

  i <- argMax(kmod)
  mlik <- logLik(chains(kmod))[i] + logPrior(chains(kmod))[i] - logp - lps2 - lpp

  permutations <- permnK(3, 5)
  kmod2 <- kmod
  zChain(kmod2) <- permuteZ(kmod, permutations[2, ])
  zz <- table(zChain(kmod)[1, ])
  zz2 <- table(zChain(kmod2)[1, ])
  ##
  ## this uses the permuted values of z directly.
  ##
  logp_theta2 <- .Call("p_theta_zpermuted", kmod2)  ## use permuted z
  logp2 <- log(mean(exp(logp_theta2)))
  ## now we have to simulate new values of z as described in Chib (We
  ## need z|y, theta* and not z|y).  To do this, we could use the z
  ## chain already simulated under the reduced model and permute these
  logp_sigma <- .Call("p_sigma2_zpermuted", kmod)
  lps1 <- log(mean(exp(logp_sigma)))

  ans <- gtools::ddirichlet(x=c(0.2, 0.3, 0.5), alpha=c(100, 150, 200))
  result <- .Call("ddirichlet_", c(0.2, 0.3, 0.5), c(100, 150, 200))
  checkEquals(log(ans), result)
##  batchExperiment(se, outdir, mcmcp.list=mcmp.list)
##  B <- getFiles(outdir, rownames(se), "batch")
##  checkTrue(is.null(readRDS(model(B))))

  ## simulate a bogus batch
##  se$plate <- factor(rep(letters[1:3], length.out=ncol(se)))
##  batchExperiment(se, outdir, mcmcp.list=mcmp.list)
##  B <- getFiles(outdir, rownames(se), "batch")
  ##  checkTrue(is.null(readRDS(model(B))))


  setClass("foo", representation(x="numeric"))
  object <- new("foo")
  object@x <- rnorm(10)

  object2 <- .Call("test_clone", object)

}

test_marginal_Moderate <- function(){
  set.seed(100)
  truth <- simulateData(N=2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.05, 0.1, 0.8))
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=500, burnin=500, thin=2)
  model <- MarginalModel(y(truth), k=3, mcmc.params=mcmcp)
  model <- startAtTrueValues(model, truth)
  model <- posteriorSimulation(model)
  checkEquals(sort(theta(model)), theta(truth), tolerance=0.15)
  if(FALSE){
    plot.ts(thetac(model), plot.type="single")
    plot.ts(sigmac(model), plot.type="single")
    plot.ts(pic(model), plot.type="single")
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
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
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=1000, burnin=100, thin=10)
  modelk <- MarginalModel(y(truth), k=3, mcmc.params=mcmcp)
  model <- posteriorSimulation(modelk)
  i <- order(theta(model))
  checkEquals(theta(model)[i], theta(truth), tolerance=0.1)
  checkEquals(colMeans(sigmac(model))[i], sigma(truth), tolerance=0.1)
  checkEquals(colMeans(pic(model))[i], p(truth), tolerance=0.15)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
    plot.ts(sigmac(model), col=1:3, plot.type="single")
  }
}



.test_cnp360 <- function(){
  library(GenomicRanges)
  se360 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp360.rds")
  r <- assays(se360)[["mean"]][1, ]
  if(FALSE) hist(r, breaks=1000, col="gray", border="gray")
  x <- computeMarginalLik(r, nchains=3)
  ## hemizygous deletion is rare
  models <- orderModels(x)
  checkTrue(k(models)[1] %in% c(2, 3))
  cn <- map(models[[1]])
  freq_hem <- mean(cn==1)
  ## get homozygous post-hoc
  cn[ r < -4 ] <- 0

  ##
  ## Computationally not feasible to explore all K models adequately
  ## with 10,000 samples
  ##
  index <- sample(seq_along(r), 2000)
  table(cn[index])
  xx <- computeMarginalLik(r[index], nchains=3)
  models2 <- orderModels(xx)
  checkTrue(k(models2[[1]]) == 2)
  df <- imputeFromSampledData(models2[[1]], r, index)
  checkTrue(!any(is.na(df$cn)))
  cn_complete <- map(models[[1]])
  cn_incomplete <- df$cn
  checkTrue(mean(cn_incomplete != cn_complete) < 0.001)

  p_complete <- cnProbability(probz(models[[1]]), 2)
  p_incomplete <- df$p
  checkTrue(mean(abs(p_complete - p_incomplete) > 0.1 ) < 0.01)
}
