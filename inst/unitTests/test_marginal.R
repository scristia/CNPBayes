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
  mp <- McmcParams(iter=500, burnin=500)
  model <- MarginalModel(data=y(truth), k=3, mcmc.params=mp)
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

  ##
  ## Explore another mode by switching the labels
  ##
  model2 <- relabel(model, zindex=3:1)
  checkEquals(zFreq(model2), zFreq(model)[3:1])
  mcmcParams(model2) <- mp
  burnin(model2) <- 0
  model2 <- posteriorSimulation(model2)
  if(FALSE){
    par(mfrow=c(2,1), las=1)
    plot.ts(thetac(model), col=1:3, plot.type="single")
    plot.ts(thetac(model2), col=1:3, plot.type="single")
  }
  pmns2 <- colMeans(thetac(model2))
  checkEquals(pmns, pmns2[3:1], tolerance=0.05)

  ##
  ## Create a model for each of the possible modes
  ##
  mmod <- useModes(model)
  burnin(mmod) <- 500
  K <- k(model)
  checkIdentical(thetac(model)[argMax(model), ], theta(mmod))
  model.list <- ModelEachMode(mmod, maxperm=5)
  zlist <- lapply(model.list, zFreq)
  ztab <- zFreq(mmod)
  ## the first model is in the same order
  permutations <- combinat::permn(1:3)
  checkTrue(all(foreach(zfreq=zlist, i=permutations, .combine="c") %do% identical(zfreq, ztab[i])))
  model <- startAtTrueValues(model, truth)
  m.y(model) <- computeMarginalProbs(model, iter=500, burnin=200)
  checkTrue(diff(range(m.y(model))) < 10 )

  ##
  ## compute marginal density for other K
  ##
  marginaly <- computeMarginalEachK(y(truth), K=1:4,
                                    mcmcp=McmcParams(iter=200, burnin=100),
                                    MAX.RANGE=10)
  checkTrue(which.max(marginaly)==3)
  bayesFactor(marginaly)
}

test_selectK_easy <- function(){
  library(GenomicRanges)
  set.seed(1000)
  means <- c(-1, 0, 1)
  sds <- c(0.1, 0.2, 0.2)
  truth <- simulateData(N=2500, p=rep(1/3, 3), theta=means, sds=sds)
  ## Also fit batch model when truth is no batch effect
  ## batch found
  outdir <- tempdir()
  se <- as(truth, "SummarizedExperiment")
  batchExperiment(se, outdir, mcmcp.list=mcmp.list)
  B <- getFiles(outdir, rownames(se), "batch")
  checkTrue(is.null(readRDS(model(B))))

  ## simulate a bogus batch
  se$plate <- factor(rep(letters[1:3], length.out=ncol(se)))
  batchExperiment(se, outdir, mcmcp.list=mcmp.list)
  B <- getFiles(outdir, rownames(se), "batch")
  checkTrue(is.null(readRDS(model(B))))
}

test_marginal_Moderate <- function(){
  set.seed(100)
  truth <- simulateData(N=2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.05, 0.1, 0.8))
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=500, burnin=100, thin=10)
  model <- MarginalModel(y(truth), k=3, mcmc.params=mcmcp)
  model <- startAtTrueValues(model, truth)
  model <- posteriorSimulation(model, mcmcp)
  checkEquals(theta(model), theta(truth), tolerance=0.1)
  if(FALSE){
    plot.ts(thetac(model), plot.type="single")
    plot.ts(sigmac(model), plot.type="single")
    plot.ts(pic(model), plot.type="single")
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  checkEquals(sigma(model), sigma(truth), tolerance=0.15)
  checkEquals(colMeans(pic(model)), p(truth), tolerance=0.2)
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

  mcmcp <- McmcParams(iter=150, burnin=25)
  modelk <- MarginalModel(y(truth), k=3, mcmc.params=mcmcp)
  modelk <- startAtTrueValues(modelk, truth)
  model <- posteriorSimulation(modelk, mcmcp)
  checkEquals(sort(theta(model)), theta(truth), tolerance=0.1)
  checkEquals(sigma(model)[order(theta(model))], sigma(truth), tolerance=0.1)
  checkEquals(p(model)[order(theta(model))], p(truth), tolerance=0.1)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
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
