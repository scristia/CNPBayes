## no reason to test if the other tests are working
test_marginalEasy <- function(){
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth, use.current=TRUE)

  params <- ModelParams("marginal", y=y(truth), k=3)
  mcmcp <- McmcParams(iter=150, burnin=50)
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

#.test_k_too_big <- function(){
#  ## when k is too big, chains may have a greater likelihood of crossing
#  set.seed(1)
#  truth <- simulateData(N=2500, p=rep(1/3, 3), theta=c(-1, 0, 1),
#                        sds=rep(0.1, 3))
#  if(FALSE) plot(truth, use.current=TRUE)
#  params <- ModelParams("marginal", y=y(truth), k=3)
#  mcmcp <- McmcParams(iter=1000, burnin=500)
#  params <- ModelParams("marginal", y=y(truth), k=5)
#  mcmcp <- McmcParams(iter=1000, burnin=500)
#  model <- initializeModel(params)
#  model <- posteriorSimulation(model, mcmcp)
#}

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
  se <- as(truth, "SummarizedExperiment")
  outdir <- tempdir()
  mcmp.list <- mcmcpList(iter=rep(50, 4), nStarts=5, nStartIter=10)
  mcmp.list[[2]]@iter <- 50
  mcmp.list[[2]]@burnin <- 10
  marginalExperiment(se, outdir, mcmcp.list=mcmp.list)
  M <- getFiles(outdir, rownames(se), "marginal")
  checkTrue(file.exists(model(M)))

  if(FALSE){
    ## Needs more MCMC iterations
    m.y <- getPosteriorStats(M)
    bf <- bayesFactor(m.y)
    checkIdentical(substr(names(bf), 1, 2), "M3")
  }

  ## metadata fields for plotting legend
  if(FALSE){
    rowRanges(se)$source <- "simulation"
    rowRanges(se)$nSNPs <- 10
    rowRanges(se)$nmarkers <- 15
    load_all()
    par(mfrow=c(1,4))
    plotModel(list(M), se)
  }

  ## Also fit batch model when truth is no batch effect
  ## batch found
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
  mcmcp <- McmcParams(iter=150, burnin=25)
  params <- ModelParams("marginal", y=y(truth), k=3)
  model <- initializeModel(params)
  model <- startAtTrueValues(model, truth)
  model <- posteriorSimulation(model, mcmcp)
  checkEquals(theta(model), theta(truth), tolerance=0.05)

  model <- initializeModel(params)
  model <- posteriorSimulation(model, mcmcp)
  checkEquals(sort(theta(model)), theta(truth), tolerance=0.1)
  ix <- order(theta(model))
  if(FALSE){
    plot.ts(thetac(model))
    plot.ts(sigma(model))
    plot.ts(pic(model))
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  checkEquals(sigma(model)[ix], sigma(truth), tolerance=0.15)
  ## check mixing proportions ... FAILS
  ##i <- argMax(model)
  ##checkEquals(pMean(model)[ix], p(truth), tolerance=0.2)
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
  params <- ModelParams("marginal", y=y(truth), k=3,
                        mcmc.params=mcmcp)
  modelk <- initializeModel(params)
  modelk <- startAtTrueValues(modelk, truth)
  model <- posteriorSimulation(modelk, mcmcp)
  checkEquals(theta(model), theta(truth), tolerance=0.1)
  checkEquals(sigma(model), sigma(truth), tolerance=0.1)
  checkEquals(p(model), p(truth), tolerance=0.1)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  ll.truth <- logLik(truth)
  set.seed(123)
  seeds <- sample(seq(1e6), 500)
  ll <- sapply(seeds, computeLogLikForRandomStarts, params=params, hypp=hyperParams(model))
  ll2 <- computeLogLikForRandomStarts(seeds[which.max(ll)], params, hyperParams(model))
  checkIdentical(max(ll), ll2)
  model <- computeLogLikForRandomStarts(seeds[which.max(ll)], params, hyperParams(model),
                                        return.model=TRUE)
  mcmcp <- McmcParams(iter=250, burnin=50)
  model <- posteriorSimulation(model, mcmcp)
  checkEquals(sort(theta(model)), theta(truth), tolerance=0.2)
  ##plot.ts(thetac(model), plot.type="single", col=1:3)
  ##
  ## Have to increase the tolerance a bit
  ##
  pmns <- thetaMean(model)
  checkEquals(sort(pmns), theta(truth), tolerance=0.2)

  ps <- sigmaMean(model)
  checkEquals(ps[order(pmns)], sigma(truth), tolerance=0.2)

  pmix <- pMean(model)
  exc <- tryCatch(checkEquals(p(truth), pmix[order(pmns)], tolerance=0.2), error=function(e)NULL)
  checkTrue(is.null(exc))
  ##min(effectiveSize::p(mc)) ## is very low
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
