.test_batchFull <- function(){
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
                             p=c(1/5, 1/3, 1-1/3-1/5))
  mcmcp <- McmcParams(iter=1000, burnin=100)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
  model <- posteriorSimulation(model, mcmcp)
  checkIdentical(batch(model), batch(truth))
  checkIdentical(y(model), y(truth))
  if(FALSE){
    op <- par(mfrow=c(1,2), las=1)
    plot(truth, use.current=TRUE)
    plot(model)
    par(op)
    op <- par(mfrow=c(1,3), las=1)
    tracePlot(model, "theta", col=1:3)
    par(op)
  }
  th <- thetaMean(model)
  j <- order(th[1,])
  ##
  ## Why would a label switch within one batch?  - in batch 'c', the
  ##   third component has label '3' while in batches 'a' and 'b' the
  ##   third component has label '1'. We assume that the proportion of
  ##   samples of a given copy number is the same between batches.
  ##   The label should carry the same meaning between batches.
  ##
  checkEquals(th[, j], theta(truth), tolerance=0.1) #!!FAILS

  ##
  ## Informative prior on mu_h (small tau2.0)
  ##
  mcmcp <- McmcParams(iter=1000, burnin=100)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params, hypp=HyperparametersBatch(k=3, tau2.0=1))
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2), las=1)
    plot(truth, use.current=TRUE)
    plot(model)
    par(op)
    op <- par(mfrow=c(1,3), las=1)
    tracePlot(model, "theta", col=1:3)
    par(op)
  }
  th <- thetaMean(model)
  j <- order(th[1,])
  checkEquals(th[, j], theta(truth), tolerance=0.1) #!!FAILS

  ##
  ## What happens when the fraction in each component truly differs
  ## between batches, and mu is small as above?
  ##
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             theta=means,
                             sds=sds,
                             p=c(1/10, 2/10, 7/10))
  zz <- z(truth)
  table(zz[batch(truth)=="c"])
  ## reverse the percentages for 3rd component and 1st component
  L <- table(batch(truth))[["c"]]
  z.batchc <- c(ceiling(rep(1, L*7/10)),
                rep(2, floor(L*2/10)),
                rep(3, ceiling(L*1/10)))
  zz[batch(truth) == "c"] <- z.batchc
  dataMean(truth) <- computeMeans(truth)
  dataPrec(truth) <- 1/computeVars(truth)
  theta(truth) <- dataMean(truth)
  sigma2(truth) <- 1/dataPrec(truth)
  ## would just be the marginal fraction across batches
  p(truth) <- as.numeric(table(z(truth))/length(y(truth)))

  ## Non-informative prior on mu_h (big tau2.0)
  mcmcp <- McmcParams(iter=100, burnin=0)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
  model <- CNPBayes:::startAtTrueValues(model, truth)
  model <- posteriorSimulation(model, mcmcp)
  checkEquals(theta(model), theta(truth), tolerance=0.1)

  ## random starts
  model <- initializeBatchModel(params)
  model <- posteriorSimulation(model, McmcParams(iter=1000, burnin=100))
  th <- matrix(thetac(model)[CNPBayes:::argMax(model), ], 3, 3)
  ## Fails even with noninformative prior
  checkEquals(th, theta(truth), tolerance=0.1)
  if(FALSE){
    op <- par(mfrow=c(1,2), las=1)
    plot(truth, use.current=TRUE)
    plot(model)
    par(op)
    op <- par(mfrow=c(1,3), las=1)
    tracePlot(model, "theta", col=1:3)
    par(op)
  }

  ## Informative prior on mu_h (small tau2.0)
  mcmcp <- McmcParams(iter=1000, burnin=100)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
  model <- posteriorSimulation(model, mcmcp)
  th <- matrix(thetac(model)[CNPBayes:::argMax(model), ], 3, 3)
  ## Fails even with noninformative prior
  checkEquals(th, theta(truth), tolerance=0.1)
}

.test_selectK_batch <- function(){
  ## Need to replace setParallelization with source code
  ## -investigate BiocParallel
  library(GenomicRanges)
  set.seed(1000)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             p=c(1/10, 1/5, 1-0.1-0.2),
                             theta=means,
                             sds=sds)
  se <- as(truth, "SummarizedExperiment")
  outdir <- tempdir()
  ##
  ## Evaluate at different K
  ##
  mcmp.list <- mcmcpList(iter=rep(100, 4))
  batchExperiment(se, outdir, mcmcp.list=mcmp.list)
  B <- getFiles(outdir, rownames(se), "batch")
  models <- readRDS(model(B))
  ##trace(plotModel, browser)
  par(mfrow=c(1,4), las=1)
  load_all()
  plotModel(list(B), se, xlim=c(-3,2), xaxt="s")
  ##
  ## 4-component has higher bic, yet the max a post. estimate has only 3 components
  ## --> inference from map would be the same
  ##
  checkIdentical(which.min(bicstat), 3L)
  if(FALSE) {
    op <- par(mfrow=c(1,2), las=1)
    CNPBayes::plot(truth, use.current=TRUE)
    CNPBayes::plot(mmodels[[3]])
    par(op)
  }
}

.test_batch_moderate <- function(){
  set.seed(100)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95,
                    -0.41, -0.4, -.395,
                    -0.1, 0, 0.05), nbatch, k, byrow=FALSE)
  sds <- matrix(0.15, nbatch, k)
  sds[,1] <- 0.3
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             p=c(1/10, 1/5, 1-0.1-0.2),
                             theta=means,
                             sds=sds)
  mcmcp <- McmcParams(iter=1000, burnin=100, nStarts=10, nStartIter=50)
  modelk <- posteriorSimulation(model, mcmcp)
  checkEquals(theta(modelk), theta(truth), tolerance=0.1)

  pmns <- thetaMean(modelk)
  j <- order(pmns[1,])
  checkEquals(pmns[, j], theta(truth), tolerance=0.03)

  ps <- sigmaMean(modelk)
  checkEquals(ps[, j], sigma(truth), tolerance=0.1)

  pmix <- pMean(modelk)
  checkEquals(pmix[j], p(truth), tolerance=0.09)

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
    ##trace(.plotBatch, browser)
    ##.plotBatch(truth, use.current=TRUE)
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
  p1 <- prop_comp1
  p2 <- 20*p1
  p3 <- 1-p1-p2
  N <- 8e3
  truth <- simulateBatchData(N,
                             theta=means,
                             sds=sds,
                             batch=rep(letters[1:3], length.out=N),
                             p=c(p1, p2, p3))
}


.test_hard1 <- function(){
  truth <- hardTruth(0.005, s=0.3)
  table(z(truth), batch(truth))
  if(FALSE) plot(truth, use.current=TRUE)
  ##
  ## Use defaults
  ##
  mcmcp <- McmcParams(iter=2000, burnin=1000,
                      nStarts=20, nStartIter=100, thin=2)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE, xlim=c(-2,1))
    CNPBayes::plot(model, xlim=c(-2,1))
    par(op)
    par(mfrow=c(1,3), las=1)
    tracePlot(model, "theta", col=1:3)
    ## we get the mixture probabilities backwards for comp 2 and 3
  }
  pmix <- colMeans(pic(model))
  ## this fails
  checkEquals(sort(pmix), p(truth), tolerance=0.1)
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
  model <- initializeBatchModel(params)
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

.test_cnp472 <- function(){
  library(GenomicRanges)
  library(oligoClasses)
  if(FALSE){
    cnp472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472.rds")
    se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
    j <- subjectHits(findOverlaps(rowData(cnp472), rowData(se.ea)))[1]
    ## just the european ancestry CN summaries
    se472 <- se.ea[j, ]
    b <- collapseBatch(se472)
    se472$batch <- b
    saveRDS(se472, file="~/Software/CNPBayes/inst/extdata/se_cnp472_EA.rds")
  }
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472_EA.rds")
  mcmcp <- McmcParams(iter=c(500, 2000, 2000, 2000),
                      burnin=c(50, 100, 100, 100),
                      thin=c(1, 2, 2, 2),
                      nStarts=20,
                      nStartIter=200)
  hypp <- CNPBayes:::HyperparametersBatch()
  hplist <- CNPBayes:::HyperParameterList(hypp, K=1:4)
  mplist <- ModelParamList(hypp, K=1:4, data=copyNumber(se472)[1, ],
                           batch=se472$batch, mcmcp=mcmcp)
  bmodlist <- foreach(hypp=hplist, param=mplist) %do% {
    initializeBatchModel(params=param, hypp=hypp)
  }
  bmodels <- foreach(k=1:4, model=bmodlist) %do% posteriorSimulation(model, mcmcp[k])
  mcmcp2 <- McmcParams(iter=1000, burnin=100)
  x <- lapply(bmodels, computeMarginalPr, mcmcp=mcmcp2)
  xx <- do.call(cbind, lapply(x, rowMeans))
  post.range <- unlist(lapply(x, posteriorRange))
  marginals <- computeMarginal(xx)
  marginals[ post.range > 5 ] <- -Inf
  ##bicstat <- sapply(bmodels, bic)
  ## this is just a batch effect. Number of components is 1
  ##checkTrue(which.min(bic(bmodels))==1)
  checkTrue(which.max(marginals) <= 3)
  if(FALSE){
    ## I think 472 is probable for 2 or 3 states. Batch model might
    ## help reduce number of false positives
    tmp=readRDS("~/Labs/ARIC/AricCNPData/data/intensity_472.rds")
    se472.2=se472[, colnames(se472)%in%colnames(tmp)]
    tmp <- tmp[, colnames(se472.2)]
    all.equal(colnames(tmp), colnames(se472.2))
    plot(assays(tmp)[["A"]][1, ], assays(tmp)[["B"]][1, ], pch=".", col="gray")
    y <- assays(se472.2)[["medr"]][1, ]
    index=which(assays(se472.2)[["medr"]] < -250)
    points(assays(tmp)[["A"]][1, index], assays(tmp)[["B"]][1, index], pch=".", col="black")
  }
  ## marginal model produces incorrect inference
}

.test_chr4_locus <- function(){
  ##
  ## Without the multiple starts, the chain for the batch model gets
  ## stuck in a local model where one of the components has zero
  ## observations
  ##
  library(GenomicRanges)
  library(oligoClasses)
  truth <- readRDS("~/Software/CNPBayes/inst/extdata/unequal_mix_model.rds")
  mcmcp <- McmcParams(iter=1000, burnin=100, nStarts=20, nStartIter=100)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  hypp <- CNPBayes:::HyperparametersBatch(k=3, tau2.0=1000)
  bmodel <- initializeBatchModel(params, hypp=hypp)
  bmodel <- posteriorSimulation(bmodel, mcmcp)

  pmns <- thetaMean(bmodel)
  j <- order(pmns[1, ])
  ps <- sigmaMean(bmodel)[, j]
  pmix <- pMean(bmodel)[j]

  checkEquals(pmns[, j], theta(truth), tolerance=0.02)
  checkEquals(pmix, as.numeric(p(truth)), tolerance=0.01)

  checkTrue(CNPBayes:::logpotentialc(bmodel)[CNPBayes:::argMax(bmodel)] > -2000)

  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(bmodel)
    par(op)

    ##trace(tracePlot, signature="BatcModel")
    par(mfrow=c(3,3))
    tracePlot(bmodel, "theta", col=1:3)
    tracePlot(bmodel, "p")
    tracePlot(bmodel, "sigma")
  }
}


.test_selectK_moderate <- function(){
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
  hypp <- CNPBayes:::HyperparametersMarginal()
  hplist <- CNPBayes:::HyperParameterList(hypp, K=1:4)
  mplist <- ModelParamList(hypp, K=1:4, data=y(truth),
                           mcmcp=mcmcp)
  modlist <- foreach(hypp=hplist, param=mplist) %do% {
    CNPBayes:::initializeModel(params=param, hypp=hypp)
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


.test_selectK_hard <- function(){
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
  hypp <- CNPBayes:::HyperparametersMarginal()
  hplist <- CNPBayes:::HyperParameterList(hypp, K=1:4)
  mplist <- ModelParamList(hypp, K=1:4, data=y(truth),
                           mcmcp=mcmcp)
  modlist <- foreach(hypp=hplist, param=mplist) %do% {
    CNPBayes:::initializeModel(params=param, hypp=hypp)
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
