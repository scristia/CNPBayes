test_batchEasy <- function(){
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

  mcmcp <- McmcParams(iter=1500, burnin=100)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=rep(letters[1:3], length.out=2500),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
  model <- posteriorSimulation(model, mcmcp)
  ##
  ## Start from true values
  ##
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  pmns <- thetaMean(model)
  j <- order(pmns[1,])
  checkEquals(pmns[, j], theta(truth), tolerance=0.03)

  ps <- sigmaMean(model)
  checkEquals(ps[, j], sigma(truth), tolerance=0.1)

  pmix <- pMean(model)
  checkEquals(pmix[j], p(truth), tolerance=0.02)

  ## initialize batch model with k=1
  params <- ModelParams("batch", y=y(truth), k=1,
                        batch=rep(letters[1:3], length.out=2500),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
}

test_unequalp_across_batch <- function(){
  library(oligoClasses)
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  p <- matrix(c(1/5, 1/3, 1-1/3-1/5,
                1/4, 1/4, 1-1/2,
                1/3, 1/5, 1-1/3-1/5), nbatch, k, byrow=TRUE)
  sds <- matrix(0.1, nbatch, k)
  ##trace(simulateBatchData, browser)
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             theta=means,
                             sds=sds,
                             p=p)

  ##
  ##
  ##
  tablez(truth)
  zz <- as.character(z(truth))
  zz[sample(which(batch(truth)=="a" & zz == 1), 150)] <- sample(c("2", "3"), 150, replace=TRUE)
  z(truth) <- factor(zz, levels=1:3)
  tablez(truth)

  zz <- as.character(z(truth))
  zz[sample(which(batch(truth)=="c" & zz == 3), 200)] <- sample(c("1", "2"), 200, replace=TRUE)
  z(truth) <- factor(zz, levels=1:3)

  truth <- simulateBatchData(N=2500, batch=batch(truth),
                             theta=means, sds=sds,
                             zz=z(truth), p=p)


  mcmcp <- McmcParams(iter=1000, burnin=1000)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=rep(letters[1:3], length.out=2500),
                        mcmc.params=mcmcp)
  bmod <- initializeBatchModel(params)
  bmodel <- posteriorSimulation(bmod, mcmcp)

  if(FALSE){
    op <- par(mfrow=c(1,2), las=1)
    plot(truth, use.current=TRUE)
    plot(bmodel)
    par(op)
  }
}


test_selectK_batch_easy <- function(){
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
  ##
  ## Evaluate at different K
  ##
  mcmcp <- McmcParams(iter=c(1000, 1000, 2000, 2000),  burnin=rep(300, 4),
                      thin=c(1, 1, 2, 2))
  mmodels <- fitMixtureModels(se, mcmcp, K=1:4)
  bicstat <- sapply(mmodels, bic)
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

test_batch_moderate <- function(){
  library(oligoClasses)
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
  se <- as(truth, "SummarizedExperiment")
  mcmcp <- McmcParams(iter=1000, burnin=1000)
  params <- ModelParams("batch", y=y(truth), k=3, batch=batch(truth),
                        mcmc.params=mcmcp)
  modelk <- initializeBatchModel(params)
  modelk <- posteriorSimulation(modelk, mcmcp)

  pmns <- thetaMean(modelk)
  j <- order(pmns[1,])
  checkEquals(pmns[, j], theta(truth), tolerance=0.03)

  ps <- sigmaMean(modelk)
  checkEquals(ps[, j], sigma(truth), tolerance=0.1)

  pmix <- pMean(modelk)
  checkEquals(pmix[j], p(truth), tolerance=0.08)

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
  mcmcp <- McmcParams(iter=2000, burnin=1000, thin=2)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
  batch(model) <- collapseBatch(model)
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE, xlim=c(-2,1))
    CNPBayes::plot(model, xlim=c(-2,1))
    par(op)
    ## we get the mixture probabilities backwards for comp 2 and 3
  }
  mc <- mcmcChains(model)
  cbind(colMeans(theta(mc)), as.numeric(theta(truth)))
  cbind(colMeans(sigma(mc)), as.numeric(sigma(truth)))
  cbind(colMeans(p(mc)), p(truth))
  ## this fails
  checkEquals(colMeans(p(mc)), p(truth), tolerance=0.1)
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

test_hard3 <- function(){
  ##
  ## Repeat above, but with smaller simulated variance (sd=0.2)
  ##
  library(oligoClasses)
  truth <- hardTruth(0.005, s=0.1)
  table(z(truth), batch(truth))
  if(FALSE) CNPBayes::plot(truth, use.current=TRUE)
  se <- as(truth, "SummarizedExperiment")
  if(FALSE) hist(oligoClasses::copyNumber(se), breaks=1000, col="gray", border="gray")
  ##
  ## Use defaults
  ##
  mcmcp <- McmcParams(iter=5000, burnin=1000, thin=5)
  mmodel <- fitMixtureModels(se, mcmcp, K=3)[[1]]
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE)
    CNPBayes::plot(mmodel)
    par(op)
    mc <- mcmcChains(mmodel)
    plot.ts(sigma(mc), col="gray")
    plot.ts(theta(mc), col="gray")
    cbind(colMeans(sigma(mc)), as.numeric(sigma(truth)))
  }
  ##mc <- mcmcChains(mmodel)
  pmns <- thetaMean(mmodel)
  j <- order(pmns[1,])
  ps <- sigmaMean(mmodel)[, j]
  pmix <- pMean(mmodel)[j]

  checkEquals(pmns[, j], theta(truth), tolerance=0.02)
  checkEquals(ps, sigma(truth), tolerance=0.15)
  checkEquals(pmix, p(truth), tolerance=0.04)
}

.test_hard4 <- function(){
  truth <- hardTruth(0.01, s=0.15)
  mcmcp <- McmcParams(iter=2000, burnin=1000, thin=2)
  table(z(truth), batch(truth))
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
  batch(model) <- collapseBatch(model)
  model <- posteriorSimulation(model, mcmcp)
  ##
  ## The standard deviation of the second component is nearly twice as
  ## high as the standard deviation of the third component
  ##
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE, xlim=c(-2,1))
    CNPBayes::plot(model, xlim=c(-2,1))
    par(op)
    mc <- mcmcChains(model)
    plot.ts(sigma(mc), col="gray")
    plot.ts(theta(mc), col="gray")
    cbind(colMeans(sigma(mc)), as.numeric(sigma(truth)))

    probs <- probz(model)
    maxprob <- apply(probs, 1, max)
    plot(maxprob~jitter(as.integer(z(truth)), amount=0.1), pch=".")
  }
  mc <- mcmcChains(model)
  checkEquals(colMeans(theta(mc)), as.numeric(theta(truth)), tolerance=0.05)
  checkEquals(colMeans(sigma(mc)), as.numeric(sigma(truth)), tolerance=0.05)
  checkEquals(colMeans(p(mc)), as.numeric(p(truth)), tolerance=0.03)
}




test_chr4_locus <- function(){
  ##
  ## Without the multiple starts, the chain for the batch model gets
  ## stuck in a local model where one of the components has zero
  ## observations
  ##
  library(GenomicRanges)
  library(oligoClasses)
  truth <- readRDS("~/Software/CNPBayes/inst/extdata/unequal_mix_model.rds")
  mcmcp <- McmcParams(iter=1000, burnin=100, nStarts=20, nStartIter=200)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  bmod <- initializeBatchModel(params)
  bmodel <- posteriorSimulation(bmod, mcmcp)

  pmns <- thetaMean(bmodel)
  j <- order(pmns[1, ])
  ps <- sigmaMean(bmodel)[, j]
  pmix <- pMean(bmodel)[j]

  checkEquals(pmns[, j], theta(truth), tolerance=0.02)
  checkEquals(pmix, as.numeric(p(truth)), tolerance=0.01)
  checkTrue(logpotential(bmodel) > -2000)

  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(bmodel)
    par(op)

    trace(tracePlot, signature="BatcModel")
    tracePlot(bmodel, "theta")
    tracePlot(bmodel, "p")
    tracePlot(bmodel, "sigma")
  }
}




test_cnp707 <- function(){
  library(GenomicRanges)
  library(oligoClasses)
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472.rds")
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  j <- subjectHits(findOverlaps(rowData(se472), rowData(se.ea)))[1]
  se472 <- se.ea[j, ]

  mcmcp <- McmcParams(iter=c(1000, 1000, 1000), burnin=c(100, 200, 200))
  bmodels <- fitMixtureModels(se472, mcmcp, K=1:3)
  bicstat <- sapply(bmodels, bic)
  ## this is just a batch effect. Number of components is 1
  checkTrue(which.min(bicstat)==1)
  ## marginal model produces incorrect inference
}

test_cnp472 <- function(){
  ##
  ## Not sure that the batch model is capable of correct inference
  ## when the number of observations in one component is very small.
  ## Not enough strength borrowed between batches.
  ##
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472.rds")
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  j <- subjectHits(findOverlaps(rowData(se472), rowData(se.ea)))[1]
  se472 <- se.ea[j, ]
  truth <- readRDS("~/Software/CNPBayes/inst/extdata/cnp472_model.rds")

  b <- collapseBatch(se472)
  mcmcp <- McmcParams(iter=10000, burnin=100, thin=10)
  params <- ModelParams("batch", y=copyNumber(se472)[1, ], k=3,
                         batch=b,
                         mcmc.params=mcmcp)
  bmod <- initializeBatchModel(params, zz=z(truth))
  ##dataMean(bmod) ## note a lot of NaN's because of batches with no observations
  bmodel <- posteriorSimulation(bmod, mcmcp)
  zz <- map(z(bmodel))
  pmns <- thetaMean(bmodel)
  n.hp <- tablez(bmodel)
  ##weightedMean <- colSums(tablez(bmodel)*pmns)/colSums(tablez(bmodel))
  j <- order(pmns[1, ])
  if(FALSE){
    thetaMean(bmodel)
    sigmaMean(bmodel)
    ## interestingly most of the incorrect inference is with respect
    ## to the middel component.
    op <- par(mfrow=c(1,2), las=1)
    plot(truth, use.current=T)
    plot(bmodel)
    par(op)

    op <- par(mfrow=c(4, 4), las=1)
    tracePlot(bmodel, "theta")
    par(op)
    tracePlot(bmodel, "p")
    tracePlot(bmodel, "sigma")

    zz <- map(bmodel)
    table(zz, z(truth))
    z.marginal <- setNames(z(truth), colnames(se472))
    z.marginal <- z.marginal[colnames(sei)]

    zz <- setNames(zz, colnames(se472))
    sei <- readRDS("~/Labs/ARIC/AricCNPData/data/intensity_472.rds")
    se472.2 <- se472[, colnames(se472) %in% colnames(sei)]
    sei <- sei[, colnames(se472.2)]
    zz <- zz[colnames(sei)]
    df <- data.frame(a=log2(assays(sei)[["A"]])[1,],
                     b=log2(assays(sei)[["B"]])[1,],
                     z=zz,
                     z.marginal=as.integer(z.marginal))
    ## acceptable
    library(lattice)
    xyplot(b~a, df, pch=20, cex=0.3, color=c("gray", "black", "royalblue")[df$z],
           panel=function(x, y, color, ..., subscripts){
             panel.xyplot(x, y, col="gray", ...)
             panel.points(x[color=="black"], y[color=="black"], col="black", ...)
             panel.points(x[color=="royalblue"], y[color=="royalblue"], col="royalblue", ...)
           })
    xyplot(b~a, df, pch=20, cex=0.3, color=c("gray", "black", "royalblue")[df$z.marginal],
           panel=function(x, y, color, ..., subscripts){
             panel.xyplot(x, y, col="gray", ...)
             panel.points(x[color=="black"], y[color=="black"], col="black", ...)
             panel.points(x[color=="royalblue"], y[color=="royalblue"], col="royalblue", ...)
           })
  }
}
