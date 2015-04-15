test_batchEasy <- function(){
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
                             p=c(1/5, 1/3, 1-1/3-1/5))
  if(FALSE){
    while(all(z(model) > 0) && i < 1000){
      theta(model) <- .Call("update_theta_batch", model)
      sigma2(model) <- .Call("update_sigma2_batch", model)
      p(model) <- .Call("update_p_batch", model)
      probs = .Call("update_multinomialPr_batch", model)
      all(rowSums(probs) == 1)
      mu(model) <- .Call("update_mu_batch", model)
      tau2(model) <- .Call("update_tau2_batch", model)
      nu.0(model) <- .Call("update_nu0_batch", model)
      sigma2.0(model) <- .Call("update_sigma20_batch", model)
      z(model) <- .Call("update_z_batch", model)
      dataMean(model) <- .Call("compute_means_batch", model)
      dataPrec(model) <- .Call("compute_prec_batch", model)
      table(z(model), batch(model))
      theta(model)
      min(probs)
      tau(model)
      i <- i+1
    }
  }
  ##model <- posteriorSimulation(model)
  yy <- y(truth)
  checkIdentical(yy, yy[order(batch(truth))])
  mcmcp <- McmcParams(iter=50, burnin=0)
  set.seed(123)
  model <- BatchModel(y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp)
  model <- startAtTrueValues(model, truth)
  checkIdentical(batch(model), batch(truth))
  checkIdentical(y(model), y(truth))
  checkIdentical(theta(model), theta(truth))
  checkIdentical(sigma(model), sigma(truth))
  checkIdentical(p(model), p(truth))
  checkIdentical(z(model), z(truth))
  iter(model) <- 50
  model2 <- CNPBayes:::.runMcmc(model)
  checkEquals(theta(model2), theta(truth), tolerance=0.1)
  if(FALSE){
    iter(model) <- 500
    model2 <- CNPBayes:::.runMcmc(model)
    checkEquals(theta(model2), theta(truth), tolerance=0.1)
  }
  ##paramUpdates(model)[c("nu.0", "sigma2.0")] <- 0L
  ##paramUpdates(model)[c("sigma2.0")] <- 0L
  if(FALSE){
    model <- .Call("mcmc_batch", model, mcmcParams(model))
    set.seed(123)
    .Call("update_sigma20_batch", model)
    set.seed(123)
    .updateSigma2.0Batch(model2)
    ## problem lies with nu.0 or sigma2.0...
    eta.0 <- 1800
    m2.0 <- 1/60
    a <- 0.5*eta.0
    b <- 0.5*eta.0*m2.0
    x <- rgamma(1000, a, rate=1/b)
    ix <- 1/x
    hist(sqrt(1/x), breaks=100)
    s <- mean(sqrt(1/x))
    theta <- rnorm(1000, -1, s)
    hist(theta, breaks=100, xlim=c(-2, 1.5))
  }
  set.seed(1)
  mcmcp <- McmcParams(iter=300, burnin=100, nStarts=20)
  model <- BatchModel(y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp)
  model <- posteriorSimulation(model)
  i <- order(theta(model)[1, ])
  checkEquals(theta(model)[, i], theta(truth), tolerance=0.1)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
    par(mfrow=c(1,3))
    tracePlot(model, "theta", col=1:3)
    plot.ts(muc(model), col=1:3)
  }

  ##
  ## check that the marginal density estimates are consistent for the
  ## K=3 model
  ##
  marginaly_k3 <- computeMarginalProbs(model, mcmcp)
  spread <- diff(range(marginaly_k3))
  checkTrue(spread < 50)
  ##
  ## Select K
  ##
  marginaly <- computeMarginalEachK2(y(truth), batch(truth), K=1:4,
                                     mcmcp=McmcParams(iter=500, burnin=100, nStarts=20),
                                     MAX.RANGE=100)
  K <- which.max(marginaly)
  checkIdentical(as.integer(K), 3L)
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
  sds[, 1] <- 0.3
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             p=c(1/10, 1/5, 1-0.1-0.2),
                             theta=means,
                             sds=sds)
  mcmcp <- McmcParams(iter=150, burnin=0)
  hypp <- HyperparametersBatch(m2.0=1/60, eta.0=1800, k=3)
  model <- BatchModel(data=y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp,
                      hypp=hypp)
  model <- startAtTrueValues(model, truth)
  model <- posteriorSimulation(model, mcmcp)
  i <- order(theta(model)[1, ])
  checkEquals(theta(model)[, i], theta(truth), tolerance=0.1)
  ## random starts
  mcmcp <- McmcParams(iter=100, burnin=100, nStarts=20)
  model <- BatchModel(data=y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp,
                      hypp=hypp)
  model <- posteriorSimulation(model)
  i <- order(theta(model)[1, ])
  checkEquals(theta(model)[, i], theta(truth), tolerance=0.1)
  checkEquals(mu(model)[i], mu(truth), tolerance=0.15)
  if(FALSE){
    zz <- as.integer(z(truth))
    ps <- c(mean(zz==1), mean(zz==2), mean(zz==3))
    modelk <- model
    plot.ts(pic(modelk), plot.type="single")
    abline(h=p(truth))
    par(mfrow=c(1,3))
    tracePlot(modelk, "theta", col=1:3)
    abline(h=theta(truth))
    plot.ts(sigmac(modelk), plot.type="single")
    abline(h=sigma(truth))
    op <- par(mfrow=c(1,2),las=1)
    ##trace(.plotBatch, browser)
    ##.plotBatch(truth, use.current=TRUE)
    CNPBayes::plot(truth, use.current=TRUE)
    CNPBayes::plot(modelk, use.current=TRUE)
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

test_hard3 <- function(){
  ##
  ## Repeat above, but with smaller simulated variance (sd=0.2)
  ##
  library(oligoClasses)
  library(GenomicRanges)
  truth <- hardTruth(0.005, s=0.1)
  table(z(truth), batch(truth))
  if(FALSE) CNPBayes::plot(truth, use.current=TRUE)
  se <- as(truth, "SummarizedExperiment")
  if(FALSE) hist(oligoClasses::copyNumber(se), breaks=1000, col="gray", border="gray")
  ##
  ## Use defaults
  ##
  mcmcp <- McmcParams(iter=100, burnin=0)
  modelk <- BatchModel(data=y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp,
                       HyperparametersBatch(k=3, m2.0=1/60, eta.0=1800))
  modelk <- startAtTrueValues(modelk, truth)
  mmodel <- posteriorSimulation(modelk, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, use.current=TRUE)
    ##trace(.plotBatch, browser)
    .plotBatch(mmodel)
    par(op)
    i <- order(theta(mmodel)[1, ])
    theta(mmodel)[, i]
    theta(truth)
    sigma(model)[, i]
    sigma(truth)
    mc <- mcmcChains(mmodel)
    par(mfrow=c(1,3))
    tracePlot(mmodel, "sigma", col=1:3)
    tracePlot(mmodel, "theta", col=1:3)
    plot.ts(sigma(mc), col="gray", )
    plot.ts(theta(mc), col="gray")
    plot.ts(mu(mc), plot.type="single", col=1:3)
    cbind(colMeans(sigma(mc)), as.numeric(sigma(truth)))
  }
  ##mc <- mcmcChains(mmodel)
  pmns <- thetaMean(mmodel)
  j <- order(pmns[1,])
  ps <- sigmaMean(mmodel)[, j]
  pmix <- pMean(mmodel)[j]
  checkEquals(pmns[, j], theta(truth), tolerance=0.04)
  checkEquals(ps, sigma(truth), tolerance=0.15)
  checkEquals(pmix, p(truth), tolerance=0.04)
  ##
  ## With multiple starts we can find the mode
  ##
  mcmcp <- McmcParams(iter=200, burnin=100, nStarts=20)
  modelk <- BatchModel(data=y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp,
                       HyperparametersBatch(k=3, m2.0=1/60, eta.0=1800, tau2.0=1000))
  mmodel <- posteriorSimulation(modelk, mcmcp)
  pmns <- thetaMean(mmodel)
  j <- order(pmns[1,])
  ps <- sigmaMean(mmodel)[, j]
  pmix <- pMean(mmodel)[j]
  checkEquals(pmns[, j], theta(truth), tolerance=0.04)
  checkEquals(ps, sigma(truth), tolerance=0.15)

  checkEquals(pmix, p(truth), tolerance=0.04)
}

.test_number_batches <- function(){
  library(foreach)
  library(GenomicRanges)
  library(devtools)
  library(oligoClasses)
  library(matrixStats)
  library(CNPAric)
  library(CNPBayes)
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  cnpids <- rownames(se.ea)
  library(HopkinsHPC)
  NW <- setParallelization(1)
  if(is.null(tryCatch(se.ea[1, ], error=function(e) return(NULL)))) stop()
  se <- se.ea
  outdir <- tempdir()
  B <- getFiles(outdir, rownames(se), model="batch")
  saved.date <- file.info(model(B))$mtime
  d <- difftime(Sys.time(), saved.date, unit="days")
  ##index <- which(d > 1.75)
  index <- 18
  ##index <- seq_len(nrow(se))
  se <- se[index, ]
  B <- B[index]
  model.files <- model(B)
  batch.files <- paste0(dirname(model(B)), "/", rownames(se), "_batch.rds")
  mcmcp <- McmcParams(iter=1000, thin=1, nStarts=20, burnin=200)
  ##bt <- readRDS(batch.files)
  ##unique(bt)
  cn <- copyNumber(se)[1, ]
  if(FALSE){
    bt <- collapseBatch(cn, se$plate)
    bt2 <- collapseBatch(cn, se$plate, THR=0.2)
    bt3 <- collapseBatch(cn, se$plate, THR=0.05)
  }
  bt4 <- collapseBatch(cn, se$plate, THR=0.005)
  mcmcp <- McmcParams(iter=50, thin=1, nStarts=1, burnin=0)
  hypp <- Hyperparameters("batch", k=3)
  bm <- BatchModel(data=cn, k=3, batch=bt4, mcmc.params=mcmcp)
  if(FALSE){
    system.time(tmp <- posteriorSimulation(bm))
    bm2 <- BatchModel(data=cn, k=3, batch=bt3, mcmc.params=mcmcp)
    system.time(tmp2 <- posteriorSimulation(bm2))
    bm3 <- BatchModel(data=cn, k=3, batch=bt4, mcmc.params=mcmcp)
    system.time(tmp3 <- posteriorSimulation(bm3))
    bm3 <- BatchModel(data=cn, k=3, batch=bt4, mcmc.params=mcmcp)
    system.time(tmp4 <- posteriorSimulation(bm3))
  }
  ##
  ## Alternative. Find modes with marginal model
  ##
  ##   - use these to initialize batch model in region of high
  ##   posterior probability
  ##
  ##
  ## what updates are the most time consuming
  ##
  u <- paramUpdates(mcmcp)
  u[1:7] <- 0L
  paramUpdates(bm) <- u
  ## updating only z
  system.time(tmp4 <- posteriorSimulation(bm)) ## 2.25
  tz <- table(z(tmp4))


  system.time(tmp <- .Call("update_z_batch", bm))
  set.seed(123)
  p <- .Call("update_z_batch", bm)

  p2 <- .Call("update_z_batch", bm)

  system.time(tmp4 <- posteriorSimulation(bm)) ## 2.25
  system.time(tmp4 <- posteriorSimulation(bm)) ## 3.13


  ## updating only theta
  u[1:8] <- 0L
  u["theta"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 18.89
  u[1:8] <- 0L
  u["sigma2"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 19
  u[1:8] <- 0L
  u["p"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 15.5
  u[1:8] <- 0L
  u["mu"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ##18
  u[1:8] <- 0L
  u["tau2"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 15.5
  u[1:8] <- 0L
  u["nu.0"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 15.5
  u[1:8] <- 0L
  u["sigma2.0"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 15.1

}
