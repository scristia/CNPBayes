test_computemeans <- function(){
  dir <- system.file("unitTests", package="CNPBayes")
  testdat <- readRDS(file.path(dir, "test_data.rds"))
  model <- BatchModel(data=testdat$y, batch=testdat$b)
  mns <- CNPBayes:::computeMeans(model)
  mns2 <- matrix(as.numeric(CNPBayes:::.computeMeansBatch(model)), 4, 2)
  checkEquals(mns, mns2)
}

test_batchEasy <- function(){
  library(oligoClasses)
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  N <- 1500
  truth <- simulateBatchData(N=N,
                             batch=rep(letters[1:3], length.out=N),
                             theta=means,
                             sds=sds,
                             p=c(1/5, 1/3, 1-1/3-1/5))
  ## The data in object 'truth' is ordered by batch automatically
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
  if(FALSE){
    model <- .Call("mcmc_batch", model, mcmcParams(model))
    set.seed(123)
    .Call("update_sigma20_batch", model)
    set.seed(123)
    .updateSigma2.0Batch(model2)
    ## problem lies with nu.0 or sigma2.0...
    eta.0 <- 1800/100
    #eta.0 <- 180
    m2.0 <- 1/(60/100)
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
  mcmcp <- McmcParams(iter=300, burnin=100, nStarts=5)
  model <- BatchModel(y(truth), batch=batch(truth), k=3, mcmc.params=mcmcp)
  model <- posteriorSimulation(model)
  i <- order(theta(model)[1, ])
  checkEquals(theta(model)[, i], theta(truth), tolerance=0.1)
  if(FALSE){
    ## plot the true posterior density (left) and the Gibbs
    ## approximation (right)
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
    op <- par(mfrow=c(1,3))
    tracePlot(model, "theta", col=1:3)
    par(op)
    ## trace plots for the cross-batch means
    plot.ts(muc(model), plot.type="single", col=1:3)
  }
}

test_kbatch <- function(){
  library(oligoClasses)
  set.seed(123)
  k <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), 3, k, byrow=FALSE)
  sds <- matrix(0.1, 3, k)
  N <- 1500
  truth <- simulateBatchData(N=N,
                             batch=rep(letters[1:3], length.out=N),
                             theta=means,
                             sds=sds,
                             p=c(1/5, 1/3, 1-1/3-1/5))
  if(FALSE){
    ##
    ## With 10k iterations and 1k iterations burnin, all models
    ## provide consistent estimates of the marginal likelihood and the
    ## K=4 model is best.  This unit test is too time consuming to be
    ## run during each package update.
    ##
    fit <- computeMarginalLik(y(truth), batch(truth), K=1:4,
                              burnin=1000,
                              T2=1000, T=10000,
                              nchains=3)
    prz <- probz(fit$models[[4]])
    cn <- map(fit$models[[4]])
    plot(r, cn, pch=20, cex=0.3)
    ## The 4th state is the map estimate in 2 individuals, so even if
    ## K = 4 has a higher marginal likelihood the copy number
    ## inference is not effected.
    trace(cnProbability, browser)
    prz <- cnProbability(prz, 4)
    plot(jitter(prz, amount=0.05), jitter(cn, amount=0.05), pch=20, cex=0.3)
    table(cn)

    pz <- cnProbability(probz(fit$models[[4]]), 4)
    r <- y(fit$models[[4]])
    plot(r, pz, pch=".")
    checkTrue(k(orderModels(fit))[1] == 3)
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
  model <- posteriorSimulation(model)
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
  mmodel <- posteriorSimulation(modelk)
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
  mmodel <- posteriorSimulation(modelk)
  pmns <- thetaMean(mmodel)
  j <- order(pmns[1,])
  ps <- sigmaMean(mmodel)[, j]
  pmix <- pMean(mmodel)[j]
  checkEquals(pmns[, j], theta(truth), tolerance=0.04)
  checkEquals(ps, sigma(truth), tolerance=0.15)
  checkEquals(pmix, p(truth), tolerance=0.04)
}

##
## Good example of heterogeneity in batch means, and some batches do
## not have the homozygous deletion because it is rare
##
test_missingcomponent <- function(){
  dir <- system.file("unitTests/data", package="CNPBayes")
  testdat <- readRDS(file.path(dir, "ea992.rds"))
  x <- BatchModel(data=testdat$cn, batch=testdat$batch, k=3)
  ##
  ## Design test to throw exception
  ##
  true_theta <- matrix(c(rep(-1.6, 7),
                         c(-0.5, -0.5, -0.5, -0.55, -0.42, -0.43, -0.43),
                         c(0, -0.02, 0, -0.01, 0.1, 0.1, 0.07)),
                       7, 3)
  truemu <- c(-1.6, -0.5, 0.1)
  truez <- z(x)
  truez[y(x) < -1] <- 1L
  truez[y(x) > -1 & y(x) <=-0.2] <- 2L
  truez[y(x) > -0.2] <- 3L
  truep <- as.numeric(table(truez)/length(truez))

  if(FALSE){
    xx <- posteriorSimulation(x)
    th <- theta(xx)
    th <- th[, order(th[1,])]
    dimnames(th) <- NULL
    checkEquals(th, true_theta, tolerance=0.5)
    checkEquals(truemu, sort(mu(xx)), tolerance=0.15)
    library(lattice)
    df <- data.frame(y=y(x), batch=batch(x))
    histogram(~y|batch, df, breaks=100, col="gray", border="gray", as.table=TRUE,
              panel=function(x, truth, ...){
                ##panel.abline(v=c(-1.5, -1, -0.5, 0, 0.5),
                ##col="gray")
                panel.histogram(x, ...)
                panel.abline(v=truth[panel.number(), ])
              }, truth=true_theta)


    zz <- as.integer(factor(z(x), levels=order(theta(x)[1,])))
    checkEquals(truez, zz, tolerance=0.1)


    checkEquals(truep, p(x), tolerance=0.3)

    ##
    ## start at truth
    ##
    z(x) <- truez
    zFreq(x) <- as.integer(table(truez))
    p(x) <- truep
    mu(x) <- truemu
    theta(x) <- true_theta
    sigma2(x)[, ] <- 0.06^2

    plot(x, use.current=TRUE)
    ##
    ## Note the zeros for some components
    ##
    table(batch(x), z(x))
    ##
    ## update mean with mu if no observations
    ##
    mns <- .Call("compute_means_batch", x)
    checkTrue(!any(is.nan(mns)))
    checkEquals(mns, true_theta, tolerance=0.1)
    dataMean(x) <- mns
    ##
    ## update variance with tau2 if no observations
    ##
    vars <- .Call("compute_vars_batch", x)
    checkTrue(!any(is.nan(vars)))
    prec <- .Call("compute_prec_batch", x)
    dataPrec(x) <- prec

    ##
    ## update theta
    ##
    th <- .Call("update_theta_batch", x)
    checkEquals(th, true_theta, tolerance=0.1)
    theta(x) <- th
    ##
    ## update mu
    ##
    mus <- .Call("update_mu_batch", x)
    checkEquals(mus, truemu, tolerance=0.1)
    mu(x) <- mus

    tau2s <- .Call("update_tau2_batch", x)
    tau2(x) <- tau2s

    set.seed(123) ;
    ps <- .Call("update_p_batch", x)
    checkEquals(ps, p(x), tolerance=0.05)

    iter(x) <- 500
    burnin(x) <- 0
    x2 <- posteriorSimulation(x)
    th <- theta(x2)
    dimnames(th) <- NULL
    checkEquals(th, true_theta, tolerance=0.1)
    checkEquals(mu(x2), truemu, tolerance=0.05)
    checkEquals(p(x2), truep, tolerance=0.05)
    if(FALSE){
      plot(x2, use.current=TRUE)
    }
  }
  ##
  ## Now test with random starting values
  ##
  ##load_all()
  xx <- BatchModel(data=testdat$cn, batch=testdat$batch, k=3)
  nStarts(xx) <- 10
  iter(xx, force=TRUE) <- 250
  burnin(xx) <- 200
  xx <- posteriorSimulation(xx)
  if(FALSE){
    plot(xx, use.current=TRUE)
    par(mfrow=c(4,3), las=1)
    tracePlot(xx, "theta")
  }
  th <- theta(xx)
  dimnames(th) <- NULL
  ix <- order(th[1, ])
  checkEquals(th, true_theta, tolerance=0.1)
  mus <- mu(xx)[ix]
  checkEquals(mus, truemu, tolerance=0.15)
  ps <- p(xx)[ix]
  checkEquals(ps, truep, tolerance=0.05)
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
  u[1:8] <- 0L
  paramUpdates(bm) <- u
  ## an entire second even if 0 updates
  system.time(tmp4 <- posteriorSimulation(bm)) ## 3.13

  ##
  ## Timings with 50 iterations, no burnin
  ##
  ## no updates
  u <- paramUpdates(mcmcp)
  u[1:8] <- 0L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ##  1.1
  system.time(tmp4 <- .Call("mcmc_batch", bm, mcmcParams(bm))) ##  1.09  0.746 when mean/prec not updated

  ## updating only z
  u <- paramUpdates(mcmcp)
  u[1:7] <- 0L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ##  1.75
  system.time(tmp4 <- posteriorSimulation(bm)) ##  1.75
  .Call("compute_means_batch", bm)
  system.time(tmp4 <- posteriorSimulation(bm)) ##  1.47
  ## updating only theta
  u[1:8] <- 0L
  u["theta"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 1.38
  u[1:8] <- 0L
  u["sigma2"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 1.38
  u[1:8] <- 0L
  u["p"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 1.1
  u[1:8] <- 0L
  u["mu"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ##1.35
  u[1:8] <- 0L
  u["tau2"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 1.1
  u[1:8] <- 0L
  u["nu.0"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 1.1
  u[1:8] <- 0L
  u["sigma2.0"] <- 1L
  paramUpdates(bm) <- u
  system.time(tmp4 <- posteriorSimulation(bm)) ## 1.1

}

.test_twocomponent_with_substantial_overlap <- function(){
}

.test_hard_fourcomp <- function(){
  set.seed(123)
  dir <- system.file("unitTests/data", package="CNPBayes")
  testdat <- readRDS(file.path(dir, "ea2.rds"))
  x <- BatchModel(data=testdat$cn, batch=testdat$batch)
  se <- as(x, "SummarizedExperiment")
  mcmcp <- McmcParams(iter=250, burnin=250, nStarts=5)
  m <- marginal(se, mcmc.params=mcmcp, maxperm=2)
  if(FALSE) plot(m[[1]], use.current=TRUE)

  b <- marginal(se, batch=testdat$batch, mcmc.params=mcmcp, maxperm=2)
  my <- rbind(summarizeMarginal(m),
              summarizeMarginal(b))
  bf <- bayesFactor(my)
  calls <- names(bf)
  ## check 4-component model is called
  checkTrue(substr(calls, 2, 2) == 4)
  if(FALSE){
    par(mfrow=c(1,3), las=1)
    plot(b[[4]], use.current=TRUE)
  }
}
