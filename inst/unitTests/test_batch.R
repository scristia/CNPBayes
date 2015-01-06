test_batchEasy <- function(){
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             .batch=rep(letters[1:3], length.out=2500),
                             .k=3,
                             .alpha=rep(1, k),
                             means=means,
                             sds=sds)

  ##
  ## Start from true values
  ##
  if(FALSE){
    model <- truth
    if(FALSE) CNPBayes::plot(truth, use.current=TRUE)
    mcmcp <- McmcParams(iter=1500, burnin=100)
    mcmcChains(model) <- McmcChains(model, mcmcp)
    model <- posteriorSimulation(model, mcmcp)
    if(FALSE){
      op <- par(mfrow=c(1,2),las=1)
      plot(truth, use.current=T)
      plot(model)
      par(op)
    }
    mc <- mcmcChains(model)
    pmns <- colMeans(theta(mc))
    checkEquals(pmns, as.numeric(theta(truth)), tolerance=0.03)

    ps <- colMeans(sigma(mc))
    checkEquals(ps, as.numeric(sigma(truth)), tolerance=0.18)

    pmix <- p(truth)
    pm_pmix <- colMeans(p(mc))
    checkEquals(pmix, pm_pmix, tolerance=0.02)
  }
  ##
  ## Start using defaults
  ##
  mcmcp <- McmcParams(iter=1000, burnin=100)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=rep(letters[1:3], length.out=2500),
                        mcmc.params=mcmcp)
  model2 <- initializeModel(params)
  ## floating point exception, core dumped
  ##trace(posteriorSimulation, browser)
  model2 <- posteriorSimulation(model2, mcmcp)
  mc <- mcmcChains(model2)
  pmns <- colMeans(theta(mc))
  checkEquals(pmns, as.numeric(theta(truth)), tolerance=0.02)

  ps <- colMeans(sigma(mc))
  checkEquals(ps, as.numeric(sigma(truth)), tolerance=0.1)

  pmix <- p(truth)
  pm_pmix <- colMeans(p(mc))
  checkEquals(pmix, pm_pmix, tolerance=0.02)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model2)
    par(op)
  }
}

test_selectK_batch_easy <- function(){
  ## Need to replace setParallelization with source code
  ## -investigate BiocParallel
  library(foreach)
  library(HopkinsHPC)
  NW <- setParallelization(7)
  if(FALSE){
    if(require(doSNOW)){
      library(doSNOW)
      cl <- makeCluster(7, type = "SOCK")
      registerDoSNOW(cl)
    }
  }
  set.seed(1000)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             .batch=rep(letters[1:3], length.out=2500),
                             .k=3,
                             .alpha=rep(1, k),
                             means=means,
                             sds=sds)

  ##
  ## Evaluate at different K
  ##
  mcmcp <- McmcParams(iter=1000, burnin=100, constrainTheta=TRUE)
  bicstat <- foreach(k = 1:7, .packages="CNPBayes", .combine="c") %dopar% {
    cat(".")
    params <- ModelParams("batch", y=y(truth), k=k,
                          batch=rep(letters[1:3], length.out=2500),
                          mcmc.params=mcmcp)
    modelk <- initializeModel(params)
    modelk <- posteriorSimulation(modelk, mcmcp)
    bic(modelk)
  }
  checkIdentical(which.min(bicstat), 3L)
}

test_selectK_batch_moderate <- function(){
  set.seed(100)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95,
                    -0.41, -0.4, -.395,
                    -0.1, 0, 0.05), nbatch, k, byrow=FALSE)
  sds <- matrix(0.15, nbatch, k)
  sds[,1] <- 0.3
  truth <- simulateBatchData(N=2500,
                             means=means,
                             sds=sds,
                             .k=3,
                             .batch=rep(letters[1:3], length.out=2500),
                             .alpha=c(100, 200, 400))
  if(FALSE) CNPBayes::plot(truth, use.current=TRUE)
  if(FALSE){
    ## unrealistic best case scenario: we start at the true values
    model <- truth
    mcmcp <- McmcParams(iter=2000, burnin=500)
    mcmcChains(model) <- McmcChains(model, mcmcp)
    model <- posteriorSimulation(model, mcmcp)
    mc <- mcmcChains(model)
    pmns <- colMeans(theta(mc))
    checkEquals(pmns, as.numeric(theta(truth)), tolerance=0.03)

    ps <- colMeans(sigma(mc))
    checkEquals(ps, as.numeric(sigma(truth)), tolerance=0.06)

    pmix <- colMeans(p(mc))
    checkEquals(pmix, p(truth), tolerance=0.08)
  }

  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
  ##
  ## Start from defaults.
  ##
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=rep(letters[1:3], length.out=2500),
                        mcmc.params=mcmcp)
  modelk <- initializeModel(params)
  modelk <- posteriorSimulation(modelk, mcmcp)
  ##
  ## This test demonstrating sensitivity to starting values.  -- chain
  ## has not yet converged.
  ##
  mc <- mcmcChains(modelk)
  pmns <- colMeans(theta(mc))
  checkEquals(pmns, as.numeric(theta(truth)), tolerance=0.03)
  if(FALSE)plot.ts(theta(mc), col="gray")
  ps <- colMeans(sigma(mc))
  checkEquals(ps, as.numeric(sigma(truth)), tolerance=0.1)
  pmix <- colMeans(p(mc))
  checkEquals(pmix, p(truth), tolerance=0.08)
}


test_selectK_batch_hard <- function(){
  if(require(doSNOW)){
    library(doSNOW)
    cl <- makeCluster(7, type = "SOCK")
    registerDoSNOW(cl)
  }
  set.seed(1234)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.9, -2, -1.85,
                    -0.45, -0.4, -0.35,
                    -0.1, 0, -0.05), nbatch, k, byrow=FALSE)
  sds <- matrix(0.3, nbatch, k)
  truth <- simulateBatchData(2500,
                             means=means,
                             sds=sds,
                             .batch=rep(letters[1:3], length.out=2500),
                             .alpha=c(5, 500, 1000))
  if(FALSE) plot(truth, use.current=TRUE)
  ##
  ## Evaluate at different K
  ##
  mcmcp <- McmcParams(iter=1000, burnin=0)
  bicstat <- foreach(k = 1:5, .packages="CNPBayes", .combine="c") %dopar% {
    cat(".")
    params <- ModelParams("batch", y=y(truth), k=3,
                          batch=rep(letters[1:3], length.out=2500),
                          mcmc.params=mcmcp)
    model <- initializeModel(params)
    model <- posteriorSimulation(model, mcmcp)

    model2 <- truth
    params <- ModelParams("batch", y=y(truth), k=3,
                          batch=rep(letters[1:3], length.out=2500),
                          mcmc.params=mcmcp)
    model2 <- posteriorSimulation(model2, mcmcp)

    mc <- mcmcChains(model)
    mc2 <- mcmcChains(model2)

    op <- par(mfrow=c(1,2), las=1)
    plot(truth, use.current=TRUE, xlim=c(-2,1))
    plot(model, xlim=c(-2,1))

    op <- par(mfrow=c(1,2), las=1)
    plot(truth, use.current=TRUE, xlim=c(-2,1))
    plot(model2, xlim=c(-2,1))

    plot.ts(theta(mc), col="gray")
    plot.ts(p(mc), col="gray")

    bic(model)

    params <- ModelParams("batch", y=y(truth), k=2,
                          batch=rep(letters[1:3], length.out=2500),
                          mcmc.params=mcmcp)
    model2 <- initializeModel(params)
    model2 <- posteriorSimulation(model2, mcmcp)
    mc <- mcmcChains(model2)
        op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=TRUE, xlim=c(-2,1))
    plot(model2, xlim=c(-2,1))

    ## see if marginal model works
    params <- ModelParams("marginal", y=y(truth), k=3)
    model3 <- initializeModel(params)

    model3 <- initializeModel(params)
    model3 <- posteriorSimulation(model3, mcmcp)
    mc <- mcmcChains(model2)

  }
  checkTrue(which.min(bicstat) == 3)

  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=TRUE, xlim=c(-2,1))
    plot(model, xlim=c(-2,1))
    par(op)
    ##
    ## variance of first component is bigger than the marginal
    ## variance of the data
    ##
    ##-- acts as an outlier component that captures extreme
    ##   observatios on either side of the modal component
    ##
    ## - constraining the variance to be less than or equal to the
    ##   marginal variance would prevent the above, forcing other
    ##   solutions
    ##
    ## Identifiability with batch.
    ##
    ##- is it reasonable to assume the variances should be about the
    ## same for a component across batches
    ##
    ##- why are the means for the first and 3rd components too small?
    ##-> 1st component has only a few observations -- lots of
    ##-##shrinkage. Not sure about last component.
    ##
    mc <- mcmcChains(model)
    ps <- colMeans(sigma(mc))
    checkEquals(ps, as.numeric(sigma(truth)), tolerance=0.1)

    pmns <- colMeans(theta(mc))
    checkEquals(pmns, as.numeric(theta(truth)), tolerance=0.15)
    pmix <- colMeans(p(mc))
    ## FALSE
    checkEquals(pmix, p(truth), tolerance=0.1)
    plot.ts(sigma2(mc), col="gray")
  }
}

test_KolmogorovSmirnov <- function(){
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.2, -1.0,
                    0, 0, 0.2,
                    0.8, 0.8, 1.0), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             .batch=rep(letters[1:3], length.out=2500),
                             .k=3,
                             .alpha=rep(1, k),
                             means=means,
                             sds=sds)
  ## - collapse using KS test for two distributions
  ks <- ksTest(truth)
  checkIdentical(ks, c(TRUE, FALSE, FALSE))
  truth2 <- collapseBatch(truth)
  checkIdentical(uniqueBatch(truth2), c("a-b", "c"))
}



test_cnp360 <- function(){
  library(GenomicRanges)
  se360 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp360.rds")
  phenodata <- readRDS("~/Labs/ARIC/aricPhenotypes/data/phenotypes.rds")
  phenodata <- phenodata[phenodata$racegrp == "W",]

  se360 <- se360[, colnames(se360) %in% phenodata$cel_id]
  r <- assays(se360)[["mean"]][1, ]
  freq <- table(se360$chemplate)
  if(FALSE) hist(r, breaks=1000, col="gray", border="gray")
  mcmcp <- McmcParams(iter=1000, burnin=500)
  params <- ModelParams("batch", y=r, k=3,
                        batch=se360$chemplate,
                        mcmc.params=mcmcp)
  modelk <- initializeModel(params)
  modelk2 <- collapseBatch(modelk, mcmcp)
  modelk3 <- posteriorSimulation(modelk2, mcmcp)


  if(FALSE){
    saveRDS(modelk3, file="~/Software/CNPBayes/inst/unitTests/model_360.rds")
    saveRDS(batch(modelk3), file="~/Software/CNPBayes/inst/extdata/batch_ks.rds")

    op <- par(mfrow=c(1,1),las=1)
    CNPBayes::plot(modelk3)
    par(op)

    ##
    ##
    ##
    index <- which(z(modelk3)==3)
    y(modelk3)[index]## large values -- could be real, or outliers. Would need BAFs to verify
    bic(modelk3)

    params <- ModelParams("batch", y=r, k=4,
                          batch=batch(modelk3),
                          mcmc.params=mcmcp)
    modelk4 <- initializeModel(params)
    modelk4 <- posteriorSimulation(modelk4, mcmcp)
    table(z(modelk4))
    bic(modelk4)

    ##
    ## call homozygous deletions post-hoc?
    ##
    params <- ModelParams("batch", y=r, k=2,
                          batch=batch(modelk3),
                          mcmc.params=mcmcp)
    modelk2 <- initializeModel(params)
    modelk2 <- posteriorSimulation(modelk2, mcmcp)
    bic(modelk2)

    ## the k = 3 model is selected by bic
    sum(r(modelk3) < -2)

    ##
    mc <- mcmcChains(model)
    ps <- colMeans(sigma(mc))
    checkEquals(ps, as.numeric(sigma(truth)), tolerance=0.1)

    pmns <- colMeans(theta(mc))
    checkEquals(pmns, as.numeric(theta(truth)), tolerance=0.15)
    pmix <- colMeans(p(mc))
    ## FALSE
    checkEquals(pmix, p(truth), tolerance=0.1)
    plot.ts(sigma2(mc), col="gray")
  }

  ##modelk <- initializeModel(params)
  ##  params <- ModelParams("marginal", y=r, k=3, batch=rep("A", length(r)),
  ##                        mcmc.params=mcmcp)
  ##  model <- initializeModel(params)
  ##  model <- posteriorSimulation(model, mcmcParams(params))
  ##  mc <- mcmcChains(model)
  ##  plot.ts(theta(mc), col="gray")
  ##  plot(model)
  ##  bic(model)
}


test_cnp472 <- function(){
  library(GenomicRanges)
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472.rds")
  phenodata <- readRDS("~/Labs/ARIC/aricPhenotypes/data/phenotypes.rds")
  phenodata <- phenodata[phenodata$racegrp == "W",]
  se472 <- se472[, colnames(se472) %in% phenodata$cel_id]
  ##
  ## Batches may differ between CNPs
  ##
  mcmcp <- McmcParams(iter=1000, burnin=500)
  r <- assays(se472)[["mean"]][1, ]
  params <- ModelParams("batch", y=r, k=3,
                        batch=se472$chemplate,
                        mcmc.params=mcmcp)
  m473k3 <- initializeModel(params)
  m473k3 <- collapseBatch(m473k3, mcmcp)
  m473k3 <- posteriorSimulation(m473k3, mcmcp)
  ## Call homozygous deletions post-hoc
  bic(m473k3)

  params <- ModelParams("batch", y=r, k=2,
                        batch=batch(m473k3),
                        mcmc.params=mcmcp)
  m472k2 <- initializeModel(params)
  m472k2 <- posteriorSimulation(m472k2, mcmcp)
  bic(m472k2)

  op <- par(mfrow=c(1,1),las=1)
  CNPBayes::plot(m472k2)
  par(op)
  ##plot(y(modelk3)~z(modelk3))
  if(FALSE){
    se_snps <- readRDS("~/Software/CNPBayes/inst/extdata/se_snps472.rds")
    se_snps <- se_snps[, colnames(se472)]
    df <- data.frame(A=as.numeric(A(se_snps)),
                     B=as.numeric(B(se_snps)),
                     z=z(modelk3))
    xyplot(B~A, df, pch=20, col=df$z, cex=0.4)
  }
}

test_cnp707 <- function(){
  library(GenomicRanges)
  se707 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp707.rds")
  phenodata <- readRDS("~/Labs/ARIC/aricPhenotypes/data/phenotypes.rds")
  phenodata <- phenodata[phenodata$racegrp == "W",]
  se707 <- se707[, colnames(se707) %in% phenodata$cel_id]
  ##
  ## Batches may differ between CNPs
  ##
  mcmcp <- McmcParams(iter=1000, burnin=500)
  r <- assays(se707)[["mean"]][1, ]
  params <- ModelParams("batch", y=r, k=3,
                        batch=se707$chemplate,
                        mcmc.params=mcmcp)
  model <- initializeModel(params)
  model <- collapseBatch(model, mcmcp)
  modelk3 <- posteriorSimulation(model, mcmcp)
  freq <- as.integer(table(z(modelk3)))
  checkIdentical(c(9718L, 0L, 0L), freq)
  ##
  ## with k = 3, we only have one component with more than 1 observation -- this is effectively k=1
  ##
  ##  params <- ModelParams("batch", y=r, k=2,
  ##                        batch=batch(modelk3),
  ##                        mcmc.params=mcmcp)
  ##  modelk2 <- initializeModel(params)
  ##  modelk2 <- posteriorSimulation(modelk2, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,1),las=1)
    CNPBayes::plot(modelk3)
    par(op)

    se_snps <- readRDS("~/Software/CNPBayes/inst/extdata/se_snps707.rds")
    library(lattice)
    se_snps <- se_snps[, colnames(se707)]
    a <- A(se_snps)
    b <- B(se_snps)
    ##cn <- assays(se_map)[["mean"]]["CNP_707", names(r)]
    cn <- matrix(z(modelk3), nrow(se_snps), ncol(se_snps), byrow=TRUE)
    snpid <- matrix(rownames(se_snps), nrow(se_snps), ncol(se_snps))
    alleles <- data.frame(A=log2(as.numeric(a)),
                          B=log2(as.numeric(b)),
                          cn=as.numeric(cn),
                          snp=as.character(snp))
    alleles$snp <- factor(alleles$snp)
    p <- xyplot(B~A|snp, alleles, pch=20, cex=0.3, layout=c(1, 3),
                panel=function(x, y, cn,  ..., subscripts){
                  cn <- cn[subscripts]
                  colors <- c("black", "gray")[cn]
                  lpoints(x[cn==2], y[cn==2], col="gray", ...)
                  lpoints(x[cn==1], y[cn==1], col="black", ...)
                }, cn=alleles$cn, auto.key=list(points=TRUE, columns=2), xlab=expression(log[2](A)),
                ylab=expression(log[2](B)))
    p
  }

}
