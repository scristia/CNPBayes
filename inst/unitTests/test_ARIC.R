.test_chr4locus <- function(){
  library(GenomicRanges)
  library(oligoClasses)
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  load("~/MyPapers/aricUricAcid/data/chr4locus.rda")
  hits <- findOverlaps(rowData(se.ea), chr4locus)
  se.ea <- se.ea[queryHits(hits), ]
  se.ea <- se.ea[rowData(se.ea)$source=="ARIC_VI", ]

  mcmcp <- McmcParams(iter=1000, burnin=1000)
  params <- ModelParams("marginal", y=copyNumber(se.ea)[1,], k=3)
  marg.model <- initializeModel(params)
  if(any(diff(theta(marg.model)) < 0.1)){
    cn <- copyNumber(se.ea)[1,]
    theta(marg.model) <- c(theta(marg.model)[1], theta(marg.model)[1]+1.5,
                           theta(marg.model)[1]+2)
  }
  ##mmfit <- mixtools:::normalmixEM(cn, arbvar = T, epsilon = 1e-06, k=3, maxit=2000)
  marg.model <- posteriorSimulation(marg.model, mcmcp)
  marg.model <- sort(marg.model)

  b.model <- BatchModel(y(marg.model), k(marg.model), batch=se.ea$plate)
  b <- collapseBatch(b.model)
  b.model <- BatchModel(y(marg.model), k(marg.model), batch=b)
  b.model <- startingValues(b.model)
  mu(b.model) <- theta(marg.model)
  dataMean(b.model) <- theta(b.model) <- matrix(theta(marg.model), nBatch(b.model), k(b.model), byrow=TRUE)
  sigma2(b.model) <- matrix(sigma2(marg.model), nBatch(b.model), k(b.model), byrow=TRUE)
  dataPrec(b.model) <- 1/sigma2(b.model)
  z(b.model) <- z(marg.model)
  sigma2.0(b.model) <- sigma2.0(marg.model)
  nu.0(b.model) <- nu.0(marg.model)
  tau2(b.model) <- sigma2(b.model)[1, ]/4
  b.model <- posteriorSimulation(b.model, mcmcp)

###plot(b.model)
  trace(.plotBatch, browser)
  .plotBatch(b.model)

  df <- data.frame(y=y(b.model), b=batch(b.model))
  library(lattice)
  histogram(~y|b, df)


  mcmcp <- McmcParams(iter=1000, burnin=1000)
  mmodel <- fitMixtureModels(se.ea[1, ], mcmcp, K=3)[[1]]
  plot(mmodel)

  mcmcp <- McmcParams(iter=2000, burnin=1000)
  mmodel2 <- posteriorSimulation(mmodel, mcmcp)

  plot(mmodel2)

  se.ea2 <- se.ea[, se.ea$plate!= "BOCKS"]
  se.ea2$plate <- batch(mmodel)
  mcmcp <- McmcParams(iter=1, burnin=0)
  mmodel2 <- fitMixtureModels(se.ea2[1, ], mcmcp, K=3)[[1]]



  mc <- mcmcChains(mmodel)
  plot.ts(theta(mc)[, 1:10])
  mcmcp <- McmcParams(iter=1000, burnin=0)
  mmodel2 <- posteriorSimulation(mmodel, mcmcp)
  ##
  ## Run the marginal model to find good starting values
  ##
  library(oligoClasses)
  mcmcp <- McmcParams(iter=1000, burnin=0)
  params <- ModelParams("marginal", y=copyNumber(se.ea)[1, ], k=3)
  mmodel2 <- initializeModel(params)
  mmodel2 <- posteriorSimulation(mmodel2, mcmcp)
  plot.ts(theta(mcmcChains(mmodel2)))
  if(FALSE) plot(mmodel2)

  mcmcp <- McmcParams(iter=1000, burnin=0)
  ##mmodels <- fitMixtureModels(se.ea[1, ], mcmcp, K=1:5)
  mmodel <- fitMixtureModels(se.ea[1, ], mcmcp, K=3)[[1]]
  plot(mmodel)



  mcmcp <- McmcParams(iter=1000, burnin=1000)
  mmodel <- posteriorSimulation(mmodel, mcmcp)
  ##bicstat <- sapply(mmodels, bic)
  model <- mmodels[[which.min(bicstat)]]
  saveRDS(model, file="cnp9.rds")
  model <- readRDS("~/Software/CNPBayes/tests/cnp9.rds")
  if(FALSE){
    CNPBayes::plot(model)
  }
  checkIdentical(which.min(bicstat), 3L)

  mcmcp <- McmcParams(iter=1000, burnin=500)
  mmodels <- fitMixtureModels(se.ea[2, ], mcmcp, K=1:5)
  bicstat <- sapply(mmodels, bic)
  model <- mmodels[[which.min(bicstat)]]
  saveRDS(model, file="cnp10.rds")
  if(FALSE){
    op <- par(mfrow=c(1,1),las=1)
    CNPBayes::plot(mmodel)
    par(op)
  }
  checkIdentical(which.min(bicstat), 3L)
}

.test_cnp360 <- function(){
  library(GenomicRanges)
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  se360 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp360.rds")
  j <- subjectHits(findOverlaps(rowData(se360), rowData(se.ea)))[1]
  se360 <- se.ea[j, ]

  phenodata <- readRDS("~/Labs/ARIC/aricPhenotypes/data/phenotypes.rds")
  phenodata <- phenodata[phenodata$racegrp == "W",]
  se360 <- se360[, colnames(se360) %in% phenodata$cel_id]


  r <- assays(se360)[["mean"]][1, ]
  freq <- table(se360$chemplate)
  if(FALSE) hist(r, breaks=1000, col="gray", border="gray")
  mcmcp <- McmcParams(iter=1000, burnin=500)

  mmodel <- fitMixtureModels()
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


.test_cnp472 <- function(){
  library(GenomicRanges)
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472.rds")
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp360.rds")
  j <- subjectHits(findOverlaps(rowData(se472), rowData(se.ea)))[1]
  se472 <- se.ea[j, ]

  ##
  ## Batches may differ between CNPs
  ##
  mcmcp <- McmcParams(iter=1000, burnin=500)
  mmodels <- fitMixtureModels(se472, mcmcp, K=1:3)

  op <- par(mfrow=c(1,1),las=1)
  CNPBayes::plot(mmodels[[2]])
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

.test_cnp707 <- function(){
  library(GenomicRanges)
  se707 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp707.rds")
  phenodata <- readRDS("~/Labs/ARIC/aricPhenotypes/data/phenotypes.rds")
  phenodata <- phenodata[phenodata$racegrp == "W",]
  se707 <- se707[, colnames(se707) %in% phenodata$cel_id]
  se707 <- SummarizedExperiment(assays=SimpleList(medr=assays(se707)[["mean"]]),
                                rowData=rowData(se707),
                                colData=colData(se707))
  colnames(colData(se707)) <- "plate"
  mcmcp <- McmcParams(iter=1000, burnin=500)
  mmodel <- fitMixtureModels(se707, mcmcp, K=3)[[1]]
  ##
  ## Batches may differ between CNPs
  ##
  freq <- as.integer(table(z(mmodel)))
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
