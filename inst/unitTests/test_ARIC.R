.test_chr4locus <- function(){
  ## defining a 'truth' for chrom 4 CNP for unit testing
  library(GenomicRanges)
  library(oligoClasses)
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  load("~/MyPapers/aricUricAcid/data/chr4locus.rda")
  hits <- findOverlaps(rowData(se.ea), chr4locus)
  se.ea <- se.ea[queryHits(hits), ]
  se.ea <- se.ea[rowData(se.ea)$source=="ARIC_VI", ]

  mcmcp <- McmcParams(iter=1000, burnin=1000)
  params <- ModelParams("marginal", y=copyNumber(se.ea)[1,], k=3)
  marg.model <- CNPBayes:::initializeModel(params)
  if(any(diff(theta(marg.model)) < 0.1)){
    cn <- copyNumber(se.ea)[1,]
    theta(marg.model) <- c(theta(marg.model)[1], theta(marg.model)[1]+1.5,
                           theta(marg.model)[1]+2)
  }
  ##mmfit <- mixtools:::normalmixEM(cn, arbvar = T, epsilon = 1e-06, k=3, maxit=2000)
  marg.model <- posteriorSimulation(marg.model, mcmcp)
  ##  saveRDS(marg.model, file="~/Software/CNPBayes/inst/unitTests/model_unequal_mix.rds")
  ##  plot(marg.model)
  b <- collapseBatch(se.ea[1,])
  params <- ModelParams("batch", y=copyNumber(se.ea)[1, ], batch=b, k=3,
                        mcmc.params=mcmcp)
  truth <- CNPBayes:::initializeModel(params)
  z(truth) <- z(marg.model)
  p(truth) <- sapply(split(z(truth), z(truth)), length)/length(z(truth))
  ylist <- split(y(truth), z(truth))
  blist <- split(truth@batch, z(truth))
  thetas <- foreach(y=ylist, b=blist, .combine="rbind") %do%{
    sapply(split(y, b), mean)
  }
  thetas <- (t(thetas))[CNPBayes:::uniqueBatch(truth), ]
  theta(truth) <- thetas
  sigmas <- foreach(y=ylist, b=blist, .combine="rbind") %do%{
    sapply(split(y, b), sd)
  }
  sigmas <- (t(sigmas))[CNPBayes:::uniqueBatch(truth), ]
  sigma2(truth) <- sigmas^2
  dataMean(truth) <- thetas
  dataPrec(truth) <- 1/sigma2(truth)
  saveRDS(truth, file="~/Software/CNPBayes/inst/extdata/unequal_mix_model.rds")
}

.test_cnp360 <- function(){
  library(GenomicRanges)
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  se360 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp360.rds")
  j <- subjectHits(findOverlaps(rowData(se360), rowData(se.ea)))[1]
  se360 <- se.ea[j, ]


  b <- collapseBatch(se360[1,])
  mcmcp <- McmcParams(iter=1000, burnin=1000)
  params <- ModelParams("batch", y=copyNumber(se360)[1, ], batch=b, k=3,
                        mcmc.params=mcmcp)
  bmodel <- CNPBayes:::initializeModel(params)
  bmodel <- posteriorSimulation(bmodel, mcmcp)
  zz <- map(bmodel)
  is_homozygous <- zz==1
  checkTrue(sum(is_homozygous) >= 2 & sum(is_homozygous) <=4)

  if(FALSE){
    params <- ModelParams("batch", y=copyNumber(se360)[1, ], batch=b, k=2,
                          mcmc.params=mcmcp)
    bmodel2 <- CNPBayes:::initializeModel(params)
    bmodel2 <- posteriorSimulation(bmodel2, mcmcp)
  }

  if(FALSE){
    op <- par(mfrow=c(1,1),las=1)
    CNPBayes::plot(bmodel)
    par(op)
  }
}


.test_cnp472 <- function(){
  library(GenomicRanges)
  library(oligoClasses)
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472.rds")
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  j <- subjectHits(findOverlaps(rowData(se472), rowData(se.ea)))[1]
  se472 <- se.ea[j, ]
  ##
  ##
  ##
  mcmcp <- McmcParams(iter=10000, burnin=200, thin=10, nStarts=10)
  mparams <- ModelParams("marginal", y=copyNumber(se472)[1, ], k=3,
                        mcmc.params=mcmcp)
  mmodel <- CNPBayes:::initializeModel(mparams)
  mmodel <- posteriorSimulation(mmodel, mcmcp)

  truth <- CNPBayes:::initializeModel(mparams)
  z(truth) <- z(mmodel)
  theta(truth) <- theta(mmodel)
  p(truth) <- sapply(split(z(truth), z(truth)), length)/length(z(truth))
  sigma2(truth) <- sigma2(mmodel)
  mu(truth) <- mu(mmodel)
  tau2(truth) <- tau2(mmodel)
  mcmcp <- McmcParams(iter=0, burnin=500)
  truth2 <- posteriorSimulation(truth, mcmcp)
  saveRDS(truth2, file="~/Software/CNPBayes/inst/extdata/cnp472_model.rds")

  sei <- readRDS("~/Labs/ARIC/AricCNPData/data/intensity_472.rds")
  table(colnames(se472) %in% colnames(sei))
  ##zz <- zz[colnames(se472) %in% colnames(sei)]
  ##zz <- colnames(se472)[colnames(se472) %in% colnames(sei)]
  se472 <- se472[, colnames(se472) %in% colnames(sei)]
  sei <- sei[, colnames(se472)]
  df <- data.frame(a=log2(assays(sei)[["A"]])[1,],
                   b=log2(assays(sei)[["B"]])[1,],
                   z=zz)
  library(lattice)
  xyplot(b~a, df, pch=20, cex=0.3)

  xyplot(b~a, df[zz==2, ], pch=20, cex=0.3)
  xyplot(b~a, df, pch=20, cex=0.3, col=zz)

  ##mcmcp <- McmcParams(iter=1000, burnin=500)
  ##mmodels <- fitMixtureModels(se472, mcmcp, K=1:3)
  ##plot(y(modelk3)~z(modelk3))
  if(FALSE){
    op <- par(mfrow=c(1,1),las=1)
    CNPBayes::plot(mmodels[[2]])
    par(op)
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
  library(oligoClasses)
  se707 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp707.rds")
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  j <- subjectHits(findOverlaps(rowData(se707), rowData(se.ea)))[1]
  se707 <- se.ea[j, ]
  b <- collapseBatch(se707[1,])


  mcmcp <- McmcParams(iter=c(500, 1000, 1000), burnin=c(100, 200, 200))
  params <- ModelParams("marginal", y=copyNumber(se707)[1, ], k=3,
                        mcmc.params=mcmcp)
  mmodels <- CNPBayes:::fitMixtureModels(y(params), mcmcp, K=1:3)
  sapply(mmodels, bic)


  mcmcp <- McmcParams(iter=1000, burnin=100)
  bparam <- ModelParams("batch", y=copyNumber(se707)[1, ], k=3, batch=b,
                        mcmc.params=mcmcp)
  bmod <- initializeBatchModel(bparam)
  bmodel <- posteriorSimulation(bmod, mcmcp)

  bparam <- ModelParams("batch", y=copyNumber(se707)[1, ], k=2, batch=b,
                        mcmc.params=mcmcp)
  bmod2 <- initializeBatchModel(bparam)
  bmodel2 <- posteriorSimulation(bmod2, mcmcp)

  bparam <- ModelParams("batch", y=copyNumber(se707)[1, ], k=1, batch=batch(bmod2),
                        mcmc.params=mcmcp)
  bmod1 <- initializeBatchModel(bparam)
  bmodel1 <- posteriorSimulation(bmod1, mcmcp)

  baffile <- list.files(datdir, pattern="BAllele", full.names=TRUE)
  lrrfile <- list.files(datdir, pattern="LogR", full.names=TRUE)
  gtfile <- list.files(datdir, pattern="Genotype", full.names=TRUE)

  bafdat <- fread(baffile, nrows=10)

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
    library(oligoClasses)
    library(crlmm)
    op <- par(mfrow=c(1,1),las=1)
    CNPBayes::plot(modelk3)
    par(op)

    se_snps <- readRDS("~/Software/CNPBayes/inst/extdata/se_snps707.rds")
    library(lattice)
    table(colnames(se707) %in% colnames(se_snps))
    se707_2 <- se707[, colnames(se707) %in% colnames(se_snps)]
    se_snps <- se_snps[, colnames(se707_2)]
    a <- log2(assays(se_snps)[["A"]])
    b <- log2(assays(se_snps)[["B"]])

    cn <- map(bmodel2)
    names(cn) <- colnames(se707)
    cn <- cn[colnames(se_snps)]
    cn <- rep(cn, each=nrow(se_snps))

    cn.marginal <- z(mmodels[[3]])
    names(cn.marginal) <- colnames(se707)
    cn.marginal <- cn.marginal[colnames(se_snps)]
    cn.marginal <- rep(cn.marginal, each=nrow(se_snps))
    ##cn <- assays(se_map)[["mean"]]["CNP_707", names(r)]
    ##cn <- matrix(z(bmodel2), nrow(se_snps), ncol(se_snps), byrow=TRUE)
    snpid <- rep(rownames(se_snps), ncol(se_snps))
    ##snpid <- matrix(rownames(se_snps), nrow(se_snps), ncol(se_snps))
    alleles <- data.frame(A=as.numeric(a),
                          B=as.numeric(b),
                          cn=as.numeric(cn),
                          cn.marginal=as.numeric(cn.marginal),
                          snp=snpid)
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
    p2 <- xyplot(B~A|snp, alleles, pch=20, cex=0.3, layout=c(1, 3),
                 panel=function(x, y, cn,  ..., subscripts){
                   cn <- cn[subscripts]
                   colors <- c("black", "gray")[cn]
                   lpoints(x[cn==2], y[cn==2], col="gray", ...)
                   lpoints(x[cn==1], y[cn==1], col="black", ...)
                 }, cn=alleles$cn.marginal, auto.key=list(points=TRUE, columns=2), xlab=expression(log[2](A)),
                 ylab=expression(log[2](B)))
  }

}
