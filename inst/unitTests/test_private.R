.test_profiling <- function(){
  ##
  ## example (from test_marginal.R)
  ##
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth, use.current=TRUE)
  params <- ModelParams("marginal", y=y(truth), k=3)
  mcmcp <- McmcParams(iter=150, burnin=0)
  ##trace(initializeModel, browser)
  model <- CNPBayes:::initializeModel(params)
  
  Rprof("m2.prof")
  model <- posteriorSimulation(model, mcmcp)
  Rprof(NULL)
  m1 <- summaryRprof("m1.prof")
  m2 <- summaryRprof("m2.prof")
  
  ##
  ## example (from test_batch.R)
  ##
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
  mcmcp <- McmcParams(iter=100, burnin=0)
  ## use the true k
  k <- 3
  ## marginal
  params <- ModelParams("marginal", y=y(truth), k=k,
                        batch=rep(letters[1:3], length.out=2500),
                        mcmc.params=mcmcp)
  model <- CNPBayes:::initializeModel(params)
  #model <- posteriorSimulation(model, mcmcp)
  updateTheta(model, constrain=FALSE)
  updateThetaCpp(model, constrain=TRUE)
  updateSigma2Cpp(model)
  ## Need to profile posteriorSimulation.
  Rprof()
  posteriorSimulation(model, mcmcp)
  
  ##
  ## Batch model
  ##
  outdir <- '~/Software'
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472_EA.rds")
  B <- getFiles(outdir, rownames(se472), "batch")
  
  batch.files <- paste0(dirname(model(B)), "/", rownames(se472), "_batch.rds")
  object <- BatchModel(copyNumber(se472)[1,], k=3, batch= batch.files)
  
  load_all()
  Rprof()
  batchExperiment(se472, outdir, test=TRUE)
  Rprof(NULL)
  summaryRprof()
  prof <- summaryRprof()
  tot <- prof$by.total
  self <- prof$by.self
  
  load_all()
  batchExperiment(se472, outdir, test=TRUE)
  outdir <- tempdir()
  Rprof("marginal.prof", interval=0.1)
  marginalExperiment(se472, outdir)
  Rprof(NULL)
  mp2 <- summaryRprof("marginal.prof")
  
  load_all()
  marginalExperiment(se472, outdir, test=TRUE)
  
  
  
  #mcmcp <- McmcParams(iter=c(500, 2000, 2000, 2000),
  #                      burnin=c(50, 100, 100, 100),
  #                      thin=c(1, 2, 2, 2),
  #                      nStarts=20,
  #                      nStartIter=200)
  #  hypp <- HyperparametersBatch()
  #  hplist <- HyperParameterList(hypp, K=1:4)
  #  mplist <- ModelParamList(hypp, K=1:4, data=copyNumber(se472)[1, ],
  #                           batch=se472$batch, mcmcp=mcmcp)
  #  bmodlist <- foreach(hypp=hplist, param=mplist) %do% {
  #    initializeBatchModel(params=param, hypp=hypp)
  #  }
  #  bmodels <- foreach(k=1:4, model=bmodlist) %do% posteriorSimulation(model, mcmcp[k])
  
  ##### Rcpp S4 testing stuff
  library(inline)
  src <- '
  S4 foo(x) ; foo.slot(".Data") = "bar" ; return(foo);
  '
  fun <- cxxfunction(signature(x="any"), src,
                     plugin="Rcpp")
  setClass( "S4ex", contains = "character",
  representation( x = "numeric" ) )
  x <- new( "S4ex", "bla", x = 10 )
  fun(x)
  str(fun(x))
}

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

.test_BatchModel <- function(){
  library(GenomicRanges)
  library(oligoClasses)
  outdir <- tempdir()
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472_EA.rds")
  B <- getFiles(outdir, rownames(se472), "batch")

  b <- factor(se472$batch, levels=unique(se472$batch))
  ix <- order(b)
  yy <- copyNumber(se472)[1, ]
  yy <- as.numeric(yy[ix])
  b <- b[ix]


  load_all()
  set.seed(134)
  mp <- ModelParams("batch",
                    y=copyNumber(se472)[1,],
                    k=3,
                    batch=se472$batch,
                    mcmc.params=McmcParams(iter=2))
  hp <- CNPBayes:::HyperparametersBatch(tau2.0=1000, k=3)
  mod <- initializeBatchModel(mp, hypp=hp)
  ##checkEquals(as.numeric(y(mod)), yy)
  ##checkEquals(batch(mod), b)
  checkEquals(round(mu(mod), 2), c(-0.28, 0.23, 0.28))
  checkEquals(round(tau2(mod), 2), c(4.16, 0.04, 71.26))
  checkIdentical(CNPBayes:::uniqueBatch(mod)[1:2], c("LAITH,WI", "AWFUL,BO"))
  checkEquals(logLikData(mod), -13799.59, tolerance=0.01)
  checkEquals(round(as.numeric(dataMean(mod)[1, ]), 2), c(-0.07, -0.06, -0.06))
  ##checkEquals(round(as.numeric(dataPrec(mod)[1, ]), 2), c(38.45, 40.10, 52.65))
  set.seed(134)
  modfit <- posteriorSimulation(mod, McmcParams(iter=10))
  th <- round(theta(modfit)[1:2, ], 2)
  dimnames(th) <- NULL
  mns <- matrix(c(-0.09, 0.05, -0.04, 0.03, -0.05, 0.07), 2, 3)
  checkEquals(mns, th)
  th <- round(sigma(modfit)[1:2, ], 2)
  dimnames(th) <- NULL
  mns <- matrix(c(0.17, 0.13, 0.16, 0.20, 0.14, 0.13), 2, 3)
  checkEquals(mns, th)
  th <- as.numeric(round(mu(modfit), 3))
  mns <- c(-0.003, -0.004, 0.000)
  checkEquals(mns, th)
  th <- as.numeric(round(tau2(modfit), 3))
  mns <- c(0.024, 0.029, 0.029)
  checkEquals(mns, th)
  checkEquals(sigma2.0(modfit), 0.1668, tolerance=0.01)
  checkEquals(nu.0(modfit), 1L)
  checkEquals(round(p(modfit), 2), c(0.23, 0.35, 0.41))
  if(FALSE){
    probs1 <- posteriorMultinomial(modfit)
    zz1 <- .updateZ(probs1)
    checkEquals(as.integer(table(zz1)), c(2244L, 3256L, 3917L))
    tabz <- CNPBayes:::tablez(modfit)
    checkEquals(as.integer(tabz[1:2, ]), c(459L, 191L, 714L, 223L, 828L, 317L))
  }
  if(FALSE){
    probs2 <- posteriorMultinomial(modfit)
    checkEquals(probs1[ix, ], probs2)
  }
  ## zz2 would not necessarily be the same as zz because a uniform 0,
  ## 1 RV (U) is simulated in .updateZ.  The relationship between U
  ## and y is not the same!
  ## This means that dataMean and dataPrec will not be exactly the
  ## same either, unless we force z to be the same as we do here.
  checkEquals(as.integer(table(zz2)), c(2244L, 3256L, 3917L), tolerance=25)
  lld <- logLikData(modfit)
  ##checkIdentical(round(lld,3), 4466.369)
  checkIdentical(round(lld,3), 5610.011)



  sds <- round(sigma(modfit)[1:2, ], 2)
  dimnames(sds) <- NULL
  sds2 <- matrix(c(0.17, 0.09, 0.18, 0.26, 0.09, 0.09), 2, 3)
  checkEquals(sds, sds2)

  zz <- CNPBayes:::tablez(modfit)[1:2, ]
  mat <- matrix(c(368L, 205L, 611L, 187L, 1022L, 339L), 2, 3)
  dimnames(zz) <- NULL
  checkEquals(zz, mat)
  llp <- logLikPhi(modfit)
  checkEquals(llp, -185.66, tolerance=0.01)

if(FALSE){
  modfit2=modfit
  modfit2@data <- modfit@data2
  modfit2@batch <- as.character(modfit2@batch2)
  modfit2@z <- modfit@z[modfit@ix]
  mns <- computeMeans(modfit)[1:2, ]
  mns2 <- computeMeans(modfit2)[1:2, ]
  checkEquals(mns, mns2)

  prec <- 1/computeVars(modfit)
  prec2 <- 1/computeVars(modfit2)
  checkEquals(prec, prec2)

  set.seed(1)
  th1 <- updateTheta(modfit)
  set.seed(1)
  th2 <- updateTheta(modfit2)
  checkEquals(th1, th2)

  set.seed(1)
  s1 <- updateSigma2(modfit)
  set.seed(1)
  s2 <- updateSigma2(modfit2)
  checkEquals(s1, s2)

  ## there is nothing random in posteriorMultinomial. Should get the
  ## same result after reordering the probabilities
  P1 <- posteriorMultinomial(modfit)
  P2 <- posteriorMultinomial(modfit2)
  checkEquals(P1[ix,1], P2[,1])

  set.seed(1)
  zz1 <- .updateZ(P1[ix, ])
  set.seed(1)
  zz2 <- .updateZ(P2)
  checkEquals(zz1, zz2)

  set.seed(1)
  mp1 <- updateMixProbs(modfit)
  set.seed(1)
  mp2 <- updateMixProbs(modfit2)
  checkEquals(mp1, mp2)

  zz <- CNPBayes:::tablez(modfit)[1:2, ]
  zz <- CNPBayes:::tablez(modfit2)[1:2, ]
}

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
  my <- rbind(CNPBayes:::summarizeMarginal(m),
              CNPBayes:::summarizeMarginal(b))
  bf <- bayesFactor(my)
  calls <- names(bf)
  ## check 4-component model is called
  checkTrue(substr(calls, 2, 2) == 4)
  if(FALSE){
    par(mfrow=c(1,3), las=1)
    plot(b[[4]], use.current=TRUE)
  }
}
