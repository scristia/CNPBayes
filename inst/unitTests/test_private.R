.test_profiling <- function(){
  ##
  ## example (from test_marginal.R)
  ##
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth)
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
  if(FALSE) plot(m[[1]])

  b <- marginal(se, batch=testdat$batch, mcmc.params=mcmcp, maxperm=2)
  my <- rbind(CNPBayes:::summarizeMarginal(m),
              CNPBayes:::summarizeMarginal(b))
  bf <- bayesFactor(my)
  calls <- names(bf)
  ## check 4-component model is called
  checkTrue(substr(calls, 2, 2) == 4)
  if(FALSE){
    par(mfrow=c(1,3), las=1)
    plot(b[[4]])
  }
}

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
    plot(truth)
    plot(model)
    par(op)
    op <- par(mfrow=c(1,3), las=1)
    tracePlot(model, "theta", col=1:3)
    par(op)
  }
  th <- CNPBayes:::thetaMean(model)
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
    plot(truth)
    plot(model)
    par(op)
    op <- par(mfrow=c(1,3), las=1)
    tracePlot(model, "theta", col=1:3)
    par(op)
  }
  th <- CNPBayes:::thetaMean(model)
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
    plot(truth)
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
  mcmcp <- McmcParams(iter=1000, burnin=100, nStarts=10)
  modelk <- posteriorSimulation(model, mcmcp)
  checkEquals(theta(modelk), theta(truth), tolerance=0.1)

  pmns <- CNPBayes:::thetaMean(modelk)
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
    CNPBayes::plot(truth)
    CNPBayes::plot(modelk)
    par(op)
    mc <- chains(modelk)
    plot.ts(sigma(mc), col="gray")
    plot.ts(theta(mc), col="gray")
    plot.ts(p(mc), col="gray")
  }
}

.test_hard1 <- function(){
  truth <- hardTruth(0.005, s=0.3)
  table(z(truth), batch(truth))
  if(FALSE) plot(truth)
  ##
  ## Use defaults
  ##
  mcmcp <- McmcParams(iter=2000, burnin=1000,
                      nStarts=20, thin=2)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  model <- initializeBatchModel(params)
  model <- posteriorSimulation(model, mcmcp)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, xlim=c(-2,1))
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
  if(FALSE) CNPBayes::plot(truth)

  ## start at true values
  mcmcp <- McmcParams(iter=1000, burnin=1000, thin=1)
  model <- truth
  model <- posteriorSimulation(model, mcmcp)
  mc <- chains(model)
  checkEquals(colMeans(theta(mc)), as.numeric(theta(truth)), tolerance=0.06)
  checkEquals(colMeans(sigma(mc)), as.numeric(sigma(truth)), tolerance=0.05)
  checkEquals(colMeans(p(mc)), as.numeric(p(truth)), tolerance=0.17)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, xlim=c(-2,1))
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
  mc <- chains(model)
  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    CNPBayes::plot(truth, xlim=c(-2,1))
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
                      nStarts=20)
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
  mcmcp <- McmcParams(iter=1000, burnin=100, nStarts=20)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=batch(truth),
                        mcmc.params=mcmcp)
  hypp <- CNPBayes:::HyperparametersBatch(k=3, tau2.0=1000)
  bmodel <- initializeBatchModel(params, hypp=hypp)
  bmodel <- posteriorSimulation(bmodel, mcmcp)

  pmns <- CNPBayes:::thetaMean(bmodel)
  j <- order(pmns[1, ])
  ps <- sigmaMean(bmodel)[, j]
  pmix <- pMean(bmodel)[j]

  checkEquals(pmns[, j], theta(truth), tolerance=0.02)
  checkEquals(pmix, as.numeric(p(truth)), tolerance=0.01)

  checkTrue(CNPBayes:::logpotentialc(bmodel)[CNPBayes:::argMax(bmodel)] > -2000)

  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth)
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
                      nStart=20)
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
    plot(truth)
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
                      nStart=20)
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

.test_selectK_nobatchEffect <- function(){
  ## Fit batch model when truth is no batch effect
  ## batch found
  library(GenomicRanges)
  set.seed(1000)
  means <- c(-1, 0, 1)
  sds <- c(0.1, 0.2, 0.2)
  truth <- simulateData(N=2500, p=rep(1/3, 3), theta=means, sds=sds)
  se <- as(truth, "SummarizedExperiment")
  library(devtools)
  load_all()
  ##m <- marginal(se, mcmc.params=McmcParams(iter=10, burnin=5, nStarts=1), K=3)

  kmod <- MarginalModel(y(truth), k=3)
  iter(kmod, force=TRUE) <- 500
  kmod <- posteriorSimulation(kmod)
  headtheta <- head(thetac(kmod))

  ##
  ## Marginal likelihood
  ##
  logp_theta <- .Call("marginal_theta", kmod)
  logp <- log(mean(exp(logp_theta)))
  checkTrue(identical(headtheta, head(thetac(kmod))))

  ##  logp_theta <- .Call("p_theta_zpermuted", kmod)
  ##  logp <- log(mean(exp(logp_theta)))
  ##  checkTrue(identical(headtheta, head(thetac(kmod))))

  ## Simulate Z under the reduced model
  kmodZ <- .Call("simulate_z_reduced1", kmod)
  checkTrue(identical(headtheta, head(thetac(kmod))))
  zz <- zChain(kmodZ)

  logp_sigma2 <- .Call("p_sigma2_zpermuted", kmodZ)
  lps2 <- log(mean(exp(logp_sigma2)))
  checkTrue(identical(headtheta, head(thetac(kmod))))

  logp_p <- .Call("p_p_unpermuted", kmod)
  lpp <- log(mean(exp(logp_p)))
  checkTrue(identical(headtheta, head(thetac(kmod))))

  i <- CNPBayes:::argMax(kmod)
  mlik <- log_lik(chains(kmod))[i] + logPrior(chains(kmod))[i] - logp - lps2 - lpp

  permutations <- permnK(3, 5)
  kmod2 <- kmod
  zChain(kmod2) <- permuteZ(kmod, permutations[2, ])
  zz <- table(zChain(kmod)[1, ])
  zz2 <- table(zChain(kmod2)[1, ])
  ##
  ## this uses the permuted values of z directly.
  ##
  logp_theta2 <- .Call("p_theta_zpermuted", kmod2)  ## use permuted z
  logp2 <- log(mean(exp(logp_theta2)))
  ## now we have to simulate new values of z as described in Chib (We
  ## need z|y, theta* and not z|y).  To do this, we could use the z
  ## chain already simulated under the reduced model and permute these
  logp_sigma <- .Call("p_sigma2_zpermuted", kmod)
  lps1 <- log(mean(exp(logp_sigma)))

  ans <- gtools::ddirichlet(x=c(0.2, 0.3, 0.5), alpha=c(100, 150, 200))
  result <- .Call("ddirichlet_", c(0.2, 0.3, 0.5), c(100, 150, 200))
  checkEquals(log(ans), result)
##  batchExperiment(se, outdir, mcmcp.list=mcmp.list)
##  B <- getFiles(outdir, rownames(se), "batch")
##  checkTrue(is.null(readRDS(model(B))))

  ## simulate a bogus batch
##  se$plate <- factor(rep(letters[1:3], length.out=ncol(se)))
##  batchExperiment(se, outdir, mcmcp.list=mcmp.list)
##  B <- getFiles(outdir, rownames(se), "batch")
  ##  checkTrue(is.null(readRDS(model(B))))


  setClass("foo", representation(x="numeric"))
  object <- new("foo")
  object@x <- rnorm(10)

  object2 <- .Call("test_clone", object)

}

.test_cnp360 <- function(){
  library(GenomicRanges)
  se360 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp360.rds")
  r <- assays(se360)[["mean"]][1, ]
  if(FALSE) hist(r, breaks=1000, col="gray", border="gray")
  x <- computeMarginalLik(r, nchains=3)
  ## hemizygous deletion is rare
  models <- orderModels(x)
  checkTrue(k(models)[1] %in% c(2, 3))
  cn <- map(models[[1]])
  freq_hem <- mean(cn==1)
  ## get homozygous post-hoc
  cn[ r < -4 ] <- 0


  probs <- seq(0, 1, 0.001)
  rgrid <- cut(r, breaks=quantile(r, probs=probs), label=FALSE,
               include.lowest=TRUE)
  rcenters <- sapply(split(r, rgrid), mean)

  ##
  ## Computationally not feasible to explore all K models adequately
  ## with 10,000 samples
  ##
  index <- sample(seq_along(r), 2000)
  table(cn[index])
  xx <- computeMarginalLik(r[index], nchains=3)
  models2 <- orderModels(xx)
  checkTrue(k(models2[[1]]) == 2)
  df <- CNPBayes:::imputeFromSampledData(models2[[1]], r, index)
  checkTrue(!any(is.na(df$cn)))
  cn_complete <- map(models[[1]])
  cn_incomplete <- df$cn
  checkTrue(mean(cn_incomplete != cn_complete) < 0.001)

  p_complete <- cnProbability(probz(models[[1]]), 2)
  p_incomplete <- df$p
  checkTrue(mean(abs(p_complete - p_incomplete) > 0.1 ) < 0.01)
}

.test_running_models_in_sequence <- function(){
  ##
  ## Start with k=4 model.  Inititialize k=3 model with modes found in
  ## the k=4 model.
  ##
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth)
  ##mp <- McmcParams(iter=500, burnin=500)

  mp <- McmcParams(iter=100, burnin=100, nStarts=1)
  model <- MarginalModel(data=y(truth), k=4, mcmc.params=mp)
  set.seed(123)
  model <- posteriorSimulation(model)
  model.nested <- NestedMarginalModel(model)
  checkTrue(k(model.nested) == k(model) -1L)
  checkEquals(sort(theta(model.nested)), c(-1, 0, 1), tolerance=0.1)
  model.nested <- posteriorSimulation(model.nested)
  checkEquals(sort(theta(model.nested)), c(-1, 0, 1), tolerance=0.1)
  model.nested <- NestedMarginalModel(model.nested)
  checkTrue(k(model.nested) ==2)
  model.nested <- NestedMarginalModel(model.nested)
  checkTrue(k(model.nested)==1)
  checkException(NestedMarginalModel(model.nested))
}

##
## Good example of heterogeneity in batch means, and some batches do
## not have the homozygous deletion because it is rare
##
.test_missingcomponent <- function(){
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

    plot(x)
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
      plot(x2)
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
    plot(xx)
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

test_computemeans <- function(){
  dir <- system.file("unitTests", package="CNPBayes")
  testdat <- readRDS(file.path(dir, "test_data.rds"))
  model <- BatchModel(data=testdat$y, batch=testdat$b)
  mns <- CNPBayes:::computeMeans(model)
  mns2 <- matrix(as.numeric(CNPBayes:::.computeMeansBatch(model)), 4, 2)
  checkEquals(mns, mns2)
}

