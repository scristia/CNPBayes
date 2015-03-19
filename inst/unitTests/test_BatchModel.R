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
  hp <- HyperparametersBatch(tau2.0=1000, k=3)
  mod <- initializeBatchModel(mp, hypp=hp)
  ##checkEquals(as.numeric(y(mod)), yy)
  ##checkEquals(batch(mod), b)
  checkEquals(round(mu(mod), 2), c(-0.28, 0.23, 0.28))
  checkEquals(round(tau2(mod), 2), c(4.16, 0.04, 71.26))
  checkIdentical(uniqueBatch(mod)[1:2], c("LAITH,WI", "AWFUL,BO"))
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
    tabz <- tablez(modfit)
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

  zz <- tablez(modfit)[1:2, ]
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

  zz <- tablez(modfit)[1:2, ]
  zz <- tablez(modfit2)[1:2, ]
}

}
