test_BatchModel <- function(){
  library(GenomicRanges)
  library(oligoClasses)
  outdir <- tempdir()
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472_EA.rds")
  B <- getFiles(outdir, rownames(se472), "batch")

  set.seed(134)
  mp <- ModelParams("batch", y=copyNumber(se472)[1,],
                    k=3, batch=se472$batch,
                    mcmc.params=McmcParams())
  hp <- HyperparametersBatch(tau2.0=1000, k=3)
  mod <- initializeBatchModel(mp, hypp=hp)
  modfit <- posteriorSimulation(mod, McmcParams(iter=2))
  lld <- logLikData(modfit)
  checkIdentical(round(lld,3), 4466.369)
}
