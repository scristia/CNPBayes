.test_posthoc <- function(){
  library(GenomicRanges)
  se.ea <- readRDS("~/Labs/ARIC/AricCNPData/data/se_medr_EA.rds")
  se360 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp360.rds")
  ##se <- se.ea[subjectHits(findOverlaps(se360, se.ea))[1]]

##  cn <- copyNumber(se)[1,]
##  cn2 <- assays(se360)[["mean"]]
##  cn <- cn[names(cn) %in% colnames(cn2)]
##  cn2 <- cn2[,names(cn)]
##  plot(cn, cn2)
  options( warn=1 )
  mcmcp <- McmcParams(iter=500, burnin=0)
  ##mmodels <- fitMixtureModels(se, mcmcp, K=1:4)
  ##bics <- sapply(mmodels, bic)

  cn <- copyNumber(se)[1, ]
  plt <- se$plate
  params <- ModelParams("batch", y=cn, k=1,
                        batch=plt,
                        mcmc.params=mcmcp)
  model <- initializeModel(params)
  message("Defining batch variable")
  model <- collapseBatch(model, mcmcp)
  options( warn=2 )

  modelk <- posteriorSimulation(model, mcmcp)


}
