context("K too big")

test_that("K too big", {
  library(SummarizedExperiment)
  se <- readRDS(system.file("extdata", "simulated_se.rds", package="CNPBayes"))
  grl <- readRDS(system.file("extdata", "grl_deletions.rds", package="CNPBayes"))
  cnv.region <- consensusCNP(grl, max.width=5e6)
  i <- subjectHits(findOverlaps(cnv.region, rowRanges(se)))
  med.summary <- matrixStats::colMedians(assays(se)[["cn"]][i, ], na.rm=TRUE)

  set.seed(1337)
  mp <- McmcParams(nStarts=100, burnin=0, iter=0)
  sb <- MarginalModel(data=med.summary, mcmc.params=mp, k=4)
  ## expect no errors even though model is over-parameterized
  sb <- posteriorSimulation(sb)
  expect_identical(k(sb), 4L)

  mb.list <- BatchModelList(data=med.summary, k=1:4,
                          batch=c(rep(1, 12), rep(2, 23)),
                          mcmc.params=mp)
  ## expect no errors even though model is over-parameterized
  mb.list <- posteriorSimulation(mb.list)
})
