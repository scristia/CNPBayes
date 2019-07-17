context("deletion pipeline")

test_that("summarize region", {
  library(SummarizedExperiment)
  ##
  ## summarize_region code chunk
  ##
  path <- system.file("extdata", package="CNPBayes")
  cnp_se <- readRDS(file.path(path, "cnp_se.rds"))
  plates <- colData(cnp_se)$Sample.Plate
  mb.subsamp <- summarize_region(cnp_se[1, ],
                                 provisional_batch=plates,
                                 THR=-1)
  expect_identical(summaries(mb.subsamp)$deletion_cutoff,
                   -1)
  if(FALSE){
    saveRDS(mb.subsamp, file="../../inst/extdata/mb_subsamp.rds")
  }
  expected <- readRDS("../../inst/extdata/mb_subsamp.rds")
  expect_equivalent(mb.subsamp, expected)
})
