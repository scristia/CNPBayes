context("Identify batches")

test_that("summarize_region", {
  skip("Unit test for batches is slow")
  library(SummarizedExperiment)
  set.seed(2463)
  ##
  ## summarize_region code chunk
  ##
  path <- system.file("extdata", package="CNPBayes")
  cnp_se <- readRDS(file.path(path, "cnp_se.rds"))
  plates <- colData(cnp_se)$plate
  mb.subsamp <- summarize_region(cnp_se[1, ],
                                 provisional_batch=plates,
                                 THR=-1)
  if(FALSE){
    saveRDS(mb.subsamp, file="../../inst/extdata/CNP_001/mb_subsamp.rds")
  }
  expect_identical(numBatch(mb.subsamp), 5L)
  expect_identical(nrow(mb.subsamp), 1010L)
})
