test_that("extdata is current", {
  library(SummarizedExperiment)
  se <- readRDS(system.file("extdata", "simulated_se.rds", package="CNPBayes"))
  se2 <- updateObject(se)
  ## if this fails, we need to save the updated object in extdata/
  expect_identical(se, se2) 
  grl <- readRDS(system.file("extdata", "grl_deletions.rds", package="CNPBayes"))
  grl2 <- updateObject(grl)
  expect_identical(grl, grl2)
})
