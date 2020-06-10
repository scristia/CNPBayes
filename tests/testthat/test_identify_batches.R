context("Identify batches")

test_that("summarize_region", {
  library(SummarizedExperiment)
  set.seed(2463)
  ##
  ## summarize_region code chunk
  ##
  path <- system.file("extdata", package="CNPBayes")
  cnp_se <- readRDS(file.path(path, "cnp_se.rds"))
  plates <- colData(cnp_se)$Sample.Plate
  ##
  ## To keep the unit test short, focus on first 20 plates
  ##
  plates1_20 <- plates[ plates %in% unique(plates)[1:20] ]
  keep <- colData(cnp_se)$Sample.Plate %in% plates1_20
  mb.subsamp <- summarize_region(cnp_se[1, keep],
                                 provisional_batch=plates1_20,
                                 THR=-1)
  ## visual check
  ## ggMixture(mb.subsamp) + xlim(c(-3.5, 1))
  expect_identical(numBatch(mb.subsamp), 5L)
  expect_identical(nrow(mb.subsamp), 1001L)
})
