test_that("cnv pipeline", {
  library(SummarizedExperiment)
  data(snp_se, package="panc.data")
  data(cnp_se, package="panc.data")
  path <- system.file("extdata", package="CNPBayes")
  mb <- readRDS(file.path(path, "mb_subsamp.rds"))
  g <- rowRanges(cnp_se)[1]
  model <- cnv_models(mb, g, snp_se)
  expect_identical(modelName(model), "MBP3")
  expect_identical(mapping(model), c("0", "2", "2"))
})
