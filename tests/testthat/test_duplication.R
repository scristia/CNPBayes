context("Duplications")

test_that("duplication_models", {
  path1 <- system.file("extdata", package="CNPBayes")
  cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))["CNP_244", ]
  snp_se <- readRDS(file.path(path1, "snp_se.rds"))
  snp_se <- subsetByOverlaps(snp_se, cnp_se)
  mb.subsamp <- readRDS(file.path(path1, "CNP_244",
                                  "mb_subsamp.rds"))
  mp <- McmcParams(burnin=100, iter=400)
  model.list <- duplication_models(mb.subsamp, snp_se, mp)
  model <- choose_model(model.list)
  expect_identical(mapping(model), c("2", "3"))
})
