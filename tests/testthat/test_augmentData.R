context("Data augmentation")
test_that("augment data", {
  extdata <- system.file("extdata", package="CNPBayes")
  hapmap <- readRDS(file.path(extdata, "hapmap.rds"))
  current <- augmentData(hapmap)
  set.seed(123)
  mp <- McmcParams(iter=1000, burnin=500, thin=1)
  models <- gibbs(model="MBP",
                  dat=current$medians,
                  batch=current$batch_index,
                  k_range=c(3, 3),
                  mp=mp,
                  max_burnin=500,
                  df=100,
                  min_effsize=200)
  ml <- marginal_lik(models[[1]])
  expect_equal(as.integer(round(ml)), 84)
})
