context("Data augmentation")
test_that("augment data", {
  extdata <- system.file("extdata", package="CNPBayes")
  hapmap <- readRDS(file.path(extdata, "hapmap.rds"))
  set.seed(123)
  current <- augmentData(hapmap)
  mp <- McmcParams(iter=250, burnin=100, thin=1)
  ## we expect a warning because we only have 250 iterations
  expect_warning(models <- gibbs(model="MBP",
                  dat=current$medians,
                  batch=current$batch_index,
                  k_range=c(3, 3),
                  mp=mp,
                  max_burnin=100,
                  df=100,
                  min_effsize=50))
  ## the Mcmc slots for iter is set incorrectly
  expect_identical(numBatch(chains(models[[1]])), 6L)
  if(FALSE) ggMixture(models[[1]])
})
