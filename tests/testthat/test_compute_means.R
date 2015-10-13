test_that("computeMeans", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mcmcp <- McmcParams(iter = 10, burnin = 10)
  model <- MarginalModel(y(truth), k = 3)
  model <- posteriorSimulation(model)
  mns <- sapply(split(y(model), z(model)), mean)
  mns <- as.numeric(mns)

  mns2 <- CNPBayes:::compute_means(model)
  expect_equal(mns, mns2)
})
