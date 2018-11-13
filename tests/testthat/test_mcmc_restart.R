context("MCMC Restart")

test_that("test_mcmc_restart", {
  set.seed(1)
  truth <- simulateData(N = 500, theta = c(-1, 0, 1), sds = c(0.1,
      0.1, 0.1), p = c(1/5, 1/6, 1 - 1/5 - 1/6))
  model <- truth
  mcmcParams(model) <- McmcParams(iter = 1, burnin = 0, nStarts=1)
  model <- cpp_mcmc(model)
  mc <- chains(model)
  ## verify that the chain value is the same as the current value slot
  expect_identical(theta(model)[1, ], theta(mc)[1, ])
})

test_that("one iter", {
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  iter(mb) <- 1L
  burnin(mb) <- 0L
  mb2 <- posteriorSimulation(mb)
  expect_true(validObject(mb2))
})
