context("MCMC Restart")

test_that("test_mcmc_restart", {
  set.seed(1)
  truth <- simulateData(N = 500, theta = c(-1, 0, 1), sds = c(0.1,
      0.1, 0.1), p = c(1/5, 1/6, 1 - 1/5 - 1/6))
  model <- truth
  mcmcParams(model, force = TRUE) <- McmcParams(iter = 1, burnin = 0, nStarts=0)
  model <- posteriorSimulation(model)
  ##tmp=replicate(20, theta(posteriorSimulation(model)))
  mc <- chains(model)
  expect_identical(theta(truth)[1, ], theta(mc)[1, ])
  expect_identical(sigma(truth)[1, ], sigma(mc)[1, ])
  model2 <- model
  mcmcParams(model2) <- McmcParams(iter = 10, burnin = 0, nStarts=0)
  model2 <- posteriorSimulation(model2)
  theta.first <- thetac(model2)[1, ]
  theta.last <- theta(model)
  expect_identical(theta.first, theta.last[1, ])
  model <- truth
  mcmcParams(model, force = TRUE) <- McmcParams(iter = 1, burnin = 1)
  model <- posteriorSimulation(model)
  mc <- chains(model)
  ## not the same as truth because of burnin
  expect_false(identical(theta(mc)[1, ], theta(truth)[1, ]))
  expect_false(identical(sigma(mc)[1, ], sigma(truth)[1, ]))
  ## verify the single iteration after burnin is stored as the first value in
  ## the chain
  expect_equal(as.numeric(theta(model)), theta(mc)[1, ])
  model <- truth
  mcmcParams(model, force = TRUE) <- McmcParams(iter = 0, burnin = 1)
  model <- posteriorSimulation(model)
  mc <- chains(model)
  expect_identical(0L, nrow(theta(mc)))
  model <- truth
  mcmcParams(model, force = TRUE) <- McmcParams(iter = 1, burnin = 0, nStarts=0)
  model <- posteriorSimulation(model)
  mcmcParams(model) <- McmcParams(iter=10, burnin=0, nStarts=0)
  model <- posteriorSimulation(model)
  mc <- chains(model)
  expect_identical(10L, nrow(theta(mc)))
  expect_identical(as.numeric(theta(truth)), theta(mc)[1, ])
  th <- theta(model)
  s <- sigma(model)
  model2 <- posteriorSimulation(model)
  mc <- chains(model2)
  expect_identical(as.numeric(th), theta(mc)[1, ])
  expect_identical(as.numeric(s), sigma(mc)[1, ])
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1, -0.8, -0.2, 0, 0.2, 0.8, 1, 1.2),
      nbatch, k, byrow = FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N = 2500, batch = rep(letters[1:3],
      length.out = 2500), p = c(1/4, 1/6, 1 - 1/4 - 1/6), theta = means,
      sds = sds)
  mcmcp <- McmcParams(iter = 1, burnin = 0, nStarts=0)
  set.seed(123)
  ##
  ## checks consistency of chain with slot for theta
  ##
  mcmcp <- McmcParams(nStarts=0, iter=1)
  modelk1 <- MultiBatchModel2(dat = y(truth),
                              hp=hpList(k=3)[["MB"]],
                              mp = mcmcp,
                              batches = batch(truth))
  modelk1 <- sortComponentLabels(modelk1)
  th1 <- as.numeric(theta(modelk1))
  ##
  ## 1 iteration and 0 starts
  ##
  modelk <- posteriorSimulation(modelk1)
  th2 <- thetac(modelk)[1, ]
  expect_identical(th2, th1)
})
