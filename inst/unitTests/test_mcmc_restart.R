test_mcmc_restart <- function(){
  set.seed(1)
  truth <- simulateData(N=500, .k=3, means=c(-1, 0, 1), sds=c(0.1, 0.1, 0.1))
  model <- truth
  ##
  ## No burnin, >=1 iterations
  ##
  ## with no burnin, the first element stored in the mcmc chains should be the initial values
  ##
  mcmcp <- McmcParams(iter=1, burnin=0)
  model <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model)
  checkIdentical(theta(mc)[1, ], theta(truth))
  checkIdentical(sigma(mc)[1, ], sigma(truth))

  ##
  ## Burnin of >=1, >=1 iterations
  ##
  ## With a burnin of 1 or more, the first element of the chain will
  ## differ from the starting values
  model <- truth
  mcmcp <- McmcParams(iter=1, burnin=1)
  model <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model)
  checkTrue(!identical(theta(mc)[1, ], theta(truth)))
  checkTrue(!identical(sigma(mc)[1, ], sigma(truth)))

  ##
  ## Burnin of >=1, 0 iterations
  ##
  ## The chains will be a length-one vector containing the last value
  ## from the burnin (And will differ from the starting values)
  model <- truth
  mcmcp <- McmcParams(iter=0, burnin=1)
  model <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model)
  checkTrue(!identical(theta(mc)[1, ], theta(truth)))
  checkTrue(!identical(sigma(mc)[1, ], sigma(truth)))

  ##
  ## Burnin =0, 0 iterations
  ##
  ## The chains will be a length-one vector containing the initial
  ## values
  model <- truth
  mcmcp <- McmcParams(iter=0, burnin=0)
  model <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model)
  checkIdentical(theta(mc)[1, ], theta(truth))
  checkIdentical(sigma(mc)[1, ], sigma(truth))

  ##
  ## Restarting
  ##
  ## Restarting a chain will resume at the last saved iteration.
  mcmcp <- McmcParams(iter=10, burnin=0)
  model <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model)
  checkIdentical(theta(mc)[1, ], theta(truth))
  checkIdentical(sigma(mc)[1, ], sigma(truth))
  ## restart
  model2 <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model2)
  checkIdentical(theta(mc)[1, ], theta(model))
  checkIdentical(sigma(mc)[1, ], sigma(model))
}
