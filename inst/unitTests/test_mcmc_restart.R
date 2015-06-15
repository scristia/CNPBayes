test_mcmc_restart <- function(){
  set.seed(1)
  truth <- simulateData(N=500, theta=c(-1, 0, 1), sds=c(0.1, 0.1, 0.1),
                        p=c(1/5, 1/6, 1-1/5-1/6))
  model <- truth
  ##
  ## No burnin, >=1 iterations
  ##
  ## with no burnin, the first element stored in the mcmc chains should be the initial values
  ##
  mcmcParams(model, force=TRUE) <- McmcParams(iter=1, burnin=0)
  model <- posteriorSimulation(model)
  mc <- chains(model)
  checkIdentical(theta(mc)[1, ], theta(truth))
  checkIdentical(sigma(mc)[1, ], sigma(truth))
  ##checkEquals(rowSums(probz(model)), rep(1.0, length(y(model))))
  ##
  ## Verify that when burnin exceeds the number of saved iterations
  ## there are no subsetting errors
  ##
  ## with no burnin, the first element stored in the mcmc chains
  ## should be the initial values
  ##
  mcmcParams(model) <- McmcParams(iter=1, burnin=10)
  ##trace(posteriorSimulation, browser)
  model <- posteriorSimulation(model)
  ##checkEquals(rowSums(probz(model)), rep(1.0, length(y(model))))
  ##
  ## Burnin of >=1, >=1 iterations
  ##
  ## With a burnin of 1 or more, the first element of the chain will
  ## differ from the starting values
  model <- truth
  mcmcParams(model, force=TRUE) <- McmcParams(iter=1, burnin=1)
  model <- posteriorSimulation(model)
  mc <- chains(model)
  checkTrue(!identical(theta(mc)[1, ], theta(truth)))
  checkTrue(!identical(sigma(mc)[1, ], sigma(truth)))
  checkEquals(theta(mc)[1, ], as.numeric(theta(model)))
  ##checkEquals(rowSums(probz(model)), rep(1.0, length(y(model))))

  model <- truth
  mcmcParams(model, force=TRUE) <- McmcParams(iter=0, burnin=1)
  model <- posteriorSimulation(model)
  mc <- chains(model)
  checkIdentical(nrow(theta(mc)), 0L)
  ##
  ## Burnin 0, 1 iteration
  ##
  model <- truth
  mcmcParams(model, force=TRUE) <- McmcParams(iter=1, burnin=0)
  model <- posteriorSimulation(model)


  ##
  ## Restarting
  ##
  ## Restarting a chain will resume at the last saved iteration.
  ##  -- a postive value for iter triggers an update of the chains
  mcmcParams(model, force=TRUE) <- McmcParams(iter=10, burnin=0)
  model <- posteriorSimulation(model)
  mc <- chains(model)
  checkIdentical(nrow(theta(mc)), 10L)
  checkIdentical(theta(mc)[1, ], as.numeric(theta(truth)))
  ## restart
  th <- theta(model)
  s <- sigma(model)
  model2 <- posteriorSimulation(model)
  mc <- chains(model2)
  checkIdentical(theta(mc)[1, ], as.numeric(th))
  checkIdentical(sigma(mc)[1, ], as.numeric(s))

  ##
  ## Batch model
  ##
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             p=c(1/4, 1/6, 1-1/4-1/6),
                             theta=means,
                             sds=sds)
  library(GenomicRanges)
  se <- as(truth, "SummarizedExperiment")
  mcmcp <- McmcParams(iter=5, burnin=0)
  library(oligoClasses)
  set.seed(123)
  modelk1 <- BatchModel(data=y(truth), k=3, mcmc.params=mcmcp,
                        batch=batch(truth))
  th1 <- as.numeric(theta(modelk1))
  modelk <- posteriorSimulation(modelk1)
  th2 <- thetac(modelk)[1,]
  checkIdentical(th1, th2)
}
