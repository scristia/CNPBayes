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
  ## Verify that when burnin exceeds the number of saved iterations
  ## there are no subsetting errors
  ##
  ## with no burnin, the first element stored in the mcmc chains should be the initial values
  ##
  mcmcp <- McmcParams(iter=1, burnin=10)
  model <- posteriorSimulation(model, mcmcp)

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
  checkIdentical(theta(mc)[2, ], as.numeric(theta(model)))

  ##
  ## Burnin of >=1, 0 iterations
  ##
  model <- truth
  mcmcp <- McmcParams(iter=0, burnin=1)
  model <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model)
  checkTrue(identical(theta(mc)[1, ], theta(truth)))
  checkTrue(nrow(theta(mc))==2)

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
  ## restart
  model2 <- posteriorSimulation(model, mcmcp)
  mc <- mcmcChains(model2)
  checkIdentical(theta(mc)[1, ], theta(model))
  checkIdentical(sigma(mc)[1, ], sigma(model))

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
                             .batch=rep(letters[1:3], length.out=2500),
                             .k=3,
                             .alpha=rep(1, k),
                             means=means,
                             sds=sds)
  se <- as(truth, "SummarizedExperiment")
  mcmcp <- McmcParams(iter=5, burnin=0)
  library(oligoClasses)
  set.seed(123)
  params <- ModelParams("batch", y=copyNumber(se)[1,], k=3,
                        batch=se$plate,
                        mcmc.params=mcmcp)
  modelk1 <- initializeModel(params)
  th1 <- as.numeric(theta(modelk1))
  modelk <- posteriorSimulation(modelk1, mcmcp)
  th2 <- thetac(modelk)[1,]
  checkIdentical(th1, th2)
  ##  set.seed(123)
##  tmp.file <- tempfile()
##  sink(tmp.file)
##  mmfit <- mixtools:::normalmixEM(y(modelk1), arbvar = FALSE, epsilon = 1e-03, k=k(modelk1), maxit=2000)
##  sink()
##  mus <- mmfit$mu
##  ##vars <- (mmfit$sigma[order(mmfit$mu)])^2
##  B <- nBatch(modelk1)
##  mus2 <- as.numeric(matrix(mus, B, k(modelk1), byrow=TRUE))
##  checkIdentical(th1, as.numeric(mus2))

##  mcmcp <- McmcParams(iter=5, burnin=0)
##  set.seed(123)
##  model <- fitMixtureModels(se, mcmcp, K=3)[[1]]
##  th3 <- thetac(model)[1,]
##  checkIdentical(th1, th3)c
}
