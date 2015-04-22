test_constructor <- function(){
  mcmc.params <- McmcParams()
  checkTrue(validObject(mcmc.params))

  hypp <- Hyperparameters()
  checkTrue(validObject(hypp))

  hypp <- Hyperparameters("batch")
  checkTrue(validObject(hypp))

  mmod <- MarginalModel()
  checkTrue(validObject(mmod))

  mc <- McmcChains()
  checkTrue(validObject(mc))

  hp <- Hyperparameters()
  checkTrue(validObject(hp))

  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  checkTrue(validObject(truth))
  ##
  ##
  mmod <- MarginalModel(data=y(truth))
  checkTrue(validObject(mmod))
  ## default n. of iterations
  checkTrue(iter(mmod) == 1000)
  checkTrue(burnin(mmod) == 100)
  checkTrue(nStarts(mmod) == 1)
  checkTrue(nrow(thetac(mmod))==1000)

  mp <- McmcParams(iter=10, thin=10, burnin=100)
  ## replacement method triggers a change in the size of the chains
  mcmcParams(mmod) <- mp
  checkTrue(nrow(thetac(mmod)) == 10)

  iter(mmod) <- 1000
  checkTrue(nrow(thetac(mmod)) == 1000)

  ##
  ## Batch model
  ##
  bmod <- BatchModel()
  checkTrue(validObject(bmod))

  k <- nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             theta=means,
                             sds=sds,
                             p=c(1/5, 1/3, 1-1/3-1/5))
  checkTrue(validObject(truth))
  ## If batch is not specified, we assume there is not batch and a
  ## MarginalModel is generated
  bmod <- BatchModel(data=y(truth), k=3)
  checkTrue(is(bmod, "MarginalModel"))
  bmod <- BatchModel(data=y(truth), batch=oligoClasses::batch(truth))
  checkTrue(is(bmod, "BatchModel"))
  checkEquals(y(truth), y(bmod))
  checkEquals(batch(truth), batch(bmod))

  ## default n. of iterations
  checkTrue(iter(bmod) == 1000)
  checkTrue(burnin(bmod) == 100)
  checkTrue(nStarts(bmod) == 1)
  checkTrue(nrow(thetac(bmod))==1000)

  iter(bmod) <- 10
  checkTrue(nrow(thetac(bmod))==10)

  model.list <- CNPBayes:::modelEachMode(bmod)
  checkTrue(all(sapply(model.list, validObject)))
}
