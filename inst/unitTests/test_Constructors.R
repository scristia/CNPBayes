test_constructor <- function(){
  mcmc.params <- McmcParams()
  checkTrue(validObject(mcmc.params))

  mmod <- MarginalModel()
  checkTrue(validObject(mmod))

  bmod <- BatchModel()
  checkTrue(validObject(bmod))

  mc <- McmcChains()
  checkTrue(validObject(mc))

  mp <- ModelParams()
  checkTrue(validObject(mp))

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
  checkTrue(nrow(thetac(mmod))==1001)

  mp <- McmcParams(iter=10, thin=10, burnin=100)
  ## replacement method triggers a change in the size of the chains
  mcmcParams(mmod) <- mp
  checkTrue(nrow(thetac(mmod)) == 11)

  ##
  ##
  ##

}
