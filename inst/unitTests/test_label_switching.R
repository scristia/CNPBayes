.test_label_switching <- function(){
  set.seed(1000)
  truth <- simulateData(N=2500,
                        theta=c(-0.4, 0, 0.3),
                        sds=c(0.1, 0.05, 0.1),
                        p=c(0.01, 0.98, 0.01))
  plot(truth, use.current=TRUE)
  ##
  ## Start with default starting values
  ##
  mcmcp <- McmcParams(iter=1000, burnin=1000, nStarts=20, nStartIter=100)
  params <- ModelParams("marginal", y=y(truth), k=3)
  m <- initializeModel(params)
  model <- posteriorSimulation(m, mcmcp)
  pmns <- thetaMean(model)
  j <- order(pmns)
  pmix <- pMean(model)
  checkEquals(pmix[j], p(truth), tolerance=0.02)

  ##
  ## Label switching in batch model
  ##
  batch <- rep(letters[1:7], c(1500, 300, 200, 250, 15, 50, 185))
  batch <- sample(batch, 2500)
  bparams <- ModelParams("batch", y=y(truth), k=3, batch=batch)
  bmodel <- initializeBatchModel(bparams)
  bmodel <- posteriorSimulation(bmodel, mcmcp)
  mapz <- map(bmodel)
  tab <- table(batch(bmodel), mapz)
  table(batch(bmodel), z(truth))


  ## Next: Once modes have been identified, rerun the chain and
  ## reorder each update with respect to distance to the modal
  ## parameters at each iteration of the mcmc.  Will this work if
  ## hyperparameters are not equivalent for the components? (priors
  ## differ for some parameters)

  ## For the marginal model, mu is a hyperparameter and has an
  ## ordering from small to big
  hyperParams(model2)
}
