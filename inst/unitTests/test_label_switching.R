.test_label_switching <- function(){
  set.seed(1000)
  truth <- simulateData(N=2500,
                        .k=3,
                        means=c(-1, 0, 1), sds=c(0.1, 0.1, 0.1))
  ##
  ## Start with default starting values
  ##
  mcmcp <- McmcParams(iter=100, burnin=200)
##  params <- ModelParams("marginal", y=y(truth), k=3)
##  model2 <- initializeModel(params)
##  model2 <- posteriorSimulation(model2, mcmcp, check_labels=TRUE)

  ## more likely to have label switching when k is too big
  params <- ModelParams("marginal", y=y(truth), k=5)
  model2 <- initializeModel(params)
  model2 <- posteriorSimulation(model2, mcmcp, check_labels=FALSE)
  model2 <- posteriorSimulation(model2, mcmcp, check_labels=TRUE)
  ##checkException(posteriorSimulation(model2, mcmcp, check_labels=TRUE))

  ## label switching in batch model
  mcmcp <- McmcParams(iter=10, burnin=0)
  params <- ModelParams("batch", y=y(truth), k=3,
                        batch=rep(letters[1:3], length.out=2500),
                        mcmc.params=mcmcp)
  bmodel <- initializeModel(params)
  bmodel <- posteriorSimulation(bmodel, mcmcp, check_labels=TRUE)

  ## Next: Once modes have been identified, rerun the chain and
  ## reorder each update with respect to distance to the modal
  ## parameters at each iteration of the mcmc.  Will this work if
  ## hyperparameters are not equivalent for the components? (priors
  ## differ for some parameters)

  ## For the marginal model, mu is a hyperparameter and has an
  ## ordering from small to big
  hyperParams(model2)
}
