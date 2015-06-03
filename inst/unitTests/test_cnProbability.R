## provide a probabilistic estimate of copy number
test_cnProbability <- function(){

  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth, use.current=TRUE)
  mp <- McmcParams(iter=200, burnin=100, nStarts=1)
  model <- MarginalModel(data=y(truth), k=3, mcmc.params=mp)
  model <- posteriorSimulation(model)

  map_model <- CNPBayes:::mapModel(model)
  checkIdentical(CNPBayes:::zChain(model)[CNPBayes:::argMax(model), ],
                 z(map_model))
  checkIdentical(thetac(model)[CNPBayes:::argMax(model), ],
                 theta(map_model))
  probs <- mapCnProbability(model)
  ## the maximum a posteriori estimate should not be too far from the
  ## posterior mean
  pz <- probz(model)
  checkEquals(head(probs), head(pz), tolerance=0.01)
}
