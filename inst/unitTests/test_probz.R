test_probz <- function(){
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  mp <- McmcParams(iter=500, burnin=500)
  set.seed(123)
  model <- MarginalModel(data=y(truth), k=3, mcmc.params=mp)
  model <- CNPBayes:::startAtTrueValues(model, truth)
  true_z <- z(truth)
  checkEquals(z(truth), z(model))

  pz <- probz(model)
  checkTrue(all(pz==0))

  iter(model, force=TRUE) <- 2
  burnin(model) <- 0
  zz <- map(posteriorSimulation(model))
  checkEquals(zz, true_z, tolerance=0.5)

  ## intentionally reorder z's
  model2 <- CNPBayes:::modelOtherModes(model, maxperm=2)[[2]]
  z2 <- zz <- z(model2)
  checkTrue(sum(zz != true_z) > 500)
  ## zz is still ordered, but MAP should be correct after running 2
  ## iterations
  mcmcParams(model2) <- mcmcParams(model)
  model2 <- posteriorSimulation(model2)
  zz <- map(model2)
  checkEquals(zz, true_z, tolerance=0.5)

  model3 <- CNPBayes:::modelOtherModes(model, maxperm=3)[[3]]
  z3 <- z(model3)
  checkTrue(sum(z3 != true_z) > 500)
  checkTrue(sum(z2 != z3) > 500)

  model3 <- posteriorSimulation(model3)
  mz3 <- map(model3)
  table(mz3, true_z)
  checkEquals(mz3, true_z, tolerance=0.5)
}
