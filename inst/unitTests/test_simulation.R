test_simulation_moderate <- function(){
##  library(devtools)
##  load_all()
  arguments <- list("sl.good" = 6.25, ## separation parameter for "good" probes
                    "sl.bad" = 0.0625, ## sep param for "bad" probes
                    "prbias" = 0.03, ## probe level bias ~ N(0,prbias)
                    "n" = 0.2, ## background noise
                    "prvar" = c(19.92985, 0.06272) ## probe variance gamma parameters (shape,scale)
                    )
  dat <- CNPBayes:::simulateProbeLevel(cnvs=1, K=4, probes=10,
                                       arguments=arguments,
                                       qual="medium")

  ## dimensions are samples x probes x cnp x components
  x <- dat[[1]]
  ## data for 15th CNP under 3-component mixture
  if(FALSE) hist(rowMeans(x[, ,1, 3]), col="gray", breaks=80)

  K <- 4
  xx <- x[, , 1, K]
  mns <- rowMeans(xx)
  pc <- prcomp(xx, center=TRUE, scale.=TRUE)$x[, 1]
  if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc

  mp <- McmcParams(iter=1e3, burnin=500, nStarts=1)
  model1 <- MarginalModel(data=pc, k=1, mcmc.params=mp)
  model2 <- MarginalModel(data=pc, k=2, mcmc.params=mp)
  model3 <- MarginalModel(data=pc, k=3, mcmc.params=mp)
  model4 <- MarginalModel(data=pc, k=4, mcmc.params=mp)

  mlist <- list(posteriorSimulation(model1),
                posteriorSimulation(model2),
                posteriorSimulation(model3),
                posteriorSimulation(model4))
  x <- computeMarginalLik2(mlist, post.iter=500)
  models <- orderModels2(x)
  checkTrue(k(models)[1] == 4)
  if(FALSE) hist(pc, breaks=100, col="gray", border="gray")
}
