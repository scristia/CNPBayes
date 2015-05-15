test_simulation <- function(){
  library(devtools)
  load_all()
  arguments <- list("sl.good" = 6.25, ## separation parameter for "good" probes
                    "sl.bad" = 0.0625, ## sep param for "bad" probes
                    "prbias" = 0.03, ## probe level bias ~ N(0,prbias)
                    "n" = 0.2, ## background noise
                    "prvar" = c(19.92985, 0.06272) ## probe variance gamma parameters (shape,scale)
                    )
  dat <- simulateProbeLevel(cnvs=1, K=4, probes=10,
                            arguments=arguments,
                            qual="easy")

  ## dimensions are samples x probes x cnp x components
  x <- dat[[1]]
  ## data for 15th CNP under 3-component mixture
  if(FALSE)
    hist(rowMeans(x[, ,1, 3]), col="gray", breaks=80)
  K <- 3
  xx <- x[, , 1, K]
  mns <- rowMeans(xx)
  pc <- prcomp(xx, center=TRUE, scale.=TRUE)$x[, 1]
  if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc
  if(FALSE)
    hist(pc, breaks=100, col="gray", border="gray")

  x <- computeMarginalLik(pc, nchains=3,
                          T=1000, T2=500,
                          burnin=200,
                          K=1:4)
  models <- orderModels(x)
  mp <- map(models[[1]])
  checkTrue(k(models)[1] >= 3)
  ##bf <- logBayesFactor(x)
  ##checkTrue(names(bf) == "M3-M2")
  ##checkTrue(bf > 100)
  if(FALSE){
    mcmcParams(model) <- mp
    model <- posteriorSimulation(model)
  }
}

test_simulation_moderate <- function(){
  library(devtools)
  load_all()
  arguments <- list("sl.good" = 6.25, ## separation parameter for "good" probes
                    "sl.bad" = 0.0625, ## sep param for "bad" probes
                    "prbias" = 0.03, ## probe level bias ~ N(0,prbias)
                    "n" = 0.2, ## background noise
                    "prvar" = c(19.92985, 0.06272) ## probe variance gamma parameters (shape,scale)
                    )
  dat <- simulateProbeLevel(cnvs=1, K=4, probes=10,
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

  x <- computeMarginalLik(pc, nchains=3,
                          T=1000, T2=500,
                          burnin=200,
                          K=1:4)
  models <- orderModels(x)
  checkTrue(k(models)[1] == 4)
  if(FALSE) hist(pc, breaks=100, col="gray", border="gray")
}
