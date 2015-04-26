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

  ## No batch is simulated -- just test the marginal model
  model <- MarginalModel(data=pc, k=3)
  if(FALSE){
    mcmcParams(model) <- mp
    model <- posteriorSimulation(model)
  }
  mp <- McmcParams(iter=200, burnin=400, nStarts=10)
  se <- as(model, "SummarizedExperiment")
  m <- marginal(se, mcmc.params=mp, maxperm=3)
  checkTrue(is(m, "MarginalModelList"))
  dlist <- CNPBayes:::modeDifference(m)
  checkEquals(dlist[[3]], matrix(0, 3, 3), tolerance=0.1)
  ##  the four mode model will likely not recover the original modes
  checkTrue(!isTRUE(all.equal(dlist[[4]], matrix(0, 3, 4), tolerance=0.1)))
  ##
  ## Additional iterations using current values to initialize chain:
  ## - *don't do anything else for k=1 and k=2.  Note that setting
  ## iter to 1 means that the chains for these objects will disappear.
  ##   -> would like to keep everything as is
  ##
  ## - additional iterations for k=3 and k=4
  ## - *in order to use current values to start the chain, must set
  ##    nStarts to 1
  ## - *must have appropriate burnin to allow modes to switch
  ## - *verify that the mode switch worked by comparing current value
  ##   of theta to values after running additional iterations
  ch3 <- chains(m[[3]])
  m2 <- m
  ## note:  updates to mcmc parameters must be in this order!
  iter(m2) <- c(0, 0, 250, 250)
  nStarts(m2) <- 1
  burnin(m2) <- c(0, 0, 250, 250)
  checkTrue(identical(chains(m[[1]]), chains(m2[[1]])))
  checkTrue(!identical(chains(m[[3]]), chains(m2[[3]])))
  m2 <- marginal(m2)
  checkTrue(identical(chains(m[[1]]), chains(m2[[1]])))
  ##
  ## summary values for K= 1 and 2 models is exactly the same
  ##
  checkTrue(identical(summary(m)[1:2, ], summary(m2)[1:2, ]))
  ##
  ## summaries for K= 3 and 4 are updated
  ##
  checkTrue(!identical(summary(m)[3:4, ], summary(m2)[3:4, ]))
  checkTrue(best(m2) == K)
  true_cn <- dat[[2]][1, , K]
  zz <- map(m2[[K]])
  checkEquals(zz, true_cn, tolerance=0.5)
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
  hist(rowMeans(x[, ,1, 3]), col="gray", breaks=80)

  K <- 3
  xx <- x[, , 1, K]
  mns <- rowMeans(xx)
  pc <- prcomp(xx, center=TRUE, scale.=TRUE)$x[, 1]
  if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc
  hist(pc, breaks=100, col="gray", border="gray")

  ## No batch is simulated -- just test the marginal model
  model <- MarginalModel(data=pc, k=3)
  mp <- McmcParams(iter=200, burnin=400, nStarts=10)
  se <- as(model, "SummarizedExperiment")
  m <- marginal(se, mcmc.params=mp, maxperm=3)
  checkTrue(best(m)==K)
}

test_simulation_difficult <- function(){

}

test_simulation_batch <- function(){

}
