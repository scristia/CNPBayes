context("Trio models")

library(testthat)

simulateTrioData <- function(){
  set.seed(123)
  ##mendelian.probs <- mendelianProb(epsilon=0)
  K <- hp@k
  
  # for gp, which is for the simulation fn only, 
  # K must have 1: 1 correspondence to states
  # this means gp and hp may differ
  # adjust as required
  
  # params same K/ states as gp for simulation to run
  # p <- c(0.24, 0.34, 0.24, 0.09)
  # theta <- c(-1.2, 0.3, 1.7, 4)
  # sigma <- c(0.2, 0.2, 0.2, 0.2)
  # params <- data.frame(cbind(p, theta, sigma))
  
  p <- c(0.24, 0.34, 0.33)
  theta <- c(-1.2, 0.3, 1.7)
  sigma <- c(0.2, 0.2, 0.2)
  params <- data.frame(cbind(p, theta, sigma))
  
  gp <- geneticParams(K=3,
                      states=0:2,
                      xi=sigma,
                      mu=theta)
  mp <- McmcParams(iter=50, burnin=5)

  nbatch <- 3
  N <- 300
  dat2 <- simulate_data_multi(params, N=N,
                              batches = rep(c(1:nbatch),
                                         length.out = 3*N),
                              error=0, GP=gp)
 
  hp <- HyperparametersTrios(k = 4)
  maplabel <- c(0,1,2,2)
  # mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), gp=gp)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  
  model <- TBM(triodata=dat2$data,
               hp=hp,
               mp=mp,
               mprob=mprob,
               maplabel=maplabel
               )
}

test_that("TBM", {
  library(tidyverse)
  model <- TBM()
  expect_is(model, "TrioBatchModel")

  model <- simulateTrioData()
  model <- posteriorSimulation(model)
  expect_is(model, "TrioBatchModel")
  if(FALSE){
    ggChains(model)[[1]]
  }
})

test_that("burnin", {
  library(tidyverse)
  model <- simulateTrioData()
  mp <- McmcParams(iter=50, burnin=5)
  # runBurnin defined for TBM in posteriorSimulation
  model <- runBurnin(model)
  
  family_member <- family_member(model)

  lookup_mprobs(model, 1L, 1L)

  tmp <- matrix(NA, 16, 4)
  k <- 1
  for(i in 1:4){
    for(j in 1:4){
      tmp[k, ] <- lookup_mprobs(model, as.integer(i), as.integer(j))
      k <- k+1
    }
  }

  update_trioPr(model)
  update_z(model)
  update_offspring(model)
  update_ztrio(model)
  runMcmc(model)



})


