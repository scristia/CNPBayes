context("Trio models")

simulateTrioData <- function(maplabel=c(0,1,2)){
  set.seed(123)
  ##mendelian.probs <- mendelianProb(epsilon=0)

  # for gp, which is for the simulation fn only, 
  # K must have 1: 1 correspondence to states
  # this means gp and hp may differ
  # adjust as required
  
  # params same K/ states as gp for simulation to run
  # p <- c(0.24, 0.34, 0.24, 0.09)
  # theta <- c(-1.2, 0.3, 1.7, 4)
  # sigma <- c(0.2, 0.2, 0.2, 0.2)
  # params <- data.frame(cbind(p, theta, sigma))
  
  # please note params must match length of gp$K
  p <- c(0.24, 0.34, 0.33)
  theta <- c(-1.2, 0.3, 1.7)
  sigma2 <- c(0.2, 0.2, 0.2)
  params <- data.frame(cbind(p, theta, sigma2))
  
  mp <- McmcParams(iter=50, burnin=5)

  nbatch <- 3
  N <- 300
  maplabel <- c(0,1,2)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  dat2 <- simulate_data_multi2(params, N=N,
                              batches = rep(c(1:nbatch),
                                         length.out = 3*N),
                              error=0, mprob, maplabel)
  hp <- HyperparametersTrios(k = 3)
  model <- TBM(triodata=dat2$data,
               hp=hp,
               mp=mp,
               mprob=mprob,
               maplabel=maplabel)
}

test_that("mprob matrix", {
  
  # the default full deletion matrix
  maplabel <- c(0,1,2) ## 0, 0, 2   2,2,2   ## k  [0,1,2]  [1, 2, 3], [2, 3, 4], [0,0,2], [2,2,2], [1, 2, 2]
  mprob.check <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  
  # deletion example
  maplabel <- c(0,1)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(0|f,m)", "p(1|f,m)", "father", "mother"))
  
  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=2,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=2,]
  expect_equal(mprob[,1], mprob.check2[,1])
  expect_equal(mprob[,2], mprob.check2[,2])
  
  # deletion example with repeat labels
  maplabel <- c(0,1,1)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(0|f,m)", "p(1|f,m)", "p(1|f,m)", "father", "mother"))
  
  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=2,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=2,]
  expect_equal(mprob[,1], mprob.check2[,1])
  expect_equal(mprob[,2], mprob.check2[,2])
  expect_equal(mprob[,3], mprob.check2[,2])
  
  # the default full deletion matrix
  maplabel <- c(2,3,4)
  mprob.check <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  
  # duplication example
  maplabel <- c(2,3)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(2|f,m)", "p(3|f,m)", "father", "mother"))
  
  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=4,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=4,]
  expect_equal(mprob[,1], mprob.check2[,1])
  expect_equal(mprob[,2], mprob.check2[,2])
  
  # duplication example with repeat labels
  maplabel <- c(3,3,4)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(3|f,m)", "p(3|f,m)", "p(4|f,m)", "father", "mother"))
  
  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=2,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=2,]
  expect_equal(mprob[,1], mprob.check2[,2])
  expect_equal(mprob[,2], mprob.check2[,2])
  expect_equal(mprob[,3], mprob.check2[,3])
  
  # the default full multi-allelic matrix
  maplabel <- c(0,1,2,3,4)
  mprob.check <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  
  # multi-allelic example
  maplabel <- c(1,2,3,4)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(1|f,m)", "p(2|f,m)", "p(3|f,m)","p(4|f,m)", "father", "mother"))
  
  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=0,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=0,]
  expect_equal(mprob[,1], mprob.check2[,2])
  expect_equal(mprob[,2], mprob.check2[,3])
  expect_equal(mprob[,3], mprob.check2[,4])
  expect_equal(mprob[,4], mprob.check2[,5])
  
  # multi-allelic example with repeat labels
  maplabel <- c(1,2,3,3)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(1|f,m)", "p(2|f,m)", "p(3|f,m)","p(3|f,m)", "father", "mother"))
  
  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=0,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=0,]
  mprob.check2 <- mprob.check2[mprob.check2$father!=4,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=4,]
  expect_equal(mprob[,1], mprob.check2[,2])
  expect_equal(mprob[,2], mprob.check2[,3])
  expect_equal(mprob[,3], mprob.check2[,4])
  expect_equal(mprob[,4], mprob.check2[,4])
})

test_that("burnin", {
  library(tidyverse)
  model <- simulateTrioData2(maplabel=c(0,1,2))

  zz <- z(model)
  m <- model@maplabel
  truth <- m [ zz ]
  ##cn.test <- z2cn(model)
  ##expect_identical(cn.test, truth)

  expect_true(validObject(model))

  mp <- McmcParams(iter=50, burnin=5)
  # runBurnin defined for TBM in posteriorSimulation
  model <- runBurnin(model)

  family_member <- family_member(model)
  expect_is(family_member, "character")

  p <- lookup_mprobs(model, 1L, 1L)
  expect_is(p, "numeric")

  update_trioPr(model)
  update_z(model)
  update_offspring(model)
  update_ztrio(model)
  runMcmc(model)
})

test_that("full example", {
  set.seed(123)
  library(tidyverse)
  model <- simulateTrioData()

  u1 <- u(model)
  model <- posteriorSimulation(model)
  u2 <- u(model)
  expect_true(!identical(u1, u2))
  
  p <- c(0.25, 0.5, 0.25)
  theta <- c(-4,-1, 2)
  sigma2 <- c(0.3, 0.3, 0.3)
  params <- data.frame(cbind(p, theta, sigma2))
  
  #p <- c(0.11, 0.26, 0.37, 0.26)
  #theta <- c(-3.5,-1.2, 0.3, 1.7)
  #sigma2 <- c(0.3, 0.3, 0.3, 0.3)
  #params <- data.frame(cbind(p, theta, sigma2))
  
  nbatch <- 3
  N <- 300
  maplabel <- c(0,1,2)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel)
  truth <- simulate_data_multi2(params, N=N,
                              batches = rep(c(1:nbatch),
                              length.out = 3*N),
                              error=0, mprob, maplabel)
  truth_sum <- component_stats(truth$data)
  
  mp <- McmcParams(iter=300, burnin=300, nStarts = 5)
  hp <- HyperparametersTrios(k = 3)
  model <- TBM(triodata=truth$data,
               hp=hp,
               mp=mp,
               mprob=mprob,
               maplabel=maplabel)
  
  model <- startAtTrueValues2(model, truth_sum, truth)
  expect_identical(truth$data$batches, batch(model))
  expect_identical(truth$params$theta, apply(theta(model), 2, mean))
  expect_equal(truth$params$sigma2, apply(sigma2(model), 2, mean), tolerance = 0.1)
  expect_identical(truth$params$p, p(model))
  expect_identical(as.integer(truth$data$copy_number), z(model))
  
  model <- posteriorSimulation(model)

  # check z are three components and that they are non-specific components
  zs <- unique(model@z)
  zs <- zs[!is.na(zs)]
  expect_equal(length(unique(zs)), length(maplabel))
  expect_equal(sort(unique(zs)), 1:hp@k)
  
  # check parameters similar
  model.theta.means <- apply(theta(model),2, mean)
  expect_equal(model.theta.means, truth$params$theta,
               scale=0.01, tolerance=1)
  model.sigma2.means <- apply(sigma2(model),2, mean)
  expect_equal(model.sigma2.means, truth$params$sigma2,
               scale=0.01, tolerance=1)
  expect_equal(p(model), truth$params$p,
               scale=0.01, tolerance=0.5)
  
  # apply maplabel conversion
  results <- z2cn(model, maplabel)
  
  # this unit test specific to maplabel c(0,1,2) - change accordingly
  expect_equal(model@z-1, results@z)
  
  expect_equal(sort(unique(results@z)), sort(unique(maplabel)))
  expect_identical(results@z, as.integer(dat2$data$copy_number))
})
