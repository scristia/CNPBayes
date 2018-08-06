library(testthat)

context("Trio models")

.test_that <- function(nm, expr) NULL

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
  p <- c(0.24, 0.43, 0.33)
  theta <- c(-1.2, 0.3, 1.7)
  sigma2 <- c(0.05, 0.05, 0.05)
  params <- data.frame(cbind(p, theta, sigma2))
  mp <- McmcParams(iter=50, burnin=5)

  nbatch <- 1
  N <- 300
  maplabel <- c(0,1,2)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0.001)
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
  mprob.check <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  # deletion example
  maplabel <- c(0,1)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(0|f,m)", "p(1|f,m)", "father", "mother"))
  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=2,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=2,]
  mprob.check3 <- mprob.check2[,1:2]
  mprob.check4 <- (mprob.check3)/(rowSums(mprob.check3))
  expect_equal(mprob[,1], mprob.check4[,1])
  expect_equal(mprob[,2], mprob.check4[,2])
  
  # deletion example with repeat labels
  maplabel <- c(0,1,1)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(0|f,m)", "p(1|f,m)", "p(1|f,m)", "father", "mother"))
  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=2,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=2,]
  mprob.check3 <- mprob.check2[,1:3]
  mprob.check3[,3] <- mprob.check3[,2]
  mprob.check4 <- (mprob.check3)/(rowSums(mprob.check3))
  expect_equal(mprob[,1], mprob.check4[,1])
  expect_equal(mprob[,2], mprob.check4[,2])
  expect_equal(mprob[,3], mprob.check4[,2])
  # the default full deletion matrix
  maplabel <- c(2,3,4)
  mprob.check <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  # duplication example
  maplabel <- c(2,3)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(2|f,m)", "p(3|f,m)", "father", "mother"))

  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=4,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=4,]
  mprob.check3 <- mprob.check2[,1:2]
  mprob.check4 <- (mprob.check3)/(rowSums(mprob.check3))
  expect_equal(mprob[,1], mprob.check4[,1])
  expect_equal(mprob[,2], mprob.check4[,2])

  # duplication example with repeat labels
  maplabel <- c(3,3,4)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(3|f,m)", "p(3|f,m)", "p(4|f,m)", "father", "mother"))

  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=2,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=2,]
  mprob.check3 <- mprob.check2[,1:3]
  mprob.check3[,1] <- mprob.check3[,2]
  mprob.check4 <- (mprob.check3)/(rowSums(mprob.check3))
  expect_equal(mprob[,1], mprob.check4[,2])
  expect_equal(mprob[,2], mprob.check4[,2])
  expect_equal(mprob[,3], mprob.check4[,3])

  # the default full multi-allelic matrix
  maplabel <- c(0,1,2,3,4)
  mprob.check <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)

  # multi-allelic example
  maplabel <- c(1,2,3,4)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  colnames.check <- colnames(mprob)
  expect_equal(colnames.check,c("p(1|f,m)", "p(2|f,m)", "p(3|f,m)","p(4|f,m)", "father", "mother"))

  mprob.check2 <- data.frame(mprob.check)
  mprob.check2 <- mprob.check2[mprob.check2$father!=0,]
  mprob.check2 <- mprob.check2[mprob.check2$mother!=0,]
  mprob.check3 <- mprob.check2[,2:5]
  mprob.check4 <- (mprob.check3)/(rowSums(mprob.check3))
  expect_equal(mprob[,1], mprob.check4[,1])
  expect_equal(mprob[,2], mprob.check4[,2])
  expect_equal(mprob[,3], mprob.check4[,3])
  expect_equal(mprob[,4], mprob.check4[,4])
})

test_that("constructor", {
  set.seed(123)
  library(dplyr)
  maplabel <- c(0, 1, 2)
  mendel <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel,error=0)
  mns <- c(-2, 0, 2)
  comp <- sample(1:3, 300, replace=TRUE)
  triodat <- tibble(log_ratio=rnorm(300, mean=mns[comp]),
                    family_member=rep(c("f", "m", "o"), 100),
                    batches=1L,
                    copy_number=comp-1)
  test <- TBM(triodata=triodat,
              mprob=mendel,
              maplabel=maplabel)
  expect_true(validObject(test))
  yy <- y(test)
  expect_identical(yy, y(test))
  expect_identical(triodat$copy_number, triodata(test)$copy_number)

  triodat2 <- triodat %>%
    mutate(family_member=rev(family_member))
  test <- TBM(triodata=triodat2,
              mprob=mendel,
              maplabel=maplabel)
  expect_true(validObject(test))
  yy <- y(test)
  expect_identical(yy, y(test))
  expect_identical(triodat2$copy_number, triodata(test)$copy_number)
})

test_that("burnin", {
  library(tidyverse)
  N = 300
  nbatch <- 1
  p <- c(0.01, 0.18, 0.81)
  theta <- c(-4,-1.5, 0.5)
  sigma2 <- c(0.1, 0.1, 0.1)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)
  hp <- HyperparametersTrios(k = 3)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0.001)
  model <- simulate_data_multi2(params, N=N,
                                         batches = rep(c(1:nbatch),
                                                       length.out = 3*N),
                                         error=0.001, mprob, maplabel)
  zz <- model$data$copy_number
  truth <- maplabel [ zz ]
  ##cn.test <- z2cn(model)
  ##expect_identical(cn.test, truth)

  expect_true(validObject(model))
  truth <- model
  mp <- McmcParams(iter=50, burnin=5)
  model <- TBM(triodata=truth$data,
               hp=hp,
               mp=mp,
               mprob=mprob,
               maplabel=maplabel)
  # runBurnin defined for TBM in posteriorSimulation
  model <- runBurnin(model)

  family_member <- family_member(model)
  expect_is(family_member, "character")

  p <- lookup_mprobs(model, 1L, 1L)
  expect_is(p, "numeric")

  update_trioPr(model)
  update_z(model)
  update_offspring(model)
  runMcmc(model)
})

test_that("posterior predictive", {
  set.seed(123)
  library(tidyverse)
  model <- simulateTrioData()
  u1 <- u(model)
  mp <- McmcParams(iter=500, burnin=200, thin=1)
  mcmcParams(model) <- mp
  model <- posteriorSimulation(model)
  p <- ggMixture(model)
  expect_is(p, "gg")
  if(FALSE){
    ggMixture(model)
    ggChains(model)
    mp <- McmcParams(iter=4000, burnin=1000, thin=1)
    mb2 <- gibbs(model="MB", dat=y(model),
                 batches=batch(model),
                 mp=mp, k_range=c(3, 3), max_burnin=8000)
    ggMixture(mb2[[1]])
    ggChains(mb2[[1]])
    str(mb2[[1]])
    expect_identical(model@data, mb2[[1]]@data)
  }
})

test_that("update_zparent and update_zchild", {
  set.seed(123)
  library(tidyverse)
  
  p <- c(0.25, 0.5, 0.25)
  theta <- c(-2, 0.3, 1.7)
  sigma2 <- c(0.3, 0.3, 0.3)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)
  
  nbatch <- 1
  N <- 300
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  truth <- simulate_data_multi2(params, N=N,
                                batches = rep(c(1:nbatch),
                                              length.out = 3*N),
                                error=0, mprob, maplabel)
  hp <- HyperparametersTrios(k = 3)
  mp <- McmcParams(iter=100, burnin=100, thin=1, nStarts=4)
  
  model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data), hp.list = hp,
                       batches=truth$data$batches,
                       mp=mp, k_range=c(3, 3), max_burnin=100)
  
  # this set should be divisible by 3 (only updating child)
  which(model[[1]]@z != update_zchild(model[[1]]))/3
  
  # this set should all not be divisible by 3 (updating parents)
  which(model[[1]]@z != update_zparents(model[[1]]))/3
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
  sigma2 <- c(0.05, 0.05, 0.05)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)

  #p <- c(0.11, 0.26, 0.37, 0.26)
  #theta <- c(-4,-1.2, 1.5, 3)
  #sigma2 <- c(0.05, 0.05, 0.05, 0.05)
  #params <- data.frame(cbind(p, theta, sigma2))
  #maplabel <- c(0,1,2,3)

  p <- c(0.25, 0.5, 0.25)
  theta <- c(-2, 0.3, 1.7)
  sigma2 <- c(0.3, 0.3, 0.3)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)

  nbatch <- 1
  N <- 300
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  truth <- simulate_data_multi2(params, N=N,
                               batches = rep(c(1:nbatch),
                                             length.out = 3*N),
                               error=0, mprob, maplabel)
  hp <- HyperparametersTrios(k = 3)
  mp <- McmcParams(iter=2000, burnin=2000, thin=1, nStarts=3)
  
  #model <- TBM(triodata=truth$data,
  #             hp=hp,
  #             mp=mp,
  #             mprob=mprob,
  #             maplabel=maplabel)
  
  model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data), hp.list = hp,
                      batches=truth$data$batches,
                      mp=mp, k_range=c(3, 3), max_burnin=4000)

  truth_sum <- component_stats(truth$data)

  index <- model@triodata$family_member=="o"
  cn <- model@triodata$copy_number

  # fix children z
  #model@z[index] <- as.integer(cn[index] + 1)

  #fix parental z. remember to reset model
  #model@z[!index] <- as.integer(cn[!index] + 1)

  mcmcParams(model) <- mp
  model <- posteriorSimulation(model)
  ggMixture(model[[1]])
  ggChains(model[[1]])

  mp2 <- McmcParams(iter=2000, burnin=4000, thin=1)
  ##mb <- MB(dat=y(model), batches=batch(model))
  mb2 <- gibbs(model="MB", dat=truth$data$log_ratio,
               batches=batch(model[[1]]),
               mp=mp, k_range=c(3, 3), max_burnin=2000)

  #model <- MultiBatchModel2(dat=y(truth), batches=batch(truth),
  #                         hp=hpList(k=3)[["MB"]])

  ggMixture(mb2[[1]])
  ggChains(mb2[[1]])

  #model <- startAtTrueValues2(model, truth_sum, truth)
  #expect_identical(truth$data$batches, batch(model))
  #expect_identical(truth$params$theta, apply(theta(model), 2, mean))
  #expect_equal(truth$params$sigma2, apply(sigma2(model), 2, mean), tolerance = 0.1)
  #expect_identical(truth$params$p, p(model))
  #expect_identical(as.integer(truth$data$copy_number), z(model))
  #model <- posteriorSimulation(model)

  # check z are three components and that they are non-specific components
  zs <- unique(model@z)
  zs <- zs[!is.na(zs)]
  expect_equal(length(unique(zs)), length(maplabel))
  expect_equal(sort(unique(zs)), 1:hp@k)

  # check parameters similar
  # model.theta.means <- apply(theta(model),2, mean)
  expect_equal(as.numeric(model@modes$theta), truth$params$theta,
               scale=0.01, tolerance=1)
  #model.sigma2.means <- apply(sigma2(model),2, mean)
  expect_equal(as.numeric(model@modes$sigma2), truth$params$sigma2,
               scale=0.01, tolerance=1)
  expect_equal(model@pi, truth$params$p,
               scale=0.01, tolerance=0.5)

  # apply maplabel conversion
  results <- z2cn(model, maplabel)
  results.mb <- z2cn(mb2[[1]], maplabel)

  # this unit test specific to maplabel c(0,1,2) - change accordingly
   expect_equal(model@z-1, results@z)
  expect_equal(mb2[[1]]@z-1, results.mb@z)

  expect_equal(sort(unique(results@z)), sort(unique(maplabel)))
  #expect_identical(results@z, as.integer(model@triodata$copy_number))
  #expect_identical(results.mb@z, as.integer(model@triodata$copy_number))
  model <- model[[1]]
  mb <- mb2[[1]]
  true.cn <- as.integer(truth$data$copy_number)[is_offspring]
  true.component <- true.cn + 1L
  truth.par <- as.integer(truth$data$copy_number)[!is_offspring] + 1L
  #expect_identical(z(mb)[!is_offspring], as.integer(truth$data$copy_number)[!is_offspring] + 1L)
  mean(z(mb)[is_offspring] == true.component)
  mean(z(model)[is_offspring] == true.component)
  mean(z(mb)[!is_offspring] == truth.par)
  mean(z(model)[!is_offspring] == truth.par)
  
  })

test_that("fix offspring", {
  set.seed(123)
  library(tidyverse)

  p <- c(0.24, 0.43, 0.33)
  theta <- c(-2, 0.3, 1.7)
  ##sigma2 <- c(0.3, 0.3, 0.3)
  sigma2 <- rep(0.05, 3)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)
  nbatch <- 1
  N <- 300
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  truth <- simulate_data_multi2(params, N=N,
                               batches = rep(c(1:nbatch),
                                             length.out = 3*N),
                               error=0, mprob, maplabel)
  hp <- HyperparametersTrios(k = 3)
  mp <- McmcParams(iter=1000, burnin=1000, thin=2)
  model <- TBM(triodata=truth$data,
               hp=hp,
               mp=mp,
               mprob=mprob,
               maplabel=maplabel)
  expect_identical(y(model), truth$data$log_ratio)
  truth_sum <- component_stats(truth$data)
  is_offspring <- model@triodata$family_member=="o"
  true.cn <- as.integer(truth$data$copy_number)[is_offspring]
  # specific to maplabel starting with 0s
  true.component <- true.cn + 1L
  model@z[is_offspring] <- true.component
  mcmcParams(model) <- mp
  ##
  up <- rep(1L, 10)
  names(up) <- c("theta", "sigma2",
                 "pi", "mu",
                 "tau2",
                 "nu.0",
                 "sigma2.0",
                 "z.parents",
                 "z.offspring",
                 "pi.parents")
  up["z.offspring"] <- 0L
  up["z.parents"] <- 1L
  mcmcParams(model) <- McmcParams(iter=500, burnin=100, param_updates=up)
  system.time(model <- posteriorSimulation(model))
  ## test that component indices for offspring have not changed
  expect_identical(z(model)[is_offspring], true.component)
  expect_identical(y(model), truth$data$log_ratio)

  mp2 <- McmcParams(iter=500, burnin=100)
  model2 <- MB(dat=truth$data$log_ratio,
               hp=hp,
               mp=mp2,
               batches=rep(1L, nrow(truth$data)))
  expect_identical(y(model2), truth$data$log_ratio)
  system.time(mb <- posteriorSimulation(model2))
  if(FALSE){
    ggMixture(mb) ## looks great
  }
  expect_identical(truth$data$log_ratio, y(model2))
  expect_identical(y(model), y(model2))
  ##
  ## the component indices from the multibatch model are nearly identical to the
  ## true component indices
  expect_true(mean(z(mb)[is_offspring] == true.component) > 0.99)
  ##
  ## Compare parental component indices that were not fixed in trio model to
  ## multi-batch estimates
  ##
  expect_true(mean(z(mb)[!is_offspring] == z(model)[!is_offspring]) > 0.99)
  ##
  ##
  expect_equal(p(mb), p(model), tolerance=0.1)
  if(FALSE){
    ggMixture(model)
    ggChains(model)
    ggMixture(mb)
    ggChains(mb)
    # apply maplabel conversion
    results <- z2cn(model, maplabel)
    results.mb <- z2cn(mb, maplabel)
    expect_identical(results@z, as.integer(model@triodata$copy_number))
    expect_identical(results.mb@z, as.integer(model@triodata$copy_number))
  }
})

test_that("fix parent", {
  set.seed(123)
  library(tidyverse)
  p <- c(0.24, 0.43, 0.33)
  theta <- c(-2, 0.3, 1.7)
  ##sigma2 <- c(0.3, 0.3, 0.3)
  sigma2 <- rep(0.05, 3)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)
  nbatch <- 1
  N <- 300
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  truth <- simulate_data_multi2(params, N=N,
                                batches = rep(c(1:nbatch),
                                              length.out = 3*N),
                                error=0, mprob, maplabel)
  hp <- HyperparametersTrios(k = 3)
  mp <- McmcParams(iter=500, burnin=100, thin=1)
  model <- TBM(triodata=truth$data,
               hp=hp,
               mp=mp,
               mprob=mprob,
               maplabel=maplabel)
  expect_identical(y(model), truth$data$log_ratio)
  truth_sum <- component_stats(truth$data)
  is_offspring <- model@triodata$family_member=="o"
  true.cn <- as.integer(truth$data$copy_number)[!is_offspring]
  
  # specific to maplabel starting with 0s
  true.component <- true.cn + 1L
  model@z[!is_offspring] <- true.component
  mcmcParams(model) <- mp
  ##
  up <- rep(1L, 10)
  names(up) <- c("theta", "sigma2",
                 "pi", "mu",
                 "tau2",
                 "nu.0",
                 "sigma2.0",
                 "z.parents",
                 "z.offspring",
                 "pi.parents")
  up["z.parents"] <- 0L
  up["z.offspring"] <- 1L
  mcmcParams(model) <- McmcParams(iter=500, burnin=100, param_updates=up)
  model <- posteriorSimulation(model)
  #performs as expected!
  expect_identical(z(model)[is_offspring], as.integer(truth$data$copy_number)[is_offspring] + 1L)

  mp2 <- McmcParams(iter=500, burnin=100)
  model2 <- MB(dat=truth$data$log_ratio,
               hp=hp,
               mp=mp2,
               batches=rep(1L, nrow(truth$data)))
  mb <- posteriorSimulation(model2)
  if(FALSE){
    ggMixture(mb) ## looks great
    ggMixture(model) ## looks great

    pz <- probz(mb)
    pz2 <- rep(NA, length(true.component))
    for(i in seq_along(true.component)){
      j <- true.component[i]
      pz2[i] <- pz[i, j]
    }
    dat <- tibble(posterior_prob=pz2, log_ratio=truth$data$log_ratio[!is_offspring],
                  truth=true.component)
    ggplot(dat, aes(truth, log_ratio)) +
      geom_jitter(width=0.1, aes(color=posterior_prob))
  }
  ##
  ## the component indices from the multibatch model are nearly identical to the
  ## true component indices
  expect_true(mean(z(mb)[!is_offspring] == true.component) > 0.99)


  ##
  ## Compare parental component indices that were not fixed in trio model to
  ## multi-batch estimates
  ##
  expect_true(mean(z(mb)[is_offspring] == z(model)[is_offspring]) > 0.99)
  ##
  ## note that the weights seem to be wrong in the trio model
  ##
  expect_equal(p(mb), p(model), tolerance=0.3)
  if(FALSE){
    ggMixture(model)
    ggChains(model)
    ggMixture(mb)
    ggChains(mb)
    
    # apply maplabel conversion
    results <- z2cn(model, maplabel)
    results.mb <- z2cn(mb, maplabel)
    expect_identical(results@z, as.integer(model@triodata$copy_number))
    expect_identical(results.mb@z, as.integer(model@triodata$copy_number))
  }
})

test_that("fix parent more difficult", {
  set.seed(123)
  library(tidyverse)
  p <- c(0.25, 0.5, 0.25)
  theta <- c(-2, 0.3, 1.7)
  sigma2 <- c(0.2, 0.2, 0.2)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)
  nbatch <- 1
  N <- 300
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  truth <- simulate_data_multi2(params, N=N,
                                batches = rep(c(1:nbatch),
                                              length.out = 3*N),
                                error=0, mprob, maplabel)
  hp <- HyperparametersTrios(k = 3)
  model <- TBM(triodata=truth$data,
               hp=hp,
               mprob=mprob,
               maplabel=maplabel)
  expect_identical(y(model), truth$data$log_ratio)
  truth_sum <- component_stats(truth$data)
  is_offspring <- model@triodata$family_member=="o"
  true.cn <- as.integer(truth$data$copy_number)[!is_offspring]
  # specific to maplabel starting with 0s
  true.component <- true.cn + 1L
  model@z[!is_offspring] <- true.component
  ##
  up <- rep(1L, 10)
  names(up) <- c("theta", "sigma2",
                 "pi", "mu",
                 "tau2",
                 "nu.0",
                 "sigma2.0",
                 "z.parents",
                 "z.offspring",
                 "pi.parents")
  up["z.parents"] <- 0L
  up["z.offspring"] <- 1L
  mcmcParams(model) <- McmcParams(iter=1000, burnin=1000, param_updates=up)
  model <- posteriorSimulation(model)
  ## test that component indices for offspring have not changed
  expect_identical(z(model)[is_offspring], as.integer(truth$data$copy_number)[is_offspring] + 1L)
  

  mp2 <- McmcParams(iter=1000, burnin=1000)
  model2 <- MB(dat=truth$data$log_ratio,
               hp=hp,
               mp=mp2,
               batches=rep(1L, nrow(truth$data)))
  model2@z[!is_offspring] <- true.component
  mb <- posteriorSimulation(model2)
  if(FALSE){
    ggMixture(mb) ## looks great
    ggMixture(model) ## looks great
  }
  
  truth.child <- as.integer(truth$data$copy_number)[is_offspring] + 1L
  expect_identical(z(mb)[is_offspring], as.integer(truth$data$copy_number)[is_offspring] + 1L)
  mean(z(mb)[!is_offspring] == true.component)
  mean(z(model)[!is_offspring] == true.component)
  mean(z(mb)[is_offspring] == truth.child)
  mean(z(model)[is_offspring] == truth.child)
  tab <- tibble(truth=true.component,
                multibatch=z(mb)[is_offspring],
                triomodel=z(model)[is_offspring])


  ##
  ## Compare parental component indices that were not fixed in trio model to
  ## multi-batch estimates
  ##
  expect_true(mean(z(mb)[is_offspring] == z(model)[is_offspring]) > 0.99)
  ##
  ##
  expect_equal(p(mb), p(model), tolerance=0.1)
  if(FALSE){
    ggMixture(model)
    ggChains(model)
    ggMixture(mb)
    ggChains(mb)
    
    # apply maplabel conversion
    results <- z2cn(model, maplabel)
    results.mb <- z2cn(mb, maplabel)
    expect_identical(results@z, as.integer(model@triodata$copy_number))
    expect_identical(results.mb@z, as.integer(model@triodata$copy_number))
  }
})


test_that("fix offspring more difficult", {
  set.seed(123)
  library(tidyverse)
  p <- c(0.25, 0.5, 0.25)
  theta <- c(-2, 0.3, 1.7)
  sigma2 <- c(0.2, 0.2, 0.2)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)
  nbatch <- 1
  N <- 300
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  truth <- simulate_data_multi2(params, N=N,
                                batches = rep(c(1:nbatch),
                                              length.out = 3*N),
                                error=0, mprob, maplabel)
  hp <- HyperparametersTrios(k = 3)
  model <- TBM(triodata=truth$data,
               hp=hp,
               mprob=mprob,
               maplabel=maplabel)
  expect_identical(y(model), truth$data$log_ratio)
  truth_sum <- component_stats(truth$data)
  is_offspring <- model@triodata$family_member=="o"
  true.cn <- as.integer(truth$data$copy_number)[is_offspring]
  # specific to maplabel starting with 0s
  true.component <- true.cn + 1L
  model@z[is_offspring] <- true.component
  ##
  up <- rep(1L, 10)
  names(up) <- c("theta", "sigma2",
                 "pi", "mu",
                 "tau2",
                 "nu.0",
                 "sigma2.0",
                 "z.parents",
                 "z.offspring",
                 "pi.parents")
  up["z.parents"] <- 1L
  up["z.offspring"] <- 0L
  mcmcParams(model) <- McmcParams(iter=1000, burnin=1000, param_updates=up)
  model <- posteriorSimulation(model)
  ## test that component indices for offspring have not changed
  expect_identical(z(model)[is_offspring], as.integer(truth$data$copy_number)[is_offspring] + 1L)
  
  
  mp2 <- McmcParams(iter=1000, burnin=1000)
  model2 <- MB(dat=truth$data$log_ratio,
               hp=hp,
               mp=mp2,
               batches=rep(1L, nrow(truth$data)))
  model2@z[is_offspring] <- true.component
  mb <- posteriorSimulation(model2)
  if(FALSE){
    ggMixture(mb) ## looks great
    ggMixture(model) ## looks great
  }
  
  truth.par <- as.integer(truth$data$copy_number)[!is_offspring] + 1L
  expect_identical(z(mb)[!is_offspring], as.integer(truth$data$copy_number)[!is_offspring] + 1L)
  mean(z(mb)[is_offspring] == true.component)
  mean(z(model)[is_offspring] == true.component)
  mean(z(mb)[!is_offspring] == truth.par)
  mean(z(model)[!is_offspring] == truth.par)
  tab <- tibble(truth=true.component,
                multibatch=z(mb)[is_offspring],
                triomodel=z(model)[is_offspring])
  
  
  ##
  ## Compare parental component indices that were not fixed in trio model to
  ## multi-batch estimates
  ##
  expect_true(mean(z(mb)[is_offspring] == z(model)[is_offspring]) > 0.99)
  ##
  ##
  expect_equal(p(mb), p(model), tolerance=0.1)
  if(FALSE){
    ggMixture(model)
    ggChains(model)
    ggMixture(mb)
    ggChains(mb)
    
    # apply maplabel conversion
    results <- z2cn(model, maplabel)
    results.mb <- z2cn(mb, maplabel)
    expect_identical(results@z, as.integer(model@triodata$copy_number))
    expect_identical(results.mb@z, as.integer(model@triodata$copy_number))
  }
})

test_that("fix multiple values", {

  p <- c(0.25, 0.5, 0.25)
  theta <- c(-2, 0.3, 1.7)
  sigma2 <- c(0.5, 0.5, 0.5)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)
  nbatch <- 1
  N <- 300
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  truth <- simulate_data_multi2(params, N=N,
                                batches = rep(c(1:nbatch),
                                              length.out = 3*N),
                                error=0, mprob, maplabel)
  hp <- HyperparametersTrios(k = 3)
  truth_sum <- component_stats(truth$data)
  is_offspring <- model@triodata$family_member=="o"
  true.cn <- as.integer(truth$data$copy_number)
  # specific to maplabel starting with 0s
  true.component <- true.cn + 1L
  
  up <- rep(1L, 10)
  mp <- McmcParams(iter=2000, burnin=1000, thin=1, param_updates=up, nStarts=4)
  mb2 <- gibbs(model="MB", dat=truth$data$log_ratio,
              batches=rep(c(1:nbatch),
                          length.out = 3*N),
              mp=mp, k_range=c(3, 3), max_burnin=8000)

  #iter1 <- model[[1]]@mcmc.params@iter
  #burnin1 <- model[[1]]@mcmc.params@burnin
  #thin1 <- model[[1]]@mcmc.params@thin
  theta.feed <- (as.matrix(mb2[[1]]@modes$theta))
  #model[[1]]@theta <- (as.matrix(mb2[[1]]@modes$theta))
  #model[[1]]@sigma2 <- (as.matrix(mb2[[1]]@modes$sigma2))
  #model@pi <- mb2[[1]]@modes$mixprob
  
  #model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data), hp.list = hp,
  #                      batches=truth$data$batches,
  #                      mp=mp, k_range=c(3, 3), max_burnin=4000)
  
  model <- gibbs_trios2(model="TBM", dat=as.tibble(truth$data), hp.list = hp,
                       batches=truth$data$batches,
                       mp=mp, k_range=c(3, 3), max_burnin=4000, theta.feed=theta.feed)
  
  up <- rep(1L, 10)
  names(up) <- c("theta", "sigma2",
                 "pi", "mu",
                 "tau2",
                 "nu.0",
                 "sigma2.0",
                 "z.parents",
                 "z.offspring",
                 "pi.parents")
  up["theta"] <- 0L
  
  mcmcParams(model[[1]]) <- McmcParams(iter=iter1, burnin=burnin1, thin=thin1, param_updates=up)
  model <- posteriorSimulation(model[[1]])
  modes(model) <- computeModes(model)
  
  #model@theta <- (as.matrix(model@modes$theta))
  model@sigma2 <- (as.matrix(model@modes$sigma2))
  #model@pi <- model@modes$mixprob
  
  up["sigma2"] <- 0L
  mcmcParams(model) <- McmcParams(iter=iter1, burnin=burnin1, thin=thin1, param_updates=up)
  model <- posteriorSimulation(model)
  modes(model) <- computeModes(model)
  
  #model@theta <- (as.matrix(model@modes$theta))
  #model@sigma2 <- (as.matrix(model@modes$sigma2))
  model@pi <- model@modes$mixprob
  
  up["pi"] <- 0L
  mcmcParams(model) <- McmcParams(iter=iter1, burnin=burnin1, thin=thin1, param_updates=up)
  model <- posteriorSimulation(model)
  
  #model@theta <- t(as.matrix(truth_sum$mean))
  #model@sigma2 <- t(as.matrix(truth_sum$sd^2))
  #model@pi <- truth_sum$p
  
  #mp2 <- McmcParams(iter=2000, burnin=2000)
  #model2 <- MB(dat=truth$data$log_ratio,
  #             hp=hp,
  #             mp=mp2,
  #             batches=rep(1L, nrow(truth$data)))
  #mb <- posteriorSimulation(model2)
  if(FALSE){
    ggMixture(mb2[[1]])
    ggMixture(model[[1]]) 
    ##ggMixture(model[[1]]) 
  }
  
  #expect_identical(results.mb@z, as.integer(model@triodata$copy_number))

  true.cn <- as.integer(truth$data$copy_number)
  true.component <- true.cn + 1L
  mean(z(mb2[[1]]) == true.component)
  mean(z(model[[1]]) == true.component)
  true.cn <- as.integer(truth$data$copy_number)[is_offspring]
  true.component <- true.cn + 1L
  truth.par <- as.integer(truth$data$copy_number)[!is_offspring] + 1L
  #expect_identical(z(mb)[!is_offspring], as.integer(truth$data$copy_number)[!is_offspring] + 1L)
  mean(z(mb2[[1]])[is_offspring] == true.component)
  mean(z(model[[1]])[is_offspring] == true.component)
  mean(z(mb2[[1]])[!is_offspring] == truth.par)
  mean(z(model[[1]])[!is_offspring] == truth.par)

  # run values 30 times - see who does better overall
  
})

test_that("gibbs implement", {
  set.seed(123)
  library(tidyverse)
  p <- c(0.24, 0.43, 0.33)
  theta <- c(-2, 0.3, 1.7)
  sigma2 <- c(0.1, 0.1, 0.1)
  params <- data.frame(cbind(p, theta, sigma2))
  hp <- HyperparametersTrios(k = 3)
  maplabel <- c(0,1,2)
  nbatch <- 1
  N <- 300
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0.001)
  truth <- simulate_data_multi2(params, N=N,
                                batches = rep(c(1:nbatch),
                                              length.out = 3*N),
                                error=0, mprob, maplabel)
  true.cn <- as.integer(truth$data$copy_number)
  true.component <- true.cn + 1L
  
  mp2 <- McmcParams(iter=4000, burnin=2000, thin=1)
  model2 <- MB(dat=truth$data$log_ratio,
               hp=hp,
               mp=mp2,
               batches=rep(1L, nrow(truth$data)))
  mb <- gibbs(model=c("SB", "SBP"), k_range=c(3, 3),
              dat=truth$data$log_ratio,
              mp=mp2, max_burnin=2000)
  
  mp <- McmcParams(iter=1000, burnin=1000, thin=1)
  mb2 <- gibbs(model="MB", dat=truth$data$log_ratio,
               batches=truth$data$batches,
               mp=mp, k_range=c(3, 3), max_burnin=2000)
  
  
  tbm1 <- gibbs_trios(model="TBM", dat=as.tibble(truth$data), hp.list = hp,
               batches=truth$data$batches, maplabel=maplabel, mprob=mprob,
               mp=mp2, k_range=c(3, 3), max_burnin=4000)
  
  if(FALSE){
    ggMixture(mb[[1]])
    ggMixture(mb[[2]])
    ggMixture(mb2[[1]])
    ggMixture(tbm1[[1]])
    ggChains(mb[[1]])
    ggChains(mb[[2]])
    ggChains(mb2[[1]])
    ggChains(tbm1[[1]])
  }
  expect_true(mean(z(mb[[1]]) == true.component) > 0.9)
  expect_true(mean(z(mb[[2]]) == true.component) > 0.9)
  expect_true(mean(z(mb2[[1]]) == true.component) > 0.9)
  expect_true(mean(z(tbm1[[1]]) == true.component) > 0.9)
})