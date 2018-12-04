library(testthat)

context("Trio models")

.test_that <- function(nm, expr) NULL

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
                    copy_number=comp-1,
                    id=paste0("trio_", rep(1:100, each=3)))
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
  model <- simulateTrioData()[["model"]]
  #model <- runBurnin(model)
  u1 <- u(model)
  model <- runBurnin(model)
  model <- posteriorSimulation(model)
  th <- theta(model)[1, ]
  true.theta <- c(-1.2, 0.3, 1.7)
  expect_equal(th, true.theta, tolerance=0.03)
  phat <- p(model)
  true.p <- c(0.24, 0.43, 0.33)
  expect_equal(phat, true.p, tolerance=0.1)
  p <- ggMixture(model)
  expect_is(p, "gg")
  if(FALSE){
    ggMixture(model)
    ggChains(model)
    mp <- McmcParams(iter=1000, burnin=1000, thin=1)
    mb2 <- gibbs(model="MB", dat=y(model),
                 batches=batch(model),
                 mp=mp, k_range=c(3, 3), max_burnin=1000)
    ggMixture(mb2[[1]])
    ggChains(mb2[[1]])
    str(mb2[[1]])
    expect_identical(model@data, mb2[[1]]@data)
  }
})

maplabel <- c(0,1,2)
mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0.001)

test_that("update_zparent and update_zchild", {
  set.seed(123)
  library(tidyverse)

  p <- c(0.25, 0.5, 0.25)
  theta <- c(-2, 0.3, 1.7)
  sigma2 <- c(0.3, 0.3, 0.3)
  params <- data.frame(cbind(p, theta, sigma2))

  nbatch <- 1
  N <- 300
  maplabel <- c(0,1,2)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0.001)

  truth <- simulate_data_multi2(params=params, N=N,
                                batches = rep(c(1:nbatch), length.out = 3*N),
                                error=0,
                                mendelian.probs=mprob,
                                maplabel=maplabel)
  hp <- HyperparametersTrios(k = 3)
  mp <- McmcParams(iter=200, burnin=100, thin=1, nStarts=4)
  expect_warning(model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data), hp.list = hp,
                                      batches=truth$data$batches,
                                      mp=mp, mprob=mprob, maplabel=maplabel,
                                      k_range=c(2, 3), max_burnin=100))
  ## MC: These are peculiar tests that are not evaluated
  # this set should be divisible by 3 (only updating child)
  which(model[[1]]@z != update_zchild(model[[1]]))/3
  # this set should all not be divisible by 3 (updating parents)
  which(model[[1]]@z != update_zparents(model[[1]]))/3
})

test_that("full example", {
  set.seed(123)
  model <- simulateTrioData()[["model"]]
  u1 <- u(model)
  model <- posteriorSimulation(model)
  u2 <- u(model)
  expect_true(!identical(u1, u2))

  p <- c(0.25, 0.5, 0.25)
  theta <- c(-4,-1, 2)
  sigma2 <- c(0.05, 0.05, 0.05)
  params <- data.frame(cbind(p, theta, sigma2))
  maplabel <- c(0,1,2)

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
  ##mp <- McmcParams(iter=2000, burnin=2000, thin=1, nStarts=3)
  mp <- McmcParams(iter=200, burnin=100, thin=1, nStarts=3)
  ## warning with small MCMC
  expect_warning(model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data), hp.list = hp,
                                      batches=truth$data$batches,
                                      mp=mp, mprob=mprob, maplabel=maplabel,
                                      k_range=c(3, 3), max_burnin=100))
  truth_sum <- component_stats(truth$data)
  index <- model[[1]]@triodata$family_member=="o"
  cn <- model[[1]]@triodata$copy_number
  # fix children z
  #model@z[index] <- as.integer(cn[index] + 1)
  #fix parental z. remember to reset model
  #model@z[!index] <- as.integer(cn[!index] + 1)
  ##mcmcParams(model) <- mp
  ##model <- posteriorSimulation(model)
  if(FALSE){
    ggMixture(model[[1]])
    ggChains(model[[1]])
  }
  ##mp2 <- McmcParams(iter=2000, burnin=2000, thin=1)
  mp2 <- McmcParams(iter=100, burnin=200, thin=1)
  expect_warning(mb2 <- gibbs(model="MB",
                              dat=truth$data$log_ratio,
                              batches=batch(model[[1]]),
                              mp=mp2,
                              k_range=c(3, 3),
                              max_burnin=100))
  if(FALSE){
    ggMixture(mb2[[1]])
    ggChains(mb2[[1]])
  }
##  #model <- startAtTrueValues2(model, truth_sum, truth)
##  #expect_identical(truth$data$batches, batch(model))
##  #expect_identical(truth$params$theta, apply(theta(model), 2, mean))
##  #expect_equal(truth$params$sigma2, apply(sigma2(model), 2, mean), tolerance = 0.1)
##  #expect_identical(truth$params$p, p(model))
##  #expect_identical(as.integer(truth$data$copy_number), z(model))
##  #model <- posteriorSimulation(model)
## check z are three components and that they are non-specific components
  model <- model[[1]]
  zs <- unique(model@z)
  zs <- zs[!is.na(zs)]
  expect_equal(length(unique(zs)), length(maplabel))
  expect_equal(sort(unique(zs)), 1:hp@k)

  # check parameters similar
  # model.theta.means <- apply(theta(model),2, mean)
  expect_equal(as.numeric(sort(model@modes$theta)), truth$params$theta,
               tolerance=1)
  #model.sigma2.means <- apply(sigma2(model),2, mean)
  expect_equal(as.numeric(model@modes$sigma2), truth$params$sigma2,
               tolerance=1.5)
  expect_equal(model@pi, truth$params$p,
               ##scale=0.01,
               tolerance=0.6)

  # apply maplabel conversion
  results <- z2cn(model, maplabel)
  results.mb <- z2cn(mb2[[1]], maplabel)

  # this unit test specific to maplabel c(0,1,2) - change accordingly
  expect_equal(model@z-1, results@z)
  expect_equal(mb2[[1]]@z-1, results.mb@z)

  expect_equal(sort(unique(results@z)), sort(unique(maplabel)))
  #expect_identical(results@z, as.integer(model@triodata$copy_number))
  #expect_identical(results.mb@z, as.integer(model@triodata$copy_number))
  mb <- mb2[[1]]
  true.cn <- as.integer(truth$data$copy_number)[index]
  true.component <- true.cn + 1L
  truth.par <- as.integer(truth$data$copy_number)[!index] + 1L
  #expect_identical(z(mb)[!is_offspring], as.integer(truth$data$copy_number)[!is_offspring] + 1L)
  mean(z(mb)[index] == true.component)
  mean(z(model)[index] == true.component)
  mean(z(mb)[!index] == truth.par)
  mean(z(model)[!index] == truth.par)
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

  mp <- McmcParams(iter=200, burnin=100, thin=1, nStarts=3)
  expect_warning(tbm1 <- gibbs_trios(model="TBM", dat=as.tibble(truth$data),
                                     batches=truth$data$batches,
                                     mp=mp, mprob=mprob, maplabel=maplabel,
                                     k_range=c(3, 3),
                                     max_burnin=100)[[1]])
  th <- theta(tbm1)[1, ]
  expect_equal(th, truth$params$theta, tolerance=0.05)
  expect_warning(mb2 <- gibbs(model="MB",
                              dat=truth$data$log_ratio,
                              batches=batch(tbm1),
                              mp=mp,
                              k_range=c(3, 3),
                              max_burnin=100))
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
  expect_true(mean(z(mb2[[1]]) == true.component) > 0.9)
  #expect_true(mean(z(mb[[2]]) == true.component) > 0.9)
  #expect_true(mean(z(mb2[[1]]) == true.component) > 0.9)
  expect_true(mean(z(tbm1) == true.component) > 0.9)
})

test_that("triodata accessor", {
  m <- simulateTrioData()[["model"]]
  expect_is(chains(m), "McmcChainsTrios")
  dat2 <- triodata(m, TRUE) 
  expect_identical(nTrios(m), nrow(triodata(m, TRUE)))

  expect_is(chains(m), "McmcChainsTrios")
  chains(m) <- chains(m)
  expect_is(chains(m), "McmcChainsTrios")

  iter(chains(m)) <- 1000L
  expect_is(chains(m), "McmcChainsTrios")

  mp <- mcmcParams(m)
  iter(mp) <- 1000L
  mcmcParams(m) <- mp
  expect_is(chains(m), "McmcChainsTrios")

  iter(m) <- 1000L
  burnin(m) <- 1000L
  expect_is(chains(m), "McmcChainsTrios")

  zf <- z(m)[is_father(m)]
  zm <- z(m)[is_mother(m)]
  zo <- z(m)[is_child(m)]
  expect_identical(zf, dat2$f)
  expect_identical(zm, dat2$m)
  expect_identical(zo, dat2$o)
})

test_that("posterior probability for mendelian inheritance", {
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
                                error=0.001, mprob,
                                maplabel)
  zz <- model$data$copy_number
  truth <- maplabel [ zz ]
  expect_true(validObject(model))
  truth <- model
  mp <- McmcParams(iter=50, burnin=5)
  model <- TBM(triodata=truth$data,
               hp=hp,
               mp=mp,
               mprob=mprob,
               maplabel=maplabel)
  m <- isMendelian(model)
  expect_identical(length(m), nrow(triodata(model, TRUE)))

  m <- isMendelian(chains(model))
  expect_identical(length(m), nTrios(model))
})
