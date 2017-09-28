context("Local tests")

hardTruth <- function(p1, s, N=1000){
  set.seed(1234)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.9, -2, -1.85,
                    -0.45, -0.4, -0.35,
                    -0.1, 0, -0.05), nbatch, k, byrow=FALSE)
  sds <- matrix(s, nbatch, k)
  ##p1 <- prop_comp1
  p2 <- 20*p1
  p3 <- 1-p1-p2
  ##N <- 2000
  truth <- simulateBatchData(N,
                             theta=means,
                             sds=sds,
                             batch=rep(letters[1:3], length.out=N),
                             p=c(p1, p2, p3))
  truth
}

.test_that <- function(nm, expr) NULL



.test_that("batch overfit galaxy", {
  set.seed(1)
  ## load data
  library(MASS)
  data(galaxies)
  ## correct 78th observation
  galaxies[78] <- 26960
  galaxies2 <- (galaxies-median(galaxies))/1000
  galaxies3 <- c(galaxies2, galaxies2 + 10)
  mp <- McmcParams(burnin=200, nStarts=5, iter=1000)
  hp <- HyperparametersMultiBatch(k=3,
                                  mu=-0.75,
                                  tau2.0=0.4,
                                  eta.0=200, ## 32 default
                                  m2.0=100) ## 0.5 default
  model.list <- gibbs(model="MB", hp.list=list(MB=hp), mp=mp,
                      k_range=c(1, 3),
                      dat=galaxies3,
                      batches=rep(1:2, each=length(galaxies)))
  ##model.list <- gibbs_batch_K(hp=hp, mp=mp, k_range=c(1, 3), dat=galaxies3,
  ##batches=rep(1:2, each=length(galaxies)))
  expect_identical(names(model.list)[1], "MB3")
})


.test_that("MultiBatchPooled", {
  set.seed(123)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
                    0, 0.05),
                  nbatch, k, byrow = FALSE)
  sds <- matrix(rep(c(0.1, 0.15, 0.2), each=3),
                nrow=nbatch, ncol=k, byrow=TRUE)
  N <- 500
  truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
                                                length.out = N),
                             p = c(1/10, 1/5, 1 - 0.1 - 0.2),
                             theta = means,
                             sds = sds)
  mp <- McmcParams(iter=1000, burnin=1000, nStarts=4, thin=1)
  hp <- HyperparametersMultiBatch(k=4,
                                  mu=-0.75,
                                  tau2.0=0.4,
                                  eta.0=32,
                                  m2.0=0.5)
  batches <- batch(truth)
  set.seed(941)
  options(warn=2, error=utils::recover)
  model <- MultiBatchPooled(dat=y(truth), mp=mp, hp=hp,
                            batches=batch(truth))
  ## fit model with k=4
  expect_warning(model <- gibbs(model="MBP", mp=McmcParams(iter=10, burnin=5, nStart=4),
                                dat=y(truth),
                                batches=batch(truth)))
 
  mp <- McmcParams(iter=1000, nStarts=4, burnin=100)
  mlist <- gibbs("MBP", hp.list=list(MBP=hp), mp=mp,
                  dat=y(truth), batches=batch(truth))
  model <- mlist[[1]]
  expect_identical(k(model), 3L)
})


.test_that("test_batch_moderate", {
  set.seed(100)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
      0, 0.05), nbatch, k, byrow = FALSE)
  sds <- matrix(0.15, nbatch, k)
  ## first component has higher variance
  sds[, 1] <- 0.3
  N <- 1000
  truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
                                                length.out = N),
                             p = c(1/10, 1/5, 1 - 0.1 - 0.2),
                             theta = means,
                             sds = sds)
  mp <- McmcParams(iter=1000, burnin=1000, thin=2, nStarts = 4)
  model.list <- gibbs(c("MB", "MBP"),
                      mp=mp,
                      dat=y(truth),
                      batches=batch(truth),
                      k_range=c(3, 3))
  expect_identical(names(model.list)[1], "MB3")
  expect_equal(theta(truth), theta(model.list[[1]]), tolerance=0.1)
})

.test_that("SingleBatchModel2", {
  set.seed(1)
  truth <- simulateData(N = 200,
                        theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  yy <- y(truth)
  s <- (yy - median(yy))/sd(yy)
  mp <- McmcParams(iter = 1000, burnin = 1000, nStarts = 1)
  x <- qInverseTau2(mn=0.5, sd=0.5)
  hp <- Hyperparameters(k=3,
                        tau2.0=0.5,
                        mu.0=0,
                        eta.0=x$eta.0,
                        m2.0=x$m2.0)
  hist(sqrt(1/rgamma(1000, 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))), breaks=250)

  summary(sqrt(1/rgamma(200, 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))))
  ##mns <- rnorm(3, 0, sqrt(1/rgamma(1, 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))))
  set.seed(123)
  m <- SingleBatchModel2(dat=y(truth),
                         mp=mp,
                         hp=hp)
  m2 <- posteriorSimulation(m)
  if(FALSE){
    ggSingleBatch(m2)
    plist <- ggSingleBatchChains(m2)
    plist[[1]]
  }
  ##
  ## takes too long
  ##
  library(purrr)
  ##mp <- McmcParams(iter = 1000, burnin = 10, nStarts = 1, thin=1)
  mp <- McmcParams(iter = 1000, burnin = 1000, nStarts = 10, thin=1)
  mod.list <- replicate(4, SingleBatchModel2(dat=y(truth), 
                                             mp=mp,
                                             hp=hp))
  mod.list2 <- map(mod.list, posteriorSimulation)
  mc.list <- mcmcList(mod.list2)
  expect_is(mc.list, "mcmc.list")
  diagnostics(mod.list2)
  model <- combineModels(mod.list2)
  diagnostics(list(model))
  ggSingleBatchChains(model)[[1]]
  ggSingleBatchChains(model)[[2]]
  ggSingleBatch(model)
})

.test_that("SB3 better than SBP3", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.45, 0),
                        sds = c(0.3, 0.1, 0.1),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mp <- McmcParams(iter = 1000, burnin = 500, thin=2, nStarts=4)
  models <- gibbs(model=c("SB", "SBP"), k_range=c(3, 3),
                  dat=y(truth),
                  mp=mp)
  sapply(models, marginal_lik)
  expect_identical(names(models), c("SB3", "SBP3"))
})


.test_that("SBP3 better than SB3", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.45, 0),
                        sds = c(0.1, 0.1, 0.1),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mp <- McmcParams(iter = 1000, burnin = 500, thin=2, nStarts=4)
  models <- gibbs(model=c("SB", "SBP"), k_range=c(3, 3),
                  dat=y(truth),
                  mp=mp)
  sapply(models, marginal_lik)
  expect_identical(names(models), c("SBP3", "SB3"))
})


.test_that("pooled", {
    set.seed(100)
    truth <- simulateData(N = 2500,
                          theta = c(-2, -0.4, 0),
                          sds = c(0.3, 0.3, 0.3),
                          p = c(0.05, 0.1, 0.8))
    mcmcp <- McmcParams(iter = 500, burnin = 500, thin = 2, nStarts=0)
    model <- SingleBatchModel(y(truth), k = 3, mcmc.params = mcmcp)
    model <- startAtTrueValues(model, truth)
    model <- posteriorSimulation(model)
    expect_equal(theta(truth), theta(model), tolerance=0.15)
    s2_pooled <- CNPBayes:::sigma2_pooled(model)
    nu0_pooled <- CNPBayes:::nu0_pooled(model)
    sigma20_pooled <- CNPBayes:::sigma2_0_pooled(model)
    s_pooled <- sqrt(s2_pooled)
    expect_equal(object=s_pooled, expected=0.3, tolerance=0.03)

    ylist <- split(y(model), z(model))
    tmp <- vector("list", length(ylist))
    for (i in 1:length(ylist)) {
        y <- ylist[[i]]
        th <- theta(model)[i]

        tmp[[i]] <- sum((y-th)^2)
    }
    r_ss <- sum(unlist(tmp))

    sigma2_n <- 0.5*(nu.0(model) * sigma2.0(model) + r_ss)
    nu_n <- length(y(model))
    set.seed(123)
    (sigma2_new <- 1/rgamma(1, 0.5*nu_n, sigma2_n))
    set.seed(123)
    (sigma2_new.cpp <- CNPBayes:::sigma2_pooled(model))
    expect_equal(sigma2_new, sigma2_new.cpp, tolerance=0.01)
})

.test_that("test_marginal_pooled2", {
    set.seed(100)
    truth <- simulateData(N = 2500, theta = c(-2, -0.5, 0),
                          sds = c(0.1, 0.1, 0.1), p = c(0.05, 0.1, 0.8))
    pooled <- SingleBatchPooled(y(truth), Hyperparameters(k=3))
    pooled <- .posteriorSimulation2(pooled)
    expect_equal(theta(pooled), theta(truth), tolerance=0.01)
    expect_equal(sigma(pooled), sigma(truth)[3], tolerance=0.01)
    if(FALSE){
      plot.ts(sigmac(pooled), col="gray", ylim=c(0, 0.3))
      abline(h=mean(sigma(truth)))
      plot.ts(thetac(pooled), col="gray", plot.type="single")
      abline(h=theta(truth))
    }
  })

.test_that("hard4", {
  library(SummarizedExperiment)
  set.seed(123)
  truth <- hardTruth(p1=0.02, s = 0.1, N=500)
  se <- as(truth, "SummarizedExperiment")
  mp <- McmcParams(iter = 1000, burnin = 100)
  model <- gibbs("MBP", dat=y(truth), batches=batch(truth),
                 mp=mp, k_range=c(3,3))
##
##
##  hp <- HyperparametersMultiBatch(k=3,
##                             mu=-0.75,
##                             tau2.0=0.4,
##                             eta.0=32,
##                             m2.0=0.5)
##
##  modelk <- MultiBatchModel2(y(truth), hp, mp=mcmcp,
##                             batches=batch(truth))
  ##  model2 <- posteriorSimulation(modelk)
  model2 <- model[[1]]
  thetas <- theta(model2)
  pmix <- p(model2)
  expect_equal(theta(truth), thetas, tolerance=0.1)
  expect_equal(sigma(truth), sigma(model2), tolerance=0.15)
  expect_equal(p(truth), pmix, tolerance=0.04)
  if(FALSE){
    hp <- HyperparametersMultiBatch(k=3,
                               mu=-0.75,
                               tau2.0=0.4,
                               eta.0=32,
                               m2.0=0.5)
    mp <- McmcParams(iter=1000L, thin=5L, burnin=1000L,
                     nStarts=4L)
    models <- gibbs_batch_K(dat=y(truth), batches=batch(truth),
                            mp=mp, hp=hp)
    map_dbl(models, marginal_lik)
    map_dbl(models, k)
  }
})

.test_that("Select K easy", {
    set.seed(1)
    means <- c(-1, 0, 1)
    sds <- c(0.1, 0.2, 0.2)
    truth <- simulateData(N = 250, p = rep(1/3, 3), theta = means,
                          sds = sds)
    mp <- McmcParams(iter=1000, burnin=200, nStarts=4, thin=1)
    mlist <- gibbs("SB", k_range=c(2, 3), mp=mp, dat=y(truth))
    expect_identical(names(mlist)[1], "SB3")
})
