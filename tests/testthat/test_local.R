context("Local tests")

.test_that <- function(nm, expr) NULL

.test_that("targeted seq", {
  set.seed(123)
  mp <- McmcParams(iter=500, burnin=1000, nStarts=4)
  extfile <- file.path(system.file("extdata", package="CNPBayes"),
                       "targeted_seq.txt")
  dat <- read.delim(extfile)[[1]]
  dat <- sample(dat, 500)
  mlist <- gibbs_K(mp=mp, dat=dat, k_range=c(2, 3))
  ##
  ## Select k=3
  ##
  expect_identical(names(mlist)[1], "SB3")
})

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
  model.list <- gibbs_batch_K(hp=hp, mp=mp, k_range=c(1, 3), dat=galaxies3,
                              batches=rep(1:2, each=length(galaxies)))
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
  expect_warning(model <- gibbs_multibatch_pooled(hp,
                                                  mp=McmcParams(iter=10, burnin=5, nStart=4),
                                                  y(truth),
                                                  batches=batch(truth)))

  mp <- McmcParams(iter=1000, nStarts=4, burnin=100)
  mlist <- gibbsMultiBatchPooled(hp=hp,
                                 mp=mp,
                                 dat=y(truth),
                                 batches=batch(truth))
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
