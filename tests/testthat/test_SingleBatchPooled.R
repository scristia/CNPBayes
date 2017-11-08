context("Pooled variance")

test_that("sigma2_pooled", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mp <- McmcParams(iter = 10, burnin = 10)
  hp <- Hyperparameters(k=3)
  model <- SingleBatchPooled(y(truth), hp, mp)
  model <- .posteriorSimulation2(model)

  set.seed(1223)
  (s2.cpp <- sigma2_pooled(model))

  squared <- (y(model)-theta(model)[z(model)])^2
  ss <- sum(squared)
  nu.n <- 0.5*(nu.0(model) + length(y(model)))
  s2.n <- 0.5*(nu.0(model)*sigma2.0(model) + ss)
  set.seed(1223)
  prec.R <- rgamma(1, nu.n, rate=s2.n)
  (s2.R <- 1/prec.R)
  expect_equal(s2.R, s2.cpp)
})

test_that("sigma2_heavy", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3/sqrt(10), 0.15/sqrt(10), 0.15/sqrt(10)),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10),
                        df=10)
#   truth <- simulateData(N = 1000, theta = c(-2, -0.2, 0),
#                         sds = c(0.3, 0.15, 0.15),
#                         p = c(0.005, 1/10, 1 - 0.005 - 1/10),
#                         df=100)
#   truth <- simulateData(N = 1000, theta = 0, sds = 0.15, p = 1, df=100)
  mp <- McmcParams(iter = 10, burnin = 10)
  mp <- McmcParams(iter = 5000, burnin = 500)
  hp <- Hyperparameters(k=3, dfr=10)
  model <- SBPt(y(truth), hp, mp)
  model <- .posteriorSimulation2(model)
x <- seq(-2.5, 2, length.out=1000)
hist(y(truth), breaks=200, freq=FALSE)
thetas <- apply(theta(chains(model)), 2, mean)
thetas <- theta(model2)
lines(x, p(model)[1]*dlocScale_t(x, dfr(model), thetas[1], sqrt(sigma2(model))))
lines(x, p(model)[2]*dlocScale_t(x, dfr(model), thetas[2], sqrt(sigma2(model))))
lines(x, p(model)[3]*dlocScale_t(x, dfr(model), thetas[3], sqrt(sigma2(model))))
  set.seed(1223)
  (s2.cpp <- sigma2_heavy(model))
#   (s2.cpp <- sigma2_pooled(model))

  squared <- (y(model)-theta(model)[z(model)])^2
  ss <- sum(squared)
  nu.n <- 0.5*(nu.0(model) + length(y(model)))
  s2.n <- 0.5*(nu.0(model)*sigma2.0(model) + ss)
  set.seed(1223)
  prec.R <- rgamma(1, nu.n, rate=s2.n)
  (s2.R <- 1/prec.R)
  expect_equal(s2.R, s2.cpp)
})


.test_that <- function(nm, expr) NULL


## no longer updating z in the way tested below
.test_that("missing component", {
  set.seed(2000)
  truth <- simulateData(N = 45, theta = c(-2, 0),
                        sds = c(0.1, 0.1),
                        p = c(1/50, 49/50))
  mp <- McmcParams(iter = 10, burnin = 0)
  hp <- Hyperparameters(k=2)
  yy <- y(truth)
  ## idea: the second observation should be re-assigned to the first component
  ## with greater probability than the other observations
  yy[1] <- -2
  yy[2] <- -0.75
  model <- SingleBatchPooled(yy, hp, mp)
  updates <- .param_updates()
  updates["theta"] <- 0L
  mp <- McmcParams(iter=5, burnin=0, param_updates=updates, nStarts=0)
  theta(model) <- c(-2, 0)
  sigma2(model) <- 0.1^2
  truez <- c(1L, rep(2L, length(yy)-1))
  z(model) <- truez
  expectedz <- truez
  expectedz[2] <- 1L
  p <- multinomialPr_pooled(model)
  model@probz <- p
  zz <- z_pooled(model)
  expect_identical(zz, expectedz)
})



.test_that("segfault", {
  ##
  ## the problem was with the update for z
  ##
  ## -- commented the code that tried to assign samples with the best chance of
  ## -- belonging to the missing component label
  load_all()
  library(testthat)
  library(purrr)
  set.seed(2000)
  truth <- simulateData(N = 45, theta = c(-2, 0),
                        sds = c(0.1, 0.1),
                        p = c(1/3, 2/3))
  mp <- McmcParams(iter = 1000, burnin = 100, nStarts=4)
  hp <- Hyperparameters(k=5,
                        mu=-0.75,
                        tau2.0=0.4,
                        eta.0=32,
                        m2.0=0.5)
  dat <- y(truth)
  set.seed(123)
  sb <- SingleBatchPooled(dat=dat, hp=hp, mp=mp)
  .posteriorSimulation2(sb)

  post <- sb
  post <- runBurnin(post)
  if(!isOrdered(post)) label_switch(post) <- TRUE
  post <- sortComponentLabels(post)
  if( iter(post) < 1 ) return(post)
  post <- runMcmc(post)

  modes(post) <- computeModes(post)
  if(isOrdered(post)){
    label_switch(post) <- FALSE
    return(post)
  }
  ## not ordered: try additional MCMC simulations
  label_switch(post) <- TRUE
  post <- sortComponentLabels(post)
  ## reset counter for posterior probabilities
  post@probz[] <- 0
  post <- runMcmc(post)
  modes(post) <- computeModes(post)
  ##mcmcParams(post) <- mp.orig
  if(isOrdered(post)){
    label_switch(post) <- FALSE
    return(post)
  }
  label_switch(post) <- TRUE
  if(params[["warnings"]]) {
    ##
    ## at this point, we've tried to run the twice after burnin and we still
    ## have mixing. Most likely, we are fitting a model with k too big
    warning("label switching: model k=", k(post))
  }
  post <- sortComponentLabels(post)


  sb2 <- .posteriorSimulation2(sb)
  .posteriorSimulation2(sb)
})


##
## Too long for unit test
##
.test_that <- function(nm, expr) NULL

.test_that("SingleBatchPooled", {
  library(purrr)
  sb <- SingleBatchPooled()
  expect_true(validObject(sb))
  set.seed(2000)
  truth <- simulateData(N = 45, theta = c(-2, 0),
                        sds = c(0.1, 0.1),
                        p = c(1/3, 2/3))
  mp <- McmcParams(iter = 1000, burnin = 1000, nStarts=4)
  hp <- Hyperparameters(k=5,
                        mu=-0.75,
                        tau2.0=0.4,
                        eta.0=32,
                        m2.0=0.5)
  model <- gibbs_singlebatch_pooled(hp, mp, dat=y(truth))
  model.list <- gibbsSingleBatchPooled(hp, mp, dat=y(truth))
  expect_identical(k(model.list[[1]]), 2L)
  expect_identical(names(model.list)[1], "SBP2")
})

test_that("Valid starts", {
  library(purrr)
  set.seed(2000)
  truth <- simulateData(N = 45, theta = c(-2, 0),
                        sds = c(0.1, 0.1),
                        p = c(1/3, 2/3))
  mp <- McmcParams(iter = 2, burnin = 0)
  hp <- Hyperparameters(k=2,
                        mu=-0.75,
                        tau2.0=0.4,
                        eta.0=32,
                        m2.0=0.5)
  expect_error(sb <- SingleBatchPooled(dat=y(truth), mp=mp, hp=hp))
  mp <- McmcParams(iter=2, burnin=100)
  sb <- SingleBatchPooled(dat=y(truth), mp=mp, hp=hp)
  expect_identical(log_lik(sb), loglik_pooled(sb))
  sb.list <- replicate(10, SingleBatchPooled(dat=y(truth),
                                             mp=mp, hp=hp))
  expect_true(!any(is.na(map_dbl(sb.list, log_lik))))
})
