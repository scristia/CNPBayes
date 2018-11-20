context("Failing tests that need to be fixed")

test_that("mcmc2 for pooled model", {
  ## with data
  set.seed(321)
  library(SummarizedExperiment); library(tidyverse)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  mb1 <- as(mb, "MultiBatchP")
  if(FALSE){
    ##
    ## random starts. the 10th one fails
    ##
    tmp <- replicate(10, MultiBatch(data=assays(mb1),
                                    iter=300,
                                    burnin=200))
    tmp.list <- lapply(tmp, posteriorSimulation)
    th.list <- lapply(tmp.list, theta)
    tmp2 <- tmp.list[[10]]
    burnin(tmp2) <- 100000L
    tmp3 <- posteriorSimulation(tmp2)
    sigma2(tmp3)[, 1] <- sigma2(tmp3)[2, 1]
    burnin(tmp3) <- 500L
    tmp4 <- posteriorSimulation(tmp3)
  }
  tmp <- posteriorSimulation(mb1)
  ##
  ## starting from scratch
  ##
  tmp <- MultiBatchP(data=assays(mb1),
                     iter=500, burnin=400)
  ##tmp.list <- lapply(tmp, posteriorSimulation)
  mb1 <- posteriorSimulation(tmp)## fails
  expect_equal(theta(mb1), theta(mb), tolerance=0.02)
  expect_identical(ncol(sigma2(chains(mb1))), k(mb))
  expect_true(validObject(mb1))
  burnin(mb1) <- 0L
  if(FALSE){
    ##plot.ts(ll)
    tmp <- as(chains(mb1), "list")
    ggplot(tmp[["sigma2"]], aes(s, sqrt(value), group=k)) +
      geom_line() +
      facet_wrap(~k)
    ggplot(tmp[["theta"]], aes(s, value, group=b)) +
      geom_line(aes(color=b)) +
      facet_wrap(~k)
    ggplot(tmp[["p"]], aes(s, value, group=k)) +
      geom_line(aes(color=k))
    ggplot(tmp[["mu"]], aes(s, value, group=k)) +
      geom_line(aes(color=k))
    ggplot(tmp[["scalars"]], aes(s, value)) +
      geom_line() +
      facet_wrap(~parameter, scales="free")
  }
})
