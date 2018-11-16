context("Failing tests that need to be fixed")

test_that("mcmc2 for pooled model", {
  ## with data
  library(SummarizedExperiment); library(tidyverse)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  ##
  ## These models should fit the data very well
  ## non-random starts, pooled model
  ## WORKS
  mb1 <- as(mb, "MultiBatchP")
  mbp <- as(mb1, "MultiBatchPooled")
  tmp <- posteriorSimulation(mbp)
  expect_equal(theta(tmp), theta(mb), tolerance=0.02)
  mb1 <- as(mbp, "MultiBatchP")
  tmp <- posteriorSimulation(mb1)
  expect_equal(theta(tmp), theta(mb), tolerance=0.02)
  ##
  ## random start
  ##
  tmp <- MultiBatchP(data=assays(mb1))
  ##theta(mb1) <- theta(tmp)
  ##sigma2(mb1) <- sigma2(tmp)
  mb1 <- posteriorSimulation(tmp)## fails
  expect_equal(theta(mb1), theta(mb), tolerance=0.02)
  expect_identical(ncol(sigma2(chains(mb1))), k(mb))
  expect_true(validObject(mb1))
  burnin(mb1) <- 0L
  if(FALSE){
    ##plot.ts(ll)
    tmp <- as(chains(mb1), "list")
    ggplot(tmp[["sigma2"]], aes(s, sqrt(value), group=k)) +
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
