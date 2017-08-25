context("Targeted seq")

test_that("targeted_seq data", {
  ##
  ## The marginal likelihood is less useful for selecting models - can we merge
  ## and then compute the marginal lik?
  ##
  set.seed(123)
  mp <- McmcParams(iter=500, burnin=1000, nStarts=25)
  extfile <- file.path(system.file("extdata", package="CNPBayes"),
                       "targeted_seq.txt")
  dat <- read.delim(extfile)[[1]]
  dat <- sample(dat, 500)
  mlist <- SingleBatchModelList(data=dat, k=2:4, mcmc.params=mp)
  expect_warning(mlist <- posteriorSimulation(mlist))
  ##expect_warning(mlist <- posteriorSimulation(mlist), "label switching: model k=4")
 ##
  ## Select k=3
  ##
  ml <- marginalLikelihood(mlist)
  select <- which.max(ml)
  expect_true(select == 2)
})
