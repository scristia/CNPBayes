context("Targeted seq")

##
## model selection takes too long
## - need to put this in CNPBayesUnitTests
##
test_that("targeted_seq data", {
  set.seed(123)
  mp <- McmcParams(iter=100, burnin=50, nStarts=4)
  extfile <- file.path(system.file("extdata", package="CNPBayes"),
                       "targeted_seq.txt")
  dat <- read.delim(extfile)[[1]]
  dat <- sample(dat, 500)
  expect_warning(
    ## warning about effective sample size
    mlist <- gibbs("SB",
                   mp=mp, dat=dat, k_range=c(2, 3),
                   max_burnin=50,
                   min_effsize=25,
                   batches=rep(1L, length(dat)))
  )
  ##  mlist <- gibbs_K(mp=mp, dat=dat, k_range=c(2, 3),
  ##                   max_burnin=50,
  ##                   min_effsize=25)
  ##
  ## Select k=3
  ##
  expect_identical(names(mlist)[1], "SB3")
})
