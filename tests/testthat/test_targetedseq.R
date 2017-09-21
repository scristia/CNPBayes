context("Targeted seq")

test_that("targeted_seq data", {
  set.seed(123)
  mp <- McmcParams(iter=500, burnin=1000, nStarts=25)
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
