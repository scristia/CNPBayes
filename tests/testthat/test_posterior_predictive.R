context("Posterior predictive distribution")

test_that("MBP", {
  set.seed(9)
  mp <- McmcParams(iter=10, burnin=5, nStarts=1, thin=1)
  mb <- MultiBatchModelExample
  mbp.tmp <- as(mb, "MultiBatchPooled")
  expect_identical(length(sigma(mbp.tmp)), 3L)
  ch <- sigma2(chains(mbp.tmp))
  expect_identical(dim(ch), c(iter(mb), 3L))
})


