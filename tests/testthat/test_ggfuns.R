context("gg functions")

test_that("ggfun", {
  set.seed(123)
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.2, 0.2, 0.2))
  tmp <- ggSingleBatch(truth)
  expect_is(tmp, "ggplot")

  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.02, 0.02, 0.02))
  ## hard to see marginal when it looks the same as the truth
  tmp <- ggSingleBatch(truth)
  expect_is(tmp, "ggplot")

  tmp <- ggSingleBatch(MarginalModelExample, bins=200)
  expect_is(tmp, "ggplot")

  bmod <- BatchModelExample
  tmp <- dnorm_poly_multibatch(bmod)
  tmp <- ggMultiBatch(bmod, bins=200)
  expect_is(tmp, "ggplot")
})
