context("gg functions")

test_that("ggfun", {
  set.seed(123)
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.2, 0.2, 0.2))
  tmp <- ggMixture(truth)
  expect_is(tmp, "ggplot")

  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.02, 0.02, 0.02))
  ## hard to see marginal when it looks the same as the truth
  tmp <- ggMixture(truth)
  expect_is(tmp, "ggplot")

  tmp <- ggMixture(SingleBatchModelExample, bins=200)
  expect_is(tmp, "ggplot")

  bmod <- MultiBatchModelExample
  tmp <- dnorm_poly_multibatch(bmod)
  tmp <- ggMixture(bmod, bins=200)
  expect_is(tmp, "ggplot")
})


test_that("multibatch plots", {
  means <- cbind(rnorm(9, -3, 0.5),
                 rnorm(9, -0.5, 0.3),
                 rnorm(9, 0, 0.3))
  means <- round(means, 1)
  sds <- matrix(0.15, nrow(means), ncol(means))
  sds[, 1] <- 0.3
  N <- 2000
  truth <- simulateBatchData(N = N,
                             batch = rep(letters[1:9], length.out = N),
                             p = c(1/10, 1/5, 1 - 0.1 - 0.2),
                             theta = means,
                             sds = sds)
  ggMixture(truth)
})
