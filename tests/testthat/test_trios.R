context("Trio models")

test_that("TBM", {
  expect_true(validObject(TBM()))

  set.seed(98765)
  ##mendelian.probs <- mendelianProb(epsilon=0)
  p <- c(0.09, 0.24, 0.34, 0.24, 0.09)
  theta <- c(-3.5,-1.2, 0.3, 1.7, 4)
  sigma <- c(0.2, 0.2, 0.2, 0.2, 0.2)
  params <- data.frame(cbind(p, theta, sigma))
  gp <- geneticParams(K=5, states=0:4, xi=c(0.2, 0.2, 0.2, 0.2, 0.2), 
                      mu=c(-3.5, -1.2, 0.3, 1.7, 4))
  dat2 <- simulate_data_multi(params, N=500, error=0, gp)


})


