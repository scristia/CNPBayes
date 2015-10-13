test_that("sigma2_pooled", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mcmcp <- McmcParams(iter = 10, burnin = 10)
  model <- CNPBayes:::SingleBatchPooledVar(y(truth), k = 3)
  model <- CNPBayes:::posteriorSimulationPooled(model, iter=1, burnin=0, thin=1)

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
