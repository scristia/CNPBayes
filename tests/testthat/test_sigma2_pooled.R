context("Pooled variance")

test_that("sigma2_pooled", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mcmcp <- McmcParams(iter = 10, burnin = 10)
  model <- SingleBatchPooledVar(y(truth), k = 3)
  model <- posteriorSimulationPooled(model, iter=1, burnin=0, thin=1)

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

test_that("missing component", {
  set.seed(2000)
  truth <- simulateData(N = 45, theta = c(-2, 0),
                        sds = c(0.1, 0.1),
                        p = c(1/50, 49/50))
  mcmcp <- McmcParams(iter = 10, burnin = 0)
  yy <- y(truth)
  ## idea: the second observation should be re-assigned to the first component with greater probability than the other observations
  yy[1] <- -2
  yy[2] <- -0.75
  model <- SingleBatchPooledVar(yy, k = 2)
  updates <- .param_updates()
  updates["theta"] <- 0L
  mp <- McmcParams(iter=5, burnin=0, param_updates=updates, nStarts=0)
  theta(model) <- c(-2, 0)
  sigma2(model) <- 0.1^2
  truez <- c(1L, rep(2L, length(yy)-1))
  z(model) <- truez
  expectedz <- truez
  expectedz[2] <- 1L
  p <- multinomialPr_pooled(model)
  model@probz <- p
  zz <- z_pooled(model)
  expect_identical(zz, expectedz)
})
