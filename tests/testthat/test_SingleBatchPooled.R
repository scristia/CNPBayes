context("Pooled variance")
test_that("sigma2_pooled", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10),
                        df=100)
  mp <- McmcParams(iter = 10, burnin = 10)
  hp <- Hyperparameters(k=3)
  model <- SingleBatchPooled(y(truth), hp, mp)
  model <- posteriorSimulation(model)
  set.seed(1223)
  squared <- (y(model)-theta(model)[z(model)])^2
  ss <- sum(squared)
  nu.n <- 0.5*(nu.0(model) + length(y(model)))
  s2.n <- 0.5*(nu.0(model)*sigma2.0(model) + ss)
  set.seed(1223)
  prec.R <- rgamma(1, nu.n, rate=s2.n)
  (s2.R <- 1/prec.R)
  expect_equal(s2.R, 0.028, tolerance=0.01)
})

test_that("gibbs", {
  ## example from convergence vignette
  set.seed(1)
  N <- 200
  n <- 81
  lrr <- c(rnorm(100, -0.5, sd=0.1), rnorm(100, 0, sd=0.1))
  mp <- McmcParams(iter=50, burnin=10, nStarts=4)
  model <- SBP(dat=lrr, mp=mp, hp=hpList(k=2)[["SBP"]])
  model2 <- posteriorSimulation(model)
  expect_equal(marginalLikelihood(model2)[[1]], 37.3, tolerance=0.1)
  expect_equal(log_prob_thetap(model2, theta(model2)), 6.8, tolerance=0.1)
  expect_warning(sbp <- gibbs_multibatch_pooled(hp=hpList(k=2)[["MBP"]],
                                                mp=mp,
                                                dat=lrr,
                                                batches=rep(1L, length(lrr)),
                                                min_effsize=20))
  expect_is(theta_multibatch_pvar(sbp), "matrix")
  expect_equal(log_prob_thetap(sbp, theta(sbp)), 7.42, tolerance=0.1)
  expect_is(marginal_theta_pooled(sbp), "numeric")
})

test_that("Valid starts", {
  library(purrr)
  set.seed(2000)
  truth <- simulateData(N = 45, theta = c(-2, 0),
                        sds = c(0.1, 0.1),
                        p = c(1/3, 2/3))
  mp <- McmcParams(iter = 2, burnin = 0)
  hp <- Hyperparameters(k=2,
                        mu=-0.75,
                        tau2.0=0.4,
                        eta.0=32,
                        m2.0=0.5)
  sb <- SingleBatchPooled(dat=y(truth), mp=mp, hp=hp)
  mp <- McmcParams(iter=2, burnin=100)
  sb <- SingleBatchPooled(dat=y(truth), mp=mp, hp=hp)
  expect_equal(log_lik(sb), -40, tolerance=0.1)
  sb.list <- replicate(10, SingleBatchPooled(dat=y(truth),
                                             mp=mp, hp=hp))
  expect_true(!any(is.na(map_dbl(sb.list, log_lik))))
})
