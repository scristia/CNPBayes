context("SingleBatchModel")

.test_that <- function(name, expr){}

test_that("test_marginal_empty_component", {
    set.seed(1)
    truth <- simulateData(N = 10, p = rep(1/3, 3), theta = c(-1,
        0, 1), sds = rep(0.1, 3))
    mp <- McmcParams(iter = 5, burnin = 5, nStarts = 1)
    model <- SingleBatchModel2(dat = y(truth),
                               hp=hpList(k = 3)[["SB"]],
                               mp = mp)
    expect_false(any(is.na(computeMeans(model))))
})


test_that("hard", {
  set.seed(1337)
  truth <- simulateData(N = 1000,
                        theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  ##expect_is(truth, "SingleBatchModel")
  expect_is(truth, "MultiBatchModel")
})

test_that("moderate", {
  set.seed(100)
  truth <- simulateData(N = 1000,
                        theta = c(-2, -0.4, 0), 
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.05, 0.15, 0.8),
                        df=100)
  expect_is(truth, "MultiBatchModel")
  expect_true(validObject(truth))
  ## verify that if we start at the true value, we remain in a region of
  ## high posterior probability after an arbitrary number of mcmc updates
  mcmcp <- McmcParams(iter = 1000, burnin = 300,
                      thin = 5, nStarts=10)
  model <- SingleBatchModel2(dat=y(truth),
                             hp=hpList(k = 3)[["SB"]],
                             mp = mcmcp)
  dfr(model) <- 100
  model <- startAtTrueValues(model, truth)
  model <- posteriorSimulation(model)
  expect_equal(theta(truth), theta(model), tolerance=0.2)
  expect_equal(sigma(truth), sigma(model), tolerance=0.15)
  expect_equal(p(truth), colMeans(pic(model)), tolerance=0.2)
})

test_that("easy", {
  set.seed(1)
  truth <- simulateData(N = 1500, p = rep(1/3, 3),
                        theta = c(-1, 0, 1),
                        sds = rep(0.1, 3),
                        df=100)
  mp <- McmcParams(iter = 100, burnin = 100)
  model <- SingleBatchModel2(dat = y(truth),
                             hp=hpList(k=3)[["SB"]],
                             mp = mp)
  u1 <- u(model)
  model <- posteriorSimulation(model)
  u2 <- u(model)
  expect_true(!identical(u1, u2))
  if (FALSE) {
    SingleBatchModelExample <- model
    save(SingleBatchModelExample, file = "data/SingleBatchModelExample.RData")
  }
  expect_equal(theta(model), theta(truth), tolerance=0.03)
  expect_equal(sigma(model), sigma(truth), tolerance=0.11)
  expect_equal(p(model), p(truth), tolerance=0.05)
  i <- argMax(model)
  expect_true(i == which.max(log_lik(chains(model))))
  expect_identical(sort(thetac(model)[i, ]), as.numeric(modes(model)[["theta"]]))
  expect_identical(sort(sigmac(model)[i, ]),
                   as.numeric(sort(sqrt(modes(model)[["sigma2"]]))))
})


test_that("cnProbability", {
  set.seed(1)
  truth <- simulateData(N = 2500, p = rep(1/3, 3),
                        theta = c(-1, 0, 1), sds = rep(0.1, 3))
  mp <- McmcParams(iter = 200, burnin = 100, nStarts = 1)
  model <- SingleBatchModel2(dat = y(truth),
                            hp=hpList(k = 3)[["SB"]], mp = mp)
  model <- posteriorSimulation(model)
  map_model <- mapModel(model)
  expect_identical(theta(map_model)[1, ],
                   thetac(model)[argMax(model), ])
  probs <- mapCnProbability(model)
  pz <- probz(model)
  expect_equal(head(pz), head(probs), tolerance=0.1)
})

context("Check computeMeans method")

test_that("computeMeans", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mp <- McmcParams(iter = 10, burnin = 10)
  model <- SingleBatchModel2(dat=y(truth), hp=hpList(k = 3)[["SB"]],
                             mp=mp)
  model <- posteriorSimulation(model)
  mns <- sapply(split(y(model), z(model)), mean)
  mns <- as.numeric(mns)

  mns2 <- computeMeans(model)
  expect_equal(mns, mns2[1, ])
})

.test_that("targeted seq", {
  set.seed(123)
  mp <- McmcParams(iter=500, burnin=1000, nStarts=4)
  extfile <- file.path(system.file("extdata", package="CNPBayes"),
                       "targeted_seq.txt")
  dat <- read.delim(extfile)[[1]]
  dat <- sample(dat, 500)
  ##
  ## The chains for the first theta (theta1) has not converged for a couple of
  ## the chains, but is in the right neighborhood. Nevertheless, the mixture
  ## model (ggMixture) looks good.
  ##
  ##mlist <- gibbs(model="SB", mp=mp, dat=dat, k_range=c(2, 3),
  ##max_burnin=10000)
  mlist <- gibbs(model="SB", mp=mp, dat=dat, k_range=c(3, 3),
                 max_burnin=5000)
  ##
  ## Select k=3
  ##
  expect_identical(names(mlist)[1], "SB3")
})
