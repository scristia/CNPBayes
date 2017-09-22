context("SingleBatchModel")

.test_that <- function(name, expr){}

test_that("test_marginal_empty_component", {
    set.seed(1)
    truth <- simulateData(N = 10, p = rep(1/3, 3), theta = c(-1,
        0, 1), sds = rep(0.1, 3))
    mp <- McmcParams(iter = 5, burnin = 5, nStarts = 1)
    model <- SingleBatchModel(data = y(truth), k = 3, mcmc.params = mp)
    expect_false(any(is.na(CNPBayes:::computeMeans(model))))
})


test_that("test_marginal_few_data", {
  expect_error(model <- SingleBatchModel(data = 0:1, k = 3))
})

test_that("hard", {
    set.seed(1337)
    truth <- simulateData(N = 1000,
                          theta = c(-2, -0.4, 0),
                          sds = c(0.3, 0.15, 0.15),
                          p = c(0.005, 1/10, 1 - 0.005 - 1/10))
    expect_is(truth, "SingleBatchModel")
    if(FALSE){
      ## mcmcp <- McmcParams(iter = 1000, burnin = 500, thin = 10,
      ##                     nStarts = 20)
      ##
      ## do enough iterations so that any label switching occurs
      ##
      mcmcp <- McmcParams(iter = 500, burnin = 200, thin = 0,
                          nStarts = 20)
      model <- SingleBatchModel(y(truth), k = 3)
      model <- posteriorSimulation(model)
      i <- order(theta(model))
      expect_identical(i, 1:3)
      expect_equal(theta(truth), theta(model), tolerance=0.15)
      expect_equal(sigma(truth), colMeans(sigmac(model)), tolerance=0.1)
      expect_equal(p(truth), colMeans(pic(model)), tolerance=0.18)
      expect_identical(numberObs(truth), 1000L)
    }
})

test_that("moderate", {
    set.seed(100)
    truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0), sds = c(0.3,
        0.15, 0.15), p = c(0.05, 0.1, 0.8))
    ## verify that if we start at the true value, we remain in a region of
    ## high posterior probability after an arbitrary number of mcmc updates
    mcmcp <- McmcParams(iter = 250, burnin = 250, thin = 2, nStarts=0)
    model <- SingleBatchModel(y(truth), k = 3, mcmc.params = mcmcp)
    model <- startAtTrueValues(model, truth)
    model <- posteriorSimulation(model)
    expect_equal(theta(truth), theta(model), tolerance=0.15)
    expect_equal(sigma(truth), sigma(model), tolerance=0.15)
    expect_equal(p(truth), colMeans(pic(model)), tolerance=0.2)
  })


test_that("easy", {
    set.seed(1)
    truth <- simulateData(N = 2500, p = rep(1/3, 3), theta = c(-1,
        0, 1), sds = rep(0.1, 3))
    mp <- McmcParams(iter = 100, burnin = 100)
    model <- SingleBatchModel2(dat = y(truth),
                               hp=hpList(k=3)[["SB"]],
                               mp = mp)
    model <- posteriorSimulation(model)
    if (FALSE) {
        SingleBatchModelExample <- model
        save(SingleBatchModelExample, file = "data/SingleBatchModelExample.RData")
    }
    expect_equal(theta(model), theta(truth), tolerance=0.03)
    expect_equal(sigma(model), sigma(truth), tolerance=0.11)
    expect_equal(p(model), p(truth), tolerance=0.05)
    i <- CNPBayes:::argMax(model)
    expect_true(i == which.max(logPrior(chains(model)) + log_lik(chains(model))))
    expect_identical(sort(CNPBayes:::thetac(model)[i, ]), modes(model)[["theta"]])
    expect_identical(sort(sigmac(model)[i, ]),
                     sort(sqrt(modes(model)[["sigma2"]])))
})

test_that("Select K easy", {
    set.seed(1)
    means <- c(-1, 0, 1)
    sds <- c(0.1, 0.2, 0.2)
    truth <- simulateData(N = 250, p = rep(1/3, 3), theta = means,
                          sds = sds)
    mp <- McmcParams(iter=1000, burnin=200, nStarts=4, thin=1)
    mlist <- gibbs("SB", k_range=c(2, 3), mp=mp, dat=y(truth))
    expect_identical(names(mlist)[1], "SB3")
})
