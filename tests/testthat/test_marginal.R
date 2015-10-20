test_that("test_constraint", {
    baf <- readRDS(system.file("extdata", "baf.rds", package = "CNPBayes"))
    set.seed(17)
    model <- MarginalModel(baf, k = 2)
    model <- posteriorSimulation(model)
})

test_that("test_marginal_empty_component", {
    set.seed(1)
    truth <- simulateData(N = 10, p = rep(1/3, 3), theta = c(-1,
        0, 1), sds = rep(0.1, 3))
    mp <- McmcParams(iter = 5, burnin = 5, nStarts = 1)
    model <- MarginalModel(data = y(truth), k = 3, mcmc.params = mp)
    expect_false(any(is.na(CNPBayes:::computeMeans(model))))
})

test_that("test_marginal_few_data", {
    expect_error(model <- MarginalModel(data = 0:1, k = 3))
})

test_that("test_marginal_hard", {
    set.seed(1337)
    truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0), sds = c(0.3,
        0.15, 0.15), p = c(0.005, 1/10, 1 - 0.005 - 1/10))
    if (FALSE)
        plot(truth)
    mcmcp <- McmcParams(iter = 1000, burnin = 500, thin = 10,
        nStarts = 20)
    model <- MarginalModel(y(truth), k = 3)
    model <- posteriorSimulation(model)
    i <- order(theta(model))
    expect_equal(theta(truth), theta(model)[i], tolerance=0.15)
    expect_equal(sigma(truth), colMeans(sigmac(model))[i], tolerance=0.1)
    expect_equal(p(truth), colMeans(pic(model))[i], tolerance=0.18)
    if (FALSE) {
        op <- par(mfrow = c(1, 2), las = 1)
        plot(truth)
        plot(model)
        par(op)
        plot.ts(sigmac(model), col = 1:3, plot.type = "single")
    }
})

test_that("test_marginal_Moderate", {
    set.seed(100)
    truth <- simulateData(N = 2500, theta = c(-2, -0.4, 0), sds = c(0.3,
        0.15, 0.15), p = c(0.05, 0.1, 0.8))
    if (FALSE)
        plot(truth)
    mcmcp <- McmcParams(iter = 500, burnin = 500, thin = 2)
    model <- MarginalModel(y(truth), k = 3, mcmc.params = mcmcp)
    model <- CNPBayes:::startAtTrueValues(model, truth)
    model <- posteriorSimulation(model)
    expect_equal(theta(truth), sort(theta(model)), tolerance=0.15)
    if (FALSE) {
        plot.ts(CNPBayes:::thetac(model), plot.type = "single")
        plot.ts(sigmac(model), plot.type = "single")
        plot.ts(pic(model), plot.type = "single")
        op <- par(mfrow = c(1, 2), las = 1)
        plot(truth)
        plot(model)
        par(op)
    }
    expect_equal(sigma(truth), sigma(model)[order(theta(model))],
                 tolerance=0.15)
    expect_equal(p(truth), colMeans(pic(model))[order(theta(model))],
                 tolerance=0.2)
  })

test_that("test_marginal_pooled", {
    set.seed(100)
    truth <- simulateData(N = 2500, theta = c(-2, -0.4, 0), sds = c(0.3,
        0.3, 0.3), p = c(0.05, 0.1, 0.8))
    if (FALSE)
        plot(truth)
    mcmcp <- McmcParams(iter = 500, burnin = 500, thin = 2)
    model <- MarginalModel(y(truth), k = 3, mcmc.params = mcmcp)
    model <- CNPBayes:::startAtTrueValues(model, truth)
    model <- posteriorSimulation(model)
    expect_equal(theta(truth), sort(theta(model)), tolerance=0.15)
    s_pooled <- sqrt(CNPBayes:::sigma2_pooled(model))
    nu0_pooled <- CNPBayes:::nu0_pooled(model)
    sigma20_pooled <- CNPBayes:::sigma2_0_pooled(model)
    expect_equal(object=s_pooled, expected=0.3, tolerance=0.02)

    ylist <- split(y(model), z(model))
    tmp <- foreach(y=ylist, th=theta(model)) %do%{
      sum((y-th)^2)
    }
    r_ss <- sum(unlist(tmp))

    sigma2_n <- 0.5*(nu.0(model) * sigma2.0(model) + r_ss)
    nu_n <- length(y(model))
    set.seed(123)
    (sigma2_new <- 1/rgamma(1, 0.5*nu_n, sigma2_n))
    set.seed(123)
    (sigma2_new.cpp <- CNPBayes:::sigma2_pooled(model))
    expect_equal(sigma2_new, sigma2_new.cpp, tolerance=0.01)

    truth <- simulateData(N = 2500, theta = c(-2, -0.5, 0),
                          sds = c(0.1, 0.1, 0.1), p = c(0.05, 0.1, 0.8))
    pooled <- CNPBayes:::SingleBatchPooledVar(data=y(truth), k=3)
    pooled <- posteriorSimulationPooled(pooled, iter=1000, burnin=0, thin=1)
    thetas <- sort(colMeans(theta(chains(pooled))[-(1:100), ]))
    expect_equal(thetas, theta(truth), tolerance=0.01)
    sigmas <- mean(tail(sigma(chains(pooled)), 900))
    expect_equal(sigmas, sigma(truth)[3], tolerance=0.01)
    if(FALSE){
      plot.ts(sigmac(pooled), col="gray", ylim=c(0, 0.3))
      abline(h=mean(sigma(truth)))
      plot.ts(thetac(pooled), col="gray", plot.type="single")
      abline(h=theta(truth))
    }
  })


test_that("test_marginalDiffK", {
    set.seed(1)
    truth <- simulateData(N = 500, p = rep(1/3, 3), theta = c(-1,
        0, 1), sds = rep(0.1, 3))
    mp <- McmcParams(iter = 5, burnin = 5, nStarts = 1)
    model <- MarginalModel(data = y(truth), k = 3, mcmc.params = mp)
    model <- posteriorSimulation(model)
    model <- posteriorSimulation(model, 2)
    expect_true(k(model) == 2)
    expect_true(length(theta(model)) == 2)
})

test_that("test_marginalEasy", {
    set.seed(1)
    truth <- simulateData(N = 2500, p = rep(1/3, 3), theta = c(-1,
        0, 1), sds = rep(0.1, 3))
    if (FALSE)
        plot(truth)
    mp <- McmcParams(iter = 5, burnin = 5, nStarts = 1)
    model <- MarginalModel(data = y(truth), k = 3, mcmc.params = mp)
    model <- CNPBayes:::startAtTrueValues(model, truth)
    model <- posteriorSimulation(model)
    if (FALSE) {
        MarginalModelExample <- model
        save(MarginalModelExample, file = "data/MarginalModelExample.RData")
    }
    if (FALSE) {
        op <- par(mfrow = c(1, 2), las = 1)
        plot(truth)
        plot(model)
        par(op)
    }
    mc <- chains(model)
    pmns <- colMeans(theta(mc))
    i <- order(pmns)
    expect_equal(pmns[i], theta(truth), tolerance=0.03)
    ps <- colMeans(sigma(mc))
    expect_equal(ps[i], sigma(truth), tolerance=0.11)
    pmix <- p(truth)
    pp <- colMeans(p(mc))
    expect_equal(pp[i], pmix, tolerance=0.05)
    i <- CNPBayes:::argMax(model)
    expect_true(i == which.max(logPrior(chains(model)) + log_lik(chains(model))))
    expect_identical(CNPBayes:::thetac(model)[i, ], modes(model)[["theta"]])
    expect_identical(sigmac(model)[i, ], sqrt(modes(model)[["sigma2"]]))
})

test_that("test_segfaultExcept", {
    baf <- readRDS(system.file("extdata", "baf.rds", package = "CNPBayes"))
    set.seed(17)
    model <- MarginalModel(baf, k = 2)
    model@.internal.constraint <- -1
    expect_error(model <- posteriorSimulation(model))
})

test_that("test_selectK_easy", {
    set.seed(1)
    means <- c(-1, 0, 1)
    sds <- c(0.1, 0.2, 0.2)
    truth <- simulateData(N = 250, p = rep(1/3, 3), theta = means,
        sds = sds)
    mp <- McmcParams(iter = 5000, burnin = 1000, nStarts = 1)
    model <- MarginalModel(data = y(truth), k = 2, mcmc.params = mp)
    mlist <- posteriorSimulation(model, k = 2:4)
    m.y <- marginalLikelihood(mlist)
    argmax <- which.max(m.y)
    expect_true(argmax == 3L)
    dm <- DensityModel(mlist[[argmax]], merge = TRUE)
    expect_identical(3L, k(dm))
})

test_that("posteriorSimulation methods", {
    set.seed(1)
    mp <- McmcParams(iter=10, burnin=0)
    model <- MarginalModel(data=rnorm(10), k=2, mcmc.params=mp)

    # normal method
    posteriorSimulation(model)

    # different component size
    posteriorSimulation(model, k=3)

    # multiple components
    post <- posteriorSimulation(model, k=2:4)
    expect_identical(length(post), 3L)
})
