test_that("test_mcmc_restart", {
    set.seed(1)
    truth <- simulateData(N = 500, theta = c(-1, 0, 1), sds = c(0.1,
        0.1, 0.1), p = c(1/5, 1/6, 1 - 1/5 - 1/6))
    model <- truth
    mcmcParams(model, force = TRUE) <- McmcParams(iter = 1, burnin = 0)
    model <- posteriorSimulation(model)
    mc <- chains(model)
    expect_identical(theta(truth), theta(mc)[1, ])
    expect_identical(sigma(truth), sigma(mc)[1, ])
    mcmcParams(model) <- McmcParams(iter = 1, burnin = 10)
    model <- posteriorSimulation(model)
    model <- truth
    mcmcParams(model, force = TRUE) <- McmcParams(iter = 1, burnin = 1)
    model <- posteriorSimulation(model)
    mc <- chains(model)
    expect_false(identical(theta(mc)[1, ], theta(truth)))
    expect_false(identical(sigma(mc)[1, ], sigma(truth)))
    expect_equal(as.numeric(theta(model)), theta(mc)[1, ])
    model <- truth
    mcmcParams(model, force = TRUE) <- McmcParams(iter = 0, burnin = 1)
    model <- posteriorSimulation(model)
    mc <- chains(model)
    expect_identical(0L, nrow(theta(mc)))
    model <- truth
    mcmcParams(model, force = TRUE) <- McmcParams(iter = 1, burnin = 0)
    model <- posteriorSimulation(model)
    mcmcParams(model, force = TRUE) <- McmcParams(iter = 10,
        burnin = 0)
    model <- posteriorSimulation(model)
    mc <- chains(model)
    expect_identical(10L, nrow(theta(mc)))
    expect_identical(as.numeric(theta(truth)), theta(mc)[1, ])
    th <- theta(model)
    s <- sigma(model)
    model2 <- posteriorSimulation(model)
    mc <- chains(model2)
    expect_identical(as.numeric(th), theta(mc)[1, ])
    expect_identical(as.numeric(s), sigma(mc)[1, ])
    set.seed(123)
    k <- 3
    nbatch <- 3
    means <- matrix(c(-1.2, -1, -0.8, -0.2, 0, 0.2, 0.8, 1, 1.2),
        nbatch, k, byrow = FALSE)
    sds <- matrix(0.1, nbatch, k)
    truth <- simulateBatchData(N = 2500, batch = rep(letters[1:3],
        length.out = 2500), p = c(1/4, 1/6, 1 - 1/4 - 1/6), theta = means,
        sds = sds)
    mcmcp <- McmcParams(iter = 5, burnin = 0)
    set.seed(123)
    modelk1 <- BatchModel(data = y(truth), k = 3, mcmc.params = mcmcp,
        batch = batch(truth))
    th1 <- as.numeric(theta(modelk1))
    modelk <- posteriorSimulation(modelk1)
    th2 <- CNPBayes:::thetac(modelk)[1, ]
    expect_identical(th2, th1)
})
