context("Marginal likelihood calibration")

test_that("galaxy model", {
    # set seed
    set.seed(1)
    # load data
    library(MASS)
    data(galaxies)
    # correct 78th observation
    galaxies[78] <- 26960
    # fit model
    mp <- McmcParams(thin=10, iter=1000, burnin=10000, nStarts=5)
    hypp <- Hyperparameters(type="marginal", k=3, m2.0=100)
    model <- SingleBatchModel(data=galaxies/1000, k=3,
                           hypp=hypp,
                           mcmc.params=mp)
    fit <- posteriorSimulation(model)
    # calculate marginal likelihood and compare to "truth"
    published.mlik <- -226.791
    m.y <- marginalLikelihood(fit, mlParams(root=1))
    m.y <- unname(m.y)
    expect_equal(m.y, published.mlik, tolerance=1, scale=1)
})
