test_that("galaxy model", {
    # set seed
    set.seed(1)

    # load data
    library(MASS)
    data(galaxies)

    # correct 78th observation
    galaxies[78] <- 26960

    # fit model
    mp <- McmcParams(thin=10, iter=1000, burnin=10000, nStarts=20)
    hypp <- Hyperparameters(type="marginal", k=3, m2.0=100)
    model <- MarginalModel(data=galaxies/1000, k=3,
                           hypp=hypp,
                           mcmc.params=mp)
    fit <- posteriorSimulation(model, k=3:4)

    # calculate marginal likelihood and compare to "truth"
    published.mlik <- -226.791  
    m.y <- marginalLikelihood(fit, 1000)
    marginal_k3 <- unname(m.y[1])
    expect_equal(marginal_k3, published.mlik, tolerance=1)

    # check that model is not overfit
#     expect_true(m.y[1] > m.y[2])
})
