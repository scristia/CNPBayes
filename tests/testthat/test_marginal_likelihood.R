test_that("overfit model", {
    set.seed(1)

    # load data
    library(MASS)
    data(galaxies)

    # correct 78th observation
    galaxies[78] <- 26960

    mp <- McmcParams(thin=1, iter=1000, burnin=1000, nStarts=1)
    hypp <- Hyperparameters(type="marginal")
    model <- MarginalModel(data=galaxies / 1000,
                           hypp=hypp,
                           mcmc.params=mp)

    fit <- posteriorSimulation(model, k=1:4)

    expect_warning(marginalLikelihood(fit))
})
