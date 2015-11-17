context("Model Selection")

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
    model <- MarginalModel(data=galaxies/1000, k=3,
                           hypp=hypp,
                           mcmc.params=mp)
    fit <- posteriorSimulation(model, k=3:4)

    # calculate marginal likelihood and compare to "truth"
    published.mlik <- -226.791  
    expect_warning(m.y <- marginalLikelihood(fit,
                                             params=list(niter=1000,
                                                         root=1,
                                                         reject.threshold=1e-50,
                                                         prop.threshold=0.5)))
    marginal_k3 <- unname(m.y[1])
    expect_equal(marginal_k3, published.mlik, tolerance=1, scale=1)

    # check that model is not overfit
    expect_true(which.max(m.y) == 1L)
})
