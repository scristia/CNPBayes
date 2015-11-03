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

test_that("batch overfit galaxy", {
    set.seed(1)

    # load data
    library(MASS)
    data(galaxies)

    # correct 78th observation
    galaxies[78] <- 26960

    mp <- McmcParams(thin=10, iter=1000, burnin=5000, nStarts=1)
    hypp <- Hyperparameters(type="batch")
    model <- BatchModel(data=c(galaxies / 1000, 
                               galaxies / 1000 + 50),
                        batch=rep(1:2, each=length(galaxies)),
                        mcmc.params=mp)

    fit <- posteriorSimulation(model, k=1:4)
    es <- lapply(fit, function(x) effectiveSize(theta(chains(x))))
    means <- sapply(es, function(x) mean(x))

    marginalLikelihood(fit)
})
