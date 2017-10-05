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
    hp <- Hyperparameters(type="marginal", k=3, m2.0=100)
    model <- gibbs_K(mp=mp, hp=hp, k_range=c(3, 3), dat=galaxies/1000)[[1]]
    ml <- marginal_lik(model)
    # calculate marginal likelihood and compare to "truth"
    published.mlik <- -226.791
    m.y <- unname(ml)
    expect_equal(m.y, published.mlik, tolerance=3, scale=1)
})
