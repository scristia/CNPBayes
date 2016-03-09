context("Check z's")

test_that("test_probz", {
    set.seed(1)
    truth <- simulateData(N = 2500, p = rep(1/3, 3), theta = c(-1, 
        0, 1), sds = rep(0.1, 3))
    mp <- McmcParams(iter = 500, burnin = 500)
    set.seed(123)
    model <- MarginalModel(data = y(truth), k = 3, mcmc.params = mp)
    model <- CNPBayes:::startAtTrueValues(model, truth)
    true_z <- z(truth)
    expect_equal(z(model), z(truth))
    pz <- probz(model)
    expect_true(all(pz == 0))
    iter(model, force = TRUE) <- 2
    burnin(model) <- 0
    zz <- map(posteriorSimulation(model))
    expect_equal(true_z, zz)
    model2 <- CNPBayes:::modelOtherModes(model, maxperm = 2)[[2]]
    z2 <- zz <- z(model2)
    expect_true(sum(zz != true_z) > 500)
    mcmcParams(model2) <- mcmcParams(model)
    model2 <- posteriorSimulation(model2)
    zz <- map(model2)
    expect_equal(true_z, zz)
    model3 <- CNPBayes:::modelOtherModes(model, maxperm = 3)[[3]]
    z3 <- z(model3)
    expect_true(sum(z3 != true_z) > 500)
    expect_true(sum(z2 != z3) > 500)
    model3 <- posteriorSimulation(model3)
    mz3 <- map(model3)
    table(mz3, true_z)
    expect_equal(true_z, mz3)
})

