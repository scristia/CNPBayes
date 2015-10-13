test_that("test_cnProbability", {
    set.seed(1)
    truth <- simulateData(N = 2500, p = rep(1/3, 3), theta = c(-1, 
        0, 1), sds = rep(0.1, 3))
    if (FALSE) 
        plot(truth)
    mp <- McmcParams(iter = 200, burnin = 100, nStarts = 1)
    model <- MarginalModel(data = y(truth), k = 3, mcmc.params = mp)
    model <- posteriorSimulation(model)
    map_model <- CNPBayes:::mapModel(model)
    expect_identical(z(map_model), CNPBayes:::zChain(model)[CNPBayes:::argMax(model), 
        ])
    expect_identical(theta(map_model), CNPBayes:::thetac(model)[CNPBayes:::argMax(model), 
        ])
    probs <- mapCnProbability(model)
    pz <- probz(model)
    expect_equal(head(pz), head(probs), tolerance=0.01)
})

