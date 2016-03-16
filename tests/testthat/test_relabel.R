context("Relabel components")

test_that("test_relabel_marginal", {
    set.seed(1)
    truth <- simulateData(N = 2500, p = rep(1/3, 3), theta = c(-1, 
        0, 1), sds = rep(0.1, 3))
    mp <- McmcParams(iter = 500, burnin = 500)
    model <- MarginalModel(data = y(truth), k = 3, mcmc.params = mp)
    model <- posteriorSimulation(model)
    zindex <- c(3, 2, 1)
    model2 <- CNPBayes:::relabel(model, zindex)
    expect_identical(as.integer(CNPBayes:::tablez(model2)[zindex]), 
        as.integer(CNPBayes:::tablez(model)))
    expect_identical(CNPBayes:::dataMean(model2)[zindex], CNPBayes:::dataMean(model))
    expect_identical(theta(model2), theta(model))
    expect_identical(mu(model2), mu(model))
    expect_identical(p(model2), p(model))
    expect_identical(sigma(model2), sigma(model))
    expect_identical(nu.0(model2), nu.0(model))
    if (FALSE) {
        burnin(model2) <- 0
        model2 <- posteriorSimulation(model2)
        head(thetac(model2))
        plot.ts(thetac(model2), plot.type = "single", col = 1:3)
        set.seed(1)
        truth <- simulateData(N = 2500, p = rep(1/3, 3), theta = c(-1, 
            0, 1), sds = rep(0.5, 3))
        mp <- McmcParams(iter = 500, burnin = 500)
        model <- MarginalModel(data = y(truth), k = 3, mcmc.params = mp)
        model <- posteriorSimulation(model)
        zindex <- c(3, 2, 1)
        model2 <- CNPBayes:::relabel(model, zindex)
        burnin(model2) <- 0
        model2 <- posteriorSimulation(model2)
        par(mfrow = c(1, 2))
        plot.ts(thetac(model), col = 1:3, plot.type = "single")
        plot.ts(thetac(model2), col = 1:3, plot.type = "single")
    }
})

