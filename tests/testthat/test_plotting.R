test_that("test_plot", {
    set.seed(1)
    truth <- simulateData(N = 2500, p = rep(1/3, 3), theta = c(-1, 
        0, 1), sds = rep(0.1, 3))
    if (FALSE) 
        plot(truth)
    mp <- McmcParams(iter = 500, burnin = 500)
    model <- MarginalModel(data = y(truth), k = 3, mcmc.params = mp)
    model <- posteriorSimulation(model)
    op <- par(mfrow = c(1, 2), las = 1)
    fig1 <- plot(truth)
    fit2 <- plot(model)
    par(op)
})

