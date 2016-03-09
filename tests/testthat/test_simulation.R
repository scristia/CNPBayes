context("Simulation")

test_that("test_simulation_moderate", {
    set.seed(42)
    arguments <- list(sl.good = 6.25, sl.bad = 0.0625, prbias = 0.03, 
        n = 0.2, prvar = c(19.92985, 0.06272))
    dat <- CNPBayes:::simulateProbeLevel(cnvs = 1, K = 4, probes = 10, 
        arguments = arguments, qual = "medium")
    x <- dat[[1]]
    if (FALSE) 
        hist(rowMeans(x[, , 1, 3]), col = "gray", breaks = 80)
    K <- 4
    xx <- x[, , 1, K]
    mns <- rowMeans(xx)
    pc <- prcomp(xx, center = TRUE, scale. = TRUE)$x[, 1]
    if (cor(pc, mns) < cor(-pc, mns)) 
        pc <- -pc
    mp <- McmcParams(iter = 1000, burnin = 500, nStarts = 1)
    model <- MarginalModel(data = pc, k = 1, mcmc.params = mp)
    mlist <- posteriorSimulation(model, k=1:4)
    m.y <- marginalLikelihood(mlist)
    expect_true(which.max(m.y) == 4L)
    if (FALSE) 
        hist(pc, breaks = 100, col = "gray", border = "gray")
})

