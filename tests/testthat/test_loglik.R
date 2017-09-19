context("Log likelihood")

test_that("test_loglik", {
    set.seed(1337)
    truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0), sds = c(0.3,
        0.15, 0.15), p = c(0.005, 1/10, 1 - 0.005 - 1/10))
    ll.truth <- log_lik(truth)
    expect_equal(computeLoglik(truth), ll.truth)
    yy <- y(truth)
    th <- theta(truth)
    sd <- sigma(truth)
    p_ <- p(truth)
    loglik <- sum(log(p_[1] * dnorm(yy, th[1], sd[1]) + p_[2] *
        dnorm(yy, th[2], sd[2]) + p_[3] * dnorm(yy, th[3], sd[3])))
    expect_equal(ll.truth, loglik)
    hp <- HyperparametersMultiBatch(k=3,
                                    mu=-0.75,
                                    tau2.0=0.4,
                                    eta.0=32,
                                    m2.0=0.5)
    mp <- McmcParams(iter=250, burnin=500, nStarts=10)
    model <- MultiBatchModel2(y(truth), hp, mp, 
                              batches=rep(1:3,
                                          length.out = length(y(truth))))
    ##    model <- MultiBatchModel(data = y(truth),
    ##                             batch = rep(letters[1:3],
    ##                                         length.out = length(y(truth))), k = 3)
    mcmcParams(model) <- mp
    model <- posteriorSimulation(model)
    ll1 <- log_lik(model)
    ll2 <- computeLoglik(model)
    expect_equal(ll2, ll1, tolerance=0.05)
    yy <- y(model)
    th <- theta(model)
    sd <- sigma(model)
    p_ <- table(model@batch, z(model))
    p_ <- p_/rowSums(p_)
    ub <- unique(model@batch)
    loglik <- rep(0, length(y(model)))
    for (b in ub) {
        this_batch <- model@batch == ub[b]
        loglik <- loglik + log(p_[b, 1] * dnorm(yy, th[b, 1],
            sd[b, 1]) + p_[b, 2] * dnorm(yy, th[b, 2], sd[b,
            2]) + p_[b, 3] * dnorm(yy, th[b, 3], sd[b, 3])) *
            this_batch
    }
    ll3 <- sum(loglik)
    expect_equal(ll3, ll2)
})
