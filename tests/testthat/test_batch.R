context("BatchModel")

test_that("test_batch_moderate", {
    set.seed(100)
    nbatch <- 3
    k <- 3
    means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
        0, 0.05), nbatch, k, byrow = FALSE)
    sds <- matrix(0.15, nbatch, k)
    sds[, 1] <- 0.3
    N <- 1000
    truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
                                                  length.out = N),
                               p = c(1/10, 1/5, 1 - 0.1 - 0.2), theta = means,
                               sds = sds)
    mcmcp <- McmcParams(iter = 1000, burnin = 500, thin = 1,
                        nStarts = 10)

    ## this parameter setting for m2.0 allows a lot of varation of the thetas
    ## between batch
    hypp <- CNPBayes:::HyperparametersBatch(m2.0 = 1/60, eta.0 = 1800,
                                            k = 3, a = 1/6, b = 180)
    model <- BatchModel(data = y(truth), batch = batch(truth),
                        k = 3, mcmc.params = mcmcp, hypp = hypp)
    model <- posteriorSimulation(model)
    i <- order(theta(model)[1, ])
    expect_equal(theta(truth), theta(model)[, i], tolerance=0.1)
    if (FALSE) {
        zz <- as.integer(z(truth))
        ps <- c(mean(zz == 1), mean(zz == 2), mean(zz == 3))
        modelk <- model
        plot.ts(pic(modelk), plot.type = "single")
        abline(h = p(truth))
        par(mfrow = c(1, 3))
        tracePlot(modelk, "theta", col = 1:3)
        abline(h = theta(truth))
        plot.ts(sigmac(modelk), plot.type = "single")
        abline(h = sigma(truth))
        op <- par(mfrow = c(1, 2), las = 1)
        CNPBayes::plot(truth)
        CNPBayes::plot(modelk)
        par(op)
        mc <- chains(modelk)
        plot.ts(sigma(mc), col = "gray")
        plot.ts(theta(mc), col = "gray")
        plot.ts(p(mc), col = "gray")
    }
})

test_that("test_batchEasy", {
    set.seed(123)
    k <- 3
    nbatch <- 3
    means <- matrix(c(-1.2, -1, -0.8, -0.2, 0, 0.2, 0.8, 1, 1.2),
        nbatch, k, byrow = FALSE)
    sds <- matrix(0.1, nbatch, k)
    N <- 1500
    truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
        length.out = N), theta = means, sds = sds, p = c(1/5,
        1/3, 1 - 1/3 - 1/5))
    yy <- y(truth)
    expect_identical(yy[order(batch(truth))], yy)
    mcmcp <- McmcParams(iter = 50, burnin = 0)
    set.seed(123)
    model <- BatchModel(y(truth), batch = batch(truth), k = 3,
        mcmc.params = mcmcp)
    model <- CNPBayes:::startAtTrueValues(model, truth)
    expect_identical(batch(truth), batch(model))
    expect_identical(y(truth), y(model))
    expect_identical(theta(truth), theta(model))
    expect_identical(sigma(truth), sigma(model))
    expect_identical(p(truth), p(model))
    expect_identical(z(truth), z(model))
    iter(model) <- 50
    if (FALSE) {
        model <- .Call("mcmc_batch", model, mcmcParams(model))
        set.seed(123)
        .Call("update_sigma20_batch", model)
        set.seed(123)
        .updateSigma2.0Batch(model2)
        eta.0 <- 1800/100
        m2.0 <- 1/(60/100)
        a <- 0.5 * eta.0
        b <- 0.5 * eta.0 * m2.0
        x <- rgamma(1000, a, rate = 1/b)
        ix <- 1/x
        hist(sqrt(1/x), breaks = 100)
        s <- mean(sqrt(1/x))
        theta <- rnorm(1000, -1, s)
        hist(theta, breaks = 100, xlim = c(-2, 1.5))
    }
    set.seed(1)
    mcmcp <- McmcParams(iter = 300, burnin = 300, nStarts = 5)
    model <- BatchModel(y(truth), batch = batch(truth), k = 3,
        mcmc.params = mcmcp)
    model <- posteriorSimulation(model)
    if (FALSE) {
        BatchModelExample <- model
        save(BatchModelExample, file = "data/BatchModelExample.RData")
    }
    i <- order(theta(model)[1, ])
    expect_equal(theta(truth), theta(model)[, i], tolerance=0.1)
    if (FALSE) {
        op <- par(mfrow = c(1, 2), las = 1)
        plot(truth)
        plot(model)
        par(op)
        op <- par(mfrow = c(1, 3))
        tracePlot(model, "theta", col = 1:3)
        par(op)
        plot.ts(muc(model), plot.type = "single", col = 1:3)
    }
})

test_that("test_hard3", {
    # embed function in test for now
    hardTruth <- function(prop_comp1=0.005, s=0.3) {
      set.seed(1234)
      k <- 3
      nbatch <- 3
      means <- matrix(c(-1.9, -2, -1.85,
                        -0.45, -0.4, -0.35,
                        -0.1, 0, -0.05), nbatch, k, byrow=FALSE)
      sds <- matrix(s, nbatch, k)
      p1 <- prop_comp1
      p2 <- 20*p1
      p3 <- 1-p1-p2
      N <- 8e3
      truth <- simulateBatchData(N,
                                 theta=means,
                                 sds=sds,
                                 batch=rep(letters[1:3], length.out=N),
                                 p=c(p1, p2, p3))
    }
    # end function
    truth <- hardTruth(0.005, s = 0.1)
    se <- as(truth, "SummarizedExperiment")
    if (FALSE)
        hist(oligoClasses::copyNumber(se), breaks = 1000, col = "gray",
            border = "gray")
    mcmcp <- McmcParams(iter = 100, burnin = 0)
    modelk <- BatchModel(data = y(truth), batch = batch(truth),
        k = 3, mcmc.params = mcmcp, CNPBayes:::HyperparametersBatch(k = 3,
            m2.0 = 1/60, eta.0 = 1800))
    modelk <- CNPBayes:::startAtTrueValues(modelk, truth)
    mmodel <- posteriorSimulation(modelk)
    if (FALSE) {
        op <- par(mfrow = c(1, 2), las = 1)
        CNPBayes::plot(truth)
        .plotBatch(mmodel)
        par(op)
        i <- order(theta(mmodel)[1, ])
        theta(mmodel)[, i]
        theta(truth)
        sigma(model)[, i]
        sigma(truth)
        mc <- chains(mmodel)
        par(mfrow = c(1, 3))
        tracePlot(mmodel, "sigma", col = 1:3)
        tracePlot(mmodel, "theta", col = 1:3)
        plot.ts(sigma(mc), col = "gray", )
        plot.ts(theta(mc), col = "gray")
        plot.ts(mu(mc), plot.type = "single", col = 1:3)
        cbind(colMeans(sigma(mc)), as.numeric(sigma(truth)))
    }
    pmns <- CNPBayes:::thetaMean(mmodel)
    j <- order(pmns[1, ])
    ps <- CNPBayes:::sigmaMean(mmodel)[, j]
    pmix <- CNPBayes:::pMean(mmodel)[j]
    expect_equal(theta(truth), pmns[, j], tolerance=0.04)
    expect_equal(sigma(truth), ps, tolerance=0.15)
    expect_equal(p(truth), pmix, tolerance=0.04)
    mcmcp <- McmcParams(iter = 200, burnin = 100, nStarts = 20)
    modelk <- BatchModel(data = y(truth), batch = batch(truth),
        k = 3, mcmc.params = mcmcp, CNPBayes:::HyperparametersBatch(k = 3,
            m2.0 = 1/60, eta.0 = 1800, tau2.0 = 1000))
    mmodel <- posteriorSimulation(modelk)
    pmns <- CNPBayes:::thetaMean(mmodel)
    j <- order(pmns[1, ])
    ps <- CNPBayes:::sigmaMean(mmodel)[, j]
    pmix <- CNPBayes:::pMean(mmodel)[j]
    expect_equal(theta(truth), pmns[, j], tolerance=0.04)
    expect_equal(sigma(truth), ps, tolerance=0.15)
    expect_equal(p(truth), pmix, tolerance=0.04)
})

test_that("test_kbatch", {
    set.seed(123)
    k <- 3
    means <- matrix(c(rnorm(5, -1, 0.1), rnorm(5, 0, 0.1), rnorm(5,
        1, 0.1)), 5, k, byrow = FALSE)
    sds <- matrix(0.1, 5, k)
    N <- 3000
    probs <- c(1/3, 1/3, 3/10, 0.02, 0.013)
    probs <- probs/sum(probs)
    batch <- sample(1:5, size = N, prob = probs, replace = TRUE)
    p <- c(1/5, 1/3)
    p <- c(p, 1 - sum(p))
    truth <- simulateBatchData(N = N, batch = batch, theta = means,
        sds = sds, p = p)
    mp <- McmcParams(iter = 1000, burnin = 250, nStarts = 5)
    kmod <- BatchModel(y(truth), batch(truth), k = 3, mcmc.params = mp)
    kmod <- posteriorSimulation(kmod)
    cn <- map(kmod)
    set.seed(1000)
    index <- sample(seq_len(N), 1000)
    kmod2 <- BatchModel(y(truth)[index], batch(truth)[index],
        k = 3, mcmc.params = mp)
    kmod2 <- posteriorSimulation(kmod2)
    yy <- setNames(y(truth), seq_along(y(truth)))
    df <- CNPBayes:::imputeFromSampledData(kmod2, yy, index)
    cn2 <- df$cn
    mean(cn != cn2)
    cn2 <- map(kmod2)
    pz <- probz(kmod2)
    pz <- mapCnProbability(kmod2)
    if (FALSE) {
        fit <- list(posteriorSimulation(kmod, k = 1), posteriorSimulation(kmod,
            k = 2), posteriorSimulation(kmod, k = 3), posteriorSimulation(kmod,
            k = 4))
        fit <- marginalLikelihood(fit)
        prz <- probz(fit$models[[4]])
        cn <- map(fit$models[[4]])
        plot(r, cn, pch = 20, cex = 0.3)
        trace(cnProbability, browser)
        prz <- cnProbability(prz, 4)
        plot(jitter(prz, amount = 0.05), jitter(cn, amount = 0.05),
            pch = 20, cex = 0.3)
        table(cn)
        pz <- cnProbability(probz(fit$models[[4]]), 4)
        r <- y(fit$models[[4]])
        plot(r, pz, pch = ".")
        expect_true(k(orderModels(fit))[1] == 3)
    }
})

test_that("test_unequal_batch_data", {
    expect_error(BatchModel(data = 1:10, batch = 1:9))
})
