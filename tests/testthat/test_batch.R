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
                             p = c(1/10, 1/5, 1 - 0.1 - 0.2),
                             theta = means,
                             sds = sds)
  mcmcp <- McmcParams(iter = 1000, burnin = 0, thin=0,
                      nStarts = 50)
  model <- BatchModel(data = y(truth), batch = batch(truth),
                      k = 3, mcmc.params = mcmcp) ##, hypp = hypp)
  model2 <- posteriorSimulation(model)

  pz <- probz(model2)
  expect_true(all(pz >= 0 | pz <= 1))
  model2 <- useModes(model2)
  expect_equal(theta(truth), theta(model2), tolerance=0.1)
  if (FALSE) {
    plist <- ggMultiBatchChains(model2)
    plist[["batch"]]
    plist3 <- ggMultiBatchChains(model4)
    plist3[["batch"]]
      zz <- as.integer(z(truth))
      ps <- c(mean(zz == 1), mean(zz == 2), mean(zz == 3))
      modelk <- model
      plot.ts(pic(modelk), plot.type = "single")
      abline(h = p(truth))
      par(mfrow = c(1, 3))
      trace(modelk, "theta", col = 1:3)
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
    ##yy <- y(truth)
    ##expect_identical(yy[order(batch(truth))], yy)
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
    ##mcmcp <- McmcParams(iter = 300, burnin = 300, nStarts = 5)
    mcmcp <- McmcParams(iter = 20, burnin = 50, nStarts = 50)
    model <- BatchModel(y(truth), batch = batch(truth), k = 3,
                        mcmc.params = mcmcp)
    model <- posteriorSimulation(model)
    expect_equal(theta(model), theta(truth),
                 scale=0.01, tolerance=1)
    if (FALSE) {
        BatchModelExample <- model
        save(BatchModelExample, file = "data/BatchModelExample.RData")
    }
})


hardTruth <- function(p1, s, N=1000){
  set.seed(1234)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.9, -2, -1.85,
                    -0.45, -0.4, -0.35,
                    -0.1, 0, -0.05), nbatch, k, byrow=FALSE)
  sds <- matrix(s, nbatch, k)
  ##p1 <- prop_comp1
  p2 <- 20*p1
  p3 <- 1-p1-p2
  ##N <- 2000
  truth <- simulateBatchData(N,
                             theta=means,
                             sds=sds,
                             batch=rep(letters[1:3], length.out=N),
                             p=c(p1, p2, p3))
  truth
}

test_that("test_stay_near_truth", {
  set.seed(123)
  library(SummarizedExperiment)
  ## embed function in test for now
  # end function
  truth <- hardTruth(p1=0.02, s = 0.1)
  se <- as(truth, "SummarizedExperiment")
  ##
  ## nStarts=0 forces chain to start at true values
  ##
  ## - these unit tests verify that the model stays in a region of high
  ## - posterior prob.
  mcmcp <- McmcParams(iter = 100, burnin = 0, nStarts=0)
  modelk <- BatchModel(data = y(truth), batch = batch(truth),
                       k = 3, mcmc.params = mcmcp,
                       HyperparametersBatch(k = 3,
                                            m2.0 = 1/60, eta.0 = 1800))
  modelk <- startAtTrueValues(modelk, truth)
  mmodel <- posteriorSimulation(modelk)
  pmns <- CNPBayes:::thetaMean(mmodel)
  ##ps <- CNPBayes:::sigmaMean(mmodel)
  s <- sigma(mmodel)
  pmix <- CNPBayes:::pMean(mmodel)
  expect_equal(theta(truth), pmns, tolerance=0.04)
  expect_equal(p(truth), pmix, tolerance=0.04)
  ## TODO: This example could be extended to focus on what to do when a batch
  ## does not have any homozygous deletions. Batch 2 has only 1 homozygous
  ## deletion
})

test_that("test_hard4", {
  library(SummarizedExperiment)
  set.seed(123)
  truth <- hardTruth(p1=0.02, s = 0.1, N=500)
  se <- as(truth, "SummarizedExperiment")
  mcmcp <- McmcParams(iter = 100, burnin = 100, nStarts = 20)
  modelk <- BatchModel(data = y(truth), batch = batch(truth),
                       k = 3, mcmc.params = mcmcp)
  model2 <- posteriorSimulation(modelk)
  thetas <- theta(model2)
  pmix <- p(model2)
  expect_equal(theta(truth), thetas, tolerance=0.1)
  expect_equal(sigma(truth), sigma(model2), tolerance=0.15)
  expect_equal(p(truth), pmix, tolerance=0.04)
})

test_that("test_kbatch", {
    set.seed(123)
    k <- 3
    means <- matrix(c(rnorm(5, -1, 0.1), rnorm(5, 0, 0.1), rnorm(5,
        1, 0.1)), 5, k, byrow = FALSE)
    sds <- matrix(0.1, 5, k)
    N <- 3000
    ## the last 2 batches are much smaller
    probs <- c(1/3, 1/3, 3/10, 0.02, 0.013)
    probs <- probs/sum(probs)
    batch <- sample(1:5, size = N, prob = probs, replace = TRUE)
    p <- c(1/5, 1/3)
    p <- c(p, 1 - sum(p))
    truth <- simulateBatchData(N = N, batch = batch, theta = means,
        sds = sds, p = p)
    mp <- McmcParams(iter = 100, burnin = 50, nStarts = 10)
    kmod <- BatchModel(y(truth), batch(truth), k = 3, mcmc.params = mp)
    kmod <- posteriorSimulation(kmod)
    cn <- map(kmod)
    set.seed(1000)
    ##ds <- downSampleEachBatch(y(kmod), nt=250, batch=batch(kmod))
    index <- sort(unique(c(sample(seq_len(N), 500), which(batch(kmod) %in% 4:5))))
    ## subsample
    mp <- McmcParams(iter = 100, burnin = 100, nStarts = 20)
    kmod2 <- BatchModel(y(kmod)[index], batch=batch(kmod)[index],
                        k=3, mcmc.params=mp)
    kmod2 <- posteriorSimulation(kmod2)
    yy <- setNames(y(truth), seq_along(y(truth)))
    df <- imputeFromSampledData(kmod2, yy, index)
    cn2 <- df$cn
    expect_true(mean(cn != cn2) < 0.01)
    cn2 <- map(kmod2)
    pz <- probz(kmod2)
    pz <- mapCnProbability(kmod2)
    tab.z <- as.integer(table(z(kmod2)))
    tab.z2 <- colSums(round(pz, 1))
    expect_equal(tab.z, tab.z2)
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

test_that("different starts", {
  set.seed(25)
  bmodel <- BatchModelExample
  bmodel.list <- replicate(10,
                           BatchModel(y(bmodel),
                                      k=3,
                                      batch=batch(bmodel)))
  model <- selectByLogLik(bmodel.list)
  logliks <- unique(sapply(bmodel.list, log_lik))
  expect_identical(length(logliks), length(bmodel.list))
  expect_equal(log_lik(model), -286.84, tolerance=0.05, scale=1)

  mmodel <- MarginalModelExample
  mmodel.list <- replicate(10,
                           MarginalModel(y(mmodel),
                                         k=3))
  logliks <- unique(sapply(mmodel.list, log_lik))
  expect_identical(length(logliks), length(mmodel.list))
  model <- selectByLogLik(mmodel.list)
  expect_equal(log_lik(model), max(logliks, na.rm=TRUE))
})
