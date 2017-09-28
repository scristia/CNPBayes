context("Constructors")

test_that("test_constructor", {
    mcmc.params <- McmcParams()
    expect_true(validObject(mcmc.params))
    hypp <- Hyperparameters()
    expect_true(validObject(hypp))
    hypp <- Hyperparameters("batch")
    expect_true(validObject(hypp))
    ##expect_warning(mmod <- SingleBatchModel())
    mmod <- SingleBatchModel2()
    expect_true(validObject(mmod))
    mc <- McmcChains()
    expect_true(validObject(mc))
    hp <- Hyperparameters()
    expect_true(validObject(hp))
    expect_true(validObject(HyperparametersMultiBatch()))
    truth <- simulateData(N = 2500, p = rep(1/3, 3), theta = c(-1,
        0, 1), sds = rep(0.1, 3))
    expect_true(validObject(truth))
    mmod <- SingleBatchModel2(dat=rnorm(100))
    mmod <- SingleBatchModel2(dat=y(truth))
    expect_true(validObject(mmod))
    iter(mmod) <- 1000L
    expect_error(iter(mmod) <- 1001)
    iter(mmod, force = TRUE) <- 1001L
    burnin(mmod) <- 1001L
    expect_true(iter(mmod) == 1001L)
    expect_true(nrow(thetac(mmod)) == 1001)
    mp <- McmcParams(iter = 10L, thin = 10L, burnin = 100L)
    mcmcParams(mmod) <- mp
    expect_true(nrow(thetac(mmod)) == 10)
    iter(mmod, force = TRUE) <- 1000L
    expect_true(nrow(thetac(mmod)) == 1000)
    bmod <- MultiBatchModel2()
    expect_true(validObject(bmod))
    k <- nbatch <- 3
    means <- matrix(c(-1.2, -1, -0.8, -0.2, 0, 0.2, 0.8, 1, 1.2),
        nbatch, k, byrow = FALSE)
    sds <- matrix(0.1, nbatch, k)
    truth <- simulateBatchData(N = 2500, batch = rep(letters[1:3],
        length.out = 2500), theta = means, sds = sds, p = c(1/5,
        1/3, 1 - 1/3 - 1/5))
    expect_true(validObject(truth))
    bmod <- MultiBatchModel2(dat=y(truth), hp=hpList(k=3)[["MB"]],
                             batches=batch(truth))
    expect_true(is(bmod, "MultiBatchModel"))
    bmod <- MultiBatchModel2(dat=y(truth), batches=batch(truth))
    expect_true(is(bmod, "MultiBatchModel"))
    expect_equal(y(bmod), y(truth))
    expect_equal(batch(bmod), batch(truth))
    expect_error(MultiBatchModel2(dat=y(truth),
                                  batches=rep(1:3, each=2), k=3))
    iter(bmod, force = TRUE) <- 10L
    expect_true(nrow(thetac(bmod)) == 10)
    iter(mmod, force = TRUE) <- 5L
    burnin(mmod) <- 0L
    mmod2 <- posteriorSimulation(mmod)
    expect_false(identical(thetac(mmod), thetac(mmod2)))
    expect_true(all(is.na(thetac(mmod))))
    burnin(mmod2) <- 2L
    expect_false(all(is.na(thetac(mmod2))))
})

test_that("construct pooled model", {
  set.seed(123)
  truth <- simulateData(N = 2500, p = rep(1/3, 3),
                        theta = c(-1, 0, 1), sds = rep(0.1, 3))
  hp <- Hyperparameters(k=5,
                        mu=-0.75,
                        tau2.0=0.4,
                        eta.0=32,
                        m2.0=0.5)
  pooled.model <- SingleBatchPooled(y(truth), hp=hp)
  expect_true(validObject(pooled.model))
})
