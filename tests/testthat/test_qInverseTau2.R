context("Check qInverseTau2 function")

test_that("test_qInverseTau2", {
    mn <- 100
    sd <- 10
    shape.rate <- CNPBayes:::.parameterizeGammaByMeanSd(mn = mn,
        sd = sd)
    a <- as.numeric(shape.rate[[1]])
    b <- as.numeric(shape.rate[[2]])
    expect_equal(mn^2/sd^2, a)
    expect_equal(mn/sd^2, b)
    tmp <- rgamma(10000, a, rate = b)
    expect_equal(mn, mean(tmp), tolerance=1)
    expect_equal(sd, sd(tmp), tolerance=0.5)
    eta.0 <- 2 * a
    m2.0 <- b/a
    x <- qgamma(seq(0, 1 - 0.001, 0.001), 0.5 * eta.0, rate = 0.5 *
                    eta.0 * m2.0)
    x <- x[is.finite(x) & x > 0]
    prec <- qInverseTau2(mn = mn, sd = sd)
    x2 <- prec$quantiles
    expect_identical(x2, x)
    if(FALSE){
      hist(sqrt(1/x2), breaks=100)

      eta.0 <- 1800
      m2.0 <- 1/60
      prec <- qInverseTau2(eta.0=1800, m2.0=1/60)
      prec <- prec$quantiles
      sds <- sqrt(1/prec)
      hist(sds, breaks=100)

      mn.prec <- (0.05)^(-2)
      hist(sqrt(1/qInverseTau2(mn=mn.prec, sd=10)$quantiles), breaks=100)
      tmp <- qInverseTau2(mn=mn.prec, sd=10)

      library(oligoClasses)
      set.seed(100)
      nbatch <- 3
      k <- 3
      means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
          0, 0.05), nbatch, k, byrow = FALSE)
      sds <- matrix(0.15, nbatch, k)
      sds[, 1] <- 0.3
      N <- 1000
      truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
          length.out = N), p = c(1/10, 1/5, 1 - 0.1 - 0.2), theta = means,
          sds = sds)
      mcmcp <- McmcParams(iter = 1000, burnin = 500, thin = 1,
          nStarts = 10)
      hypp <- CNPBayes:::HyperparametersBatch(m2.0 = 1/60, eta.0 = 1800,
          k = 3, a = 1/6, b = 180)
      model <- BatchModel(data = y(truth), batch = batch(truth),
          k = 3, mcmc.params = mcmcp, hypp = hypp)
      model <- posteriorSimulation(model)
    }
})
