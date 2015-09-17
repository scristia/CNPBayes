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
    expect_equal(mn, mean(tmp))
    expect_equal(sd, sd(tmp))
    eta.0 <- 2 * a
    m2.0 <- b/a
    x <- qgamma(seq(0, 1 - 0.001, 0.001), 0.5 * eta.0, rate = 0.5 * 
        eta.0 * m2.0)
    prec <- qInverseTau2(mn = mn, sd = sd)
    x2 <- prec$quantiles
    expect_identical(x2, x)
})

