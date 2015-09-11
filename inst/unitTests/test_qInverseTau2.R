test_qInverseTau2 <- function(){
  mn <- 100
  sd <- 10
  shape.rate <- CNPBayes:::.parameterizeGammaByMeanSd(mn=mn, sd=sd)
  a <- as.numeric(shape.rate[[1]])
  b <- as.numeric(shape.rate[[2]])
  checkEquals(a, mn^2/sd^2)
  checkEquals(b, mn/sd^2)
  tmp <- rgamma(10e3, a, rate=b)
  checkEquals(mean(tmp), mn, tolerance=1)
  checkEquals(sd(tmp), sd, tolerance=0.5)

  eta.0 <- 2*a
  m2.0 <- b/a
  x <- qgamma(seq(0, 1-0.001, 0.001), 0.5*eta.0, rate=0.5*eta.0*m2.0)
  prec <- qInverseTau2(mn=mn, sd=sd)
  x2 <- prec$quantiles
  checkIdentical(x, x2)
}
