context("K too big")

test_that("K too big", {
  set.seed(1337)
  mp <- McmcParams(nStarts=100, burnin=0, iter=0)
  file <- system.file("extdata", "DeletionModelExample.rds", package="CNPBayes")
  sb <- readRDS(file)
  mcmcParams(sb) <- mp
  ## expect no errors even though model is over-parameterized
  sb <- posteriorSimulation(sb)
  expect_identical(k(sb), 4L)

  mb.list <- BatchModelList(data=y(sb), k=1:4,
                            batch=c(rep(1, 12), rep(2, 23)),
                            mcmc.params=mp)
  ## expect no errors even though model is over-parameterized
  mb.list <- posteriorSimulation(mb.list)
})
