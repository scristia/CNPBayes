context("Posterior predictive distribution")

test_that("MBP", {
  set.seed(9)
  mp <- McmcParams(iter=10, burnin=5, nStarts=1, thin=1)
  mb <- MultiBatchModelExample
  mbp.tmp <- as(mb, "MultiBatchPooled")
  expect_identical(length(sigma(mbp.tmp)), 3L)
  ch <- sigma2(chains(mbp.tmp))
  expect_identical(dim(ch), c(iter(mb), 3L))

  hp <- hpList(k=3)[["MBP"]]
  mbp <- MultiBatchPooled(dat=y(mb),
                          hp=hp,
                          mp=mp,
                          batches=batch(mb))
  ch <- sigma2(chains(mbp))
  expect_identical(dim(ch), c(iter(mbp), 3L))
  mbp <- posteriorSimulation(mbp)
  mbp2 <- as(mbp, "MultiBatchP")
  fl <- flags(mbp2)
  expect_true(fl$label_switch)
})
