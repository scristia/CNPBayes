context("Posterior predictive distribution")

test_that("SB", {
  library(dplyr)
  library(ggplot2)
  set.seed(149)
  tmp <- rst(1000, mean=0, sigma=1)
  expect_equal(mean(tmp), 0, scale=1, tolerance=0.05)
  expect_equal(sd(tmp), 1, scale=1, tolerance=0.05)
  model <- SingleBatchModelExample
  mp <- McmcParams(iter=50, burnin=10)
  mcmcParams(model) <- mp
  model <- posteriorSimulation(model)
  tab <- posteriorPredictive(model)
  expect_is(tab, "tbl_df")
  if(FALSE){
    ggPredictive(model, tab)
  }
})

test_that("MB", {
  bmodel <- MultiBatchModelExample
  mp <- McmcParams(iter=10, burnin=5)
  mcmcParams(bmodel) <- mp
  bmodel <- posteriorSimulation(bmodel)
  tab <- posteriorPredictive(bmodel)
  expect_is(tab, "tbl_df")
  expect_identical(colnames(tab), c("y", "batch", "component"))
  if(FALSE){
    ggPredictive(bmodel, tab)
  }
})

test_that("SBP", {
  set.seed(134)
  mp <- McmcParams(iter=10, burnin=5, nStarts=1, thin=1)
  sb <- SingleBatchModelExample
  hp <- hpList(k=3)[["SBP"]]
  sbp <- SingleBatchPooled(dat=sample(y(sb)),
                           hp=hp,
                           mp=mp)
  sbp2 <- posteriorSimulation(sbp)
  expect_warning(tab <- posteriorPredictive(sbp2))
  expect_is(tab, "tbl_df")
  if(FALSE) ggPredictive(sbp2, tab)
})

.test_that <- function(nm, expr) NULL

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
  expect_warning(mbp <- posteriorSimulation(mbp))
  tab <- posteriorPredictive(mbp)
  expect_is(tab, "tbl_df")
  if(FALSE) ggPredictive(mbp, tab)
})
