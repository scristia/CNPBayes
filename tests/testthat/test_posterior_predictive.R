context("Posterior predictive distribution")

test_that("SB", {
  library(dplyr)
  library(ggplot2)
  set.seed(149)
  model <- SingleBatchModelExample
  mp <- McmcParams(iter=500, burnin=50)
  mcmcParams(model) <- mp
  model <- posteriorSimulation(model)
  tab <- posteriorPredictive(model) %>%
    mutate(y=round(y, 2))
  expect_equal(tab[1, ], tibble(y=0.96, component=3L))
  if(FALSE){
    ggPredictive(model, tab)
  }
})

test_that("MB", {
  bmodel <- MultiBatchModelExample
  mp <- McmcParams(iter=500, burnin=150)
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
  mp <- McmcParams(iter=1000, burnin=500, nStarts=4, thin=1)
  sb <- SingleBatchModelExample
  hp <- hpList(k=3)[["SBP"]]
  sbp <- SingleBatchPooled(dat=sample(y(sb)),
                           hp=hp,
                           mp=mp)
  expect_equal(sigma(sbp), 0.1, tolerance=0.05)
  sbp2 <- posteriorSimulation(sbp)
  tab <- posteriorPredictive(sbp2)
  if(FALSE) ggPredictive(sbp2, tab)
})

.test_that <- function(nm, expr) NULL

test_that("MBP", {
  set.seed(9)
  mp <- McmcParams(iter=1000, burnin=250, nStarts=4, thin=2)
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
  tab <- posteriorPredictive(mbp)
  if(FALSE) ggPredictive(mbp, tab)
})
