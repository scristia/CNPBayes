context("Posterior predictive distribution")

test_that("posteriorPredictive", {
  library(dplyr)
  library(ggplot2)
  set.seed(149)
  model <- SingleBatchModelExample
  mp <- McmcParams(iter=500, burnin=50)
  mcmcParams(model) <- mp
  model <- posteriorSimulation(model)
  tab <- posteriorPredictive(model)
  if(FALSE){
    ggPredictive(model, tab)
  }

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
