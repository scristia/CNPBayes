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
  mp <- McmcParams(iter=500, burnin=150, nStarts=20)
  mcmcParams(bmodel) <- mp
  bmodel <- posteriorSimulation(bmodel)a
  tab <- posteriorPredictive(bmodel)
  ##ggplot(tab, aes(y, fill=factor(component))) +
  ##geom_density(adjust=1/5)
  ##
  ##ggplot(filter(tab, component==1), aes(y)) + geom_density(adjust=1/5)
  ##ggplot(filter(tab, component==2), aes(y)) + geom_density(adjust=1/5)
  expect_is(tab, "tibble")
  expect_identical(colnames(tab), c("y", "batch", "component"))
  if(FALSE){
    ggPredictive(bmodel, tab)
  }
})
