context("gg functions")

test_that("ggfun", {
  set.seed(123)
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.2, 0.2, 0.2))
  mcmcParams(truth) <- McmcParams(iter=500, burnin=200, thin=1, nStarts=1)
  expect_warning(model <- posteriorSimulation(truth))

  pred <- predictive(model)
  pred2 <- longFormatKB(pred, K=k(model), B=numBatch(model)) %>%
    set_colnames(c("s", "y", "b", "component")) %>%
    mutate(model=modelName(model))

  ##predictive.old <- posteriorPredictive(model)
  tmp <- ggMixture(model, bins=50)
  expect_is(tmp, "ggplot")
})
