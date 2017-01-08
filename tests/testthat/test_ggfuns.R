context("gg functions")

test_that("ggfun", {
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.2, 0.2, 0.2))
  ggSingleBatch(truth)


  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.02, 0.02, 0.02))
  ## hard to see marginal when it looks the same as the truth
  ggSingleBatch(truth)

  ggSingleBatch(MarginalModelExample, bins=200)

  BatchModelExample, bins=200

})
