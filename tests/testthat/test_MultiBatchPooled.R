context("MultiBatch pooled variances")

## constructor
test_that("MultiBatchPooled", {
  model <- MultiBatchModel2()
  expect_is(model, "MultiBatchModel")
  model <- MultiBatchPooled()
  expect_is(model, "MultiBatchPooled")
  expect_equivalent(sigma(model), numeric())
})
