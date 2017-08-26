context("Order component labels")

test_that("ordered_thetas", {
  expect_true(isOrdered(SingleBatchModelExample))
  expect_false(isOrdered(MultiBatchModelExample))

  bmodel <- sortComponentLabels(MultiBatchModelExample)
  expect_true(isOrdered(bmodel))

  mmodel <- SingleBatchModelExample
  theta(mmodel) <- rev(theta(mmodel))
  expect_false(isOrdered(mmodel))
})
