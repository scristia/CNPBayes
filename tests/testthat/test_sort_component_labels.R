context("Order component labels")

test_that("ordered_thetas", {
  expect_true(isOrdered(SingleBatchModelExample))
  expect_false(isOrdered(BatchModelExample))

  bmodel <- sortComponentLabels(BatchModelExample)
  expect_true(isOrdered(bmodel))

  mmodel <- SingleBatchModelExample
  theta(mmodel) <- rev(theta(mmodel))
  expect_false(isOrdered(mmodel))
})
