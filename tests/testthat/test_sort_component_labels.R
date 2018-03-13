context("Order component labels")

test_that("ordered_thetas", {
  expect_false(isOrdered(MultiBatchModelExample))
  bmodel <- sortComponentLabels(MultiBatchModelExample)
  expect_true(isOrdered(bmodel))
})
