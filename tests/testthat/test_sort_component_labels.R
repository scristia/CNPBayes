context("Order component labels")

test_that("ordered_thetas", {
  expect_true(isOrdered(MultiBatchModelExample))
  theta(MultiBatchModelExample)[1, 3] <- -5
  expect_false(isOrdered(MultiBatchModelExample))
  bmodel <- sortComponentLabels(MultiBatchModelExample)
  expect_true(isOrdered(bmodel))
})
