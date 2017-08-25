context("Deprecated")

test_that("Deprecated", {
  expect_warning(DensityModel(SingleBatchModelExample))
})
