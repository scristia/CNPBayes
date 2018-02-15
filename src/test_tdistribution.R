context("Simulating from t distribution")

test_that("rst", {
  set.seed(2)
  tmp1 <- rst(1)
  set.seed(2)
  tmp2 <- rlocScale_t(1, mu=0, sigma=1, df=100)
  expect_identical(tmp1, tmp2)
})
