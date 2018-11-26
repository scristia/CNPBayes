context("Simulating from t distribution")

test_that("rst", {
  set.seed(2)
  u <- rchisq(1, 100)
  tmp1a <- rst(1)
  set.seed(2)
  tmp1 <- rst(1, u=u)
  ## tmp1 and tmp1a are not identical -- not sure why
  set.seed(2)
  df <- 100
  tmp2 <- rlocScale_t(1, mu=0, sigma=1, df=df, u)
  expect_identical(tmp1, tmp2)
})
