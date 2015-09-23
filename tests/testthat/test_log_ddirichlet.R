test_that("log_ddirichlet", {
  p <- c(0.81, 0.049, 0.141)
  alpha <- c(73, 4, 8)
  x.cpp <- log_ddirichlet_(p, alpha)
  x.gtools <- gtools::ddirichlet(p, alpha)
  expect_equal(exp(x.cpp), x.gtools)
})
