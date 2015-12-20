context("Batches")

test_that("Small batches", {
    model <- BatchModel(data=rnorm(1000),
                        batch=c(rep(1L, 990),
                                rep(2L, 10)))
})
