context("Disjoin method in IRanges")

test_that("disjoin", {
  ir <- IRanges(c(1, 1, 4, 10), c(6, 3, 8, 10))
  ## works fine
  disjoin(ir)  # IRanges(c(1, 4, 7, 10), c(3, 6, 8, 10))
  g <- GRanges(c("chr1", "chr5"), IRanges(c(1, 3), c(2, 4)))
  ## does not work
  disjoin(g)
})
