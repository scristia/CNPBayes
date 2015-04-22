test_SummarizedExperiment <- function(){
  library(GenomicRanges) ## S4Vectors alone is insufficient.  I think
  ## GenomicRanges/IRanges adds Rle methods
  library(S4Vectors)
  Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2))
  gr <- GRanges(seqnames=Rle(paste("chr", c(1, 2, 1, 3), sep=""), c(1, 3, 2, 4)),
                ranges=IRanges(1:10, width=10:1, names=letters[1:10]),
                strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
                score=1:10,
                GC=seq(1, 0, length=10))
  ##disjoin(gr)
  cn <- matrix(rnorm(100), 10, 10)
  colnames(cn) <- letters[1:10]
  rr <- GRanges(rep("chr1", 10), IRanges(1:10, 2:11))
  se <- SummarizedExperiment(assays=cn, rowRanges=rr)
  validObject(se[1:5, ])
}
