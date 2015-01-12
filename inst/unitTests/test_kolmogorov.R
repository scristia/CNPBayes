test_KolmogorovSmirnov <- function(){
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.2, -1.0,
                    0, 0, 0.2,
                    0.8, 0.8, 1.0), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)

  truth <- simulateBatchData(N=2500,
                             .batch=rep(letters[1:3], length.out=2500),
                             .k=3,
                             .alpha=rep(1, k),
                             means=means,
                             sds=sds)
  b <- collapseBatch(truth)
  checkIdentical(unique(b), c("a,b", "c"))
}
