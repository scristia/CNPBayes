test_KolmogorovSmirnov <- function(){
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.2, -1.0,
                    0, 0, 0.2,
                    0.8, 0.8, 1.0), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)

  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             theta=means,
                             sds=sds,
                             p=c(1/3,1/3, 1/3))
  b <- collapseBatch(truth)
  checkIdentical(unique(b), c("a,b", "c"))

  b2 <- collapseBatch(y(truth), as.character(oligoClasses::batch(truth)))
  checkIdentical(b, b2)
  tmpfile <- tempfile()
  saveBatch(truth, batch.file=tmpfile)
  checkTrue(file.exists(tmpfile))

  m <- matrix(y(truth), nrow=1)
  colnames(m) <- paste0("s", seq_len(ncol(m)))
  se <- SummarizedExperiment(assays=SimpleList(medr=m))
  se$plate <- oligoClasses::batch(truth)
  tmpfile <- tempfile()
  saveBatch(se, batch.file=tmpfile)
  checkTrue(file.exists(tmpfile))
}
