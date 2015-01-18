.plotBatchP <- function(object, use.current=FALSE, show.batch=TRUE, ...){
  L <- length(y(object))
  hist(object)
  xx <- seq(min(observed(object)), max(observed(object)),  length.out=10e3)
  if(!use.current){
    thetas <- thetaMean(object)
    sds <- sigmaMean(object)
    sds <- orderCols(thetas, sds)
    P <- pMean(object)
    P <- orderCols(thetas, P)
    thetas <- orderCols(thetas, thetas)
  } else {
    thetas <- theta(object)
    sds <- sigma(object)
    P <- p(object)
  }
  marginal <- matrix(NA, length(xx), nBatch(object)*k(object))
  cols <- brewer.pal(max(k(object), 3),  "Set1")
  B <- batch(object)
  marginal.prob <- matrix(NA, length(xx), k(object))
  batchPr <- table(batch(object))/length(y(object))
  m <- 1
  for(j in seq_len(k(object))){
    for(b in uniqueBatch(object)){
      p.x <- dnorm(xx, mean=thetas[b, j], sd=sds[b, j])
      p.x <- batchPr[b] * P[b, j] * p.x
      if(show.batch) lines(xx, p.x, col=cols[j], lwd=2)
      marginal[, m] <- p.x
      m <- m+1
    }
  }
  marginal.cum.prob <- rowSums(marginal)
  limits <- list(range(y(object), na.rm=TRUE), range(marginal.cum.prob, na.rm=TRUE))
  lines(xx, marginal.cum.prob, col="black", lwd=2)
}
