setAs("MixtureModel", "SummarizedExperiment", function(from, to){
  cnmat <- matrix(y(from), 1, length(y(from)))
  cnmat <- oligoClasses::integerMatrix(cnmat, 1000)
  message("making something up for rowRanges...")
  rr <- GRanges(rep("chr1", nrow(cnmat)),
                IRanges(seq_len(nrow(cnmat)), width=1L))
  names(rr) <- paste0("CNP", seq_len(nrow(cnmat)))
  se <- SummarizedExperiment(assays=SimpleList(medr=cnmat),
                             rowRanges=rr,
                             colData=DataFrame(plate=batch(from)))
  se
})

#' @rdname collapseBatch-method
#' @aliases collapseBatch,SummarizedExperiment-method
setMethod("collapseBatch", "SummarizedExperiment", function(object, plate, THR=0.1){
  plate <- as.character(object$plate)
  collapseBatch(copyNumber(object)[1, ], plate, THR=THR)
})

#' @rdname collapseBatch-method
#' @aliases collapseBatch,numeric-method
setMethod("collapseBatch", "numeric", function(object, plate, THR=0.1){
  N <- choose(length(unique(plate)), 2)
  cond2 <- TRUE
  while(N > 1 && cond2){
    cat('.')
    B <- plate
    plate <- .collapseBatch(object, plate, THR=THR)
    cond2 <- !identical(B, plate)
    N <- choose(length(unique(plate)), 2)
  }
  makeUnique(plate)
})


ksTest <- function(object){
  B <- batch(object)
  uB <- uniqueBatch(object)
  yy <- y(object)
  ks <- matrix(NA, choose(nBatch(object), 2), 4)
  i <- 1
  for(j in seq_along(uB)){
    for(k in seq_along(uB)){
      if(k <= j) next()
      b1 <- uB[j]
      b2 <- uB[k]
      stat <- suppressWarnings(ks.test(yy[B==b1], yy[B==b2]))
      ks[i, ] <- c(b1, b2, stat$statistic, stat$p.value)
      i <- i+1
    }
  }
  colnames(ks) <- c("batch1", "batch2", "stat", "pval")
  ks
}

.collapseBatch <- function(yy, B, THR=0.1){
  uB <- unique(B)
  ## One plate can pair with many other plates.
  for(j in seq_along(uB)){
    for(k in seq_along(uB)){
      if(k <= j) next() ## next k
      b1 <- uB[j]
      b2 <- uB[k]
      stat <- suppressWarnings(ks.test(yy[B==b1], yy[B==b2]))
      if(stat$p.value < THR) next()
      b <- paste(b1, b2, sep=",")
      B[B %in% b1 | B %in% b2] <- b
      ## once we've defined a new batch, return the new batch to the
      ## calling function
      return(B)
    }
  }
  B
}


#' @export
saveBatch <- function(se, batch.file, THR=0.1){
  if(file.exists(batch.file)){
    bt <- readRDS(batch.file)
    return(bt)
  }
  bt <- collapseBatch(se, THR=THR)
  saveRDS(bt, file=batch.file)
  bt
}




#' @export
setMethod("copyNumber", "SummarizedExperiment", function(object, ...){
  assays(object)[["medr"]]/1000
})
