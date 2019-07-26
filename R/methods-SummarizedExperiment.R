integerMatrix <- function(x, scale=100) {
        if(!is(x, "matrix")) stop("argument x must be a matrix")
        dms <- dimnames(x)
        if(scale != 1){
                xx <- as.integer(x*scale)
        } else xx <- as.integer(x)
        x <- matrix(xx, nrow(x), ncol(x))
        dimnames(x) <- dms
        return(x)
}

setAs("MixtureModel", "SummarizedExperiment", function(from, to){
  cnmat <- matrix(y(from), 1, length(y(from)))
  cnmat <- integerMatrix(cnmat, 1000)
  message("making something up for rowRanges...")
  rr <- GRanges(Rle("chr1", nrow(cnmat)),
                IRanges(seq_len(nrow(cnmat)), width=1L))
  rr <- GRanges(rep("chr1", nrow(cnmat)),
                IRanges(seq_len(nrow(cnmat)), width=1L))
  names(rr) <- paste0("CNP", seq_len(nrow(cnmat)))
  se <- SummarizedExperiment(assays=SimpleList(medr=cnmat),
                             rowRanges=rr,
                             colData=DataFrame(plate=batch(from)))
  se
})


setMethod("collapseBatch", "SummarizedExperiment", function(object, provisional_batch, THR=0.1){
  batch <- as.character(object$provisional_batch)
  collapseBatch(assays(object)[["medr"]][1, ]/1000)
})


setMethod("collapseBatch", "numeric", function(object, provisional_batch, THR=0.1, nchar=8){
  N <- choose(length(unique(provisional_batch)), 2)
  if(N == 1){
    B <- provisional_batch
    provisional_batch <- .combineBatches(object, provisional_batch, THR=THR)
  }
  cond2 <- TRUE
  while(N > 1 && cond2){
    B <- provisional_batch
    provisional_batch <- .combineBatches(object, provisional_batch, THR=THR)
    cond2 <- !identical(B, provisional_batch)
    N <- choose(length(unique(provisional_batch)), 2)
  }
  makeUnique(provisional_batch, nchar)
})


##
## Combine the most similar batches first.
##
.combineBatches <- function(yy, B, THR=0.1){
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

saveBatch <- function(se, batch.file, THR=0.1){
  if(file.exists(batch.file)){
    bt <- readRDS(batch.file)
    return(bt)
  }
  bt <- collapseBatch(se, THR=THR)
  saveRDS(bt, file=batch.file)
  bt
}
