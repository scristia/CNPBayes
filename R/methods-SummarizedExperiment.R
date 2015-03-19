setAs("MixtureModel", "SummarizedExperiment", function(from, to){
  cnmat <- matrix(y(from), 1, length(y(from)))
  cnmat <- oligoClasses::integerMatrix(cnmat, 1000)
  message("making something up for rowRanges...")
  rr <- GRanges(rep("chr1", nrow(cnmat)),
                IRanges(seq_len(nrow(cnmat)), width=1L))
  names(rr) <- paste0("CNP", seq_len(nrow(cnmat)))
  se <- SummarizedExperiment(assays=SimpleList(medr=cnmat),
                             rowData=rr,
                             colData=DataFrame(plate=batch(from)))
  se
})

setMethod("collapseBatch", "numeric", function(object, plate){
  N <- choose(length(unique(plate)), 2)
  cond2 <- TRUE
  while(N > 1 && cond2){
    cat('.')
    B <- plate
    plate <- .collapseBatch(object, plate)
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

.collapseBatch <- function(yy, B){
  uB <- unique(B)
  ## One plate can pair with many other plates.
  for(j in seq_along(uB)){
    for(k in seq_along(uB)){
      if(k <= j) next() ## next k
      b1 <- uB[j]
      b2 <- uB[k]
      stat <- suppressWarnings(ks.test(yy[B==b1], yy[B==b2]))
      if(stat$p.value < 0.1) next()
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
saveBatch <- function(se, batch.file){
  if(file.exists(batch.file)){
    bt <- readRDS(batch.file)
    return(bt)
  }
  bt <- collapseBatch(se)
  saveRDS(bt, file=batch.file)
  bt
}


setMethod("collapseBatch", "SummarizedExperiment", function(object, plate){
  plate <- as.character(object$plate)
  collapseBatch(copyNumber(object)[1, ], plate)
##  N <- choose(length(unique(object$plate)), 2)
##  cond2 <- TRUE
##  while(N > 1 && cond2){
##    cat('.')
##    B <- object$plate
##    object$plate <- .collapseBatch(copyNumber(object)[1,], object$plate)
##    cond2 <- !identical(B, object$plate)
##    N <- choose(length(unique(object$plate)), 2)
##  }
##  makeUnique(object$plate)
})


#' @export
setMethod("copyNumber", "SummarizedExperiment", function(object, ...){
  assays(object)[["medr"]]/1000
})

setMethod("fitMixtureModels", "SummarizedExperiment", function(object, mcmcp, K=1:5, batch){
  cn <- copyNumber(object)[1, ]
  if(missing(batch)){
    object$plate <- collapseBatch(object)
  } else object$plate <- batch
  message("Fitting ", length(K), " mixture models")
  j <- NULL
  fit <- foreach(j = seq_along(K)) %do% {
    cat(".")
    kk <- K[j]
    mp <- mcmcp[j]
    params <- ModelParams("batch", y=cn, k=kk,
                          batch=object$plate,
                          mcmc.params=mp)
    modelk <- initializeBatchModel(params)
    modelk <- suppressMessages(posteriorSimulation(modelk, mp))
    modelk
  }
  fit
})


updateModel <- function(object, mcmcp, zz){
  params <- ModelParams("batch", y=y(object), k=k(object),
                        batch=batch(object),
                        mcmc.params=mcmcp)
  initializeBatchModel(params, zz=factor(zz, levels=seq_len(k(object))))
}


## K and batch ignored
setMethod("fitMixtureModels", "BatchModel", function(object, mcmcp, K, batch){
  posteriorSimulation(object, mcmcp)
})



## Allows restarting a chain
##setMethod("fitMixtureModels", "BatchModel", function(object, mcmcp, K=1:5){
##  ##
##  ## K is ignored
##  ##
##  params <- ModelParams("batch", y=y(object),
##                        k=2,
##                        batch=batch(object),
##                        mcmc.params=mcmcp)
##  model <- initializeModel(params)
##  message("Defining batch variable")
##  model <- collapseBatch(model, mcmcp)
##  message("Fitting ", length(K), " mixture models")
##  fit <- foreach(j = seq_along(K)) %do% {
##    cat(".")
##    kk <- K[j]
##    params <- ModelParams("batch", y=cn, k=kk,
##                          batch=batch(model),
##                          mcmc.params=mcmcp)
##    modelk <- initializeModel(params)
##    modelk <- posteriorSimulation(modelk, mcmcp)
##    modelk
##  }
##  fit
##})
