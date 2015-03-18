setMethod("tracePlot", "BatchModel", function(object, name, ...){
  ilist <- foreach(j=1:nBatch(object)) %do% seq(j, nBatch(object)*k(object), nBatch(object))
  uB <- uniqueBatch(object)
  if(name=="theta"){
    ##op <- par(mfrow=c(3, 3), las=1)
    foreach(k=1:nBatch(object)) %do% {
      plot.ts(thetac(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k], ...)
    }
    ##par(op)
  }
  if(name=="sigma"){
    ##op <- par(mfrow=c(nBatch(object)/3, 3), las=1)
    foreach(k=1:nBatch(object)) %do% {
      plot.ts(sigmac(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k],...)
    }
    ##par(op)
  }
  if(name=="p"){
    ##op <- par(mfrow=c(1, k(object)), las=1)
    ##foreach(k=1:nBatch(object)) %do% {
    plot.ts(pic(object),  ...)
    ##plot.ts(pic(object), col="gray", ...)
    ##par(op)
  }
  if(name=="mu"){
    ##op <- par(mfrow=c(1, k(object)), las=1)
    plot.ts(muc(object),  ...)
    ##par(op)
  }
  if(name=="tau"){
    ##op <- par(mfrow=c(1, k(object)), las=1)
    plot.ts(tauc(object),  ...)
    ##par(op)
  }
})


setMethod("plot", "BatchModel", function(x, y, use.current=FALSE, show.batch=TRUE, ...){
  .plotBatch(x, use.current, show.batch, ...)
})


##
##
## use empirical, batch=specific mixing probabilities
##
##
.plotBatch <- function(object, use.current=FALSE, show.batch=TRUE, ...){
  L <- length(y(object))
  hist(object, ...)
  xx <- seq(min(observed(object)), max(observed(object)),  length.out=100)
  if(!use.current){
    thetas <- thetaMean(object)
    sds <- sigmaMean(object)
    ##P <- pMean(object)  does not give the batch specific probabilities
    mapz <- factor(map(object), levels=seq_len(k(object)))
    tabz <- table(batch(object), mapz)
    tabz <- tabz/rowSums(tabz)
    tabz <- tabz[uniqueBatch(object), , drop=FALSE]
    P <- tabz
    ##P <- mapz/rowSums(mapz)
  } else {
    thetas <- theta(object)
    sds <- sigma(object)
    P <- p(object)
    P <- matrix(P, nBatch(object), k(object), byrow=TRUE)
    rownames(P) <- uniqueBatch(object)
  }
  marginal <- matrix(NA, length(xx), nBatch(object)*k(object))
  cols <- brewer.pal(max(k(object), 3),  "Set1")
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

#' @export
plot.PosteriorFiles <- function(x, y, bayes.factor, m.y, ...){
  browser()
  best.model <- substr(names(bayes.factor), 1, 2)
  se <- y
  object <- x
  if(length(object) > 1){
    message("More than one CNP in model file; only plotting first CNP")
    object <- object[1]
  }
  model <- readRDS(model(object)[1])
  if(isMarginalModel(object)){
    names(model) <- paste0("M", 1:4)
  } else names(model) <- paste0("B", 1:4)
  for(j in 1:4){
    idm <- names(model)[j]
    plot(model[[j]], ...)
    if(idm==best.model){
      rect(xleft=-5, xright=5, ybottom=0, ytop=50, col="lightblue", border="lightblue")
      ##hist(y(model[[j]]), breaks=length(y(model[[j]]))/10,
      ##col="gray", border="gray", add=TRUE)
      plot(model[[j]], ..., add=TRUE)
    }
    lim <- par("usr")
    m <- round(m.y[idm], 0)
    lgnd <- bquote(log(m[y])==.(m))
    legend("topleft",
           legend=as.expression(lgnd), bg="white",
           bty="n", title=names(model)[j],
           title.col="blue")
    if(j==1 && isMarginalModel(object)){
      pos <- c(chromosome(se)[i],
               paste0("start: ", prettyNum(start(se)[i], big.mark=",")),
               paste0("end: ", prettyNum(end(se)[i], big.mark=",")),
               paste0(rowData(se)$nSNPs_affy6[i], " / ",
                      rowData(se)$nmarkers_affy6[i]),
               paste0("source: ", rowData(se)$source[i]))
      legend("left", legend=pos, bg="gray90", box.col="gray90")
    }
    if(idm==best.model){
      lgnd <- bquote(.(names(bayes.factor))==.(bayes.factor))
      legend("topright", legend=as.expression(lgnd), bg="white", bty="n")
    }
  }
}

#' @export
plotModel <- function(model.list, se, ...){
  options(digits=3)
  options(scipen=5)
  m.y <- unlist(lapply(model.list, getPosteriorStats))
  bf <- bayesFactor(m.y)
  model1 <- model.list[[1]]
  plot.PosteriorFiles(model.list[[1]], se,
                      bayes.factor=bf, m.y=m.y,
                      ...,
                      ##xlim=xlim,
                      xaxt="n",
                      breaks=150, main="",
                      yaxt="n",
                      use.current=TRUE)
  if(length(model.list) == 1) return()
  plot.PosteriorFiles(model.list[[2]], se,
                      bayes.factor=bf, m.y=m.y,
                      ...,
                      ##xlim=xlim,
                      xaxt="n",
                      breaks=150, main="",
                      yaxt="n",
                      use.current=TRUE)
  mtext(rownames(se), 3, line=0, outer=TRUE)
  ##    plot.PosteriorFiles(B[i], se.ea,
  ##                        bayes.factor=bf, m.y=m.y, xlim=xlim,
  ##                        breaks=150, main="",
  ##                        yaxt="n", use.current=TRUE)
}
