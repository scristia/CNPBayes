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


setMethod("hist", "MixtureModel", function(x, ...){
  op <- par(las=1)
  yy <- y(x)
  ##if(rsample > 0) yy <- sample(yy, rsample)
  args <- list(...)
  if(!"breaks" %in% names(args)){
    L <- length(yy)
    hist(yy, breaks=L/10,
         col="gray"
       , border="gray",  xlab="y",
         freq=FALSE, ...)
  } else {
    hist(yy, col="gray", border="gray", xlab="y",
         freq=FALSE, ...)
  }
  par(op)
})




##
##
## use empirical, batch=specific mixing probabilities
##
##
.plotBatch <- function(object, use.current=FALSE, show.batch=TRUE, ...){
  L <- length(y(object))
  hist(object, ...)
  xx <- seq(min(observed(object)), max(observed(object)),  length.out=250)
  if(!use.current){
    thetas <- thetaMean(object)
    sds <- sigmaMean(object)
    pz <- probz(object)
    mapz <- apply(pz, 1, which.max)
    mapz <- factor(mapz, levels=seq_len(k(object)))
    tabz <- table(batch(object), mapz)
    tabz <- tabz/rowSums(tabz)
    P <- tabz
  } else {
    thetas <- theta(object)
    sds <- sigma(object)
    P <- p(object)
    P <- matrix(P, nBatch(object), k(object), byrow=TRUE)
    rownames(P) <- uniqueBatch(object)
  }
  marginal <- matrix(NA, length(xx), nBatch(object)*k(object))
  cols <- brewer.pal(max(k(object), 3),  "Set1")
  ## compute marginal density for each batch
  ## weight the 'overall' curve by the batch frequencies
  batchPr <- table(batch(object))/length(y(object))
  marginal <- drawEachComponent(xx, uniqueBatch(object), thetas, sds, P, batchPr, cols)
  marginal.cum.prob <- rowSums(marginal)
  limits <- list(range(y(object), na.rm=TRUE), range(marginal.cum.prob, na.rm=TRUE))
  lines(xx, marginal.cum.prob, col="black", lwd=2)
}

drawdens <- function(x, mean, sd, p1, p2, col){
  p.x <- p1*p2*dnorm(x, mean, sd)
  lines(x, p.x, col=col, lwd=2)
  p.x
}

drawEachComponent <- function(x, batches, thetas, sds, P, batchPr, cols){
  mlist <- list()
  for(j in seq_len(ncol(thetas))){
    mlist[[j]] <- drawEachBatch(x, batches, thetas[, j], sds[, j], P[, j], batchPr, cols[j])
  }
  marginal <- do.call(cbind, mlist)
}

drawEachBatch <- function(x, batches, thetas, sds, p1, p2, col){
  marginal <- matrix(NA, length(x), length(batches))
  for(b in seq_along(batches)){
    marginal[, b] <- drawdens(x, thetas[b], sds[b], p1[b], p2[b], col)
  }
  marginal
}




#' @export
plot.PosteriorFiles <- function(x, y, bayes.factor, m.y, ...){
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
      pos <- c(chromosome(se),
               paste0("start: ", prettyNum(start(se), big.mark=",")),
               paste0("end: ", prettyNum(end(se), big.mark=",")),
               ##paste0(rowData(se)$nSNPs_affy6[i], " / ",
               paste0(rowData(se)$nSNPs, " / ",
                      ##rowData(se)$nmarkers_affy6[i]),
                      rowData(se)$nmarkers),
               paste0("source: ", rowData(se)$source))
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
  ##model.list is a list of PosteriorFiles
  options(digits=3)
  if(nrow(se) > 1) message("Using first row of SummarizedExperiment")
  options(scipen=5)
  m.y <- unlist(lapply(model.list, getPosteriorStats))
  bf <- bayesFactor(m.y)
  model1 <- model.list[[1]]
  if(length(model.list) > 1) xaxt <- "n" else xaxt="s"
  plot.PosteriorFiles(model.list[[1]], se[1, ],
                      bayes.factor=bf, m.y=m.y,
                      ...,
                      ##xlim=xlim,
                      xaxt=xaxt,
                      breaks=150, main="",
                      yaxt="n",
                      use.current=TRUE)
  if(length(model.list) == 1) return()
  plot.PosteriorFiles(model.list[[2]], se[1, ],
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
