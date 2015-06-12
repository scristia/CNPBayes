#' @rdname tracePlot-method
#' @aliases tracePlot,BatchModel-method
setMethod("tracePlot", "BatchModel", function(object, name, ...){
  j <- NULL
  ilist <- foreach(j=1:nBatch(object)) %do% seq(j, nBatch(object)*k(object), nBatch(object))
  uB <- uniqueBatch(object)
  if(name=="theta"){
    foreach(k=1:nBatch(object)) %do% {
      plot.ts(thetac(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k], ...)
    }
  }
  if(name=="sigma"){
    foreach(k=1:nBatch(object)) %do% {
      plot.ts(sigmac(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k],...)
    }
  }
  if(name=="p"){
    plot.ts(pic(object),  ...)
  }
  if(name=="mu"){
    plot.ts(muc(object),  ...)
  }
  if(name=="tau"){
    plot.ts(tauc(object),  ...)
  }
})

##
## use empirical, batch=specific mixing probabilities
##
##
## plot.PosteriorFiles <- function(x, y, bayes.factor, m.y, ...){
##   best.model <- substr(names(bayes.factor), 1, 2)
##   se <- y
##   object <- x
##   model <- readRDS(model(object)[1])
##   if(isMarginalModel(object)){
##     names(model) <- paste0("M", 1:4)
##   } else names(model) <- paste0("B", 1:4)
##   for(j in 1:4){
##     idm <- names(model)[j]
##     plot(model[[j]], ...)
##     if(idm==best.model){
##       rect(xleft=-5, xright=5, ybottom=0, ytop=50, col="lightblue", border="lightblue")
##       ##hist(y(model[[j]]), breaks=length(y(model[[j]]))/10,
##       ##col="gray", border="gray", add=TRUE)
##       plot(model[[j]], ..., add=TRUE)
##     }
##     lim <- par("usr")
##     m <- round(m.y[idm], 0)
##     lgnd <- bquote(log(m[y])==.(m))
##     legend("topleft",
##            legend=as.expression(lgnd), bg="white",
##            bty="n", title=names(model)[j],
##            title.col="blue")
##     rr <- tryCatch(rowRanges(se), error=function(e) NULL)
##     if(is.null(rr)) rr <- rowData(se)
##     if(j==1 && isMarginalModel(object)){
##       pos <- c(chromosome(se),
##                paste0("start: ", prettyNum(start(se), big.mark=",")),
##                paste0("end: ", prettyNum(end(se), big.mark=",")),
##                ##paste0(rowRanges(se)$nSNPs_affy6[i], " / ",
##                paste0(rr$nSNPs, " / ",
##                       ##rowRanges(se)$nmarkers_affy6[i]),
##                       rr$nmarkers),
##                paste0("source: ", rr$source))
##       legend("left", legend=pos, bg="gray90", box.col="gray90")
##     }
##     if(idm==best.model){
##       lgnd <- bquote(.(names(bayes.factor))==.(bayes.factor))
##       legend("topright", legend=as.expression(lgnd), bg="white", bty="n")
##     }
##     if(!isMarginalModel(object)){
##       at <- pretty(c(-3,1), n=6)
##       axis(1, at=at, labels=at, cex.axis=0.9)
##     }
##   }
## }

summarizeMarginal <- function(model.list){
  xlist <- lapply(model.list, m.y)
  if(isMarginalModel(model.list[[1]])){
    nms <- paste0("M", seq_along(model.list))
  } else nms <- paste0("B", seq_along(model.list))
  mns <- sapply(xlist, mean)
  rg <- sapply(xlist, function(x) diff(range(x)))
  df <- data.frame(mean=mns, range=rg)
  rownames(df) <- nms
  df
}
