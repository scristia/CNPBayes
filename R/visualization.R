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

# summarizeMarginal <- function(model.list){
#   xlist <- lapply(model.list, m.y)
#   if(isMarginalModel(model.list[[1]])){
#     nms <- paste0("M", seq_along(model.list))
#   } else nms <- paste0("B", seq_along(model.list))
#   mns <- sapply(xlist, mean)
#   rg <- sapply(xlist, function(x) diff(range(x)))
#   df <- data.frame(mean=mns, range=rg)
#   rownames(df) <- nms
#   df
# }
