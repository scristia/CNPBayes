#' @rdname tracePlot-method
#' @aliases tracePlot,BatchModel-method
setMethod("tracePlot", "BatchModel", function(object, name, ...){
  ilist <- lapply(1:nBatch(object),
                  function(j, object) {
                      seq(j, 
                          nBatch(object) * k(object), 
                          nBatch(object))
                  }, object=object)
  uB <- uniqueBatch(object)
  if(name=="theta"){
    for(k in 1:nBatch(object)) {
      plot.ts(thetac(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k], ...)
    }
  }
  if(name=="sigma"){
    for(k in 1:nBatch(object)) {
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
