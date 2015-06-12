ModelList <- function(model_list=list(),
                      names=character(),
                      data=numeric()){
  ##browser()
  if(length(model_list) > 0){
    ##data <- y(model_list[[1]])
    model_list <- stripData(model_list)
  }
  new("ModelList", model_list=model_list, names=names, data=data)
}

MarginalModelList <- function(model_list=list(),
                              names=character(),
                              data=numeric()){
  if(length(model_list) > 0){
    data <- y(model_list[[1]])
    model_list <- stripData(model_list)
  }
  new("MarginalModelList", model_list=model_list, names=names,
      data=data)
}


maxperm <- function(object) object@maxperm
modeIndex <- function(object) object@mode_index
modeList <- function(object) object@mode_list
modeDifference <- function(object) object@mode_list

BatchModelList <- function(model_list=list(), names=character(), data=numeric()){
  if(length(model_list) > 0){
    data <- y(model_list[[1]])
    model_list <- stripData(model_list)
  }
  new("BatchModelList", model_list=model_list, names=names, data=data)
}

elementType <- function(object) object@elementType

setValidity("MarginalModelList", function(object){
  ## todo; validity methods
  msg <- TRUE
  msg
})

setValidity("BatchModelList", function(object){
  msg <- TRUE
  msg
})

setMethod("modelList", "ModelList", function(object) {
  setNames(object@model_list, names(object))
})

setMethod("names", "ModelList", function(x) x@names)
setMethod("length", "ModelList", function(x) length(modelList(x)))

#' @rdname k-method
#' @aliases k,ModelList-method
setMethod("k", "ModelList", function(object) sapply(object@model_list, k))

setMethod("show", "ModelList", function(object){
  cat("Object of class ", class(object), "\n")
  cat(" K: ", paste(k(object), collapse=", "), "\n")
  cat(" see summary(), modelList()\n")
})

#' Extract a particular Marginal or BatchModel but still as class ModelList.
#' @aliases [,ModelList-method
#' @docType methods
#' @rdname extract-methods
setMethod("[", "ModelList", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@model_list <- modelList(x)[i]
    x@names <- names(x)[i]
  }
  x
})

#' @rdname y-method
#' @aliases y,ModelList-method
setMethod("y", "ModelList", function(object)  object@data )

#' Extract a particular Marginal or BatchModel. 
#' @aliases [[,ModelList-method
#' @docType methods
#' @rdname extract-methods
setMethod("[[", "ModelList", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    yy <- y(x)
    x <- modelList(x)[[i]]
    x@data <- yy
  }
  x
})

#' Replace a MarginalModel or BatchModel in a ModelList.
#' @aliases [[<-,ModelList-method
#' @docType method
#' @rdname extract-methods
setReplaceMethod("[[", c("ModelList"), # "MixtureModel"),
                 function(x, i, j, ..., drop=FALSE, value){
  if(!missing(i)){
    x@model_list[[i]] <- value
  }
  x
})

setReplaceMethod("modelList", c("ModelList", "MixtureModel"),
                 function(object, value){
                   object@model_list <- value
                   object
                 })

setReplaceMethod("mcmcParams", "ModelList", function(object, value){
  for(i in seq_along(object)){
    mcmcParams(object[[i]]) <- value
  }
  object
})

#' @rdname iter-method
#' @aliases iter,ModelList-method
setMethod("iter", "ModelList", function(object) sapply(modelList(object), iter))

#' @rdname burnin-method
#' @aliases burnin,ModelList-method
setMethod("burnin", "ModelList", function(object) sapply(modelList(object), burnin))

#' @rdname nStarts-method
#' @aliases nStarts,ModelList-method
setMethod("nStarts", "ModelList", function(object) sapply(modelList(object), nStarts))

#' @rdname iter-method
#' @aliases iter<-,ModelList-method
setReplaceMethod("iter", "ModelList", function(object, value){
  if(length(value) == length(object)){
    for(i in seq_along(object)){
      iter(object[[i]]) <- value[i]
    }
  } else {
    for(i in seq_along(object)){
      iter(object[[i]]) <- value[1]
    }
  }
  object
})

#' @rdname burnin-method
#' @aliases burnin<-,ModelList-method
setReplaceMethod("burnin", "ModelList", function(object, value){
  if(length(value) == length(object)){
    for(i in seq_along(object)){
      burnin(object[[i]]) <- value[i]
    }
  } else {
    for(i in seq_along(object)){
      burnin(object[[i]]) <- value[1]
    }
  }
  object
})

#' @rdname nStarts-method
#' @aliases nStarts<-,ModelList-method
setReplaceMethod("nStarts", "ModelList", function(object, value){
  if(length(value) == length(object)){
    for(i in seq_along(object)){
      nStarts(object[[i]]) <- value[i]
    }
  } else {
    for(i in seq_along(object)){
      nStarts(object[[i]]) <- value[1]
    }
  }
  object
})

setMethod("sapply", "ModelList", function(X, FUN, ...,
                                          simplify=TRUE,
                                          USE.NAMES=TRUE){
  FUN <- match.fun(FUN)
  answer <- lapply(X, FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
    names(answer) <- X
  if (!identical(simplify, FALSE) && length(answer))
    simplify2array(answer, higher = (simplify == "array"))
  else answer
})

setMethod("lapply", "ModelList", function(X, FUN, ...){
  results <- vector("list", length(X))
  for(i in seq_along(X)){
    results[[i]] <- FUN(X[[i]], ...)
  }
  results
})
