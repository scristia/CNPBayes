MarginalModelList <- function(model_list=list(), names=character(),
                              maxperm=3L, mode_index=list(),
                              mode_list=list()){
  if(length(model_list) > 0){
    type <- unique(sapply(model_list, class))
  } else type="MarginalModel"
  new("MarginalModelList", model_list=model_list, elementType=type,
      names=names, maxperm=maxperm, mode_index=mode_index,
      mode_list=mode_list)
}


maxperm <- function(object) object@maxperm
modeIndex <- function(object) object@mode_index
modeList <- function(object) object@mode_list
modeDifference <- function(object) object@mode_list

BatchModelList <- function(model_list=list(), names=character(),
                           maxperm=3L, mode_index=list(),
                           mode_list=list()){
  new("BatchModelList", model_list=model_list, elementType="BatchModel",
      names=names, maxperm=maxperm, mode_index=mode_index,
      mode_list=mode_list)
}

elementType <- function(object) object@elementType

setValidity("MarginalModelList", function(object){
  msg <- TRUE
  type <- elementType(object)
  if(type != "MarginalModel") {
    msg <- "Elements of MarginalModelList must be of class 'MarginalModel'"
    return(msg)
  }
  msg
})

setValidity("BatchModelList", function(object){
  msg <- TRUE
  type <- elementType(object)
  if(type != "BatchModel") {
    msg <- "Elements of BatchModelList must be of class 'BatchModel'"
    return(msg)
  }
  msg
})

setMethod("modelList", "ModelList", function(object) object@model_list)
setMethod("names", "ModelList", function(x) x@names)
setMethod("length", "ModelList", function(x) length(modelList(x)))
setMethod("k", "ModelList", function(object) sapply(object@model_list, k))

setMethod("show", "ModelList", function(object){
  cat("Object of class 'MarginalModelList'\n")
  if(length(object) == 0){
    cat("  K           : 0\n")
    cat("  top two     : NA \n")
    cat("  bayes factor: NA \n")
    return()
  }
  cat(" K               : ", paste(k(object), collapse=", "), "\n")
  bf <- bayesFactor(summary(object))
  cat(" top two         : ", names(bf), "\n")
  cat(" log bayes factor: ", round(bf, 1), "\n")
  cat(" see summary(), modelList()\n")
})

setMethod("m.y", "ModelList", function(object){
  mylist <- setNames(lapply(modelList(object), m.y), names(object))
  mylist
})

## returns another ModelList
setMethod("[", "ModelList", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@model_list <- modelList(x)[i]
    x@names <- names(x)[i]
  }
  x
})

## returns MarginalModel or BatchModel
setMethod("[[", "ModelList", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x <- modelList(x)[[i]]
  }
  x
})

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


setMethod("summary", "ModelList", function(object, ...){
  my <- m.y(object)
  mns <- sapply(my, mean)
  rg <- sapply(my, function(x) diff(range(x)))
  df <- data.frame(k=k(object),
                   mean=mns,
                   range=rg,
                   max_delta_mode=sapply(modeDifference(object), max))
  rownames(df) <- names(object)
  df
})

setReplaceMethod("mcmcParams", "ModelList", function(object, value){
  for(i in seq_along(object)){
    mcmcParams(object[[i]]) <- value
  }
  object
})

setMethod("iter", "ModelList", function(object) sapply(modelList(object), iter))
setMethod("burnin", "ModelList", function(object) sapply(modelList(object), burnin))
setMethod("nStarts", "ModelList", function(object) sapply(modelList(object), nStarts))


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

setMethod("best", "ModelList", function(object){
  bf <- bayesFactor(summary(object))
  ##paste0("k = ", substr(names(bf), 2, 2))
  as.integer(substr(names(bf), 2, 2))
})

setMethod("lapply", "ModelList", function(X, FUN, ...){
  results <- vector("list", length(X))
  for(i in seq_along(X)){
    results[[i]] <- FUN(X[[i]], ...)
  }
  results
})
