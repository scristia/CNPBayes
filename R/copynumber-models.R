#' @include AllClasses.R
NULL

setGeneric("mapping", function(object) standardGeneric("mapping"))
setGeneric("mapping<-", function(object, value) standardGeneric("mapping<-"))

setMethod("mapping", "SingleBatchCopyNumber", function(object){
  object@mapping
})

setMethod("mapping", "MultiBatchCopyNumber", function(object){
  object@mapping
})

setReplaceMethod("mapping", c("SingleBatchCopyNumber", "numeric"),
                 function(object, value){
                   object@mapping <- value
                   object
                 })

setReplaceMethod("mapping", c("MultiBatchCopyNumber", "numeric"),
                 function(object, value){
                   object@mapping <- value
                   object
                 })

setGeneric("numberStates", function(model) standardGeneric("numberStates"))

setMethod("numberStates", "SingleBatchCopyNumber", function(model){
  length(unique(mapping(model)))
})

setMethod("numberStates", "MultiBatchCopyNumber", function(model){
  length(unique(mapping(model)))
})

setGeneric("probCopyNumber", function(model) standardGeneric("probCopyNumber"))

manyToOneMapping <- function(model){
  comp <- seq_len(k(model))
  map <- mapping(model)
  !identical(comp, map)
}

.prob_copynumber <- function(model){
  if(!manyToOneMapping(model)){
    return(probz(model))
  }
  S <- numberStates(model)
  N <- numberObs(model)
  result <- matrix(NA, N, S)
  if(S == 1) {
    result[] <- 1
    return(result)
  }
  ##
  ## If we've reached this point, there is a many-to-one mapping of components
  ## to copy number states.
  ##
  map <- mapping(model)
  pz <- probz(model)

  map.list <- split(seq_len(k(model)), map)
  for(i in seq_along(map.list)){
    j <- map.list[[i]]
    p <- pz[, j, drop=FALSE]
    if(ncol(p) > 1){
      result[, i] <- rowSums(p)
    } else result[, i] <- as.numeric(p)
  }
  result
}

setMethod("probCopyNumber", "SingleBatchCopyNumber", function(model){
  .prob_copynumber(model)
})

setMethod("probCopyNumber", "MultiBatchCopyNumber", function(model){
  .prob_copynumber(model)
})

setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))


.remap <- function(z, map){
  for(i in seq_along(map)){
    z[z == i] <- map[i]
  }
  z
}

.relabel_z <- function(object){
  if(!manyToOneMapping(object)) return(z(object))
  zz <- z(object)
  map <- mapping(object)
  if(numberStates(object) == 1){
    zz[zz != map[1]] <- map[1]
    return(zz)
  }
  ##
  ## multiple states and many-to-one mapping
  ##
  zz <- .remap(zz, map)
  zz
}

setMethod("copyNumber", "SingleBatchCopyNumber", function(object){
  .relabel_z(object)
})

setMethod("copyNumber", "MultiBatchCopyNumber", function(object){
  .relabel_z(object)
})


#' Parameters for merging two components
#'
#' @param threshold numeric value in [0, 0.5]. For a given observation (sample),
#'   mixture component probabilities > threshold and less than 1-threshold are
#'   combined.
#' @param proportion.subjects numeric value in [0, 1]. Two components are
#'   combined if the fraction of subjects with component probabilities in the
#'   range [threshold, 1-threshold] exceeds this value.
mapParams <- function(threshold=0.1, proportion.subjects=0.5){
  list(threshold=threshold,
       proportion.subjects=proportion.subjects)
}


.indices_to_merge <- function(model, index){
  if(length(index) >= 2) return(index)
  if(length(index) == 1){
    ## merge first component with second
    if(index == 1) return(1:2)
    ## merge last component with 2nd to last component
    if(index == k(model)) return(c(index-1, index))
    ##
    ## index is neither first nor last
    ##  -> merge to closest neighbor
    dist1 <- theta(model)[index] - theta(model)[index-1]
    dist2 <- theta(model)[index+1] - theta(model)[index]
    if(dist1 < dist2){
      index <- c(index-1, index)
    } else{
      index <- c(index, index+1)
    }
  }
  index
}

mapComponents <- function(model, params=mapParams()){
  p <- probz(model)
  K <- mapping(model)
  threshold <- params[["threshold"]]
  uncertain.component <- p > threshold & p <= (1-threshold)
  ##
  ## what fraction of subjects have low posterior probabilities
  ##
  frac.uncertain <- colMeans(uncertain.component)
  cutoff <- params[["proportion.subjects"]]
  if(all(frac.uncertain < cutoff)){
    return(K)
  }
  index <- which(frac.uncertain >= cutoff)
  ##if(length(index) < 2) stop("expect index vector to be >= 2")
  ## assume that components are ordered and successive indices should be merged
  index <- .indices_to_merge(model, index)
  for(i in index){
    index2 <- head(index, 2)
    K[index] <- min(K[index])
    index <- index[-i]
    if(length(index) == 0) break()
  }
  K
}

SingleBatchCopyNumber <- function(model){
  sb.model <- as(model, "SingleBatchCopyNumber")
  mapping(sb.model) <- seq_len(k(model))
  sb.model
}

MultiBatchCopyNumber <- function(model){
  mb.model <- as(model, "MultiBatchCopyNumber")
  mapping(mb.model) <- seq_len(k(model))
  mb.model
}


