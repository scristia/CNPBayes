#' @include AllClasses.R
NULL

#' Mapping of mixture component indices to copy number states
#'
#' @param object a \code{SingleBatchCopyNumber} or \code{MultiBatchCopyNumber} instance
#' @export
#' @examples
#' cn.model <- SingleBatchCopyNumber(MarginalModelExample)
#' map <- mapComponents(cn.model)
#' mapping(cn.model) <- map
#' @rdname CopyNumber-methods
setGeneric("mapping", function(object) standardGeneric("mapping"))

#' @export
#' @rdname CopyNumber-methods
setGeneric("mapping<-", function(object, value) standardGeneric("mapping<-"))

#' @aliases mapping,SingleBatchCopyNumber-method
#' @rdname CopyNumber-methods
setMethod("mapping", "SingleBatchCopyNumber", function(object){
  object@mapping
})

#' @aliases mapping,MultiBatchCopyNumber-method
#' @rdname CopyNumber-methods
setMethod("mapping", "MultiBatchCopyNumber", function(object){
  object@mapping
})

#' @aliases mapping,MultiBatchCopyNumber-method
#' @rdname CopyNumber-methods
setMethod("mapping", "MultiBatchCopyNumberPooled", function(object){
  object@mapping
})

#' @param value a numeric vector mapping component indices to copy number state indices
#' @aliases mapping,SingleBatchCopyNumber,numeric-method
#' @rdname CopyNumber-methods
setReplaceMethod("mapping", c("SingleBatchCopyNumber", "numeric"),
                 function(object, value){
                   object@mapping <- value
                   object
                 })

#' @aliases mapping,MultiBatchCopyNumber,numeric-method
#' @rdname CopyNumber-methods
setReplaceMethod("mapping", c("MultiBatchCopyNumber", "numeric"),
                 function(object, value){
                   object@mapping <- value
                   object
                 })

#' @aliases mapping,MultiBatchCopyNumberPooled,numeric-method
#' @rdname CopyNumber-methods
setReplaceMethod("mapping", c("MultiBatchCopyNumberPooled", "numeric"),
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

setMethod("numberStates", "MultiBatchCopyNumberPooled", function(model){
  length(unique(mapping(model)))
})

#' Posterior probabilities for copy number states
#'
#' In contrast to posterior probabilities for mixture components, this function
#' returns posterior probabilities for distinct copy number states.
#' @param model a \code{SingleBatchCopyNumber} or \code{MultiBatchCopyNumber} instance
#' @rdname CopyNumber-methods
#' @export
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

#' @aliases probCopyNumber,SingleBatchCopyNumber-method
#' @rdname CopyNumber-methods
setMethod("probCopyNumber", "SingleBatchCopyNumber", function(model){
  .prob_copynumber(model)
})

#' @aliases probCopyNumber,MultiBatchCopyNumber-method
#' @rdname CopyNumber-methods
setMethod("probCopyNumber", "MultiBatchCopyNumber", function(model){
  .prob_copynumber(model)
})

#' @aliases probCopyNumber,MultiBatchCopyNumberPooled-method
#' @rdname CopyNumber-methods
setMethod("probCopyNumber", "MultiBatchCopyNumberPooled", function(model){
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


#' Parameters for mapping mixture components to distinct copy number states
#'
#' @param threshold numeric value in [0, 0.5]. For a given observation (sample),
#'   mixture component probabilities > threshold and less than 1-threshold are
#'   combined.
#' 
#' @param proportion.subjects numeric value in [0, 1]. Two components are
#'   combined if the fraction of subjects with component probabilities in the
#'   range [threshold, 1-threshold] exceeds this value.
#' @param outlier.variance.ratio if the ratio of the component variance to the
#'   median variance of the other component exceeds the value of this argument,
#'   the component is considered to correspond to outliers. 
#' @param max.homozygous length-2 numeric vector of cutoffs used for
#'   establishing a homozygous deletion component. The first element is the
#'   cutoff for the mean log R ratios when there are 2 or more states. The
#'   second element is the cutoff for the mean log R ratio when there are 3 or
#'   more states.
#' @param min.foldchange a length-one numeric vector. When there are 3 or more
#'   states, we compute the ratio of the distance between the means of the
#'   sorted components 1 and 2 (the 2 components with lowest means) and the
#'   distance between the means of components 2 and 3. If the ratio (i) exceeds
#'   the value specified by this parameter, (ii) there are 3 or more states, and
#'   (iii) the first component has a mean less than \code{max_homozygous[2]}, we
#'   infer that the first component is a homozygous deletion.
#' @export
#' @examples
#' mapParams()
#' @seealso \code{\link{mapComponents}}
mapParams <- function(threshold=0.1, proportion.subjects=0.5,
                      outlier.variance.ratio=5,
                      max.homozygous=c(-1.5, -0.5), ## cutoff for 2 state and 3 state models
                      min.foldchange=1.5){
  list(threshold=threshold,
       proportion.subjects=proportion.subjects,
       outlier.variance.ratio=outlier.variance.ratio,
       max.homozygous=setNames(max.homozygous, c("2state", "3state")),
       min.foldchange=min.foldchange)
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

isOutlier <- function(model, params=mapParams()){
  vars <- sigma2(model)
  var.ratio <- vars/median(vars)
  var.ratio > params[["outlier.variance.ratio"]]
}


#' Map mixture components to distinct copy number states
#'
#' @param params a list of mapping parameters
#' @examples
#' mm <- MarginalModelExample
#' cn.model <- SingleBatchCopyNumber(mm)
#' mapping(cn.model) <- mapComponents(cn.model)
#' mapping(cn.model)
#' \dontrun{
#'  ggSingleBatch(mm)
#' }
#' ## Batch model
#' bmodel <- BatchModelExample
#' bmodel <- MultiBatchCopyNumber(bmodel)
#' mapping(bmodel) <- mapComponents(bmodel)
#' mapping(bmodel)
#' @export
mapComponents <- function(model, params=mapParams()){
  p <- probz(model)
  K <- mapping(model)
  threshold <- params[["threshold"]]
  ##
  ## the denominator should not include observations classified with probability near 1 to
  ##
  select <- rowSums(p > 0.99) == 0
  p <- p[select, , drop=FALSE]
  if(nrow(p) == 0) return(K)
  frac.uncertain <- colMeans(p >= threshold & p <= (1-threshold))
  ##
  ## what fraction of subjects have low posterior probabilities
  ##
  ##frac.uncertain <- n.uncertain/as.numeric(table(z(model)))
  cutoff <- params[["proportion.subjects"]]
  if(all(frac.uncertain < cutoff)){
    return(K)
  }
  vars <- sigma2(model)
  if(!isMarginalModel(model)){
    vars <- colMeans(vars)
  }
  var.ratio <- vars/median(vars)
  index <- which(frac.uncertain >= cutoff)
  var.ratio <- var.ratio[index]
  outlier.component <- index[var.ratio > params[["outlier.variance.ratio"]]]
  if(length(outlier.component) > 1) stop("multiple outlier components")
  ## assume that components are ordered and successive indices should be merged
  index <- .indices_to_merge(model, index)
  for(i in index){
    index2 <- head(index, 2)
    K[index] <- min(K[index])
    index <- index[-i]
    if(length(index) == 0) break()
  }
  if(length(outlier.component) > 0){
    ## this prevents collapsing the outlier component with other states
    K[outlier.component] <- outlier.component
  }
  K
}

#' Constructor for copy number models (single-batch)
#'
#' Used for mapping mixture components to distict copy number states
#'
#' @return a \code{SingleBatchCopyNumber} instance
#' @examples
#' SingleBatchCopyNumber(SingleBatchModelExample)
#' @export
#' @rdname CopyNumber-methods
SingleBatchCopyNumber <- function(model){
  sb.model <- as(model, "SingleBatchCopyNumber")
  mapping(sb.model) <- seq_len(k(model))
  sb.model
}

#' @return a \code{MultiBatchCopyNumber} instance
#' @export
#' @rdname CopyNumber-methods
MultiBatchCopyNumber <- function(model){
  mb.model <- as(model, "MultiBatchCopyNumber")
  mapping(mb.model) <- seq_len(k(model))
  mb.model
}

#' @return a \code{MultiBatchCopyNumber} instance
#' @export
#' @rdname CopyNumber-methods
MultiBatchCopyNumberPooled <- function(model){
  mb.model <- as(model, "MultiBatchCopyNumberPooled")
  mapping(mb.model) <- seq_len(k(model))
  mb.model
}

#' Constructs a CopyNumberModel from SB, SBP, MB, or MBP models
#'
#' @param model a SB, SBP, MB, or MBP model
#' @export
#' @examples
#' sb <- SingleBatchModelExample
#' cn.model <- CopyNumberModel(sb)
#' @rdname CopyNumber-methods
setGeneric("CopyNumberModel", function(model) standardGeneric("CopyNumberModel"))

#' @rdname CopyNumber-methods
#' @aliases CopyNumberModel,SingleBatchModel-method
setMethod("CopyNumberModel", "SingleBatchModel", function(model){
  model <- SingleBatchCopyNumber(model)
  mapping(model) <- mapComponents(model)
  model
})

#' @rdname CopyNumber-methods
#' @aliases CopyNumberModel,MultiBatchModel-method
setMethod("CopyNumberModel", "MultiBatchModel", function(model){
  model <- MultiBatchCopyNumber(model)
  mapping(model) <- mapComponents(model)
  model
})

#' @rdname CopyNumber-methods
#' @aliases CopyNumberModel,MultiBatchPooled-method
setMethod("CopyNumberModel", "MultiBatchPooled", function(model){
  model <- MultiBatchCopyNumberPooled(model)
  mapping(model) <- mapComponents(model)
  model
})

#' Map mixture components to copy number states
#'
#' @examples
#' cn.model <- SingleBatchCopyNumber(MarginalModelExample)
#' mapping(cn.model) <- mapComponents(cn.model)
#' mapCopyNumber(cn.model)
#' @param params  a list of parameters for mapping component indices to copy number
#'
#' @return a factor vector of length k
#'
#' @export
#' @rdname CopyNumber-methods
mapCopyNumber <- function(model,
                          params=mapParams()) {
  ##  labels <- c("Homozygous Deletion",
  ##              "Hemizygous Deletion",
  ##              "Diploid",
  ##              "Addition")
  map <- mapping(model)
  map <- as.integer(factor(map))  ## 1, 1, 3 turned into 1, 1, 2
  many.to.one <- manyToOneMapping(model)
  K <- k(model)
  is.outlier <- isOutlier(model)
  if(K >= 4 && !many.to.one){
    map <- map - 1
    if(any(is.outlier)){
      map[is.outlier] <- "outlier"
    }
    map <- factor(map)
    return(map)
  }
  ##
  ## check for homozygous deletion
  ##
  max.homozygous <- params[["max.homozygous"]]
  thetas <- theta(model)
  thetas <- tapply(thetas, map, mean)
  homo_deletion <- (thetas[1] <= max.homozygous[["2state"]])
  ##
  ## Above is too stringent of a criteria in many situations. Since the
  ## difference of means comparing homozygous to hemizygous is several fold
  ## larger than the difference between hemizygous and diploid, we can allow a
  ## less stringent cutoff for the mean LRR
  ##
  min.foldchange <- params[["min.foldchange"]]
  if(length(thetas) >= 3 && !homo_deletion){
    deltas <- diff(thetas)
    fold.change.of.deltas <- deltas[1]/deltas[2]
    homo_deletion <- (thetas[1] <= max.homozygous[["3state"]]) &&
      (fold.change.of.deltas >= min.foldchange)
  }
  if (homo_deletion) {
    map <- map - 1
    ## check outlier
    map[is.outlier] <- "outlier"
    map <- factor(map)
    return(map)
  }
  ##
  ## First component is not homozygous
  ##
  taken <- which.min(abs(thetas-0))
  if(taken >= 3){
    ## this should not occur if there is not a homozygous component. However,
    ## its possible that the merging step failed to merge to hemizygous
    ## deletion components.
    ##
    ## Allow diploid to be third comp. if it is the modal component. Otherwise,
    ## diploid component should be the 2nd component if there are no homozygous
    ## deletions.
    modal.comp <- which.max(table(z(model)))
    if(modal.comp >= 3){
      taken <- modal.comp
    } else{
      taken <- 2
    }
  }
  if(taken==2){
    ## nothing to do. second component is diploid
    map[is.outlier] <- "outlier"
    map <- factor(map)
    return(map)
  }
  if(taken==3){
    map <- map - 1
    map[is.outlier] <- "outlier"
    map <- factor(map)
    return(map)
  }
  if(taken > 3) stop("check labeling")
  ## taken is 1
  map <- map + 1
  map[is.outlier] <- "outlier"
  map <- factor(map)
  map
}
