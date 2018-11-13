#' @include AllClasses.R
NULL

#' @aliases mapping,MultiBatchCopyNumber-method
#' @rdname mapping
setMethod("mapping", "MultiBatchCopyNumber", function(object){
  object@mapping
})

#' @aliases mapping,MultiBatchCopyNumber-method
#' @rdname mapping
setMethod("mapping", "MultiBatchCopyNumberPooled", function(object){
  object@mapping
})

#' @aliases mapping,MultiBatchCopyNumber,numeric-method
#' @rdname mapping
setReplaceMethod("mapping", c("MultiBatchCopyNumber", "character"),
                 function(object, value){
                   object@mapping <- value
                   object
                 })

#' @aliases mapping,MultiBatchCopyNumberPooled,numeric-method
#' @rdname mapping
setReplaceMethod("mapping", c("MultiBatchCopyNumberPooled", "character"),
                 function(object, value){
                   object@mapping <- value
                   object
                 })

setMethod("numberStates", "MultiBatchCopyNumber", function(model){
  length(unique(mapping(model)))
})

setMethod("numberStates", "MultiBatchCopyNumberPooled", function(model){
  length(unique(mapping(model)))
})

manyToOneMapping <- function(model){
  copynumber <- NULL
  tab <- tibble(comp=seq_len(k(model)),
                copynumber=mapping(model)) %>%
    group_by(copynumber) %>%
    summarize(n=n())
  ##!identical(comp, map)
  any(tab$n > 1)
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

#' @aliases probCopyNumber,MultiBatchCopyNumber-method
#' @rdname probCopyNumber
setMethod("probCopyNumber", "MultiBatchCopyNumber", function(model){
  .prob_copynumber(model)
})

#' @aliases probCopyNumber,MultiBatchCopyNumberPooled-method
#' @rdname probCopyNumber
setMethod("probCopyNumber", "MultiBatchCopyNumberPooled", function(model){
  .prob_copynumber(model)
})

#' @aliases copyNumber,MultiBatchCopyNumber-method
#' @rdname copyNumber
setMethod("copyNumber", "MultiBatchCopyNumber", function(object){
  component_labels <- mapping(object)
  zz <- map_z(object)
  cn <- component_labels[zz]
  cn
})

#' @aliases copyNumber,MultiBatchCopyNumberPooled-method
#' @rdname copyNumber
setMethod("copyNumber", "MultiBatchCopyNumberPooled", function(object){
  component_labels <- mapping(object)
  zz <- map_z(object)
  cn <- component_labels[zz]
  cn
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
#' @examples
#' mapParams()
#' @export
mapParams <- function(threshold=0.1,
                      proportion.subjects=0.5,
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

isPooled <- function(model){
  class(model) %in% c("SingleBatchPooled", "MultiBatchPooled",
                      "SingleBatchCopyNumberPooled",
                      "MultiBatchCopyNumberPooled")
}
#' @export
#' @rdname CopyNumber-methods
SingleBatchCopyNumber <- function(model){
  sb.model <- as(model, "SingleBatchCopyNumber")
  mapping(sb.model) <- as.character(seq_len(k(model)))
  sb.model
}

#' @return a \code{MultiBatchCopyNumber} instance
#' @export
#' @rdname CopyNumber-methods
MultiBatchCopyNumber <- function(model){
  mb.model <- as(model, "MultiBatchCopyNumber")
  mapping(mb.model) <- as.character(seq_len(k(model)))
  mb.model
}

#' @return a \code{MultiBatchCopyNumber} instance
#' @export
#' @rdname CopyNumber-methods
MultiBatchCopyNumberPooled <- function(model){
  mb.model <- as(model, "MultiBatchCopyNumberPooled")
  mapping(mb.model) <- as.character(seq_len(k(model)))
  mb.model
}

#' @rdname CopyNumber-methods
#' @aliases CopyNumberModel,SingleBatchModel-method
setMethod("CopyNumberModel", "SingleBatchModel",
          function(model, params=mapParams()){
            model.sb <- SingleBatchCopyNumber(model)
            ##model.sb <- mapComponentsxo(model.sb, params)
            model
})

#' @rdname CopyNumber-methods
#' @aliases CopyNumberModel,MultiBatchModel-method
setMethod("CopyNumberModel", "MultiBatchModel", function(model, params=mapParams()){
  model <- MultiBatchCopyNumber(model)
  ##model.cn <- mapComponents(model, params)
  model
})

#' @rdname CopyNumber-methods
#' @aliases CopyNumberModel,MultiBatchPooled-method
setMethod("CopyNumberModel", "MultiBatchPooled", function(model, params=mapParams()){
  model <- MultiBatchCopyNumberPooled(model)
  ##model.cn <- mapComponents(model, params)
  model
})

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
