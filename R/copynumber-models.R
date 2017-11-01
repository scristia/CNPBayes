#' @include AllClasses.R
NULL

#' @aliases mapping,SingleBatchCopyNumber-method
#' @rdname mapping
setMethod("mapping", "SingleBatchCopyNumber", function(object){
  object@mapping
})

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

#' @aliases mapping,SingleBatchCopyNumber,numeric-method
#' @rdname mapping
setReplaceMethod("mapping", c("SingleBatchCopyNumber", "numeric"),
                 function(object, value){
                   object@mapping <- value
                   object
                 })

#' @aliases mapping,MultiBatchCopyNumber,numeric-method
#' @rdname mapping
setReplaceMethod("mapping", c("MultiBatchCopyNumber", "numeric"),
                 function(object, value){
                   object@mapping <- value
                   object
                 })

#' @aliases mapping,MultiBatchCopyNumberPooled,numeric-method
#' @rdname mapping
setReplaceMethod("mapping", c("MultiBatchCopyNumberPooled", "numeric"),
                 function(object, value){
                   object@mapping <- value
                   object
                 })

setMethod("numberStates", "SingleBatchCopyNumber", function(model){
  length(unique(mapping(model)))
})

setMethod("numberStates", "MultiBatchCopyNumber", function(model){
  length(unique(mapping(model)))
})

setMethod("numberStates", "MultiBatchCopyNumberPooled", function(model){
  length(unique(mapping(model)))
})

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
#' @rdname probCopyNumber
setMethod("probCopyNumber", "SingleBatchCopyNumber", function(model){
  .prob_copynumber(model)
})

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

.remap <- function(z, map){
  for(i in seq_along(map)){
    z[z == i] <- map[i]
  }
  z
}

.relabel_z <- function(object){
  if(!manyToOneMapping(object)) return(z(object))
  ##zz <- map_z(object)
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

#' @aliases copyNumber,SingleBatchCopyNumber-method
#' @rdname copyNumber
setMethod("copyNumber", "SingleBatchCopyNumber", function(object){
  cn <- .relabel_z(object)
  tab <- tibble(y=y(object), cn=cn, z=z(object),
                theta=round(theta(object)[z(object)], 3))
  thetas <- matrix(theta(object), nrow=length(y(object)),
                   ncol=length(theta(object)), byrow=TRUE)
  ## Define z* (a candidate for a new z) to be the nearest mode
  z.candidate <- abs(y(object)-thetas) %>% as.tibble %>%
    apply(1, which.min)
  tab <- tibble(y=y(object),
                theta=theta(object)[z(object)],
                theta.star=theta(object)[z.candidate])
  condition1 <- with(tab,
                     theta <= y & y <= theta.star)
  condition2 <- with(tab,
                     theta.star <= y & y <= theta)
  zstar <- ifelse(condition1 | condition2, z(object), z.candidate)
  object@z <- zstar
  cn <- .relabel_z(object)
##  ## only necessary for models with unequal variances
##  ##tab %>% filter(z != zstar)
##  cn.range <- tab %>%
##    ##group_by(cn) %>%
##    group_by(z) %>%
##    summarize(miny=min(y),
##              maxy=max(y),
##              theta=unique(theta),
##              cn=unique(cn),
##              nearest_mode=) 
##  K <- max(cn)
##  k <- 1
##  while(k < K){
##    cn[ tab$y < cn.range$maxy[k] & cn == k+1 ] <- k
##    k <- k+1
##  }
  cn
})

#' @aliases copyNumber,MultiBatchCopyNumber-method
#' @rdname copyNumber
setMethod("copyNumber", "MultiBatchCopyNumber", function(object){
  .relabel_z(object)
})

#' @aliases copyNumber,MultiBatchCopyNumberPooled-method
#' @rdname copyNumber
setMethod("copyNumber", "MultiBatchCopyNumberPooled", function(object){
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
#' @seealso \code{\link{mapComponents}}
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

.merge_components <- function(model, j){
  x <- component.left <- component.right <- NULL
  dens <- dnorm_poly(model) %>% filter(component!="SB") %>%
    as.tibble %>%
    arrange(component, x) %>%
    filter(component %in% c(j, j+1))
  stats <- dens %>%
    group_by(component) %>%
    summarize(n=n(),
              max.y=max(y),
              x.at.maxy=x[y==max(y)][1])
  ##
  ## where does difference in component 1 density and component 2 density change sign
  ##
  ##d1 <- filter(dens, cnstate==state1 & !duplicated(x))
  d1 <- filter(dens, component == j & !duplicated(x))
  d2 <- filter(dens, component == j+1) %>%
    filter(!duplicated(x))
  keep <- d1$y <= stats$max.y[1] &
    d1$x >= stats$x.at.maxy[1] &
    d2$y <= stats$max.y[2] &
    d2$x <= stats$x.at.maxy[2]
  d1 <- d1[keep, ]
  d2 <- d2[keep, ]
  ## number of quantiles where difference in density is negative and number of
  ## quantiles where difference in density is positive
  probs <- tibble(component.left=d1$y, component.right=d2$y) %>%
    mutate(sign=sign(component.left - component.right)) %>%
    group_by(sign) %>%
    summarize(n=n())
  needs.merge <- ifelse(any(probs$n <= 6) || nrow(probs) < 2, TRUE, FALSE)
  ##
  ## sometimes we have multiple hemizygous deletion components that should be merged
  ##  -- allow a greater separation between components that still merges
  if(!needs.merge){
    ## only allow separable modes to be merged if the right-most component is
    ## not merged with the diploid component
    both.hemizygous <- all(stats$x.at.maxy > -1 & stats$x.at.maxy < -0.2)
    needs.merge <- ifelse((any(probs$n <= 15) || nrow(probs) < 2)
                          && both.hemizygous, TRUE, FALSE)

    both.homozygous <- all( stats$x.at.maxy < -0.9 )
    if(both.homozygous) needs.merge <- TRUE
  }
  needs.merge
}

.merge_mb <- function(model, j){
  x <- NULL
  if(class(model) == "MultiBatchCopyNumberPooled"){
    dens <- dnorm_poly_multibatch_pooled(model) %>% filter(component!="overall") %>%
      arrange(component, x) %>%
      filter(component %in% c(j, j+1))
  } else {
    dens <- dnorm_poly_multibatch(model) %>% filter(component!="overall") %>%
      arrange(component, x) %>%
      filter(component %in% c(j, j+1)) 
  }
  stats <- dens %>%
    group_by(batch, component) %>%
    summarize(n=n(),
              max.y=max(y),
              x.at.maxy=x[y==max(y)][1]) %>%
    filter(batch != "overall")

  ubatch <- levels(dens$batch)
  ubatch <- ubatch[ubatch != "overall"]
  merge.var <- rep(NA, length(ubatch))
  component.left <- component.right <- NULL
  for(i in seq_along(ubatch)){
    d1 <- filter(dens, component == j & batch==ubatch[i]) %>%
      filter(y > 0)
    d2 <- filter(dens, component == (j+1) & batch==ubatch[i]) %>%
      filter(y > 0)
    d1 <- filter(d1, x >= min(d2$x))
    stats2 <- filter(stats, batch == ubatch[i])
    keep <- d1$y <= stats2$max.y[1] &
      d1$x >= stats2$x.at.maxy[1] &
      d2$y <= stats2$max.y[2] &
      d2$x <= stats2$x.at.maxy[2]
    d1 <- d1[keep, ]
    d2 <- d2[keep, ]
    ## number of quantiles where difference in density is negative and number of
    ## quantiles where difference in density is positive
    probs <- tibble(component.left=d1$y, component.right=d2$y) %>%
      mutate(sign=sign(component.left - component.right)) %>%
      group_by(sign) %>%
      summarize(n=n())
    merge.var[i] <- ifelse(any(probs$n <= 6) || nrow(probs) < 2, TRUE, FALSE)
    if(!merge.var[i]){
      is.hemizygous <- all(stats$x.at.maxy > -1 & stats$x.at.maxy < -0.2)
      merge.var[i] <- ifelse((any(probs$n <= 15) || nrow(probs) < 2)
                             && is.hemizygous, TRUE, FALSE)
    }
  }
  ## weight by size of batch
  needs.merge <- sum(merge.var * batchElements(model)) / length(y(model)) > 0.5
  needs.merge
}


setMethod("mergeComponents", "SingleBatchCopyNumber", function(model, j){
  .merge_components(model, j)
})

setMethod("mergeComponents", "MultiBatchCopyNumberPooled", function(model, j){
  .merge_mb(model, j)
})

setMethod("mergeComponents", "MultiBatchCopyNumber", function(model, j){
  .merge_mb(model, j)
})


#' Map mixture components to distinct copy number states
#'
#'
#'
#' @param params a list of mapping parameters
#' @param model a SB, MB, SBP, or MBP model
#' @examples
#' ## Batch model
#' bmodel <- MultiBatchModelExample
#' cn.model <- CopyNumberModel(bmodel)
#' mapping(cn.model)
#' \dontrun{
#'   ggMixture(cn.model)
#' }
#' @seealso \code{\link{CopyNumberModel}} \code{\link{mapParams}}
#' @export
mapComponents <- function(model, params=mapParams()){
  if(length(modes(model)) > 0){
    model <- useModes(model)
  }
  K <- seq_len(k(model) - 1)
  cn.states <- seq_len(k(model))
  counter <- 1
  for(j in K){
    is.merge <- mergeComponents(model, j)
    if(!is.merge){
      ## only increment counter when model does not change
      counter <- counter + 1
      next()
    }
    ## model changes
    cn.states[j+1] <- counter
    cn.states[ seq_along(cn.states) > (j + 1)] <- cn.states[ seq_along(cn.states) > (j + 1)] - 1
    mapping(model) <- cn.states
  }
  model
}


##  p <- probz(model)
##  p <- p/rowSums(p)
##  K <- mapping(model)
##  threshold <- params[["threshold"]]
##  ##
##  ## the denominator should not include observations classified with probability near 1 to
##  ##
##  select <- rowSums(p > 0.99) == 0
##  p <- p[select, , drop=FALSE]
##  if(nrow(p) == 0) return(K)
##  zz <- map_z(model)[select]
##  Ns <- table(factor(zz, levels=K))
##  ##n.uncertain <- colSums(p >= threshold & p <= (1-threshold))
##  ##n.uncertain <- n.uncertain
##  ##frac.uncertain <- n.uncertain/(Ns + 1)
##  frac.uncertain <- colMeans(p >= threshold & p <= (1-threshold))
##  ##
##  ## what fraction of subjects have low posterior probabilities
##  ##
##  cutoff <- params[["proportion.subjects"]]
##  if(all(frac.uncertain < cutoff)){
##    return(K)
##  }
##  vars <- sigma2(model)
##  if(!isSB(model) && !isPooled(model)){
##    vars <- colMeans(vars)
##  }
##  is_pooled <- length(vars) == 1
##  not_pooled <- !is_pooled
##  index <- which(frac.uncertain >= cutoff)
##  if(not_pooled){
##    ## only relevant for non-pooled models
##    var.ratio <- vars/median(vars)
##    var.ratio <- var.ratio[index]
##    outlier.component <- index[var.ratio > params[["outlier.variance.ratio"]]]
##    if(length(outlier.component) > 1) stop("multiple outlier components")
##  } else var.ratio <- 1
##  ## assume that components are ordered and successive indices should be merged
##  index <- .indices_to_merge(model, index)
##  for(i in index){
##    index2 <- head(index, 2)
##    K[index] <- min(K[index])
##    index <- index[-i]
##    if(length(index) == 0) break()
##  }
##  if(not_pooled){
##    if(length(outlier.component) > 0){
##      ## this prevents collapsing the outlier component with other states
##      K[outlier.component] <- outlier.component
##    }
##  }
##  K
##}

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

#' @rdname CopyNumber-methods
#' @aliases CopyNumberModel,SingleBatchModel-method
setMethod("CopyNumberModel", "SingleBatchModel",
          function(model, params=mapParams()){
            model.sb <- SingleBatchCopyNumber(model)
            model.sb <- mapComponents(model.sb, params)
            model.sb
})

#' @rdname CopyNumber-methods
#' @aliases CopyNumberModel,MultiBatchModel-method
setMethod("CopyNumberModel", "MultiBatchModel", function(model, params=mapParams()){
  model <- MultiBatchCopyNumber(model)
  model.cn <- mapComponents(model, params)
  model.cn
})

#' @rdname CopyNumber-methods
#' @aliases CopyNumberModel,MultiBatchPooled-method
setMethod("CopyNumberModel", "MultiBatchPooled", function(model, params=mapParams()){
  model <- MultiBatchCopyNumberPooled(model)
  model.cn <- mapComponents(model, params)
  model.cn
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
