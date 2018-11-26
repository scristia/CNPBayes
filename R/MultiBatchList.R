#' @include MultiBatch.R
NULL

setClass("MultiBatchList", representation(data="tbl_df",
                                          down_sample="integer",
                                          specs="tbl_df",
                                          parameters="list",
                                          chains="list",
                                          current_values="list",
                                          summaries="list",
                                          flags="list"))

modelSpecs <- function(models, K, data, down_sample) {
  if(missing(models)){
    models <- c("SB", "SBP", "MB", "MBP")
    if(missing(K)) K <- 1:4
    models <- lapply(models, paste0, 1:4) %>%
      unlist
  }
  if(missing(data)) data <- modelData()
  model.list <- vector("list", length(models))
  for(i in seq_along(models)){
    model.list[[i]] <- model_spec(models[i], data, down_sample)
  }
  tab <- do.call(rbind, model.list) %>%
    filter(! model %in% c("SBP1", "MBP1"))
  ##  if("SB1" %in% tab$model && "MB1" %in% tab$model){
  ##    tab <- filter(tab, model != "MB1")
  ##  }
  if(all(tab$number_batches == 1)){
    ## keep only single batch models
    tab <- filter(tab, substr(model, 1, 2) == "SB")
  }
  tab
}

setValidity("MultiBatchList", function(object){
  msg <- TRUE
  if(any(diff(batch(object)) < 0)){
    msg <- "not in batch order"
    return(msg)
  }
  msg
})

setMethod("current_values", "MultiBatchList", function(object){
  ## a list
  object@current_values
})

setMethod("specs", "MultiBatchList", function(object){
  object@specs
})

##
## Subsetting
##
setMethod("[[", c("MultiBatchList", "numeric"), function(x, i){
  ## return MultiBatch instance
  model <- specs(x)$model[i]
  hp <- hyperParams(x)
  k(hp) <- specs(x)$k[i]
  parameters(x)[["hp"]] <- hp
  if(specs(x)$number_batches[i] == 1){
    assays(x)$batch <- 1L
  }
  nm <- substr(model, 1, 3)
  if(nm == "SBP" || nm == "MBP"){
    mb <- MultiBatchP(model=model,
                      data=assays(x),
                      down_sample=down_sample(x),
                      specs=specs(x)[i, ],
                      parameters=parameters(x),
                      current_values=current_values(x)[[i]],
                      chains=chains(x)[[i]],
                      summaries=summaries(x)[[i]],
                      flags=flags(x)[[i]])
    return(mb)
  }
  mb <- MultiBatch(model=model,
                   data=assays(x),
                   down_sample=down_sample(x),
                   specs=specs(x)[i, ],
                   parameters=parameters(x),
                   current_values=current_values(x)[[i]],
                   chains=chains(x)[[i]],
                   summaries=summaries(x)[[i]],
                   flags=flags(x)[[i]])
  mb
})

setMethod("[[", c("MultiBatchList", "character"), function(x, i){
  j <- match(i, names(x))
  x[[j]]
})

setMethod("[", "MultiBatchList", function(x, i, j, ...){
  ## return MultiBatchList 
  MultiBatchList(data=assays(x),
                 down_sample=down_sample(x),
                 specs=specs(x)[i, ],
                 parameters=parameters(x),
                 current_values=current_values(x)[i],
                 summaries=summaries(x)[i],
                 chains=chains(x)[i],
                 flags=flags(x)[i])
})



setReplaceMethod("[[", "MultiBatchList", function(x, i, value){
  ## return MultiBatchList instance
  down_sample(x) <- down_sample(value)
  current_values(x)[[i]] <- current_values(value)
  summaries(x)[[i]] <- summaries(value)
  flags(x)[[i]] <- flags(value)
  chains(x)[[i]] <- chains(value)
  x
})

## value is a MultiBatchList
setReplaceMethod("[", "MultiBatchList", function(x, i, j, value){
  if(length(value) != length(i)) stop("Length of replacement and i must be the same")
  current_values(x)[i] <- current_values(value)
  summaries(x)[i] <- summaries(value)
  flags(x)[i] <- flags(value)
  chains(x)[i] <- chains(value)
  x
})


setMethod("setModes", "MultiBatchList", function(object){
  modal.ordinates <- modes(object)
  current_values(object) <- modal.ordinates
  object
})

setMethod("modes", "MultiBatchList", function(object){
  summary.list <- summaries(object)
  mode.list <- lapply(summary.list, "[[", "modes")
  mode.list
})

setReplaceMethod("current_values", c("MultiBatchList", "list"),
                 function(object, value){
                   object@current_values <- value
                   object
})

setMethod("nrow", "MultiBatchList", function(x) nrow(assays(x)))

setMethod("parameters", "MultiBatchList", function(object){
  object@parameters
})

setReplaceMethod("parameters", c("MultiBatchList", "list"), function(object, value){
  object@parameters <- value
  object
})

setMethod("mcmcParams", "MultiBatchList", function(object){
  parameters(object)[["mp"]]
})

setMethod("chains", "MultiBatchList", function(object){
  object@chains
})

setReplaceMethod("chains", c("MultiBatchList", "list"), function(object, value){
  object@chains <- value
  object
})

setMethod("numBatch", "MultiBatchList", function(object) as.integer(specs(object)$number_batches[[1]]))

setReplaceMethod("mcmcParams", c("MultiBatchList", "McmcParams"), function(object, value){
  mp <- value
  S <- iter(mp)
  B <- numBatch(object)
  K <- k(object)
  if(iter(object) != S){
    if(S > iter(object)){
      parameters(object)[["mp"]] <- value
      ## create a new chain
      chains(object) <- listChains2(specs(object), parameters(object))
    } else {
      parameters(object)[["mp"]] <- value
      index <- seq_len(S)
      chains(object) <- lapply(chains(object), "[", index)
    }
    return(object)
  }
  ## if we've got to this point, it must be safe to update mcmc.params
  ## (i.e., size of chains is not effected)
  parameters(object)[["mp"]] <- value
  object
})

setMethod("hyperParams", "MultiBatchList", function(object){
  parameters(object)[["hp"]]
})

setReplaceMethod("hyperParams", c("MultiBatchList", "Hyperparameters"),
                 function(object, value){
                   parameters(object)[["hp"]] <- value
                   object
                 })

setMethod("burnin", "MultiBatchList", function(object){
  burnin(mcmcParams(object))
})

setReplaceMethod("burnin", c("MultiBatchList", "numeric"),
                 function(object, value){
                   burnin(mcmcParams(object)) <- as.numeric(value)
                   object
})

setMethod("iter", "MultiBatchList", function(object){
  iter(mcmcParams(object))
})

setMethod("k", "MultiBatchList", function(object) specs(object)$k)

setReplaceMethod("iter", c("MultiBatchList", "numeric"),
                 function(object, value){
                   mp <- mcmcParams(object)
                   iter(mp) <- value
                   mcmcParams(object) <- mp
                   object
                 })

setMethod("thin", "MultiBatchList", function(object){
  thin(mcmcParams(object))
})

setReplaceMethod("thin", c("MultiBatchList", "numeric"), function(object, value){
  thin(mcmcParams(object)) <- as.integer(value)
  object
})

setMethod("nStarts", "MultiBatchList", function(object) nStarts(mcmcParams(object)))

setReplaceMethod("nStarts", c("MultiBatchList", "numeric"),
                 function(object, value){
                   nStarts(object) <- as.integer(value)
                   object
                 })

setReplaceMethod("nStarts", c("MultiBatchList", "integer"),
                 function(object, value){
                   nStarts(mcmcParams(object)) <- value
                   object
                 })

setMethod("flags", "MultiBatchList", function(object) object@flags)

setReplaceMethod("flags", c("MultiBatchList", "list"), function(object, value){
  object@flags <- value
  object
})

setMethod("assays", "MultiBatchList", function(x, ..., withDimnames) x@data)

setReplaceMethod("assays", c("MultiBatchList", "tbl_df"), function(x, value){
  x@data <- value
  x
})



setReplaceMethod("downSampledData", c("MultiBatchList", "tbl_df"), function(x, value){
  x@data <- value
  x
})

##
## Data accessors
##
setMethod("batch", "MultiBatchList", function(object){
  downSampledData(object)[["batch"]]
})

setReplaceMethod("oned", c("MultiBatchList", "numeric"), function(object, value){
  downSampledData(object)[["oned"]] <- value
  object
})

setMethod("oned", "MultiBatchList", function(object){
  downSampledData(object)[["oned"]]
})

setMethod("down_sample", "MultiBatchList", function(object) object@down_sample)

setMethod("downSampledData", "MultiBatchList", function(object){
  assays(object)[down_sample(object), ]
})

setReplaceMethod("down_sample", "MultiBatchList", function(object, value){
  object@down_sample <- value
  object
})

setMethod("dfr", "MultiBatchList", function(object){
  df <- dfr(hyperParams(object))
  df
})


##
## Summaries
##
setMethod("summaries", "MultiBatchList", function(object){
  object@summaries
})

setReplaceMethod("summaries", c("MultiBatchList", "list"), function(object, value){
  object@summaries <- value
  object
})


setMethod("marginal_lik", "MultiBatchList", function(object){
  sum.list <- summaries(object)
  sapply(sum.list, "[[", "marginal_lik")
})

setMethod("show", "MultiBatchList", function(object){
  ##callNextMethod()
  cls <- class(object)
  n.mcmc <- iter(object) * thin(object)
  saved.mcmc <- iter(object)
  if(nrow(object) > 0){
    b <- table(batch(object)) %>%
      range
    b.range <- paste(b, collapse="-")
  } else b.range <- "0"
  cat(paste0("An object of class ", cls), "\n")
  cat("     n. models           :", nrow(specs(object)), "\n")
  cat("     n. obs              :", nrow(assays(object)), "\n")
  cat("     n. sampled          :", nrow(downSampledData(object)), "\n")
  cat("     n. batches          :", nBatch(object), "\n")
  cat("     k range             :", paste0(min(k(object)), "-", max(k(object))), "\n")
  cat("     nobs/batch          :", b.range, "\n")
  cat("     saved mcmc          :", saved.mcmc, "\n")
})

listChains2 <- function(model_specs, parameters){
  S <- iter(parameters[["mp"]])
  B <- model_specs$number_batches
  K <- model_specs$k
  mc.list <- vector("list", nrow(model_specs))
  for(i in seq_along(mc.list)){
    nm <- substr(model_specs$model[i], 1, 3)
    if(nm == "SBP" || nm == "MBP"){
      mc.list[[i]] <- initialize_mcmcP(K[i], S, B[i])
    } else {
      mc.list[[i]] <- initialize_mcmc(K[i], S, B[i])
    }
  }
  names(mc.list) <- model_specs$model
  mc.list
}

listValues <- function(model_specs, data, hp){
  values.list <- vector("list", nrow(model_specs))
  for(i in seq_along(values.list)){
    nm <- substr(model_specs$model[i], 1, 3)
    if(nm == "SBP" || nm == "MBP"){
      values.list[[i]] <- modelValuesP(model_specs[i, ], data, hp)
    } else {
      values.list[[i]] <- modelValues2(model_specs[i, ], data, hp)
    }
  }
  names(values.list) <- model_specs$model
  values.list
}

hyperParamList <- function(k){
  lapply(k, function(x) Hyperparameters(k=x))
}

listSummaries <- function(model_specs){
  sum.list <- vector("list", nrow(model_specs))
  for(i in seq_along(sum.list)){
    nm <- substr(model_specs$model[i], 1, 3)
    if(nm == "SBP" || nm == "MBP"){
      sum.list[[i]] <- modelSummariesP(model_specs[i, ])
    } else {
      sum.list[[i]] <- modelSummaries(model_specs[i, ])
    }
  }
  names(sum.list) <- model_specs$model
  sum.list
}

listFlags <- function(model_specs){
  flag.list <- replicate(nrow(model_specs), modelFlags(), simplify=FALSE)
  names(flag.list) <- model_specs$model
  flag.list
}


MultiBatchList <- function(models,
                           data=modelData(),
                           down_sample=seq_len(nrow(data)),
                           K,
                           specs=modelSpecs(models, K, data,
                                            down_sample),
                           num_models=nrow(specs),
                           burnin=1000L,
                           iter=1000L,
                           thin=1L,
                           nStarts=4L,
                           max_burnin=32000L,
                           hp=Hyperparameters(),
                           mp=McmcParams(iter=iter,
                                         thin=thin,
                                         nStarts=nStarts,
                                         burnin=burnin,
                                         max_burnin=max_burnin),
                           parameters=modelParameters(hp=hp, mp=mp),
                           chains=listChains2(specs, parameters),
                           current_values,
                           summaries=listSummaries(specs),
                           flags=listFlags(specs)){
  if(missing(current_values)){
    current_values <- listValues(specs,
                                 data[down_sample, , drop=FALSE],
                                 hp)
  }
  new("MultiBatchList",
      data=data,
      down_sample=down_sample,
      specs=specs,
      parameters=parameters,
      chains=chains,
      current_values=current_values,
      summaries=summaries,
      flags=flags)
}



setMethod("downSampleModel", "MultiBatchList", function(object, N=1000, i){
  if(!missing(N)){
    if(N >= nrow(assays(object))){
      return(object)
    }
  }
  if(missing(i)){
    i <- sort(sample(seq_len(nrow(object)), N, replace=TRUE))
  }
  b <- assays(object)$batch[i]
  current.vals <- current_values(object)
  replaceVec <- function(x, nm, i){
    x[[nm]] <- x[[nm]][i]
    x
  }
  replaceMat <- function(x, nm, i){
    x[[nm]] <- x[[nm]][i, , drop=FALSE]
    x
  }
  current.vals <- lapply(current.vals, replaceVec, "u", i)
  current.vals <- lapply(current.vals, replaceVec, "z", i)
  current.vals <- lapply(current.vals, replaceMat, "probz", i)
  current_values(object) <- current.vals
  down_sample(object) <- i
  object@specs$number_sampled <- length(i)
  object
})

setMethod("upSampleModel", "MultiBatchList", function(object){
  for(i in seq_along(object)){
    mb <- object[[i]]
    up <- upSampleModel(mb)
    object[[i]] <- up
  }
  object
})

setMethod("length", "MultiBatchList", function(x) {
  length(chains(x))
})

##
## Coersion
##
setAs("MultiBatch", "MultiBatchList", function(from){
  MultiBatchList(models=modelName(from),
                 data=assays(from),
                 down_sample=down_sample(from),
                 specs=specs(from),
                 parameters=parameters(from),
                 current_values=list(current_values(from)),
                 summaries=list(summaries(from)),
                 chains=list(chains(from)),
                 flags=list(flags(from)))
})

setAs("MultiBatchList", "MultiBatch", function(from){
  from[[1]]
})

extractFromModelList <- function(from, FUN){
  models <- sapply(from, modelName)
  nmodels <- elementNROWS(models)
  cv <- vector("list", sum(nmodels))
  k <- 1
  for(i in seq_along(models)){
    m <- from[[i]]
    if(length(models[[i]]) == 1){
      cv[[k]] <- FUN(m)
      k <- k+1
      next()
    }
    for(j in seq_along(m)){
      cv[[k]] <- FUN(m[[j]])
      k <- k+1
    }
  }
  names(cv) <- unlist(models)
  cv
}

setAs("list", "MultiBatchList", function(from){
  it <- sapply(from, iter)
  if(length(unique(it)) > 1){
    stop("models can not be combined")
  }
  models <- sapply(from, modelName)
  ##
  ## pass data from MultiBatch model, if any
  ##
  ix <- grep("MB", models)[1]
  dat <- assays(from[[ix]])
  specs <- modelSpecs(models,
                      sapply(from, k),
                      data=dat,
                      down_sample=down_sample(from[[1]]))
  MultiBatchList(models=sapply(from, modelName),
                 data=dat,
                 down_sample=down_sample(from[[1]]),
                 parameters=parameters(from[[1]]),
                 current_values=lapply(from, current_values),
                 summaries=lapply(from, summaries),
                 chains=lapply(from, chains),
                 flags=lapply(from, flags))
})

fitModelK <- function(model.list){
  sb <- model.list[[1]]
  mod.list <- model.list[-1]
  mod.list2 <- vector("list", length(mod.list))
  ##
  ## only fit multibatch models for given k if the
  ## corresponding single-batch model converges
  ##
  sb2 <- mcmc2( sb )
  if( convergence(sb2) ){
    for(j in seq_along(mod.list)){
      tmp <- mcmc2(mod.list[[j]], guide=sb2)
      mod.list2[[j]] <- tmp
    }
    mod.list3 <- c(sb2, mod.list2)
    names(mod.list3) <- sapply(mod.list3, modelName)
    converged <- sapply(mod.list3, convergence)
    mod.list3 <- mod.list3[ converged ]
    ix <- order(sapply(mod.list3, marginal_lik), decreasing=TRUE)
    mod.list3 <- mod.list3[ix]
    if(length(mod.list3) >= 1){
      result <- as(mod.list3, "MultiBatchList")
    }
  } else {
    result <- sb2
  }
  result
}

setMethod("mcmc2", "MultiBatchList", function(object, guide){
  mlist <- listModelsByDecreasingK(object)
  k4 <- fitModelK(mlist[[1]])
  k3 <- fitModelK(mlist[[2]])
  k2 <- fitModelK(mlist[[3]])
  k1 <- fitModelK(mlist[[4]])
  list(k4, k3, k2, k1)
})

##setMethod("sapply", "MultiBatchList",
##          function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE){
##
##          })

setMethod("names", "MultiBatchList", function(x){
  specs(x)$model
})

setMethod("max_burnin", "MultiBatchList", function(object){
  max_burnin(mcmcParams(object))
})

setMethod("marginal_lik", "MultiBatchList", function(object){
  ml <- sapply(summaries(object), "[[", "marginal_lik")
  names(ml) <- names(object)
  ml
})


setMethod("lapply", "MultiBatchList", function(X, FUN, ...){
  result <- vector("list", length(X))
  for(i in seq_along(X)){
    result[[i]] <- FUN(X[[i]], ...)
  }
  result
})

setMethod("sapply", "MultiBatchList", function(X, FUN, ..., simplify=TRUE, USE.NAMES=TRUE){
  result <- lapply(X, FUN, ...)
  if(simplify){
    result <- unlist(result)
  }
  result
})

##fitSingleBatch <- function(model.list){
##  nms <- names(model.list)
##  sb.model <- paste0("SB", k(model.list)[1])
##  ix <- match(sb.model, nms)
##  SB <- models[[sb.model]]
##  SB2 <- mcmc2(SB)
##  SB2
##}

setMethod("modelName", "MultiBatchList", function(object){
  specs(object)$model
})

setReplaceMethod("max_burnin", "MultiBatchList", function(object, value){
  mp <- mcmcParams(object)
  max_burnin(mp) <- as.integer(value)
  mcmcParams(object) <- mp
  object
})

#' Use values from SingleBatch model to simulate reasonable starting values for the SingleBatch-pooled and MultiBatch models of the same number of components
#'
setMethod("singleBatchGuided", c("MultiBatchList", "MultiBatch"),
          function(x, guide){
            modes(guide) <- computeModes(guide)
            if(any(k(x) != k(guide))){
              stop("models 'x' and 'guide' must have same number of components")
            }
            ns <- nStarts(x)
            mod.list <- vector("list", length(x))
            for(j in seq_along(x)){
              mod.list[[j]] <- replicate(ns, singleBatchGuided(x[[j]], guide))
            }
            names(mod.list) <- modelName(x)
            mod.list
          })

