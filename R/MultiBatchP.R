#' @include MultiBatch.R
NULL

setValidity("MultiBatchP", function(object){
  msg <- TRUE
  if(nrow(dataMean(object) > 0)){
    vals <- dataMean(object)[, 1]
    if(!is(vals, "numeric")){
      msg <- "data.mean is not numeric"
      return(msg)
    }
  }
  if(nrow(dataPrec(object)) > 0){
    vals <- dataPrec(object)[, 1]
    if(!is(vals, "numeric")){
      msg <- "data.prec is not numeric"
      return(msg)
    }
  }
  if(!identical(ncol(sigma2(object)), 1L)){
    msg <- "sigma2 matrix should have a single column"
    return(msg)
  }
  if(!identical(ncol(sigma2(chains(object))), nBatch(object))){
    msg <- "chains for sigma2 does not have the correct dimension"
    return(msg)
  }
  msg
})

mcmc_chainsP <- function(specs, parameters){
  B <- specs$number_batches
  K <- specs$k
  S <- iter(parameters$mp)
  N <- nStarts(parameters$mp)
  initialize_mcmcP(K, S, B)
}

modelValuesP <- function(specs, data, hp){
  vals <- modelValues2(specs, data, hp)
  s2 <- vals[["sigma2"]]
  if(is.null(dim(s2))){
    return(vals)
  }
  s2 <- apply(s2, 1, "mean")
  s2 <- matrix(s2, nrow=specs$number_batches, ncol=1)
  vals[["sigma2"]] <- s2
  vals
}

modelSummariesP <- function(specs){
  B <- specs$number_batches
  K <- specs$k
  data.mean <- matrix(as.numeric(NA), nrow=B, ncol=K)
  data.prec <- matrix(as.numeric(NA), nrow=B, ncol=1)
  zfreq <- integer(K)
  marginal_lik <- as.numeric(NA)
  modes <- list()
  mapping <- seq_len(K)
  list(data.mean=data.mean,
       data.prec=data.prec,
       zfreq=zfreq,
       marginal_lik=marginal_lik,
       modes=modes,
       mapping=mapping)
}

##
## Constructor
##
MultiBatchP <- function(model="MBP3",
                        data=modelData(),
                        specs=model_spec(model, data),
                        iter=1000L,
                        burnin=500L,
                        thin=1L,
                        nStarts=4L,
                        hp=Hyperparameters(k=specs$k),
                        mp=McmcParams(iter=iter, thin=thin,
                                      burnin=burnin,
                                      nStarts=nStarts),
                        parameters=modelParameters(mp=mp, hp=hp),
                        chains=mcmc_chainsP(specs, parameters),
                        current_values=modelValuesP(specs, data, hp),
                        summaries=modelSummariesP(specs),
                        flags=modelFlags()){
  ##
  ## When there are multiple batches in data, but the model specification is one of SB[X]
  ##
  is_SB <- substr(model, 1, 2) == "SB"
  if(nrow(data) > 0 && is_SB){
    data$batch <- 1L
  }
  model <- new("MultiBatchP",
               data=data,
               specs=specs,
               parameters=parameters,
               chains=chains,
               current_values=current_values,
               summaries=summaries,
               flags=flags)
  s <- summaries(model)
  s$modes <- current_values(model)
  summaries(model) <- s
  model
}

setMethod("computePrec", "MultiBatchP", function(object){
  z(object) <- map_z(object)
  tib <- assays(object) %>%
    mutate(z=z(object)) %>%
    group_by(batch) %>%
    summarize(prec=1/var(oned))
  matrix(tib$prec, nrow=nrow(tib))
})

## setMethod("downSampleModel", "MultiBatchP", function(object, N=1000, i){
##   if(!missing(N)){
##     if(N >= nrow(assays(object))){
##       return(object)
##     }
##   }
##   ## by sorting, batches are guaranteed to be ordered
##   if(missing(i)){
##     i <- sort(sample(seq_len(nrow(object)), N, replace=TRUE))
##   }
##   b <- assays(object)$batch[i]
##   current.vals <- current_values(object)
##   current.vals[["u"]] <- current.vals[["u"]][i]
##   current.vals[["z"]] <- current.vals[["z"]][i]
##   current.vals[["probz"]] <- current.vals[["probz"]][i, , drop=FALSE]
##   mb <- MultiBatchP(model=modelName(object),
##                     data=assays(object),
##                     ##down_sample=i,
##                     parameters=parameters(object),
##                     current_values=current.vals,
##                     chains=mcmc_chainsP( specs(object), parameters(object) ))
##   dataMean(mb) <- computeMeans(mb)
##   dataPrec(mb) <- computePrec(mb)
##   zFreq(mb) <- as.integer(table(z(mb)))
##   mb
## })

setAs("MultiBatch", "MultiBatchP", function(from){
  vals <- current_values(from)
  vals[["sigma2"]] <- matrix(rowMeans(vals[["sigma2"]]),
                             numBatch(from),
                             1)
  s <- summaries(from)
  m <- s$modes
  if(length(m) > 0){
    m[["sigma2"]] <- matrix(rowMeans(m[["sigma2"]]),
                            numBatch(from),
                            1)
    s$modes <- m
  }
  if(!all(is.na(s$data.prec))){
    s$data.prec <- matrix(rowMeans(s[["data.prec"]]),
                          numBatch(from),
                          1)
  } else {
    s$data.prec <- 1/m[["sigma2"]]
  }
  if(all(is.na(s$data.mean))){
    s$data.mean <- theta(from)
  }
  ch <- chains(from)
  sigma2(ch) <- matrix(as.numeric(NA),
                       iter(from),
                       numBatch(from))
  if(numBatch(from) > 1){
    specs(from)$model <- paste0("MBP", k(from))
  }else specs(from)$model <- paste0("SBP", k(from))
  model <- new("MultiBatchP",
               data=assays(from),
               specs=specs(from),
               parameters=parameters(from),
               chains=ch,
               current_values=vals,
               summaries=s,
               flags=flags(from))
})

setAs("MultiBatchModel", "MultiBatchP", function(from){
  mb <- as(from, "MultiBatch")
  mbp <- as(mb, "MultiBatchP")
  mbp
})

setAs("MultiBatchP", "MultiBatchPooled", function(from){
  flag1 <- as.integer(flags(from)[[".internal.constraint"]])
  flag2 <- as.integer(flags(from)[[".internal.counter"]])
  be <- as.integer(table(batch(from)))
  names(be) <- unique(batch(from))
  dat <- assays(from)
  th <- theta(from)
  KB <- nrow(th) * ncol(th)
  pred <- numeric(KB)
  ##pred <- predictive(from)
  zs <- integer(KB)
  obj <- new("MultiBatchPooled",
             k=k(from),
             hyperparams=hyperParams(from),
             theta=theta(from),
             sigma2=sigma2(from)[, 1],
             mu=mu(from),
             tau2=tau2(from),
             nu.0=nu.0(from),
             sigma2.0=sigma2.0(from),
             pi=p(from),
             data=dat$oned,
             data.mean=dataMean(from),
             data.prec=dataPrec(from)[, 1],
             predictive=pred,
             zstar=zs,
             z=z(from),
             u=u(from),
             zfreq=zFreq(from),
             probz=probz(from),
             logprior=logPrior(from),
             loglik=log_lik(from),
             mcmc.chains=chains(from),
             mcmc.params=mcmcParams(from),
             batch=batch(from),
             batchElements=be,
             label_switch=label_switch(from),
             marginal_lik=marginal_lik(from),
             .internal.constraint=flag1,
             .internal.counter=flag2)
  m <- modes(from)
  if(length(m) == 0){
    m <- computeModes(from)
  }
  ix <- match(c("nu.0", "p"), names(m))
  if(length(ix) == 2){
    names(m)[ix] <- c("nu0", "mixprob")
  }
  m$zfreq <- table(map_z(obj))
  modes(obj) <- m
  obj
})

setAs("MultiBatchPooled", "MultiBatchP", function(from){
  values <- extractValues(from)
  values[["sigma2"]] <- matrix(values[["sigma2"]], nBatch(from), 1)
  flags <- extractFlags(from)
  data <- extractData(from)
  params <- extractParameters(from)
  summaries <- extractSummaries(from)
  summaries[["data.prec"]] <- matrix(summaries[["data.prec"]],
                                     numBatch(from), 1)
  specs <- model_spec(modelName(from), data)
  modal.ordinates <- modes(from)
  mb <- MultiBatchP(data=data,
                    specs=specs,
                    parameters=params,
                    chains=chains(from),
                    current_values=values,
                    summaries=summaries,
                    flags=flags)
  if(length(modal.ordinates) > 0 ){
    ix <- match(c("nu0", "mixprob"), names(modal.ordinates))
    names(modal.ordinates)[ix] <- c("nu.0", "p")
    modal.ordinates$z <- map_z(from)
    modal.ordinates$probz <- probz(from)
    modal.ordinates$u <- u(from)
    modal.ordinates$sigma2 <- matrix(modes(from)[["sigma2"]],
                                     numBatch(from), 1)
    m <- modal.ordinates[names(current_values(mb))]
    modes(mb) <- m
  }
  mb
})


setAs("MultiBatchP", "list", function(from){
  ns <- nStarts(from)
  ##
  ## This initializes a list of models each with starting values simulated independently from the hyperparameters
  ##
  mp <- mcmcParams(from)
  nStarts(mp) <- 1L
  ## the starts will be independent by not transferring current_values
  mb.list <- replicate(ns, MultiBatchP(model=modelName(from),
                                       data=assays(from),
                                       mp=mp,
                                       hp=hyperParams(from),
                                       chains=chains(from),
                                       summaries=summaries(from),
                                       current_values=current_values(from)))
  mb.list
})


setMethod("compute_marginal_lik", "MultiBatchP", function(object, params){
  if(missing(params)){
    params <- mlParams(root=1/2,
                       reject.threshold=exp(-100),
                       prop.threshold=0.5,
                       prop.effective.size=0)
  }
  mbm <- as(object, "MultiBatchPooled")
  ml <- tryCatch(marginalLikelihood(mbm, params), warning=function(w) NULL, error=function(e) NULL)
  if(!is.null(ml)){
    summaries(object)[["marginal_lik"]] <- ml
    message("     marginal likelihood: ", round(ml, 2))
  } else {
    message("Unable to compute marginal likelihood")
  }
  object
})

setMethod("computeModes", "MultiBatchP", function(object){
  modes <- callNextMethod(object)
  if(iter(object) > 0){
    i <- argMax(object)[1]
    mc <- chains(object)
    B <- specs(object)$number_batches
    sigma2max <- matrix(sigma2(mc)[i, ], B, 1)
    modes[["sigma2"]] <- sigma2max
  }
  modes
})

setReplaceMethod("mcmcParams", c("MultiBatchP", "McmcParams"), function(object, value){
  it <- iter(object)
  if(it != iter(value)){
    if(iter(value) > iter(object)){
      parameters(object)[["mp"]] <- value
      ## create a new chain
      ch <- mcmc_chainsP(specs(object), parameters(object))
    } else {
      parameters(object)[["mp"]] <- value
      index <- seq_len(iter(value))
      ch <- chains(object)[index, ]
    }
    chains(object) <- ch
    return(object)
  }
  ## if we've got to this point, it must be safe to update mcmc.params
  ## (i.e., size of chains is not effected)
  parameters(object)[["mp"]] <- value
  object
})
