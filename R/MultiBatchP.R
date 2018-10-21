setClass("MultiBatchP", contains="MultiBatch")


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
  data.mean <- matrix(nrow=B, ncol=K)
  data.prec <- matrix(nrow=B, ncol=1)
  zfreq <- integer(K)
  marginal_lik <- as.numeric(NA)
  modes <- list()
  list(data.mean=data.mean,
       data.prec=data.prec,
       zfreq=zfreq,
       marginal_lik=marginal_lik,
       modes=modes)
}

##
## Constructor
##
MultiBatchP <- function(model="MBP3",
                        data=modelData(),
                        ## by default, assume no downsampling
                        down_sample=seq_len(nrow(data)),
                        specs=model_spec(model, data, down_sample),
                        iter=1000L,
                        thin=1L,
                        nStarts=4L,
                        hp=Hyperparameters(k=specs$k),
                        mp=McmcParams(iter=iter, thin=thin,
                                      nStarts=nStarts),
                        parameters=modelParameters(mp=mp, hp=hp),
                        chains=mcmc_chainsP(specs, parameters),
                        current_values=modelValuesP(specs, data[down_sample, ], hp),
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
               down_sample=down_sample,
               specs=specs,
               parameters=parameters,
               chains=chains,
               current_values=current_values,
               summaries=summaries,
               flags=flags)
  model
}

setMethod("computePrec", "MultiBatchP", function(object){
  z(object) <- map_z(object)
  tib <- downSampledData(object) %>%
    mutate(z=z(object)) %>%
    group_by(batch) %>%
    summarize(prec=1/var(oned))
  matrix(tib$prec, nrow=nrow(tib))
})

setMethod("downSampleModel", "MultiBatchP", function(object, N=1000, i){
  if(!missing(N)){
    if(N >= nrow(assays(object))){
      return(object)
    }
  }
  ## by sorting, batches are guaranteed to be ordered
  if(missing(i)){
    i <- sort(sample(seq_len(nrow(object)), N, replace=TRUE))
  }
  b <- assays(object)$batch[i]
  current.vals <- values(object)
  current.vals[["u"]] <- current.vals[["u"]][i]
  current.vals[["z"]] <- current.vals[["z"]][i]
  current.vals[["probz"]] <- current.vals[["probz"]][i, , drop=FALSE]
  mb <- MultiBatchP(model=modelName(object),
                    data=assays(object),
                    down_sample=i,
                    parameters=parameters(object),
                    current_values=current.vals,
                    chains=mcmc_chains( specs(object), parameters(object) ))
  dataMean(mb) <- computeMeans(mb)
  dataPrec(mb) <- computePrec(mb)
  zFreq(mb) <- as.integer(table(z(mb)))
  mb
})


setAs("MultiBatchP", "MultiBatchPooled", function(from){
  flag1 <- as.integer(flags(from)[[".internal.constraint"]])
  flag2 <- as.integer(flags(from)[[".internal.counter"]])
  be <- as.integer(table(batch(from)))
  names(be) <- unique(batch(from))
  dat <- downSampledData(from)
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
             data.prec=dataPrec(from),
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
})

setAs("MultiBatchPooled", "MultiBatchP", function(from){
  values <- extractValues(from)
  values[["sigma2"]] <- matrix(values[["sigma2"]], nBatch(from), 1)
  flags <- extractFlags(from)
  data <- extractData(from)
  params <- extractParameters(from)
  summaries <- extractSummaries(from)
  summaries[["data.prec"]] <- matrix(summaries[["data.prec"]],
                                     nBatch(from), 1)
  down_sample <- seq_len(nrow(data))
  specs <- model_spec(modelName(from), data, down_sample)
  modal.ordinates <- modes(from)
  mb <- MultiBatchP(data=data,
                    ## By default, assume no downsampling
                    down_sample=down_sample,
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
    m <- modal.ordinates[names(values(mb))]
    modes(mb) <- m
  }
  mb
})
