#' @include AllGenerics.R
NULL

#' An object for running MCMC simulations.
#'
#' BatchModel and MarginalModel both inherit from this class.
#' @slot data a tibble with one-dimensional summaries (oned), id, and batch
#' @slot parameters list of parameters
#' @slot chains object of class McmcChains
#' @slot current_values current value of each chain
#' @slot summaries list of empirical data and model summaries
#' @slot flags list of model flags
#' @slot down_sample integer vector for downsampling the full data
#' @export
setClass("MultiBatch", representation(data="tbl_df",
                                      down_sample="integer",
                                      specs="tbl_df",
                                      parameters="list",
                                      chains="McmcChains",
                                      current_values="list",
                                      summaries="list",
                                      flags="list"))

setValidity("MultiBatch", function(object){
  msg <- TRUE
  if(length(u(object)) != length(down_sample(object))){
    msg <- "Incorrect length of u-vector"
    return(msg)
  }
  nb <- length(unique(batch(object)))
  if(nb != specs(object)$number_batches ){
    msg <- "Number of batches in model specs differs from number of batches in data"
    return(msg)
  }
  nr <- nrow(theta(object))
  if(nr != nb){
    msg <- "Number of batches in current values differs from number of batches in model specs"
    return(msg)
  }
  nr <- nrow(dataMean(object))
  if(nr != nb){
    msg <- "Number of batches in model summaries differs from number of batches in model specs"
    return(msg)
  }
  S <- iter(object)
  th <- theta(chains(object))
  if( S != nrow(th) || ncol(th) != nr * k(object) ) {
    msg <- "Dimension of chain parameters is not consistent with model specs"
    return(msg)
  }
  B <- batch(object)
  is.sorted <- all(diff(B) >= 0)
  if(!is.sorted) {
    msg <- "down sample index must be in batch order"
    return(msg)
  }
  msg
})


model_spec <- function(model, data, down_sample) {
  if(missing(down_sample)) down_sample <- seq_len(nrow(data))
  models <- c("SB", "SBP", "MB", "MBP")
  K <- 1:5
  avail.models <- lapply(models, paste0, K) %>%
    unlist
  if(missing(model)) model <- "MB3" else model <- model[1]
  if(!model %in% avail.models) stop("Model not recognized")
  number_batches <- length(unique(data$batch))
  k <- substr(model, nchar(model), nchar(model)) %>%
    as.integer
  is_SB <- substr(model, 1, 2) == "SB"
  if(nrow(data) > 0){
    number_batches <- ifelse(is_SB, 1, number_batches)
  }
  number_obs <- nrow(data)
  tab <- tibble(model=model,
                k=k,
                number_batches=number_batches,
                number_obs=number_obs,
                number_sampled=length(down_sample))
  tab
}

listChains1 <- function(model_specs, parameters){
  S <- iter(parameters[["mp"]])
  num.chains <- nStarts(parameters[["mp"]])
  B <- model_specs$number_batches
  K <- model_specs$k
  mc.list <- replicate(num.chains, initialize_mcmc(K, S, B))
  mc.list
}

mcmc_chains <- function(specs, parameters){
  B <- specs$number_batches
  K <- specs$k
  S <- iter(parameters$mp)
  N <- nStarts(parameters$mp)
  initialize_mcmc(K, S, B)
}

##
## Constructor
##
MultiBatch <- function(model="MB3",
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
                       chains=mcmc_chains(specs, parameters),
                       current_values=modelValues2(specs, data[down_sample, ], hp),
                       summaries=modelSummaries(specs),
                       flags=modelFlags()){

  model <- new("MultiBatch",
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

setMethod("show", "MultiBatch", function(object){
  ##callNextMethod()
  cls <- class(object)
  n.mcmc <- iter(object) * thin(object)
  saved.mcmc <- iter(object)
  cat(paste0("An object of class ", cls), "\n")
  cat("     n. obs              :", nrow(assays(object)), "\n")
  cat("     n. sampled          :", nrow(downSampledData(object)), "\n")
  cat("     n. batches          :", nBatch(object), "\n")
  cat("     k                   :", k(object), "\n")
  cat("     nobs/batch          :", table(batch(object)), "\n")
  cat("     saved mcmc          :", saved.mcmc, "\n")
  cat("     log lik (s)         :", round(log_lik(object), 1), "\n")
  cat("     log prior (s)       :", round(logPrior(object), 1), "\n")
  cat("     log marginal lik (s):", round(marginal_lik(object), 1), "\n")
})

setClass("MultiBatchPooled2", contains="MultiBatch")

setClass("GenotypeModel", contains="MultiBatch",
         representation(mapping="character"))

setClass("GenotypePooled", contains="MultiBatch",
         representation(mapping="character"))


setGeneric("specs", function(object) standardGeneric("specs"))
setGeneric("specs<-", function(object, value) standardGeneric("specs<-"))

setMethod("specs", "MultiBatch", function(object) object@specs)

harmonizeDimensions <- function(object){
  L <- length(unique(batch(object)))
  spec <- specs(object)
  if(L != spec$number_batches){
    spec$number_batches <- L
    object@specs <- spec
  }
  nr <- nrow(theta(object))
  if(L != nr){
    values(object) <- modelValues2( spec, downSampledData(object), hyperParams(object) )
  }
  ncols1 <- k( object ) * L
  ncols2 <- ncol(theta(chains(object)))
  if( ncols1 != ncols2 ){
    chains(object) <- mcmc_chains( spec, parameters(object) )
  }
  if( nrow(dataMean(object)) != L){
    summaries(object) <- modelSummaries( spec )
  }
  object
}

setReplaceMethod("specs", "MultiBatch", function(object, value){
  object@specs <- value
  object <- harmonizeDimensions(object)
  object
})

setReplaceMethod("chains", c("MultiBatch", "McmcChains"), function(object, value){
  object@chains <- value
  object
})

modelParameters <- function(hp=Hyperparameters(),
                            mp=McmcParams()){
  list(hp=hp, mp=mp)
}

extractParameters <- function(old){
  list(hp=hyperParams(old),
       mp=mcmcParams(old))
}

modelValues2 <- function(specs, data, hp){
  n.sampled <- specs$number_sampled
  B <- specs$number_batches
  K <- specs$k

  mu <- sort(rnorm(k(hp), mu.0(hp), sqrt(tau2.0(hp))))
  tau2 <- 1/rgamma(k(hp), 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
  p <- rdirichlet(1, alpha(hp))[1, ]
  sim_theta <- function(mu, tau, B) sort(rnorm(B, mu, tau))
  . <- NULL
  theta <- map2(mu, sqrt(tau2), sim_theta, B) %>%
    do.call(cbind, .) %>%
    apply(., 1, sort) %>%
    t
  if(K == 1) theta <- t(theta)
  nu.0 <- 3.5
  sigma2.0 <- 0.25
  sigma2 <- 1/rgamma(k(hp) * B, 0.5 * nu.0, 0.5 * nu.0 * sigma2.0) %>%
    matrix(B, K)
  u <- rchisq(n.sampled, hp@dfr)
  z <- sample(seq_len(k(hp)), prob=p, replace=TRUE, size=n.sampled)
  logprior <- numeric()
  loglik <- numeric()
  probz <- matrix(nrow=n.sampled, ncol=K)
  list(theta=theta,
       sigma2=sigma2,
       nu.0=nu.0,
       sigma2.0=sigma2.0,
       p=p,
       mu=mu,
       tau2=tau2,
       u=u,
       z=z,
       logprior=logprior,
       loglik=loglik,
       probz=probz)
}

modelValues <- function(theta,
                        sigma2,
                        nu.0,
                        sigma2.0,
                        p,
                        mu,
                        tau2,
                        u,
                        z,
                        logprior,
                        loglik,
                        probz){
  list(theta=theta,
       sigma2=sigma2,
       nu.0=nu.0,
       sigma2.0=sigma2.0,
       p=p,
       mu=mu,
       tau2=tau2,
       u=u,
       z=z,
       logprior=logprior,
       loglik=loglik,
       probz=probz)
}

extractValues <- function(old){
  modelValues(theta=theta(old),
              sigma2=sigma2(old),
              nu.0=nu.0(old),
              sigma2.0=sigma2.0(old),
              p=p(old),
              mu=mu(old),
              tau2=tau2(old),
              u=u(old),
              z=z(old),
              logprior=logPrior(old),
              loglik=log_lik(old),
              probz=probz(old))
}

modelSummaries <- function(specs){
  B <- specs$number_batches
  K <- specs$k
  data.mean <- matrix(nrow=B, ncol=K)
  data.prec <- matrix(nrow=B, ncol=K)
  zfreq <- integer(K)
  marginal_lik <- as.numeric(NA)
  modes <- list()
  list(data.mean=data.mean,
       data.prec=data.prec,
       zfreq=zfreq,
       marginal_lik=marginal_lik,
       modes=modes)
}

extractSummaries <- function(old){
  list(data.mean=dataMean(old),
       data.prec=dataPrec(old),
       zfreq=zFreq(old),
       marginal_lik=marginal_lik(old),
       modes=modes(old))
}

modelFlags <- function(.internal.constraint=5e-4,
                       .internal.counter=0L,
                       label_switch=FALSE){
  list(.internal.constraint=.internal.constraint,
       .internal.counter=.internal.counter,
       label_switch=label_switch)
}

extractFlags <- function(old){
  list(.internal.constraint=old@.internal.constraint,
       .internal.counter=old@.internal.counter,
       label_switch=label_switch(old))
}

modelData <- function(id=character(),
                      oned=numeric(),
                      batch=integer()){
  tibble(id=id,
         oned=oned,
         batch=batch)
}

extractData <- function(old){
  tibble(id=as.character(seq_along(y(old))),
         oned=y(old),
         batch=batch(old))
}


##
## Coersion
##
setAs("MultiBatchModel", "MultiBatch", function(from){
  values <- extractValues(from)
  flags <- extractFlags(from)
  data <- extractData(from)
  params <- extractParameters(from)
  summaries <- extractSummaries(from)
  down_sample <- seq_len(nrow(data))
  specs <- model_spec(modelName(from), data, down_sample)
  mb <- MultiBatch(data=data,
                   ## By default, assume no downsampling
                   down_sample=down_sample,
                   specs=specs,
                   parameters=params,
                   chains=chains(from),
                   current_values=values,
                   summaries=summaries,
                   flags=flags)
  mb
})

setAs("MultiBatch", "MultiBatchModel", function(from){
  flag1 <- as.integer(flags(from)[[".internal.constraint"]])
  flag2 <- as.integer(flags(from)[[".internal.counter"]])
  be <- as.integer(table(batch(from)))
  names(be) <- unique(batch(from))
  dat <- downSampledData(from)
  obj <- new("MultiBatchModel",
             k=k(from),
             hyperparams=hyperParams(from),
             theta=theta(from),
             sigma2=sigma2(from),
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

setAs("MultiBatch", "list", function(from){
  ns <- nStarts(from)
  mb.list <- replicate(ns, as(from, "MultiBatchModel"))
  mb.list <- lapply(mb.list, function(x) {nStarts(x) <- 1; return(x)})
  mb.list
})

##setAs("list", "MultiBatch", function(from){
##  mb.list
##})




##
## Accessors
##
setMethod("values", "MultiBatch", function(x, ...){
  x@current_values
})

setMethod("theta", "MultiBatch", function(object){
  th <- values(object)[["theta"]]
  th
})



setReplaceMethod("values", c("MultiBatch", "list"),
                 function(object, value){
  object@current_values <- value
  object
})

setReplaceMethod("theta", c("MultiBatch", "matrix"),
                 function(object, value){
                   values(object)[["theta"]] <- value
                   object
                 })

setMethod("sigma2", "MultiBatch", function(object){
  values(object)[["sigma2"]]
})

setReplaceMethod("sigma2", c("MultiBatch", "matrix"),
                 function(object, value){
                   values(object)[["sigma2"]] <- value
                   object
                 })

setMethod("sigma", "MultiBatch", function(object){
  sqrt(sigma2(object))
})

setReplaceMethod("sigma", c("MultiBatch", "matrix"),
                 function(object, value){
                   sigma2(object) <- value^2
                   object
                 })

setMethod("p", "MultiBatch", function(object){
  values(object)[["p"]]
})

setReplaceMethod("p", c("MultiBatch", "numeric"),
                 function(object, value){
                   values(object)[["p"]] <- value
                   object
                 })

setMethod("nu.0", "MultiBatch", function(object){
  values(object)[["nu.0"]]
})

setReplaceMethod("nu.0", c("MultiBatch", "numeric"),
                 function(object, value){
                   values(object)[["nu.0"]] <- value
                   object
                 })

setMethod("sigma2.0", "MultiBatch", function(object){
  values(object)[["sigma2.0"]]
})

setReplaceMethod("sigma2.0", c("MultiBatch", "numeric"),
                 function(object, value){
                   values(object)[["sigma2.0"]] <- value
                   object
                 })

setMethod("mu", "MultiBatch", function(object){
  values(object)[["mu"]]
})

setReplaceMethod("mu", c("MultiBatch", "numeric"),
                 function(object, value){
                   values(object)[["mu"]] <- value
                   object
                 })

setMethod("tau2", "MultiBatch", function(object){
  values(object)[["tau2"]]
})

setReplaceMethod("tau2", c("MultiBatch", "numeric"),
                 function(object, value){
                   values(object)[["tau2"]] <- value
                   object
                 })

setMethod("log_lik", "MultiBatch", function(object){
  values(object)[["loglik"]]
})

setReplaceMethod("log_lik", c("MultiBatch", "numeric"),
                 function(object, value){
                   values(object)[["loglik"]] <- value
                   object
                 })

setMethod("logPrior", "MultiBatch", function(object){
  values(object)[["logprior"]]
})

setReplaceMethod("logPrior", c("MultiBatch", "numeric"),
                 function(object, value){
                   values(object)[["logprior"]] <- value
                   object
                 })

setMethod("probz", "MultiBatch", function(object){
  values(object)[["probz"]]
})

zProb <- function(object){
  probz(object)/(iter(object)-1)
}

setReplaceMethod("probz", c("MultiBatch", "matrix"),
                 function(object, value){
                   values(object)[["probz"]] <- value
                   object
                 })

setMethod("z", "MultiBatch", function(object){
  values(object)[["z"]]
})

setReplaceMethod("z", c("MultiBatch", "numeric"),
                 function(object, value){
                   values(object)[["z"]] <- value
                   object
                 })


setMethod("u", "MultiBatch", function(object){
  values(object)[["u"]]
})

setMethod("nrow", "MultiBatch", function(x) nrow(assays(x)))

setReplaceMethod("u", c("MultiBatch", "numeric"),
                 function(object, value){
                   values(object)[["u"]] <- value
                   object
                 })

setMethod("parameters", "MultiBatch", function(object){
  object@parameters
})

setReplaceMethod("parameters", c("MultiBatch", "list"), function(object, value){
  object@parameters <- value
  object
})

setMethod("mcmcParams", "MultiBatch", function(object){
  parameters(object)[["mp"]]
})

setMethod("chains", "MultiBatch", function(object){
  object@chains
})

setReplaceMethod("chains", c("MultiBatch", "list"), function(object, value){
  object@chains <- value
  object
})

setReplaceMethod("mcmcParams", c("MultiBatch", "McmcParams"), function(object, value){
  it <- iter(object)
  if(it != iter(value)){
    if(iter(value) > iter(object)){
      parameters(object)[["mp"]] <- value
      ## create a new chain
      mcmc_chains <- chains(object)
    } else {
      parameters(object)[["mp"]] <- value
      index <- seq_len(iter(value))
      mcmc_chains <- chains(object)[index, ]
    }
    chains(object) <- mcmc_chains
    return(object)
  }
  ## if we've got to this point, it must be safe to update mcmc.params
  ## (i.e., size of chains is not effected)
  parameters(object)[["mp"]] <- value
  object
})

setMethod("hyperParams", "MultiBatch", function(object){
  parameters(object)[["hp"]]
})

setReplaceMethod("hyperParams", c("MultiBatch", "Hyperparameters"),
                 function(object, value){
                   parameters(object)[["hp"]] <- value
                   object
                 })

setMethod("burnin", "MultiBatch", function(object){
  burnin(mcmcParams(object))
})

setReplaceMethod("burnin", c("MultiBatch", "numeric"),
                 function(object, value){
                   burnin(mcmcParams(object)) <- as.numeric(value)
                   object
})

setMethod("iter", "MultiBatch", function(object){
  iter(mcmcParams(object))
})

setMethod("k", "MultiBatch", function(object) k(hyperParams(object)))

setReplaceMethod("iter", c("MultiBatch", "numeric"),
                 function(object, value){
                   mp <- mcmcParams(object)
                   iter(mp) <- value
                   mcmcParams(object) <- mp
                   object
                 })

setMethod("thin", "MultiBatch", function(object){
  thin(mcmcParams(object))
})

setReplaceMethod("thin", c("MultiBatch", "numeric"), function(object, value){
  thin(mcmcParams(object)) <- as.integer(value)
  object
})

setMethod("nStarts", "MultiBatch", function(object) nStarts(mcmcParams(object)))

setReplaceMethod("nStarts", c("MultiBatch", "numeric"),
                 function(object, value){
                   nStarts(object) <- as.integer(value)
                   object
                 })

setReplaceMethod("nStarts", c("MultiBatch", "integer"),
                 function(object, value){
                   nStarts(mcmcParams(object)) <- value
                   object
                 })

setMethod("flags", "MultiBatch", function(object) object@flags)

setReplaceMethod("flags", c("MultiBatch", "list"), function(object, value){
  object@flags <- value
  object
})

setMethod("label_switch", "MultiBatch", function(object){
  flags(object)[["label_switch"]]
})

setReplaceMethod("label_switch", c("MultiBatch", "logical"),
                 function(object, value){
                   flags(object)[["label_switch"]] <- value
                   object
                 })


setMethod("assays", "MultiBatch", function(x, ..., withDimnames) x@data)

setReplaceMethod("assays", c("MultiBatch", "tbl_df"), function(x, value){
  x@data <- value
  x <- harmonizeDimensions(x)
  x
})

setGeneric("downSampledData<-", function(x, value) standardGeneric("downSampledData<-"))

setReplaceMethod("downSampledData", c("MultiBatch", "tbl_df"), function(x, value){
  x@data <- value
  x
})

##
## Data accessors
##
setMethod("batch", "MultiBatch", function(object){
  downSampledData(object)[["batch"]]
})

setReplaceMethod("oned", c("MultiBatch", "numeric"), function(object, value){
  downSampledData(object)[["oned"]] <- value
  object
})

setMethod("oned", "MultiBatch", function(object){
  downSampledData(object)[["oned"]]
})

setMethod("zFreq", "MultiBatch", function(object){
  summaries(object)[["zfreq"]]
})

setReplaceMethod("zFreq", "MultiBatch", function(object, value){
  summaries(object)[["zfreq"]] <- value
  object
})



setReplaceMethod("batch", c("MultiBatch", "numeric"), function(object, value){
  downSampledData(object)[["batch"]] <- as.integer(value)
  L <- length(unique(batch(object)))
  if( L != specs(object)$number_batches ){
    spec(object)$number_batches <- L
    chains(object) <- mcmc_chains( specs(object), parameters(object) )
  }
  object
})

##
## Summaries
##
setMethod("summaries", "MultiBatch", function(object){
  object@summaries
})

setReplaceMethod("summaries", c("MultiBatch", "list"), function(object, value){
  object@summaries <- value
  object
})

setMethod("dataMean", "MultiBatch", function(object){
  summaries(object)[["data.mean"]]
})

setReplaceMethod("dataMean", c("MultiBatch", "matrix"), function(object, value){
  summaries(object)[["data.mean"]] <- value
  object
})

setMethod("dataPrec", "MultiBatch", function(object){
  summaries(object)[["data.prec"]]
})

setReplaceMethod("dataPrec", c("MultiBatch", "matrix"), function(object, value){
  summaries(object)[["data.prec"]] <- value
  object
})

setMethod("marginal_lik", "MultiBatch", function(object){
  summaries(object)[["marginal_lik"]]
})


setReplaceMethod("marginal_lik", c("MultiBatch", "numeric"), function(object, value){
  summaries(object)[["marginal_lik"]] <- value
  object
})

setMethod("tablez", "MultiBatch", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})

setMethod("showSigmas", "MultiBatch", function(object){
  sigmas <- round(sqrt(sigma2(object)), 2)
  sigmas <- c("\n", paste0(t(cbind(sigmas, "\n")), collapse="\t"))
  sigmas <- paste0("\t", sigmas[2])
  sigmas <- paste0("\n", sigmas[1])
  sigmas
})

setMethod("showMeans", "MultiBatch", function(object){
  thetas <- round(theta(object), 2)
  mns <- c("\n", paste0(t(cbind(thetas, "\n")), collapse="\t"))
  mns <- paste0("\t", mns[2])
  mns <- paste0("\n", mns[1])
  mns
})

##
## Summary statistics
##
setMethod("computeModes", "MultiBatch", function(object){
  i <- argMax(object)[1]
  mc <- chains(object)
  B <- specs(object)$number_batches
  K <- k(object)
  thetamax <- matrix(theta(mc)[i, ], B, K)
  sigma2max <- matrix(sigma2(mc)[i, ], B, K)
  pmax <- p(mc)[i, ]
  mumax <- mu(mc)[i, ]
  tau2max <- tau2(mc)[i,]
  ##
  ## We do not store u.  Just use current u.
  ##
  currentu <- u(object)
  ## We do not store z.  Use map_z
  zz <- map_z(object)
  pz <- probz(object)
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i],
                p=pmax,
                mu=mumax,
                tau2=tau2max,
                u=currentu,
                z=zz,
                logprior=logPrior(mc)[i],
                loglik=log_lik(mc)[i],
                probz=pz)
})

setReplaceMethod("modes", "MultiBatch", function(object, value){
  summaries(object)[["modes"]] <- value
  object
})

setMethod("computeMeans", "MultiBatch", function(object){
  z(object) <- map_z(object)
  tib <- downSampledData(object) %>%
    mutate(z=z(object)) %>%
    group_by(batch, z) %>%
    summarize(mean=mean(oned))
  nr <- length(unique(tib$batch))
  nc <- length(unique(tib$z))
  m <- spread(tib, z, mean) %>%
    ungroup %>%
    select(-batch) %>%
    as.matrix
  dimnames(m) <- NULL
  m
})

setMethod("computePrec", "MultiBatch", function(object){
  z(object) <- map_z(object)
  tib <- downSampledData(object) %>%
    mutate(z=z(object)) %>%
    group_by(batch, z) %>%
    summarize(prec=1/var(oned))
  nr <- length(unique(tib$batch))
  nc <- length(unique(tib$z))
  m <- spread(tib, z, prec) %>%
    ungroup %>%
    select(-batch) %>%
    as.matrix
  dimnames(m) <- NULL
  m
})

summarizeModel <- function(object){
  stats <- list(modes=computeModes(object),
                data.mean=computeMeans(object),
                data.prec=computePrec(object),
                zfreq=as.integer(table(z(object))),
                marginal_lik=marginal_lik(object))
}

collectFlags <- function(model.list){
  getConstraint <- function(model) model@.internal.constraint
  getCounter <- function(model) model@.internal.counter
  nlabel_swap <- sum(map_lgl(model.list, label_switch))
  n.internal.constraint <- sum(map_dbl(model.list, getConstraint))
  n.internal.counter <- sum(map_dbl(model.list, getCounter))
  flags <- list(label_switch=nlabel_swap > 0,
                .internal.constraint=n.internal.constraint,
                .internal.counter=n.internal.counter)
}

setGeneric("setModes", function(object) standardGeneric("setModes"))

setMethod("setModes", "MultiBatch", function(object){
  modal.ordinates <- modes(object)
})

setMethod("modes", "MultiBatch", function(object) summaries(object)[["modes"]])

combineChains <- function(model.list){
  ch.list <- map(model.list, chains)
  th <- map(ch.list, theta) %>% do.call(rbind, .)
  s2 <- map(ch.list, sigma2) %>% do.call(rbind, .)
  ll <- map(ch.list, log_lik) %>% unlist
  pp <- map(ch.list, p) %>% do.call(rbind, .)
  n0 <- map(ch.list, nu.0) %>% unlist
  s2.0 <- map(ch.list, sigma2.0) %>% unlist
  logp <- map(ch.list, logPrior) %>% unlist
  .mu <- map(ch.list, mu) %>% do.call(rbind, .)
  .tau2 <- map(ch.list, tau2) %>% do.call(rbind, .)
  zfreq <- map(ch.list, zFreq) %>% do.call(rbind, .)
  mc <- new("McmcChains",
            theta=th,
            sigma2=s2,
            pi=pp,
            mu=.mu,
            tau2=.tau2,
            nu.0=n0,
            sigma2.0=s2.0,
            zfreq=zfreq,
            logprior=logp,
            loglik=ll)
  mc
}

## update the current values with the posterior means across all chains
combineModels <- function(model.list){
  mc <- combineChains(model.list)
  hp <- hyperParams(model.list[[1]])
  mp <- mcmcParams(model.list[[1]])
  nStarts(mp) <- length(model.list)
  tib <- tibble(oned=y(model.list[[1]])) %>%
    mutate(id=seq_along(oned),
           id=as.character(id),
           batch=batch(model.list[[1]]))
  param.list <- list(mp=mp, hp=hp)
  mb <- MultiBatch(modelName(model.list[[1]]),
                   data=tib,
                   parameters=param.list,
                   chains=mc)
  summaries(mb) <- summarizeModel(mb)
  flags(mb) <- collectFlags(model.list)
  ##
  ## Set current values to the modal ordinates
  ##  (except for u which is not stored)
  ##  ( for z, we use the map estimate)
  values(mb) <- computeModes(mb)
  mb
}


##
## Markov chain Monte Carlo
##
setGeneric("mcmc2", function(object) standardGeneric("mcmc2"))


setMethod("isOrdered", "MultiBatch", function(object){
  thetas <- theta(object)
  checkOrder <- function(theta) identical(order(theta), seq_along(theta))
  is_ordered <- apply(thetas, 1, checkOrder)
  all(is_ordered)
})

setMethod("sortComponentLabels", "MultiBatch", function(model){
  is_ordered <- isOrdered(model)
  if(is_ordered) return(model)
  ## thetas are not all ordered
  thetas <- theta(model)
  s2s <- sigma2(model)
  K <- k(model)
  ix <- order(thetas[1, ])
  B <- specs(model)$number_batches
  . <- NULL
  tab <- tibble(z_orig=z(model),
                z=z(model),
                batch=batch(model)) %>%
    mutate(index=seq_len(nrow(.)))
  z_relabel <- NULL
  for(i in seq_len(B)){
    ix.next <- order(thetas[i, ])
    thetas[i, ] <- thetas[i, ix.next]
    s2s[i, ] <- s2s[i, ix]
    index <- which(tab$batch == i)
    tab2 <- tab[index, ] %>%
      mutate(z_relabel=factor(z, levels=ix.next)) %>%
      mutate(z_relabel=as.integer(z_relabel))
    tab$z[index] <- tab2$z_relabel
  }
  ps <- p(model)[ix]
  mu(model) <- mu(model)[ix]
  tau2(model) <- tau2(model)[ix]
  sigma2(model) <- s2s
  theta(model) <- thetas
  p(model) <- ps
  z(model) <- tab$z
  model
})

setMethod("posteriorSimulation", "list", function(object, k){
  ind.starts <- object
  for(i in seq_along(ind.starts)){
    ind.starts[[i]] <- posteriorSimulation(ind.starts[[i]])
  }
  no_label_swap <- !map_lgl(ind.starts, label_switch)
  ind.starts <- ind.starts[ no_label_swap ]
  ind.starts <- ind.starts[ selectModels(ind.starts) ]
  mlist <- mcmcList(ind.starts)
  neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
  if(is.null(neff)){
    neff <- 0
  } else {
    ## remove parameters that were not updated
    neff <- neff[ neff > 0 ]
  }
  hp <- hyperParams(ind.starts[[1]])
  mp <- mcmcParams(ind.starts[[1]])
  gr <- gelman_rubin(mlist, hp)
  message("     Gelman-Rubin: ", round(gr$mpsrf, 2))
  message("     eff size (median): ", round(median(neff), 1))
  message("     eff size (mean): ", round(mean(neff), 1))
  if((mean(neff) > min_effsize(mp)) && r$mpsrf < min_GR(mp)) {
    convergence.flag <- FALSE
  } else {
    convergence.flag <- TRUE
  }
  mb <- combineModels(ind.starts)
  summaries(mb)[["GR"]] <- gr
  summaries(mb)[["effsize"]] <- mean(neff)
  flags(mb)[["convergence"]] <- convergence.flag
  ##
  ## Compute summary statistics
  ##
  mb
})

continueMcmc <- function(mp){
  burnin(mp) <- as.integer(burnin(mp) * 2)
  mp@thin <- as.integer(thin(mp) + 2)
  nStarts(mp) <- nStarts(mp) + 1
  mp
}

setMethod("mcmc2", "MultiBatch", function(object){
  mb <- object
  mp <- mcmcParams(mb)
  if(iter(mp) < 500)
    warning("Very few Monte Carlo simulations specified")
  maxb <- max_burnin(mp)
  while(burnin(mp) < maxb && thin(mp) < 100){
    message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
    mcmcParams(object) <- mp
    mod.list <- as(mb, "list")
    mb <- posteriorSimulation(mod.list)
    if(! flags(mb)[["convergence"]]) break()
    mp <- continueMcmc(mp)
  }
  if(! flags(mb)[["convergence"]]) {
    mb <- compute_marginal_lik(mb)
  }
  mb
})

setGeneric("down_sample", function(object) standardGeneric("down_sample"))
setGeneric("down_sample<-", function(object, value) standardGeneric("down_sample<-"))

setGeneric("downSampledData", function(object) standardGeneric("downSampledData"))

setMethod("down_sample", "MultiBatch", function(object) object@down_sample)

setMethod("downSampledData", "MultiBatch", function(object){
  assays(object)[down_sample(object), ]
})

setReplaceMethod("down_sample", "MultiBatch", function(object, value){
  object@down_sample <- value
  object
})

downSampleModel <- function(object, N=1000){
  if(N >= nrow(assays(object))){
    return(object)
  }
  ## by sorting, batches are guaranteed to be ordered
  ix <- sort(sample(seq_len(nrow(object)), N, replace=TRUE))
  b <- assays(object)$batch[ix]
  current.vals <- values(object)
  current.vals[["u"]] <- current.vals[["u"]][ix]
  current.vals[["z"]] <- current.vals[["z"]][ix]
  current.vals[["probz"]] <- current.vals[["probz"]][ix, , drop=FALSE]
  mb <- MultiBatch(model=modelName(object),
                   data=assays(object),
                   down_sample=ix,
                   parameters=parameters(object),
                   current_values=current.vals,
                   chains=mcmc_chains( specs(object), parameters(object) ))
  dataMean(mb) <- computeMeans(mb)
  dataPrec(mb) <- computePrec(mb)
  zFreq(mb) <- as.integer(table(z(mb)))
  mb
}
