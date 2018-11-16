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
  if(iter(object) != iter(chains(object))){
    msg <- "Number of iterations not the same between MultiBatch model and chains"
    return(msg)
  }
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
##  if(nrow(assays(object)) != length(down_sample(object))){
##    msg <- "down sample index must be the same length as the number of rows in assay"
##    return(msg)
##  }
  if(!identical(dim(zstar(object)), dim(predictive(object)))){
    msg <- "z* and predictive matrices in MCMC chains should be the same dimension"
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
                number_batches=as.integer(number_batches),
                number_obs=as.integer(number_obs),
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
  ##
  ## When there are multiple batches in data, but the model specification is one of SB[X]
  ##
  is_SB <- substr(model, 1, 2) == "SB"
  if(nrow(data) > 0 && is_SB){
    data$batch <- 1L
  }
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
  ml <- marginal_lik(object)
  include_ml <- length(ml) > 0 && !is.na(ml)
  ##cat(paste0("An object of class ", cls), "\n")
  cat("Model name:", modelName(object), "\n")
  cat("   n. obs              :", nrow(assays(object)), "\n")
  cat("   n. sampled          :", nrow(downSampledData(object)), "\n")
  cat("   n. batches          :", nBatch(object), "\n")
  cat("   k                   :", k(object), "\n")
  cat("   nobs/batch          :", table(batch(object)), "\n")
  cat("   saved mcmc          :", saved.mcmc, "\n")
  cat("   log lik (s)         :", round(log_lik(object), 1), "\n")
  ##cat("     log prior (s)       :", round(logPrior(object), 1), "\n")
  if(include_ml)
    cat("   log marginal lik    :", round(ml, 1), "\n")
})

setMethod("numBatch", "MultiBatch", function(object) as.integer(specs(object)$number_batches[[1]]))

setClass("MultiBatchPooled2", contains="MultiBatch")

setClass("GenotypeModel", contains="MultiBatch",
         representation(mapping="character"))

setClass("GenotypePooled", contains="MultiBatch",
         representation(mapping="character"))



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
    current_values(object) <- modelValues2( spec, downSampledData(object), hyperParams(object) )
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
  if(nrow(data) == 0){
    K <- specs$k
    vals <- list(theta=matrix(nrow=0, ncol=K),
                 sigma2=matrix(nrow=0, ncol=K),
                 nu.0=numeric(),
                 sigma2.0=numeric(),
                 p=numeric(K),
                 mu=numeric(K),
                 tau2=numeric(K),
                 u=numeric(),
                 z=numeric(),
                 logprior=numeric(),
                 loglik=numeric(),
                 probz=matrix(nrow=0, ncol=K))
    return(vals)
  }
  n.sampled <- specs$number_sampled
  B <- specs$number_batches
  K <- specs$k
  alpha(hp) <- rep(1, K)
  mu <- sort(rnorm(K, mu.0(hp), 3))
  tau2 <- 1/rgamma(K, 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
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
  sigma2 <- 1/rgamma(K * B, 0.5 * nu.0, 0.5 * nu.0 * sigma2.0) %>%
    matrix(B, K)
  u <- rchisq(n.sampled, dfr(hp))
  z <- sample(seq_len(K), prob=p, replace=TRUE, size=n.sampled)
  logprior <- numeric()
  loglik <- numeric()
  probz <- matrix(0L, nrow=n.sampled, ncol=K)
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
  s.list <- list(data.mean=dataMean(old),
                 data.prec=dataPrec(old),
                 zfreq=zFreq(old),
                 marginal_lik=marginal_lik(old),
                 modes=modes(old))
  x <- s.list[["data.mean"]]
  if(!is.numeric(x[, 1])){
    x <- matrix(as.numeric(x),
                nrow(x),
                ncol(x))
    ix <- order(x[1, ], decreasing=FALSE)
    x <- x[, ix, drop=FALSE]
    s.list[["data.mean"]] <- x
    x <- s.list[["data.prec"]]
    x <- matrix(as.numeric(x),
                nrow(x),
                ncol(x))
    if(ncol(x) > 1) x <- x[, ix, drop=FALSE]
    s.list[["data.prec"]] <- x
  }
  s.list
}

modelFlags <- function(.internal.constraint=5e-4,
                       .internal.counter=0L,
                       label_switch=FALSE,
                       warn=FALSE){
  list(.internal.constraint=.internal.constraint,
       .internal.counter=.internal.counter,
       label_switch=label_switch,
       warn=warn)
}

extractFlags <- function(old){
  list(.internal.constraint=old@.internal.constraint,
       .internal.counter=old@.internal.counter,
       label_switch=label_switch(old),
       warn=FALSE)
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
  modal.ordinates <- modes(from)
  mb <- MultiBatch(data=data,
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
    m <- modal.ordinates[names(current_values(mb))]
    modes(mb) <- m
  }
  mb
})

setAs("MultiBatch", "MultiBatchModel", function(from){
  flag1 <- as.integer(flags(from)[[".internal.constraint"]])
  flag2 <- as.integer(flags(from)[[".internal.counter"]])
  be <- as.integer(table(batch(from)))
  names(be) <- unique(batch(from))
  dat <- downSampledData(from)
  KB <- prod(dim(theta(from)))
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
             predictive=numeric(KB),
             zstar=integer(KB),
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
  modal.ordinates <- modes(from)
  if(length(modal.ordinates) > 0 ){
    ix <- match(c("nu.0", "p"), names(modal.ordinates))
    names(modal.ordinates)[ix] <- c("nu0", "mixprob")
    modal.ordinates$z <- map_z(from)
    modal.ordinates$probz <- probz(from)
    modal.ordinates$u <- u(from)
    modes(obj) <- modal.ordinates
  }
  obj
})

setAs("MultiBatch", "list", function(from){
  ns <- nStarts(from)
  ##
  ## This initializes a list of models each with starting values simulated independently from the hyperparameters
  ##
  mb.list <- replicate(ns, MultiBatch(model=modelName(from),
                                      data=assays(from),
                                      down_sample=down_sample(from),
                                      mp=mcmcParams(from),
                                      hp=hyperParams(from),
                                      chains=chains(from)))
  mb.list <- lapply(mb.list, as, "MultiBatchModel")
  mb.list <- lapply(mb.list, function(x) {nStarts(x) <- 1; return(x)})
  mb.list
})

##setAs("list", "MultiBatch", function(from){
##  mb.list
##})




##
## Accessors
##
setMethod("current_values", "MultiBatch", function(object){
  object@current_values
})

setMethod("theta", "MultiBatch", function(object){
  th <- current_values(object)[["theta"]]
  th
})



setReplaceMethod("current_values", c("MultiBatch", "list"),
                 function(object, value){
  object@current_values <- value
  object
})

setReplaceMethod("theta", c("MultiBatch", "matrix"),
                 function(object, value){
                   current_values(object)[["theta"]] <- value
                   object
                 })

setMethod("sigma2", "MultiBatch", function(object){
  current_values(object)[["sigma2"]]
})

setReplaceMethod("sigma2", c("MultiBatch", "matrix"),
                 function(object, value){
                   current_values(object)[["sigma2"]] <- value
                   object
                 })

setMethod("sigma_", "MultiBatch", function(object){
  sqrt(sigma2(object))
})

setReplaceMethod("sigma_", c("MultiBatch", "matrix"),
                 function(object, value){
                   sigma2(object) <- value^2
                   object
                 })

setMethod("p", "MultiBatch", function(object){
  current_values(object)[["p"]]
})

setReplaceMethod("p", c("MultiBatch", "numeric"),
                 function(object, value){
                   current_values(object)[["p"]] <- value
                   object
                 })

setMethod("nu.0", "MultiBatch", function(object){
  current_values(object)[["nu.0"]]
})

setReplaceMethod("nu.0", c("MultiBatch", "numeric"),
                 function(object, value){
                   current_values(object)[["nu.0"]] <- value
                   object
                 })

setMethod("sigma2.0", "MultiBatch", function(object){
  current_values(object)[["sigma2.0"]]
})

setReplaceMethod("sigma2.0", c("MultiBatch", "numeric"),
                 function(object, value){
                   current_values(object)[["sigma2.0"]] <- value
                   object
                 })

setMethod("mu", "MultiBatch", function(object){
  current_values(object)[["mu"]]
})

setReplaceMethod("mu", c("MultiBatch", "numeric"),
                 function(object, value){
                   current_values(object)[["mu"]] <- value
                   object
                 })

setMethod("tau2", "MultiBatch", function(object){
  current_values(object)[["tau2"]]
})

setReplaceMethod("tau2", c("MultiBatch", "numeric"),
                 function(object, value){
                   current_values(object)[["tau2"]] <- value
                   object
                 })

setMethod("log_lik", "MultiBatch", function(object){
  current_values(object)[["loglik"]]
})

setReplaceMethod("log_lik", c("MultiBatch", "numeric"),
                 function(object, value){
                   current_values(object)[["loglik"]] <- value
                   object
                 })

setMethod("logPrior", "MultiBatch", function(object){
  current_values(object)[["logprior"]]
})

setReplaceMethod("logPrior", c("MultiBatch", "numeric"),
                 function(object, value){
                   current_values(object)[["logprior"]] <- value
                   object
                 })

setMethod("probz", "MultiBatch", function(object){
  current_values(object)[["probz"]]
})

zProb <- function(object){
  probz(object)/(iter(object)-1)
}

setReplaceMethod("probz", c("MultiBatch", "matrix"),
                 function(object, value){
                   current_values(object)[["probz"]] <- value
                   object
                 })

setMethod("z", "MultiBatch", function(object){
  current_values(object)[["z"]]
})

setReplaceMethod("z", c("MultiBatch", "numeric"),
                 function(object, value){
                   current_values(object)[["z"]] <- value
                   object
                 })


setMethod("u", "MultiBatch", function(object){
  current_values(object)[["u"]]
})

setMethod("nrow", "MultiBatch", function(x) nrow(assays(x)))

setReplaceMethod("u", c("MultiBatch", "numeric"),
                 function(object, value){
                   current_values(object)[["u"]] <- value
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
      ch <- mcmc_chains(specs(object), parameters(object))
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
                   iter(chains(object)) <- as.integer(value)
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

setMethod("down_sample", "MultiBatch", function(object) object@down_sample)

setReplaceMethod("down_sample", "MultiBatch", function(object, value){
  object@down_sample <- value
  object
})

setMethod("downSampledData", "MultiBatch", function(object){
  assays(object)[down_sample(object), ]
})

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
                nu.0=nu.0(mc)[i],
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

setMethod("setModes", "MultiBatch", function(object){
  modal.ordinates <- modes(object)
  current_values(object) <- modal.ordinates
  object
})

setMethod("useModes", "MultiBatch", function(object) setModes(object))


setMethod("modes", "MultiBatch", function(object) summaries(object)[["modes"]])

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
  n.internal.counter <- map(model.list, getCounter) %>%
    unlist %>%
    sum
  flags <- list(label_switch=nlabel_swap > 0,
                .internal.constraint=n.internal.constraint,
                .internal.counter=n.internal.counter)
}

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
  pred <- map(ch.list, predictive) %>% do.call(rbind, .)
  zz <- map(ch.list, zstar) %>% do.call(rbind, .)
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
            loglik=ll,
            predictive=pred,
            zstar=zz,
            iter=nrow(th),
            k=k(model.list[[1]]),
            B=numBatch(model.list[[1]]))
  mc
}

## update the current values with the posterior means across all chains

combineModels <- function(model.list){
  mc <- combineChains(model.list)
  pz.list <- lapply(model.list, probz)
  pz <- Reduce("+", pz.list) %>%
    "/"(rowSums(.))
  hp <- hyperParams(model.list[[1]])
  mp <- mcmcParams(model.list[[1]])
  ##nStarts(mp) <- length(model.list)
  nStarts(mp) <- 1L
  iter(mp) <- iter(mp) * length(model.list)
  nStarts(mp) <- 1
  tib <- tibble(oned=y(model.list[[1]])) %>%
    mutate(id=seq_along(oned),
           id=as.character(id),
           batch=batch(model.list[[1]]))
  param.list <- list(mp=mp, hp=hp)
  if(is(model.list[[1]], "MultiBatchPooled")){
    mb <- MultiBatchP(modelName(model.list[[1]]),
                      data=tib,
                      parameters=param.list,
                      chains=mc)
  } else {
    mb <- MultiBatch(modelName(model.list[[1]]),
                     data=tib,
                     parameters=param.list,
                     chains=mc)
  }
  probz(mb) <- pz
  summaries(mb) <- summarizeModel(mb)
  current_values(mb)[["probz"]] <- pz
  flags(mb) <- collectFlags(model.list)
  ##
  ## Set current values to the modal ordinates
  ##  (except for u which is not stored)
  ##  ( for z, we use the map estimate)
  ##current_values(mb) <- computeModes(mb)
  current_values(mb) <- summaries(mb)[["modes"]]
  mb
}


##
## Markov chain Monte Carlo
##
continueMcmc <- function(mp){
  burnin(mp) <- as.integer(burnin(mp) * 2)
  mp@thin <- as.integer(thin(mp) + 2)
  nStarts(mp) <- nStarts(mp) + 1
  mp
}

setFlags <- function(mb.list){
  mb1 <- mb.list[[1]]
  mp <- mcmcParams(mb1)
  hp <- hyperParams(mb1)
  mcmc_list <- mcmcList(mb.list)
  r <- gelman_rubin(mcmc_list, hp)
  mb <- combineModels(mb.list)
  tmp <- tryCatch(validObject(mb), error=function(e) NULL)
  if(is.null(tmp)) browser()
  flags(mb)[["fails_GR"]] <- r$mpsrf > min_GR(mp)
  neff <- tryCatch(effectiveSize(mcmc_list), error=function(e) NULL)
  if(is.null(neff)){
    neff <- 0
  }else {
    neff <- neff[ neff > 0 ]
  }
  flags(mb)[["small_effsize"]] <- mean(neff) < min_effsize(mp)
  mb
}



setMethod("convergence", "MultiBatch", function(object){
  !flags(object)$label_switch && !flags(object)[["fails_GR"]] && !flags(object)[["small_effsize"]]
})

setMethod("convergence", "MultiBatchList", function(object){
  sapply(object, convergence)
})

setMethod("max_burnin", "MultiBatch", function(object) {
  max_burnin(mcmcParams(object))
})

##
## burnin 100 iterations and choose the top nStart models by log lik
##
startingValues2 <- function(object){
  object2 <- object
  ns <- nStarts(object)
  burnin(object2) <- 100L
  iter(object2) <- 0L
  nStarts(object2) <- ns*2
  obj.list <- as(object2, "list")
  obj.list2 <- posteriorSimulation(obj.list)
  ll <- sapply(obj.list2, log_lik)
  obj.list2 <- obj.list2[ is.finite(ll) ]
  ll <- ll[ is.finite(ll) ]
  ix <- order(ll, decreasing=TRUE)
  if(length(ix) >= ns) {
    ix <- ix[ seq_len(ns) ] ## top ns models
    obj.list3 <- obj.list2[ix]
    obj.list3 <- lapply(obj.list3, function(x, i) {
      iter(x) <- i
      x
    }, i=iter(object))
    obj.list3 <- lapply(obj.list3, function(x, i) {
      burnin(x) <- i
      x
    }, i=burnin(object))
    return(obj.list3)
  }
  stop("problem identifying starting values")
}

setMethod("mcmc2", "MultiBatch", function(object, guide){
  mb <- object
  mp <- mcmcParams(mb)
  K <- specs(object)$k
  if(iter(mp) < 500)
    if(flags(object)$warn) warning("Very few Monte Carlo simulations specified")
  maxb <- max(max_burnin(mp), burnin(mp))
  while(burnin(mp) <= maxb && thin(mp) < 100){
    message("  k: ", K, ", burnin: ", burnin(mp), ", thin: ", thin(mp))
    mcmcParams(mb) <- mp
    ##
    ## Convert to list of MultiBatchModels with independent starting values
    ##
    if(missing(guide)){
      mb.list <- startingValues2(mb)
    } else {
      mb.list <- replicate(nStarts(mb), singleBatchGuided(mb, guide))
      if(class(mb.list[[1]]) == "MultiBatchP"){
        mb.list <- lapply(mb.list, as, "MultiBatchPooled")        
      } else {
        mb.list <- lapply(mb.list, as, "MultiBatchModel")
      }
    }
    ##
    ## Run posterior simulation on each
    ##
    mb.list <- posteriorSimulation(mb.list)
    mb <- setFlags(mb.list)
    ## if no flags, move on
    if( convergence(mb) ) break()
    mp <- continueMcmc(mp)
  }
  if( convergence(mb) ) {
    mb <- setModes(mb)
    mb <- compute_marginal_lik(mb)
  }
  stopifnot(validObject(mb))
  mb
})

setMethod("compute_marginal_lik", "MultiBatch", function(object, params){
  if(missing(params)){
    params <- mlParams(root=1/2,
                       reject.threshold=exp(-100),
                       prop.threshold=0.5,
                       prop.effective.size=0)
  }
  mbm <- as(object, "MultiBatchModel")
  ml <- tryCatch(marginalLikelihood(mbm, params), warning=function(w) NULL)
  if(!is.null(ml)){
    summaries(object)[["marginal_lik"]] <- ml
    message("     marginal likelihood: ", round(ml, 2))
  } else {
    ##warning("Unable to compute marginal likelihood")
    message("Unable to compute marginal likelihood")
  }
  object
})

##setMethod("marginalLikelihood", "MultiBatch", function(model, params){
##  mb <- as(model, "MultiBatchModel")
##  ml <- marginalLikelihood(mb, params)
##  ml
##})

setMethod("downSampleModel", "MultiBatch", function(object, N=1000, i){
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
  current.vals <- current_values(object)
  current.vals[["u"]] <- current.vals[["u"]][i]
  current.vals[["z"]] <- current.vals[["z"]][i]
  current.vals[["probz"]] <- current.vals[["probz"]][i, , drop=FALSE]
  mb <- MultiBatch(model=modelName(object),
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

setMethod("probability_z", "MultiBatch", function(object){
  thetas <- theta(object)
  sigmas <- sigma(object)
  p.comp <- p(object)
  df <- dfr(object)
  K <- seq_len(k(object))
  B <- specs(object)$number_batches %>%
                   seq
  N <- specs(object)$number_obs
  pz <- matrix(NA, N, max(K))
  dat <- assays(object)
  batches <- dat$batch
  for(b in B){
    j <- which(batches == b)
    m <- thetas[b, ]
    ss <- sigmas[b, ]
    yy <- dat$oned[j]
    for(k in K){
      temp <- p.comp[k] * dst(yy, df=df, mu=m[k], sigma=ss[k])
      pz[j, k] <- temp
    }
  }
  pz2 <- pz/rowSums(pz)
  ## the probz slot expects a frequency
  freq <- pz2 * (iter(object) - 1)
  freq2 <- matrix(as.integer(freq), nrow=nrow(freq), ncol=ncol(freq))
  freq2
})

setMethod("dfr", "MultiBatch", function(object){
  df <- dfr(hyperParams(object))
  df
})

setMethod("upsample_z", "MultiBatch", function(object){
  down_sample(object) <- seq_len(specs(object)$number_obs)
  object@specs$number_sampled <- specs(object)$number_obs
  current_values(object)[["u"]] <- rchisq(specs(object)$number_obs, dfr(object))
  mbm <- as(object, "MultiBatchModel")
  zz <- update_z(mbm)
  zz
})

setMethod("upSampleModel", "MultiBatch", function(object){
  ds <- down_sample(object)
  N <- specs(object)$number_obs
  us <- seq_len(N)
  if(identical(ds, us)) return(object)
  object <- setModes(object)
  current_values(object)[["u_up"]] <- rchisq(N, dfr(object))
  current_values(object)[["probz_up"]] <- probability_z(object)
  current_values(object)[["z_up"]] <- upsample_z(object)
  object
})

##singleBatchGuided <- function(model, sb){
##  modes(sb) <- computeModes(sb)
##  sb <- setModes(sb)
##  sp <- specs(model)
##  vals <- current_values(sb)
##}

setMethod("modelName", "MultiBatch", function(object) specs(object)$model)

setReplaceMethod("max_burnin", "MultiBatch", function(object, value){
  mp <- mcmcParams(object)
  max_burnin(mp) <- as.integer(value)
  mcmcParams(object) <- mp
  object
})

setMethod("predictive", "MultiBatch", function(object) predictive(chains(object)))
setMethod("zstar", "MultiBatch", function(object) zstar(chains(object)))

setMethod("singleBatchGuided", c("MultiBatch", "MultiBatch"), function(x, guide){
  stopifnot(k(x) == k(guide))
  means <- theta(guide)[1, ]
  sds <- sqrt(sigma2(guide))[1, ]
  ##
  ## number of means to simulate depends on the model
  ##
  mu(x) <- mu(guide)
  tau2(x) <- tau2(guide)
  B <- numBatch(x)
  K <- k(x)
  th <- theta(x)
  if(FALSE){
    ## Is prior to informative for these not to give reasonable values of theta?
    ##
    ## -- tau seems much too big -- is prior driving tau to larger values
    ## -- simulated values of theta too disperse
    for(j in seq_len(K)){
      th[, j] <- rnorm(B, mu(guide)[j], tau(guide)[j])
    }
  }
  for(j in seq_len(K)){
    th[, j] <- rnorm(B, theta(guide)[, j], sds[j]/2)
  }
  theta(x) <- th
  nu.0(x) <- nu.0(guide)
  sigma2.0(x) <- sigma2.0(guide)
  ## 1/sigma2 ~gamma
  ## sigma2 is invgamma
  NC <- ncol(sigma2(x))
  if(NC == 1){
    w <- as.numeric(table(z(guide)))
    sigma2(x) <- matrix((sum((w * sds)/sum(w)))^2, B, 1)
  } else{
    sigma2(x) <- matrix(sds/2, B, K, byrow=TRUE)
  }
  nStarts(x) <- 1L
  ##
  ## shouldn't have to initialize z since z is the first update of the gibbs sampler (and its update would be conditional on the above values)
  x
})

listModelsByDecreasingK <- function(object){
  N <- nrow(object)
  object2 <- augmentData2(object)
  sp <- specs(object2)
  object2.list <- split(object2, sp$k)
  object2.list <- rev(object2.list)
  object2.list
}
