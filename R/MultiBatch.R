#' @include AllGenerics.R
NULL

setValidity("MultiBatch", function(object){
  msg <- TRUE
  if(iter(object) != iter(chains(object))){
    msg <- "Number of iterations not the same between MultiBatch model and chains"
    return(msg)
  }
  if(length(u(object)) != nrow(assays(object))){
    msg <- "Incorrect length of u-vector"
    return(msg)
  }
  nb <- length(unique(batch(object)))
  if(nb != specs(object)$number_batches[1] ){
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
    msg <- "sample index must be in batch order"
    return(msg)
  }
  if(!identical(dim(zstar(object)), dim(predictive(object)))){
    msg <- "z* and predictive matrices in MCMC chains should be the same dimension"
    return(msg)
  }
  if(!is.matrix(p(object))){
    msg <- "mixture probabilities should be a matrix with dimensions B x K"
    return(msg)
  }
  if(nrow(p(object)) != numBatch(object)){
    msg <- "matrix of mixture probabilities should have dimensions B x K"
    return(msg)
  }
  if(ncol(p(chains(object))) != numBatch(object) * k(object)){
    msg <- "matrix of mixture probabilities in McmcChains should have dimensions B x K"
    return(msg)
  }
  msg
})


setValidity("CnList", function(object){
  msg <- TRUE
  if(iter(object) != iter(chains(object))){
    msg <- "Number of iterations not the same between MultiBatch model and chains"
    return(msg)
  }
  if(length(u(object)) != nrow(assays(object))){
    msg <- "Incorrect length of u-vector"
    return(msg)
  }
  nb <- length(unique(batch(object)))
  if(nb != specs(object)$number_batches[1] ){
    msg <- "Number of batches in model specs differs from number of batches in data"
    return(msg)
  }
  nr <- nrow(theta(object))
  if(nr != nb){
    msg <- "Number of batches in current values differs from number of batche in model specs"
    return(msg)
  }
  ##  nr <- nrow(dataMean(object))
  ##  if(nr != nb){
  ##    msg <- "Number of batches in model summaries differs from number of batches in model specs"
  ##    return(msg)
  ##  }
  S <- iter(object)
  th <- theta(chains(object))
  if( S != nrow(th) || ncol(th) != nr * k(object) ) {
    msg <- "Dimension of chain parameters is not consistent with model specs"
    return(msg)
  }
  B <- batch(object)
  is.sorted <- all(diff(B) >= 0)
  if(!is.sorted) {
    msg <- "sample index must be in batch order"
    return(msg)
  }
  if(!identical(dim(zstar(object)), dim(predictive(object)))){
    msg <- "z* and predictive matrices in MCMC chains should be the same dimension"
    return(msg)
  }
  msg
})

model_spec <- function(model, data) {
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
                number_obs=as.integer(number_obs))
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
#' @export
MultiBatch <- function(model="MB3",
                       data=modelData(),
                       ## by default, assume no downsampling
                       specs=model_spec(model, data),
                       iter=1000L,
                       burnin=200L,
                       thin=1L,
                       nStarts=4L,
                       max_burnin=burnin,
                       hp=Hyperparameters(k=specs$k),
                       mp=McmcParams(iter=iter, thin=thin,
                                     burnin=burnin,
                                     nStarts=nStarts,
                                     max_burnin=max_burnin),
                       parameters=modelParameters(mp=mp, hp=hp),
                       chains=mcmc_chains(specs, parameters),
                       current_values=modelValues2(specs, data, hp),
                       summaries=modelSummaries(specs),
                       flags=modelFlags()){
  is_SB <- substr(model, 1, 2) == "SB"
  if(nrow(data) > 0 && is_SB){
    data$batch <- 1L
  }
  if(nrow(data) > 0){
    if("batch" %in% colnames(data)){
      data <- data[order(data$batch), , drop=FALSE]
    }
  }
  if(!"batch" %in% colnames(data)) data$batch <- 1L
  model <- new("MultiBatch",
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
  cat("   n. batches          :", nBatch(object), "\n")
  cat("   k                   :", k(object), "\n")
  cat("   nobs/batch          :", table(batch(object)), "\n")
  cat("   saved mcmc          :", saved.mcmc, "\n")
  cat("   log lik (s)         :", round(log_lik(object), 1), "\n")
  ##cat("     log prior (s)       :", round(logPrior(object), 1), "\n")
  if(include_ml)
    cat("   log marginal lik    :", round(ml, 1), "\n")
})

setMethod("length", "CnList", function(x)  nrow(specs(x)))

setMethod("show", "CnList", function(object){
  ##callNextMethod()
  cls <- class(object)
  L <- length(object)
  cat(L, "candidate genotypes for model", modelName(object)[1], "\n")
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
     current_values(object) <- modelValues2( spec, assays(object),
                                            hyperParams(object) )
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
  ##object <- harmonizeDimensions(object)
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
                 p=matrix(nrow=0, ncol=K),
                 mu=numeric(K),
                 tau2=numeric(K),
                 u=numeric(),
                 z=numeric(),
                 logprior=numeric(),
                 loglik=numeric(),
                 probz=matrix(nrow=0, ncol=K))
    return(vals)
  }
  n.sampled <- specs$number_obs
  B <- specs$number_batches
  K <- specs$k
  alpha(hp) <- rep(1, K)
  mu <- sort(rnorm(K, mu.0(hp), 2))
  tau2 <- 1/rgamma(K, 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
  p <- rdirichlet(1, alpha(hp))[1, ]
  p <- matrix(p, B, K, byrow=TRUE)
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
  pmeans <- colMeans(p)
  z <- sample(seq_len(K), prob=pmeans, replace=TRUE, size=n.sampled)
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
  mapping <- seq_len(K)
  list(data.mean=data.mean,
       data.prec=data.prec,
       zfreq=zfreq,
       marginal_lik=marginal_lik,
       modes=modes,
       mapping=mapping)
}

setMethod("mapping", "MultiBatch", function(object){
  summaries(object)$mapping
})


setReplaceMethod("mapping", "MultiBatch", function(object, value){
  summaries(object)$mapping <- value
  object
})

setMethod("copyNumber", "MultiBatch", function(object){
  component_labels <- mapping(object)
  zz <- map_z(object)
  cn <- component_labels[zz]
  cn
})

extractSummaries <- function(old){
  s.list <- list(data.mean=dataMean(old),
                 data.prec=dataPrec(old),
                 zfreq=zFreq(old),
                 logprior=logPrior(old),
                 loglik=log_lik(old),
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
                       small_effsize=NA,
                       fails_GR=NA,
                       warn=FALSE){
  ##
  ## fails_GR is only set when multiple chains are run
  ##
  list(.internal.constraint=.internal.constraint,
       .internal.counter=.internal.counter,
       label_switch=label_switch,
       small_effsize=NA,
       fails_GR=NA,
       warn=warn)
}

extractFlags <- function(old){
  list(.internal.constraint=old@.internal.constraint,
       .internal.counter=old@.internal.counter,
       label_switch=label_switch(old),
       ## these slots are not available in old
       small_effsize=NA,
       fails_GR=NA,
       warn=FALSE)
}

modelData <- function(id=character(),
                      oned=numeric(),
                      batch=integer(),
                      is_simulated=logical()){
  tibble(id=id,
         oned=oned,
         batch=batch,
         is_simulated=is_simulated)
}


setMethod("isSimulated", "MultiBatch", function(object){
  if('is_simulated' %in% colnames(assays(object))){
    is_sim <- assays(object)$is_simulated
  } else {
    is_sim <- rep(FALSE, nrow(object))
  }
  is_sim
})

setMethod("isSimulated", "MultiBatchList", function(object){
  if('is_simulated' %in% colnames(assays(object))){
    is_sim <- assays(object)$is_simulated
  } else {
    is_sim <- rep(FALSE, nrow(object))
  }
  is_sim
})

setMethod("isSimulated", "MixtureModel", function(object){
  rep(FALSE, length(y(object)))
})


extractData <- function(old){
  tibble(id=as.character(seq_along(y(old))),
         oned=y(old),
         batch=batch(old),
         is_simulated=FALSE)
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
  specs <- model_spec(modelName(from), data)
  modal.ordinates <- modes(from)
  mb <- MultiBatch(data=data,
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
  dat <- assays(from)
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
  mp <- mcmcParams(from)
  nStarts(mp) <- 1L
  mb.list <- replicate(ns, MultiBatch(model=modelName(from),
                                      data=assays(from),
                                      mp=mp,
                                      hp=hyperParams(from),
                                      chains=chains(from),
                                      summaries=summaries(from),
                                      current_values=current_values(from)))
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

setMethod("current_values2", "MultiBatch", function(object){
  cv <- object@current_values
  cv <- cv[ !names(cv) %in% c("u", "z", "probz") ]
  cv
})

setReplaceMethod("current_values2", "MultiBatch", function(object, value){
  object@current_values[ names(value) ] <- value
  object
})

setMethod("summaries2", "MultiBatch", function(object){
  x <- object@summaries
  m <- x$modes
  m <- m[ !names(x) %in% c("u", "z", "probz")]
  x$modes <- m
  x
})

setReplaceMethod("summaries2", c("MultiBatch", "list"), function(object, value){
  x <- object@summaries
  value_modes <- value$modes
  value_modes <- value_modes[ !names(value_modes) %in% c("u", "z", "probz") ]
  x$modes[ names(value_modes) ] <- value_modes
  nms <- c("data.mean", "data.prec", "marginal_lik",
           "zfreq", "mapping")
  x[nms] <- value[ nms ]
  x <- x[ names(value) ]
  object@summaries <- x
  object
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

setReplaceMethod("p", c("MultiBatch", "matrix"),
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

setReplaceMethod("assays", "MultiBatch", function(x, value){
  x@data <- value
  ##x <- harmonizeDimensions(x)
  x
})

##
## Data accessors
##
setMethod("batch", "MultiBatch", function(object){
  assays(object)[["batch"]]
})

setReplaceMethod("oned", c("MultiBatch", "numeric"), function(object, value){
  assays(object)[["oned"]] <- value
  object
})

setMethod("oned", "MultiBatch", function(object){
  assays(object)[["oned"]]
})

setMethod("zFreq", "MultiBatch", function(object){
  summaries(object)[["zfreq"]]
})

setReplaceMethod("zFreq", "MultiBatch", function(object, value){
  summaries(object)[["zfreq"]] <- value
  object
})



setReplaceMethod("batch", c("MultiBatch", "numeric"), function(object, value){
  assays(object)[["batch"]] <- as.integer(value)
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
  if(iter(object) == 0){
    modes <- list(theta=theta(object),
                  sigma2=sigma2(object),
                  nu.0=nu.0(object),
                  sigma2.0=sigma2.0(object),
                  p=p(object),
                  mu=mu(object),
                  tau2=tau2(object),
                  u=u(object),
                  z=z(object),
                  logprior=logPrior(object),
                  loglik=log_lik(object),
                  probz=probz(object))
    return(modes)
  }
  i <- argMax(object)[1]
  mc <- chains(object)
  B <- specs(object)$number_batches
  K <- k(object)
  thetamax <- matrix(theta(mc)[i, ], B, K)
  sigma2max <- matrix(sigma2(mc)[i, ], B, K)
  pmax <- matrix(p(mc)[i, ], B, K)
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
  tib <- assays(object) %>%
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

setGeneric("computeMixProbs", function(object) standardGeneric("computeMixProbs"))
setMethod("computeMixProbs", "MultiBatch", function(object){
  z(object) <- map_z(object)
  tib2 <- assays(object) %>%
    group_by(batch) %>%
    summarize(N=n())
  ## frequency of copy number states in each batch
  tib <- assays(object) %>%
    mutate(z=z(object)) %>%
    group_by(batch, z) %>%
    summarize(n=n()) %>%
    left_join(tib2, by="batch")
  tib3 <- expand.grid(unique(batch(object)), seq_len(k(object))) %>%
    as_tibble() %>%
    set_colnames(c("batch", "z")) %>%
    left_join(tib, by=c("batch", "z"))
  p <- matrix(tib3$n/tib3$N, nrow(tib2), ncol=k(object),
              byrow=TRUE)
  dimnames(p) <- NULL
  p
})

setMethod("computePrec", "MultiBatch", function(object){
  z(object) <- map_z(object)
  tib <- assays(object) %>%
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
                marginal_lik=marginal_lik(object),
                mapping=seq_len(k(object)))
}

collectFlags <- function(model.list){
  getConstraint <- function(model) flags(model)$.internal.constraint
  getCounter <- function(model) flags(model)$.internal.counter
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

#' @export
combineModels <- function(model.list){
  mc <- combineChains(model.list)
  pz.list <- lapply(model.list, probz)
  pz <- Reduce("+", pz.list) %>%
    "/"(rowSums(.))
  hp <- hyperParams(model.list[[1]])
  mp <- mcmcParams(model.list[[1]])
  iter(mp) <- iter(mp) * length(model.list)
  nStarts(mp) <- 1
  mb <- model.list[[1]]
  param.list <- list(mp=mp, hp=hp)
  if(is(mb, "MultiBatchP")){
    mb <- MultiBatchP(modelName(model.list[[1]]),
                      data=assays(mb),
                      parameters=param.list,
                      chains=mc)
  } else {
    mb <- MultiBatch(modelName(model.list[[1]]),
                     data=assays(mb),
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
  r <- tryCatch(gelman_rubin(mcmc_list, hp), error=function(e) NULL)
  mb <- combineModels(mb.list)
  if(is.null(r)){
    flags(mb)[["fails_GR"]]
  }
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
  ml <- mcmcList(object)
  effsize <- ml %>%
    effectiveSize %>%
    median
  gelman_rub <- gelman_rubin(ml)$mpsrf
  flags(object)[["small_effsize"]] <- effsize < min_effsize(object)
  flags(object)[["fails_GR"]] <- gelman_rub > min_GR(object)
  smallcount <- flags(object)$.internal.counter > 10
  fl <- c(label_switch=flags(object)[["label_switch"]],
          small_effsize=flags(object)[["small_effsize"]],
          fails_GR=flags(object)[["fails_GR"]],
          high_internal_counter=smallcount)
  any_flags <- any(fl)
  !any_flags
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

setMethod("mcmc2", c("MultiBatch", "missing"), function(object, guide){
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
    mb.list <- as(mb, "list")
    ##
    ## Run posterior simulation on each
    ##
    mb.list <- posteriorSimulation(mb.list)
    mb <- setFlags(mb.list)
    assays(mb) <- assays(object)
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


setMethod("mcmc2", c("MultiBatch", "MultiBatch"), function(object, guide){
  mb <- object
  mp <- mcmcParams(mb)
  K <- specs(object)$k
  if(iter(mp) < 500)
    if(flags(object)$warn) warning("Very few Monte Carlo simulations specified")
  maxb <- max(max_burnin(mp), burnin(mp))
  while(burnin(mp) <= maxb && thin(mp) < 100){
    message("  k: ", K, ", burnin: ", burnin(mp), ", thin: ", thin(mp))
    mcmcParams(mb) <- mp
    mb.list <- replicate(nStarts(mb), singleBatchGuided(mb, guide))
    mb.list <- posteriorSimulation(mb.list)
    mb <- setFlags(mb.list)
    assays(mb) <- assays(object)
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

##setMethod("downSampleModel", "MultiBatch", function(object, N=1000, i){
##  if(!missing(N)){
##    if(N >= nrow(assays(object))){
##      return(object)
##    }
##  }
##  ## by sorting, batches are guaranteed to be ordered
##  nr <- nrow(assays(object))
##  if(missing(i)){
##    i <- sort(sample(seq_len(nr), N, replace=FALSE))
##  }
##  down_sample <- rep(FALSE, nr)
##  down_sample[i] <- TRUE
##  assays(object)$down_sample <- down_sample
##
##  b <- assays(object)$batch[i]
##  current.vals <- current_values(object)
##  current.vals[["u"]] <- current.vals[["u"]][i]
##  current.vals[["z"]] <- current.vals[["z"]][i]
##  current.vals[["probz"]] <- current.vals[["probz"]][i, , drop=FALSE]
##  mb <- MultiBatch(model=modelName(object),
##                   data=assays(object),
##                   parameters=parameters(object),
##                   current_values=current.vals,
##                   chains=mcmc_chains( specs(object), parameters(object) ))
##  dataMean(mb) <- computeMeans(mb)
##  dataPrec(mb) <- computePrec(mb)
##  zFreq(mb) <- as.integer(table(z(mb)))
##  mb
##})

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
  ##down_sample(object) <- rep(TRUE, specs(object)$number_obs)
  object@specs$number_sampled <- specs(object)$number_obs
  current_values(object)[["u"]] <- rchisq(specs(object)$number_obs, dfr(object))
  mbm <- as(object, "MultiBatchModel")
  zz <- update_z(mbm)
  zz
})

setMethod("upSampleModel", "MultiBatch",
          function(downsampled.model, data.full){
            dat.ds <- assays(downsampled.model)
            ##dat.full <- assays(full.model)
            full.model <- arrange(data.full, batch) %>%
              MultiBatchList(data=.) %>%
              "[["(modelName(downsampled.model))
            m <- modes(downsampled.model)
            m[["u"]] <- modes(full.model)[["u"]]
            m[["z"]] <- modes(full.model)[["z"]]
            m[["probz"]] <- modes(full.model)[["probz"]]
            modes(full.model) <- m

            cv <- current_values(downsampled.model)
            cv[["u"]] <- current_values(full.model)[["u"]]
            cv[["z"]] <- current_values(full.model)[["z"]]
            cv[["probz"]] <- current_values(full.model)[["probz"]]
            current_values(full.model) <- cv
            burnin(full.model) <- 0L
            ##full.model <- posteriorSimulation(full.model)
            full.model
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

setMethod("singleBatchGuided", c("MultiBatchP", "MultiBatch"), function(x, guide){
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
  w <- as.numeric(table(z(guide)))
  sigma2(x) <- matrix((sum((w * sds)/sum(w)))^2, B, 1)
  nStarts(x) <- 1L
  x <- as(x, "MultiBatchPooled")
  m <- modes(x)
  m[["mixprob"]] <- matrix(modes(guide)[["p"]][1, ], B, k(x), byrow=TRUE)
  m[["theta"]] <- theta(x)
  m[["sigma2"]] <- sigma2(x)
  modes(x) <- m
  ##
  ## shouldn't have to initialize z since z is the first update of the gibbs sampler (and its update would be conditional on the above values)
  x
})

setMethod("singleBatchGuided", c("MultiBatchP", "MultiBatchP"), function(x, guide){
  stopifnot(k(x) == k(guide))
  ##means <- theta(guide)[1, ]
  means <- colMeans(theta(guide))
  sds <- median(sigma(guide)[, 1])
  mu(x) <- mu(guide)
  tau2(x) <- tau2(guide)
  B <- numBatch(x)
  K <- k(x)
  th <- t(replicate(B, rnorm(k(x), means, sds)))
  theta(x) <- th
  nu.0(x) <- nu.0(guide)
  sigma2.0(x) <- sigma2.0(guide)
  sigma2(x)[, 1] <- 2*sigma2(guide) ## start at more diffuse value
  nStarts(x) <- 1L
  x <- as(x, "MultiBatchPooled")
  x
})

setMethod("singleBatchGuided", c("MultiBatch", "MultiBatchP"), function(x, guide){
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
  x <- as(x, "MultiBatchModel")
  ##
  ## shouldn't have to initialize z since z is the first update of the gibbs sampler (and its update would be conditional on the above values)
  x
})


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
  x <- as(x, "MultiBatchModel")
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


.find_surrogates <- function(x, B, THR=0.1){
  if(length(unique(B))==1) return(list(pval=1, batches=B))
  B2 <- B
  dat <- tibble(x=x, batch=B) %>%
    group_by(batch) %>%
    summarize(n=n()) %>%
    arrange(-n)
  uB <- unique(dat$batch)
  ##uB <- unique(B)
  ## One plate can pair with many other plates.
  stat <- matrix(NA, length(uB), length(uB))
  comparison <- stat
  pval <- stat
  for(j in seq_along(uB)){
    for(k in seq_along(uB)){
      if(k <= j) next() ## next k
      ## get rid of duplicate values
      x <- x + runif(length(x), -1e-10, 1e-10)
      b1 <- uB[j]
      b2 <- uB[k]
      ## edits
      tmp <- ks.test(x[B==b1], x[B==b2])
      stat[j, k] <- tmp$statistic
      bb <- c(b1, b2) %>%
        sort %>%
        paste(collapse="::")
      comparison[j, k] <- bb
      pval[j, k] <- tmp$p.value
    }
  }
  ## Combine the most similar plates according to the KS test statistic (smallest test statistic)
  stat2 <- as.numeric(stat)
  comp2 <- as.character(comparison)
  pval2 <- as.numeric(pval)
  ##min.index <- which.min(stat2)
  max.index <- which.max(pval2)
  sim.batches <- strsplit(comp2[max.index], "::")[[1]]
  max.pval <- pval2[max.index]
  if(max.pval < THR) return(list(pval=max.pval, batches=B))
  comp3 <- gsub("::", ",", comp2[max.index])
  B2[ B %in% sim.batches ] <- comp3
  result <- list(pval=max.pval, batches=B2)
  result
}

find_surrogates <- function(dat, THR=0.1, min_oned=-1){
  ## do not define batches based on homozygous deletion
  dat2 <- filter(dat, oned > min_oned)
  current <- dat2$provisional_batch
  oned <- dat2$oned
  latest <- NULL
  while(!identical(current, latest)){
    if(is.null(latest)) latest <- current
    current <- latest
    latest <- .find_surrogates(oned, current, THR)$batches
  }
  result <- tibble(provisional_batch=dat2$provisional_batch,
                   batch=latest) %>%
    group_by(provisional_batch) %>%
    summarize(batch=unique(batch))
  if("batch" %in% colnames(dat)){
    dat <- select(dat, -batch)
  }
  dat3 <- dat %>%
    left_join(result, by="provisional_batch")
  dat3
}

#' Estimate batch from any sample-level surrogate variables that capture aspects of sample processing, such as the PCR experiment (e.g., the 96 well chemistry plate), laboratory, DNA source, or DNA extraction method.
#'
#' In high-throughput assays, low-level summaries of copy number at
#' copy number polymorphic loci (e.g., the mean log R ratio for each
#' sample, or a principal-component derived summary) often differ
#' between groups of samples due to technical sources of variation
#' such as reagents, technician, or laboratory.  Technical (as opposed
#' to biological) differences between groups of samples are referred
#' to as batch effects.  A useful surrogate for batch is the chemistry
#' plate on which the samples were hybridized. In large studies, a
#' Bayesian hierarchical mixture model with plate-specific means and
#' variances is computationally prohibitive.  However, chemistry
#' plates processed at similar times may be qualitatively similar in
#' terms of the distribution of the copy number summary statistic.
#' Further, we have observed that some copy number polymorphic loci
#' exhibit very little evidence of a batch effect, while other loci
#' are more prone to technical variation.  We suggest combining plates
#' that are qualitatively similar in terms of the Kolmogorov-Smirnov
#' two-sample test of the distribution and to implement this test
#' independently for each candidate copy number polymophism identified
#' in a study.  The \code{collapseBatch} function is a wrapper to the
#' \code{ks.test} implemented in the \code{stats} package that
#' compares all pairwise combinations of plates.  The \code{ks.test}
#' is performed recursively on the batch variables defined for a given
#' CNP until no batches can be combined. For smaller values of THR, plates are more likely to be judged as similar and combined.
#' @export
#' @param object a MultiBatch instance
#' @param THR scalar for the p-value cutoff from the K-S test.  Two batches with p-value > THR will be combined to define a single new batch
#'
#' @details All pairwise comparisons of batches are performed.  The two most similar batches are combined if the p-value exceeds THR.  The process is repeated recursively until no two batches can be combined.
#' @return MultiBatch object
setMethod("findSurrogates", "MultiBatch",
          function(object, THR=0.1, min_oned=-1){
  dat <- assays(object) %>%
    select(c(id, provisional_batch, oned))
  message("Putting rows of data in batch order")
  dat2 <- find_surrogates(dat, THR, min_oned) %>%
    mutate(batch=factor(batch, levels=unique(batch)),
           batch_labels=as.character(batch),
           batch=as.integer(batch)) %>%
    filter(!duplicated(id)) %>%
    arrange(batch) %>%
    select(c(provisional_batch, batch, batch_labels))  %>%
    filter(!duplicated(provisional_batch))
  if(any(is.na(dat2$batch))){
    ## randomly assign to one of available batches
    nna <- sum(is.na(dat2$batch))
    ub <- unique(dat2$batch)
    ub <- ub[!is.na(ub)]
    dat2$batch[is.na(dat2$batch)] <- sample(ub, nna, replace=TRUE)
  }
  ##
  ## There is a many to one mapping from provisional_batch to batch
  ## Since each sample belongs to a single plate, samples in the downsampled data will only belong to a single batch
  full.data <- assays(object)
  full.data2 <- full.data %>%
    select(-batch) %>%
    left_join(dat2, by="provisional_batch") %>%
    ##filter(!duplicated(id)) %>%
    arrange(batch)
  assays(object) <- full.data2
  L <- length(unique(full.data2$batch))
  if( L == specs(object)$number_batches ) return(object)
  spec <- specs(object)
  spec$number_batches <- L
  specs(object) <- spec
  current_values(object) <- modelValues2(specs(object),
                                         assays(object),
                                         parameters(object)[["hp"]])
  s <- modelSummaries(specs(object))
  s$data.mean <- computeMeans(object)
  s$data.prec <- computePrec(object)
  summaries(object) <- s
  chains(object) <- initialize_mcmc(k(object),
                                    iter(object),
                                    numBatch(object))
  object
})

.candidate_mapping <- function(model){
  if(k(model) == 5){
    candidate_models <- list(c(0, 1, 2, 3, 4),
                             c(0, 1, 2, 3, 3),
                             c(0, 1, 2, 2, 2),
                             c(1, 2, 2, 2, 3),
                             c(2, 2, 2, 2, 2),
                             c(2, 2, 3, 3, 4)) %>%
      lapply(as.character)
  }
  if(k(model) == 4){
    candidate_models <- list(c(0, 1, 2, 2),
                             c(0, 2, 2, 2),
                             c(2, 2, 2, 2),
                             c(0, 1, 1, 2),
                             c(0, 1, 2, 3),
                             c(1, 2, 3, 3)) %>%
      lapply(as.character)
  }
  if(k(model) == 3){
    candidate_models <- list(c(0, 1, 2),
                             ## hemizygous component can not be
                             ## distinguished
                             c(0, 2, 2),
                             ##c(1, 1, 2),
                             c(1, 2, 2),
                             c(2, 2, 2),
                             c(1, 2, 3),
                             c(2, 3, 4),
                             c(2, 3, 3)) %>%
      lapply(as.character)
  }
  if(k(model) == 2){
    candidate_models <- list(c(0, 2),
                             c(1, 2),
                             c(2, 3),
                             c(2, 2)) %>%
      lapply(as.character)
  }
  if(k(model) == 1){
    candidate_models <- list("2")
  }
  tibble(model=modelName(model),
         cn.model=sapply(candidate_models,
         paste, collapse=","))
}

setMethod("[[", c("CnList", "numeric"), function(x, i){
  spec <- specs(x)[i, ]
  model <- spec$model
  nm <- substr(model, 1, 3)
  if(nm == "SBP" || nm == "MBP"){
    mb <- MultiBatchP(model=model,
                      data=assays(x),
                      specs=spec,
                      parameters=parameters(x),
                      current_values=current_values(x),
                      chains=chains(x),
                      summaries=summaries(x),
                      flags=flags(x))
  } else {
    mb <- MultiBatch(model=model,
                     data=assays(x),
                     specs=spec,
                     parameters=parameters(x),
                     current_values=current_values(x),
                     chains=chains(x),
                     summaries=summaries(x),
                     flags=flags(x))
  }
  cn.map <- strsplit(specs(mb)$cn.model, ",")[[1]]
  mapping(mb) <- cn.map
  mb
})

isAugmented <- function(model){
  ix <- grep("augment_", id(model))
  if(length(ix) == 0){
    return(rep(FALSE, nrow(model)))
  }
  seq_len(nrow(model)) %in% ix
}

#' @export
CnList <- function(mb){
  mb <- mb[ !isAugmented(mb) ]
  cn.specs <- .candidate_mapping(mb) %>%
    mutate(k=k(mb),
           number_batches=numBatch(mb),
           number_obs=nrow(assays(mb)))
  clist <- new("CnList",
               data=assays(mb),
               specs=cn.specs,
               parameters=parameters(mb),
               chains=chains(mb),
               current_values=current_values(mb),
               summaries=summaries(mb),
               flags=flags(mb))
  clist
}

setMethod("probCopyNumber", "MultiBatch", function(model){
  .prob_copynumber(model)
})

setMethod("baf_loglik", "CnList", function(object, snpdat){
  clist <- object
  blik <- sapply(clist, modelProb, snpdata=snpdat)
  sp <- select(specs(clist), c("model", "cn.model", "k")) %>%
    mutate(baf_loglik=blik)
  ix <- order(sp$baf_loglik, decreasing=TRUE)
  sp[ix, ]
})

setMethod("lapply", "CnList", function(X, FUN, ...){
  result <- vector("list", length(X))
  for(i in seq_along(X)){
    result[[i]] <- FUN(X[[i]], ...)
  }
  result
})

setMethod("sapply", "CnList",
          function(X, FUN, ..., simplify=TRUE, USE.NAMES=TRUE){
  result <- lapply(X, FUN, ...)
  if(simplify){
    result <- unlist(result)
  }
  result
})

setMethod("numberStates", "MultiBatch", function(model){
  length(unique(mapping(model)))
})

setMethod("numberObs", "MultiBatch", function(model) {
  specs(model)$number_obs[1]
})

#' @export
setGeneric("id", function(object) standardGeneric("id"))
setMethod("id", "MultiBatch", function(object) assays(object)$id)
setMethod("id", "MultiBatchList", function(object) assays(object)$id)

#' @export
id2 <- function(object){
  id(object)[!isSimulated(object)]
}

setMethod("[", c("MultiBatch", "numeric"), function(x, i, j, ..., drop=FALSE){
  nbatch1 <- numBatch(x)
  x@data <- x@data[i, , drop=FALSE]
  ubatch <- unique(x@data$batch)
  cv <- current_values(x)
  cv$probz <- cv$probz[i, , drop=FALSE]
  cv$u <- cv$u[i]
  cv$z <- cv$z[i]
  cv$theta <- cv$theta[ubatch, , drop=FALSE]
  cv$sigma2 <- cv$sigma2[ubatch, , drop=FALSE]
  cv$p <- cv$p[ubatch, , drop=FALSE]
  x@current_values <- cv
  specs(x)$number_batches <- length(ubatch)
  specs(x)$number_obs <- nrow(x@data)
  nbatch2 <- numBatch(x)
  L2 <- length(unique(batch(x)))
  ##sp <- specs(x)
  ##sp$number_batches <- L2
  ##sp$number_obs <- length(i)
  ##specs(x) <- sp
  ##if( L == L2 ) return(x)
  current_values(x)[["theta"]] <- computeMeans(x)
  current_values(x)[["sigma2"]] <- 1/computePrec(x)
  current_values(x)[["p"]] <- computeMixProbs(x)
  summaries(x)[["data.mean"]] <- theta(x)
  summaries(x)[["data.prec"]] <- 1/sigma2(x)
  if(L2 == nbatch1) {
    ## leave chains alone and return
    return(x)
  }
  if(substr(modelName(x), 1, 3) == "MBP"){
    chains(x) <- initialize_mcmcP(k(x), iter(x), numBatch(x))
  } else {
    chains(x) <- initialize_mcmc(k(x), iter(x), numBatch(x))
  }
  assays(x)$batch <- as.integer(factor(batch(x)))
  x
})

setMethod("[", c("MultiBatch", "logical"), function(x, i, j, ..., drop=FALSE){
  i <- which(i)
  x <- x[i]
  x
})


setMethod("toSingleBatch", "MultiBatchP", function(object){
  mbp <- object
  B <- numBatch(mbp)
  if(B == 1) return(mbp)
  sp <- specs(mbp)
  sp$number_batches <- 1L
  model <- gsub("MBP", "SBP", modelName(mbp))
  dat <- assays(mbp)
  dat$batch <- 1L
  mbp2 <- MultiBatchP(model,
                      specs=sp,
                      data=dat,
                      parameters=parameters(mbp))
  mbp2
})

setMethod("toSingleBatch", "MultiBatch", function(object){
  mb <- object
  B <- numBatch(mb)
  if(B == 1) return(mb)
  sp <- specs(mb)
  sp$number_batches <- 1L
  model <- gsub("MB", "SB", modelName(mb))
  dat <- assays(mb)
  dat$batch <- 1L
  mb2 <- MultiBatch(model,
                    specs=sp,
                    data=dat,
                    parameters=parameters(mb))
  mb2
})

gr <- function(x){
  stat <- mcmcList(x) %>%
    gelman_rubin
  stat$mpsrf
}

#' @export
is_high_mpsrf <- function(m){
  mpsrf <- tryCatch(gr(m), error=function(e) NULL)
  if(is.null(mpsrf) || mpsrf > 1.25){
    return(TRUE)
  }
  FALSE
}

#' @export
setGeneric("fails_gr<-", function(object, value) standardGeneric("fails_gr<-"))

setReplaceMethod("fails_gr", "MultiBatch", function(object, value){
  flags(object)[["fails_GR"]] <- value
  object
})

#' @export
reset <- function(from, to){
  sp <- specs(to)
  nsim <- sum(isSimulated(to))
  sp$number_obs <- nrow(from) + nsim
  specs(from) <- sp
  is_pooled <- substr(modelName(to), 1, 3) %in% c("SBP", "MBP") 
  if(is_pooled){
    from <- as(from, "MultiBatchP")
  }
  if(nsim > 0){
    dat <- assays(from) %>%
      mutate(is_simulated=FALSE)
    simdat <- assays(to)[isSimulated(to), ] %>%
      select(colnames(dat))
    dat <- bind_rows(dat, simdat) %>%
      arrange(batch)
    assays(from) <- dat
  } else {
    assays(from)$is_simulated <- FALSE
  }
  B <- numBatch(to)
  if(B == 1){
    batch(from) <- 1L
  }
  current_values2(from) <- current_values2(to)
  K <- k(to)
  df <- dfr(to)
  current_values(from)$u <- rchisq(nrow(from), df)
  current_values(from)$probz <- matrix(0, nrow(from), K)
  summaries2(from) <- summaries2(to)
  modes(from)$u <- current_values(from)$u
  modes(from)$probz <- current_values(from)$probz
  parameters(from) <- parameters(to)
  if(!is_pooled){
    chains(from) <- mcmc_chains(sp, parameters(to))
  } else {
    chains(from) <- mcmc_chainsP(sp, parameters(to))
  }
  from
}

#' Downsampling will not work well if the range of the downsampled data is not similar to the range of the original data
#'
#' @export
sample2 <- function(mb, N){
  ix <- sort(sample(seq_len(nrow(mb)), N, replace=TRUE))
  r <- oned(mb)[ix]
  minr <- min(r)
  if(minr > -2 && min(oned(mb)) < -2){
    ix2 <- which(oned(mb) < -2)
    ix <- sort(c(ix, ix2))
  }
  ix
}

setMethod("min_effsize", "MultiBatch", function(object){
  min_effsize(mcmcParams(object))
})

setMethod("min_GR", "MultiBatch", function(object){
  min_GR(mcmcParams(object))
})

#' @export
incrementK <- function(object){
  model_name <- modelName(object)
  K <- substr(model_name, nchar(model_name), nchar(model_name))
  Kplus <- (as.integer(K) + 1)  %>%
    as.character
  model_name <- gsub(K, Kplus, model_name)
  model_name
}

#' @export
genotypeModel <- function(model, snpdat){
  keep <- !duplicated(id(model)) & !isSimulated(model)
  gmodel <- model[ keep ]
  snpdat2 <- snpdat[, id(gmodel) ]
  clist <- CnList(gmodel)
  (stats <- baf_loglik(clist, snpdat2))
  mapping(gmodel) <- strsplit(stats$cn.model[1], ",")[[1]]
  gmodel
}

#' @export
genotypeData <- function(gmodel, snpdat, min_probz=0.9){
  snpdat <- snpdat[, id(gmodel)]
  maxpz <- probz(gmodel) %>%
    "/"(rowSums(.)) %>%
    rowMax
  bafdat <- assays(snpdat)[["baf"]] %>%
    as_tibble() %>%
    mutate(rsid=rownames(snpdat)) %>%
    gather("id", "BAF", -rsid)
  cndat <- tibble(id=id(gmodel),
                  batch=batch(gmodel),
                  oned=oned(gmodel),
                  pz=maxpz,
                  z=map_z(gmodel)) %>%
    mutate(cn=mapping(gmodel)[z]) %>%
    mutate(cn=factor(cn))
  bafdat <- left_join(bafdat, cndat, by="id")
  bafdat2 <- filter(bafdat, pz > min_probz)
  bafdat2
}

homozygousdel_mean <- function(object, THR=-1) {
  mn <- mean(oned(object)[ oned(object) < THR])
  if(is.na(mn)) mn <- THR - 1
  mn
}

isPooledVar <- function(object) ncol(sigma2(object))==1

meanSdHomDel <- function(object, THR){
  list(homozygousdel_mean(object, THR),
       sd(oned(object)[ oned(object) < THR]))
}

sdModel2_3 <- function(mod_2.3){
  if(isPooledVar(mod_2.3)){
    sigma2_ <- sigma2(mod_2.3)
  } else {
    sigma2_ <- cbind(sigma2(mod_2.3)[1, 2],
                     sigma2(mod_2.3))
  }
  sigma2_
}

.augment_homozygous <- function(mb, mean_sd, THR=-1, phat=0.01){
  if(all(is.na(mean_sd[[2]]))) mean_sd[[2]] <- 0.1
  freq.hd <- assays(mb) %>%
    filter(!is_simulated) %>%
    group_by(batch) %>%
    summarize(N=n(),
              n=sum(likely_hd)) %>%
    filter(n/N < 0.05)
  if(nrow(freq.hd) == 0){
    obsdat <- assays(mb)
    return(obsdat)
  }
  loc.scale <- tibble(theta=mean_sd[[1]],
                      sigma2=mean_sd[[2]]^2,
                      phat=phat,
                      batch=seq_len(nrow(theta(mb))))
  loc.scale <- left_join(freq.hd, loc.scale, by="batch") %>%
    select(-N)
  start.index <- length(grep("augment", id(mb))) + 1
  any_simulated <- any(isSimulated(mb))
  imp.hd <- impute(mb, loc.scale, start.index=start.index)
  if(any_simulated){
    simulated <- bind_rows(imp.hd,
                           filter(assays(mb), is_simulated))
  } else simulated <- imp.hd
  obsdat <- assays(mb) %>%
    filter(is_simulated==FALSE)
  simdat <- bind_rows(obsdat, simulated) %>%
    arrange(batch) %>%
    mutate(homozygousdel_mean=mean_sd[[1]])
  simdat
}

augment_homozygous <- function(mb.subsamp, THR){
  mean_sd <- meanSdHomDel(mb.subsamp, -1)
  rare_homozygous <- sum(oned(mb.subsamp) < THR) < 5
  ##expect_false(rare_homozygous)
  if(rare_homozygous){
    simdat <- augment_homozygous(mb.subsamp, mean_sd, THR)
  } else {
    simdat <- assays(mb.subsamp) %>%
      arrange(batch) %>%
      mutate(homozygousdel_mean=mean_sd[[1]])
  }
  simdat
}

ok_hemizygous <- function(sb3){
  varratio <- max(sigma2(sb3))/min(sigma2(sb3))
  p(sb3)[2] < 0.1 || varratio > 100
}

.warmup <- function(tib, model){
  mbl <- replicate(10, MultiBatchList(data=tib)[[model]])
  for(j in seq_along(mbl)){
    cat(".")
    mb <- mbl[[j]]
    iter(mb) <- 0
    burnin(mb) <- 100
    mb <- tryCatch(posteriorSimulation(mb),
                   warning=function(w) NULL)
    if(is.null(mb)) next()
    mbl[[j]] <- mb
  }
  mbl
}

warmup <- function(tib, model1, model2=NULL){
  ##
  ##
  mbl <- .warmup(tib, model1)
  ml <- sapply(mbl, log_lik)
  if(is(ml, "list")){
    mbl <- mbl[ lengths(ml) > 0 ]
    ml <- ml[ lengths(ml) > 0 ]
    ml <- unlist(ml)
  }
  if(is.null(model2)){
    model <- mbl[[which.max(ml)]]
    return(model)
  }
  if(!is.null(model2)){
    mbl2 <- .warmup(tib, model2)
    ml2 <- sapply(mbl2, log_lik)
    if(is(ml2, "list")){
      mbl2 <- mbl2[ lengths(ml2) > 0 ]
      ml2 <- ml2[ lengths(ml2) > 0 ]
      ml2 <- unlist(ml2)
    }
    if(max(ml2, na.rm=TRUE) > max(ml, na.rm=TRUE) + 50){
      model <- mbl2[[which.max(ml2)]]
    } else {
      model <- mbl[[which.max(ml)]]
    }
  }
  return(model)
}

stop_early <- function(model){
  pz <- probz(model)
  maxprob <- rowMax(pz)
  pz <- probz(model) %>%
    "/"(rowSums(.)) %>%
    rowMax
  mean_maxp <- mean(pz > 0.99)
  mean_maxp > 0.995
}

revertToMultiBatch <- function(sb){
  adat <- assays(sb)
  batch_levels <- unique(assays(sb)$batch_labels)
  adat$batch <- as.integer(factor(assays(sb)$batch_labels,
                                  levels=batch_levels))
  model_name <- gsub("SB", "MB", modelName(sb))
  mb <- MultiBatch(model_name, data=adat)
  mb
}


augment_hemizygous <- function(mod_2.3, sb3){
  is_pooledvar <- TRUE
  ##
  ## This variable is needed later on
  ##
  batch_labels <- assays(mod_2.3)$batch_labels
  batch_labels <- batch_labels %>%
    factor(., levels=unique(.)) %>%
    as.integer(.) %>%
    unique(.) %>%
    sort()
  ## normalize probabilities
  pz <- probz(mod_2.3) %>%
    "/"(rowSums(.)) %>%
    rowMax
  ##
  ## batches with high posterior probabilities
  ##
  ubatch <- unique(batch(mod_2.3)[ pz >= 0.95 ])
  ##
  ## Find number of samples assigned to the hemizygous deletion
  ## component with high probability.
  ##
  ## Check if any of the batches have fewer than 10 subjects
  ## with high posterior probability
  ##
  tmp <- tibble(batch=batch(mod_2.3),
                z=map_z(mod_2.3),
                pz=pz) %>%
    group_by(batch) %>%
    summarize(N=n(),
              n=sum(z==1 & pz > 0.9)) %>%
    filter(n < 10)
  ##
  ##
  ##
  is_dropped <- !batch_labels %in% ubatch |
    batch_labels %in% tmp$batch
  condition3 <- any(is_dropped)
  if(condition3){
    ##
    ## possible hemizygous deletion missing in
    ##   one component (e.g., CNP023)
    ##
    ## There are batches with fewer than 10 samples
    ## assigned to hemiz. del component with probability > 0.9
    ##
    dropped_batches <- uniqueBatch(mod_2.3)[ is_dropped ]
    if(modelName(sb3) == "SBP3"){
      hemvar <- sigma2(sb3)[, 1]
    } else {
      hemvar <- sigma2(sb3)[, 2]
    }
    message("Find mean and variance of hemizygous deletion component")
    loc.scale.hem <- tibble(theta=theta(sb3)[2],
                            sigma2=hemvar,
                            phat=p(sb3)[2],
                            batch=dropped_batches,
                            theta.diploid=theta(mod_2.3)[dropped_batches, 2]) %>%
      mutate(delta=theta.diploid-theta)
    message("Augment data with additional hemizygous deletions")
    impdat <- impute(mod_2.3, loc.scale.hem, 1)
    obsdat <- assays(mod_2.3) %>%
      mutate(is_simulated=FALSE)
    simdat <- bind_rows(obsdat,
                        impdat) %>%
      arrange(batch)
  }
  mb <- MultiBatchList(data=simdat)[[modelName(mod_2.3)]]
  mcmcParams(mb) <- mcmcParams(mod_2.3)
  theta(mb) <- theta(mod_2.3)
  theta(mb)[is_dropped, 1] <- loc.scale.hem$theta
  sigma2(mb) <- matrix(pmin(sigma2(mod_2.3)[, 1], hemvar), ncol=1)
  mb
}

batchLevels <- function(mb){
  levs <- unique(assays(mb)$batch_labels)
  levels <- unique(as.integer(factor(assays(mb)$batch_labels, levels=levs)))
  sort(levels [ !is.na(levels) ])
}

.mcmcWithHomDel <- function(simdat, mod_2.3){
  is_pooledvar <- ncol(sigma2(mod_2.3))==1
  hdmean <- simdat$homozygousdel_mean[1]
  batch_labels <- batchLevels(mod_2.3)
  ##
  ## Rerun 3 batch model with augmented data
  ##
  mbl <- MultiBatchList(data=simdat)
  model <- incrementK(mod_2.3)
  mod_1.3 <- mbl[[ model ]]
  theta(mod_1.3) <- cbind(hdmean, theta(mod_2.3))
  if(is_pooledvar){
    sigma2(mod_1.3) <- sigma2(mod_2.3)
  } else {
    sigma2(mod_1.3) <- cbind(sigma2(sb3)[1], sigma2(mod_2.3))
  }
  burnin(mod_1.3) <- 200
  iter(mod_1.3) <- 0
  mod_1.3 <- tryCatch(posteriorSimulation(mod_1.3),
                      error=function(e) NULL)
  bad_start <- is.null(mod_1.3) || is.nan(log_lik(mod_1.3))
  if(bad_start){
    fdat <- filter(assays(mod_1.3), oned > -1)  %>%
      filter(is_simulated=FALSE)
    mod_1.3 <- warmup(fdat, model)
  }
  if(is.null(mod_1.3)) return(NULL)
  internal.count <- flags(mod_1.3)$.internal.counter
  any_dropped <- TRUE
  if(internal.count < 100){
    iter(mod_1.3) <- 1000
    mod_1.3 <- posteriorSimulation(mod_1.3)
  }
  mod_1.3
}

mcmcWithHomDel <- function(mb, sb, restricted_model){
  hdmean <- homozygousdel_mean(sb)
  mn_sd <- list(hdmean, sdModel2_3(restricted_model))
  ## If at least one batch has fewer than 5% subjects with
  ## homozygous deletion, augment the data for homozygous deletions
  mb.observed <- mb[ !isSimulated(mb) ]
  mb1 <- filter(assays(restricted_model),
                isSimulated(restricted_model)) %>%
    bind_rows(assays(mb.observed)) %>%
    arrange(batch) %>%
    MultiBatchList(data=.) %>%
    "[["("MB3")
  simdat <- .augment_homozygous(mb1, mn_sd, -1,
                               phat=max(p(sb)[1], 0.05))
  mod_1.3 <- .mcmcWithHomDel(simdat, restricted_model)
  mod_1.3
}


mcmcHomDelOnly <- function(simdat, restricted_model, model){
  mbl <- MultiBatchList(data=simdat)
  mod_1.3 <- mbl[[ model ]]
  hdmean <- homozygousdel_mean(mod_1.3)
  theta(mod_1.3) <- cbind(hdmean, theta(restricted_model))
  is_pooledvar <- ncol(sigma2(mod_1.3)) == 1
  ##expect_false(is_pooledvar)
  if(!is_pooledvar){
    multibatchvar <- sigma2(restricted_model)
    s2 <- replicate(k(mod_1.3)-1, multibatchvar, simplify=FALSE) %>%
      do.call(cbind, .)
    singlebatchvar <- sigma2(sb3)[, 1]
    foldchange_singlebatchvar <-
      singlebatchvar/median(multibatchvar[, 1])
    if(foldchange_singlebatchvar > 5) {
      singlebatchvar <- median(multibatchvar[, 1])*5
    }
    s2 <- cbind(singlebatchvar, s2)
  } else {
    s2 <- sigma2(restricted_model)
  }
  sigma2(mod_1.3) <- s2
  burnin(mod_1.3) <- 200
  iter(mod_1.3) <- 1000
  mod_1.3 <- mcmc_homozygous(mod_1.3)
  mod_1.3
}

ok_model <- function(mod_1.3, restricted_model){
  if(is.null(mod_1.3)){
    return(FALSE)
  }
  batch_labels <- batchLevels(mod_1.3)
  internal.count <- flags(mod_1.3)$.internal.counter
  pz <- probz(mod_1.3) %>%
    "/"(rowSums(.)) %>%
    rowMax
  ubatch <- batch(mod_1.3)[ pz >= 0.95 & z(mod_1.3) == 2] %>%
    unique
  any_dropped <- any(!batch_labels %in% ubatch)
  ## Check that variance estimates are comparable to restricted_model
  varratio <- sigma2(mod_1.3)/sigma2(restricted_model)
  bad_pooled_variance <- any(varratio > 4) ||
    internal.count >= 100 ||
    any_dropped
  !bad_pooled_variance
}

summarize_region <- function(se, provisional_batch, THR=-1){
  ##
  ## Flag homozygous deletions
  ##
  message("Flagging apparent homozygous deletions")
  dat <- tibble(id=colnames(se),
                oned=assays(se)[["MEDIAN"]][1, ],
                provisional_batch=provisional_batch) %>%
    mutate(likely_hd = oned < THR)
  dat.nohd <- filter(dat, !likely_hd)
  ##
  ## Group chemistry plates, excluding homozygous deletions
  ##
  ix <- sample(seq_len(nrow(dat.nohd)), 1000, replace=TRUE)
  message("Downsampling non-homozygous deletions and identify batches from surrogates")
  mb.subsamp <- dat.nohd[ix, ] %>%
    bind_rows(filter(dat, likely_hd)) %>%
    mutate(is_simulated=FALSE) %>%
    MultiBatch("MB3", data=.) %>%
    findSurrogates(0.001, THR)

  ##print(a)
  batches <- assays(mb.subsamp) %>%
    group_by(provisional_batch) %>%
    summarize(batch=unique(batch))
  pr.batch <- assays(mb.subsamp)$provisional_batch
  stopifnot(all(pr.batch %in% dat$provisional_batch))
  dat <- left_join(dat, batches, by="provisional_batch")
  ##
  ## We need the number in `hdmean` later. Where to keep it?
  ##
  hdmean <- median(dat$oned[dat$likely_hd])
  ##expect_equal(hdmean, -3.887, tolerance=0.001)
  ##
  message("Check batches")
  ##
  batchfreq <- assays(mb.subsamp) %>%
    group_by(batch) %>%
    summarize(n=n())
  if(any(batchfreq$n < 50)){
    batchfreq <- filter(batchfreq, n < 50)
    adat <- assays(mb.subsamp) %>%
      filter(!batch %in% batchfreq$batch)
    bdat <- filter(dat, batch %in% batchfreq$batch)
    adat2 <- bind_rows(adat, bdat)
    mb.subsamp <- MultiBatch("MB2", data=adat2)
  }
  mb.subsamp
}

genotype_model <- function(mod_1.3, snpdat){
  pz <- probz(mod_1.3) %>%
    "/"(rowSums(.))
  max_zprob <- rowMax(pz)
  is_high_conf <- max_zprob > 0.95
  ##mean(is_high_conf)
  gmodel <- mod_1.3[ !isSimulated(mod_1.3) &  is_high_conf ]
  keep <- !duplicated(id(gmodel))
  gmodel <- gmodel[ keep ]
  snpdat2 <- snpdat[, id(gmodel) ]
  ##identical(colnames(snpdat), id(final.model))
  clist <- CnList(gmodel)
  stats <- baf_loglik(clist, snpdat2)
  mapping(gmodel) <- strsplit(stats$cn.model[1], ",")[[1]]
  mapping(mod_1.3) <- mapping(gmodel)
  mod_1.3
}

join_baf_oned <- function(mod_1.3, snpdat){
  xlimit <- range(oned(mod_1.3))
  if(diff(xlimit) < 4){
    xlimit <- c(-3, 1)
  }
  maxpz <- probz(mod_1.3) %>%
    "/"(rowSums(.)) %>%
    rowMax
  bafdat <- assays(snpdat)[["baf"]] %>%
    as_tibble() %>%
    mutate(rsid=rownames(snpdat)) %>%
    gather("id", "BAF", -rsid)
  ##
  ## tibble of copy number probabilities
  ##
  gmodel <- mod_1.3
  cndat <- tibble(id=id(gmodel),
                  batch=batch(gmodel),
                  oned=oned(gmodel),
                  pz=maxpz,
                  z=map_z(gmodel)) %>%
    mutate(cn=mapping(gmodel)[z]) %>%
    mutate(cn=factor(cn))
  bafdat <- left_join(bafdat, cndat, by="id")
  bafdat2 <- filter(bafdat, pz > 0.9)
  bafdat2
}

component_labels <- function(model){
  cnlabels <- paste0(seq_len(k(model)),
                     "%->%",
                     mapping(model))
  labs <- as.expression(parse(text=cnlabels))
  labs
}

list_mixture_plots <- function(mod_1.3, bafdat2){
  A <- ggMixture(mod_1.3) +
    xlab(expression(paste("Median ", log[2], " R ratio"))) +
    ylab("Density\n")
  ## predictive densities excluding simulated data
  A2 <- ggMixture(mod_1.3[ !isSimulated(mod_1.3) ]) +
    xlab(expression(paste("Median ", log[2], " R ratio"))) +
    ylab("Density\n")
  ##
  ## PLot BAFs
  ##
  labs <- component_labels(mod_1.3)
  xlab <- expression(paste("\n", "Mixture component"%->%"Copy number"))
  legtitle <- "Mixture\ncomponent\nprobability"
  B <- ggplot(bafdat2, aes(factor(z), BAF)) +
    geom_hline(yintercept=c(0, 1/3, 0.5, 2/3, 1), color="gray95") +
    geom_jitter(aes(color=pz), width=0.1, size=0.3) +
    scale_y_continuous(expand=c(0, 0.05)) +
    scale_x_discrete(breaks=seq_len(k(gmodel)),
                     labels=labs) +
    theme(panel.background=element_rect(fill="white", color="gray30"),
          legend.key=element_rect(fill="white")) +
    xlab(xlab) +
    ylab("BAF\n") +
    guides(color=guide_legend(title=legtitle))
  list(augmented=A, ## observed + augmented data
       observed=A2, ## observed data only
       baf=B)
}

mixture_layout <- function(figure_list, augmented=TRUE){
  if(augmented) {
    A <- figure_list[["augmented"]]
  } else A <- figure_list[["observed"]]
  B <- figure_list[["baf"]]
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1, 2,
                                           widths=c(0.6, 0.4))))
  pushViewport(viewport(layout.pos.row=1,
                        layout.pos.col=1))
  pushViewport(viewport(width=unit(0.96, "npc"),
                        height=unit(0.9, "npc")))
  print(A, newpage=FALSE)
  popViewport(2)
  pushViewport(viewport(layout.pos.row=1,
                        layout.pos.col=2))
  pushViewport(viewport(width=unit(0.96, "npc"),
                        height=unit(0.6, "npc")))
  print(B, newpage=FALSE)
}
