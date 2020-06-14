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
          if(!"is_simulated" %in% colnames(data)){
              data$is_simulated <- FALSE
          }
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
  if(length(obj@.internal.counter)==0)
    obj@.internal.counter=0L
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


setMethod("assays", "MultiBatch", function(x, withDimnames, ...) x@data)

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
      set_colnames(c("batch", "z")) %>%
      as_tibble() %>%
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
    ##message("  k: ", K, ", burnin: ", burnin(mp), ", thin: ", thin(mp))
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
    if(length(unique(B)) == 1) return(list(pval=1, batches=B))
    B2 <- B
    dat <- tibble(x=x, batch=B) %>%
        group_by(batch) %>%
        summarize(n=n()) %>%
        arrange(-n)
    uB <- unique(dat$batch)
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

find_surrogates <- function(dat, THR=0.1){
  ## do not define batches based on homozygous deletion
    ##dat2 <- filter(dat, oned > min_oned)
    dat2 <- filter(dat, !likely_deletion)
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
          function(object, THR=0.1){
    dat <- assays(object) %>%
        select(c(id, provisional_batch, oned, likely_deletion))
    ##message("Putting rows of data in batch order")
    ##dat2 <- find_surrogates(dat, THR, min_oned) %>%
    dat2 <- find_surrogates(dat, THR) %>%    
        mutate(batch=factor(batch, levels=unique(batch)),
               batch_labels=as.character(batch),
               batch=as.integer(batch)) %>%
        filter(!duplicated(id)) %>%
        arrange(batch) %>%
        select(c(provisional_batch, batch, batch_labels, likely_deletion))  %>%
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
    batch_mapping <- select(dat2, c(provisional_batch, batch))
    ## Remove the previous batch assignment that was arbitrary
    full.data <- assays(object) %>%
        select(-batch)
    ## update batch assignment by joining
    full.data2 <- full.data %>%
        left_join(batch_mapping, by="provisional_batch") %>%
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

setMethod("findSurrogates", "tbl_df",
          function(object, THR=0.1){
    dat <- object %>%
        select(c(id, provisional_batch, oned, likely_deletion))
    ##message("Putting rows of data in batch order")
    ##dat2 <- find_surrogates(dat, THR, min_oned) %>%
    dat2 <- find_surrogates(dat, THR) %>%    
        mutate(batch=factor(batch, levels=unique(batch)),
               batch_labels=as.character(batch),
               batch=as.integer(batch)) %>%
        filter(!duplicated(id)) %>%
        arrange(batch) %>%
        select(c(provisional_batch, batch, batch_labels, likely_deletion))  %>%
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
    batch_mapping <- select(dat2, c(provisional_batch, batch))
    ## Remove the previous batch assignment that was arbitrary
    full.data <- object %>%
        select(-batch)
    ## update batch assignment by joining
    full.data2 <- full.data %>%
        left_join(batch_mapping, by="provisional_batch") %>%
        arrange(batch) 
    full.data2
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
  ##means <- computeMeans(x)
  ##precs <- computePrec(x)
  ##ps <- computeMixProbs(x)
  ##  if(any(is.na(ps))) ps <- p(x)
  ##  if(any(is.na(means))) means <- theta(x)
  ##  if(any(is.na(precs))) precs <- 1/sigma2(x)
  current_values(x)[["theta"]] <- theta(x)
  current_values(x)[["sigma2"]] <- sigma2(x)
  current_values(x)[["p"]] <- p(x)
  summaries(x)[["data.mean"]] <- theta(x)
  summaries(x)[["data.prec"]] <- 1/sigma2(x)
  if(L2 == nbatch1) {
    ## leave chains alone and return
    return(x)
  }
  ## can we keep the chains for batches still in the model
  B <- matrix(seq_len(nbatch1), nbatch1, k(x)) %>%
    as.numeric
  keep <- B %in% ubatch
  ch <- chains(x)
  ch@theta <- theta(ch)[, keep, drop=FALSE]
  ch@pi <- p(ch)[, keep, drop=FALSE]
  ch@predictive <- predictive(ch)[, keep, drop=FALSE]
  ch@zstar <- ch@zstar[, keep, drop=FALSE]
  if(substr(modelName(x), 1, 3) == "MBP"){
    ##chains(x) <- initialize_mcmcP(k(x), iter(x), numBatch(x))
  } else {
    ##chains2 <- initialize_mcmc(k(x), iter(x), numBatch(x))
    ch@sigma2 <- sigma2(ch)[, keep, drop=FALSE]
  }
  x@chains <- ch
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
    rowMaxs
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

homozygousdel_var <- function(object, THR=-1) {
    v <- var(oned(object)[ oned(object) < THR])
    if(is.na(v)) v <- sqrt(0.3)
    v
}

isPooledVar <- function(object) ncol(sigma2(object))==1

meanSdHomDel <- function(object, THR){
  list(homozygousdel_mean(object, THR),
       sd(oned(object)[ oned(object) < THR]))
}

sdRestrictedModel <- function(restricted){
  if(isPooledVar(restricted)){
    sigma2_ <- sigma2(restricted)
  } else {
    sigma2_ <- cbind(sigma2(restricted)[1, 2],
                     sigma2(restricted))
  }
  sqrt(sigma2_)
}

.augment_homozygous <- function(mb, mean_sd, phat=0.01){
    if(all(is.na(mean_sd[[2]]))) mean_sd[[2]] <- 0.1
    freq.hd <- assays(mb) %>%
        filter(!is_simulated) %>%
        group_by(batch) %>%
        summarize(N=n(),
                  n=sum(likely_deletion)) %>%
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

deletion_midpoint <- function(mb){
    vals <- assays(mb) %>%
        group_by(likely_deletion) %>%
        summarize(min_oned=min(oned, na.rm=TRUE),
                  max_oned=max(oned, na.rm=TRUE))
    midpoint <- mean(c(vals$max_oned[vals$likely_deletion],
                       vals$min_oned[!vals$likely_deletion]))
    
    midpoint
}

augment_homozygous <- function(mb.subsamp){
    ##THR <- summaries(mb.subsamp)$deletion_cutoff
    THR <- deletion_midpoint(mb.subsamp)
    mean_sd <- meanSdHomDel(mb.subsamp, THR)
    rare_homozygous <- sum(oned(mb.subsamp) < THR) < 5
    ##expect_false(rare_homozygous)
    if(rare_homozygous){
        simdat <- .augment_homozygous(mb.subsamp, mean_sd)
    } else {
        simdat <- assays(mb.subsamp) %>%
            arrange(batch) %>%
            mutate(homozygousdel_mean=mean_sd[[1]])
    }
    simdat
}

ok_hemizygous <- function(sb){
  varratio <- max(sigma2(sb))/min(sigma2(sb))
  !(p(sb)[2] < 0.1 || varratio > 100)
}

.warmup <- function(tib, model, Nrep=10, .burnin=100){
    ##if(model=="MBP2") browser()
    mbl <- replicate(Nrep, MultiBatchList(data=tib)[[model]])
    for(j in seq_along(mbl)){
        cat(".")
        mb <- mbl[[j]]
        iter(mb) <- 0
        burnin(mb) <- .burnin
        mb <- tryCatch(posteriorSimulation(mb),
                       warning=function(w) NULL)
        if(is.null(mb)) next()
        mbl[[j]] <- mb
    }
    mbl
}

warmup <- function(tib, model1, model2=NULL, model2.penalty=50,
                   Nrep=10, .burnin=100){
    ##
    ##
    message("Warmup with ", model1)
    mbl <- .warmup(tib, model1, Nrep=Nrep, .burnin=.burnin)
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
        message("Warmup with ", model2)        
        mbl2 <- .warmup(tib, model2, Nrep=Nrep, .burnin=.burnin)
        ml2 <- sapply(mbl2, log_lik)
      if(is(ml2, "list")){
          mbl2 <- mbl2[ lengths(ml2) > 0 ]
          ml2 <- ml2[ lengths(ml2) > 0 ]
          ml2 <- unlist(ml2)
      }
        if(all(is.na(ml2))) return(model)
        if(max(ml2, na.rm=TRUE) > max(ml, na.rm=TRUE) + model2.penalty){
            model <- mbl2[[which.max(ml2)]]
        } else {
            model <- mbl[[which.max(ml)]]
        }
    }
    return(model)
}

stop_early <- function(model, min_prob=0.99, prop_greater=0.995){
  pz <- probz(model)
  maxprob <- rowMaxs(pz)
  pz <- probz(model) %>%
    "/"(rowSums(.)) %>%
    rowMaxs
  mean_maxp <- mean(pz > min_prob)
  mean_maxp > prop_greater
}

revertToMultiBatch <- function(sb){
  adat <- assays(sb)
  batch_levels <- unique(assays(sb)$batch_labels)
  adat$batch <- as.integer(factor(assays(sb)$batch_labels,
                                  levels=batch_levels))
  model_name <- gsub("SB", "MB", modelName(sb))
  mb <- MultiBatch(model_name, data=adat)
  mcmcParams(mb) <- mcmcParams(sb)
  mb
}

modal_theta <-  function(model){
  p_ <- p(model)
  th_ <- rep(NA, nrow(p_))
  TH <- theta(model)
  for(i in seq_along(th_)){
    j <- which.max(p_[i, ])
    th_[i] <- TH[i, j]
  }
  th_
}

## rare_component: of restricted model
## diploid_component: of SB model
augment_rarecomponent <- function(restricted,
                                  sb,
                                  rare_component_restricted=1,
                                  rare_component_sb=2,
                                  diploid_component_sb=3,
                                  use_restricted_theta=FALSE){
    batch_labels <- batchLevels(restricted)
    ## normalize probabilities
    pz <- probz(restricted) %>%
        "/"(rowSums(.)) %>%
        rowMaxs
    ##
    ## batches with high posterior probabilities
    ##
    ubatch <- unique(batch(restricted)[ pz >= 0.95 ])
    ##
    ## Find number of samples assigned to the rare component with high
    ## probability.
    ##
    ## Check if any of the batches have fewer than 10 subjects
    ## with high posterior probability
    ##
    tmp <- tibble(batch=batch(restricted),
                  z=map_z(restricted),
                  pz=pz) %>%
        group_by(batch) %>%
        summarize(N=n(),
                  n=sum(z==rare_component_restricted & pz > 0.9)) %>%
        filter(n < 10)
    is_dropped <- !batch_labels %in% ubatch |
        batch_labels %in% tmp$batch
    any_dropped <- any(is_dropped)
    if(!any_dropped) return(restricted)
    ##
    ## There are batches with fewer than 10 samples
    ## assigned to mixture component with probability > 0.9
    ##
    dropped_batches <- uniqueBatch(restricted)[ is_dropped ]
    if(ncol(sigma2(sb)) == 1){
        component_var <- sigma2(sb)[, 1]
    } else {
        ##
        ## estimate variance from diploid mixture component
        ##
        component_var <- sigma2(sb)[, diploid_component_sb]
    }
    msg <- "Find mean and variance of rare component"
    ##message(msg)
    i <- dropped_batches
    ##        if(use_restricted_theta){
    ##            j <- rare_component_restricted
    ##            theta_ <- median(theta(restricted)[, j])
    ##            diploid_component_restricted <-  2
    ##            theta_diploid <- modal_theta(restricted)
    ##            delta <- theta_diploid - median(theta_diploid)
    ##            theta_ <- (theta_ + delta)[dropped_batches]
    ##            p_ <- pmax(median(p(restricted)[, j]), 0.05)
    ##        } else {
    j <- rare_component_sb
    theta_ <- theta(sb)[j]
    p_ <- pmax(p(sb)[j], 0.05)
    loc.scale <- tibble(theta=theta_,
                        sigma2=component_var,
                        phat=p_,
                        batch=dropped_batches)
    ##message("Augment data with additional hemizygous deletions")
    start_index <- length(grep("augment", id(restricted))) + 1
    impdat <- impute(restricted, loc.scale, start.index=start_index)
    ## removes all simulated and then adds back only the
    ## data simulated above
    simdat <- bind_rows(assays(restricted),
                        impdat) %>%
        arrange(batch)
    mb <- MultiBatchList(data=simdat)[[modelName(restricted)]]
    mcmcParams(mb) <- mcmcParams(restricted)
    theta(mb) <- theta(restricted)

    i_ <- is_dropped
    j_ <- rare_component_restricted
    theta(mb)[i_, j_] <- loc.scale$theta
    if(ncol(sigma2(restricted)) == 1){
        j <- 1
    } else{
        j <- rare_component_restricted
    }
    ##component_var only defined in conditional statement above
    sigma2(mb) <- matrix(pmin(sigma2(restricted)[, j],
                              component_var),
                         ncol=1)
    mb
}

augment_rareduplication <- function(sb3,
                                    mod_2.4,
                                    full_data,
                                    THR){
    densities <- compute_density(mod_2.4, THR)
    modes <- round(compute_modes(densities), 3)
    loc.scale <- tibble(theta=theta(sb3)[3],
                        sigma2=sigma2(sb3),
                        phat=max(p(sb3)[1, 3], 0.05),
                        batch=seq_len(numBatch(mod_2.4)))
    ## shift the batches according to location of mode
    loc.scale$theta <- loc.scale$theta + modes
    start.index <- length(grep("augment", id(mod_2.4))) + 1
    imp.dup <- impute(mod_2.4, loc.scale,
                      start.index=start.index) %>%
        mutate(likely_deletion=FALSE)
    obsdat <- full_data %>%
        mutate(is_simulated=FALSE)
    simdat <- bind_rows(obsdat, imp.dup) %>%
        arrange(batch)
    simdat
}

augment_rarehemdel <- function(sb3,
                               mod_2.4,
                               full_data){
    loc.scale <- tibble(theta=theta(sb3)[1],
                        sigma2=sigma2(sb3),
                        phat=max(p(sb3)[1, 1], 0.05),
                        batch=seq_len(numBatch(mod_2.4)))
    ##densities <- compute_density(mod_2.4)
    ##modes <- round(compute_modes(densities), 3)
    ## shift the batches according to location of mode
    ##loc.scale$theta <- loc.scale$theta + modes
    start.index <- length(grep("augment", id(mod_2.4))) + 1
    imp.dup <- impute(mod_2.4,
                      loc.scale,
                      start.index=start.index) %>%
        mutate(likely_deletion=FALSE)
    obsdat <- full_data %>%
        mutate(is_simulated=FALSE)
    simdat <- bind_rows(filter(assays(mod_2.4), is_simulated),
                        imp.dup)
    dat <- bind_rows(obsdat, simdat) %>%
        arrange(batch)
    dat
}

augment_rarehomdel <- function(restricted, sb4, mb.subsamp){
    p_ <- cbind(p(sb4)[1, 1], p(restricted)) %>%
        "/"(rowSums(.))
    dat <- assays(mb.subsamp)
    hdmean <- median(dat$oned[dat$likely_deletion])
    hdvar <- var(dat$oned[dat$likely_deletion])
    if(is.na(hdmean)) hdmean <- -4
    theta_ <- cbind(hdmean, theta(restricted))
    is_pooledvar <- ncol(sigma(restricted)) == 1
    if(is_pooledvar){
        sigma2_ <- sigma2(restricted)
    } else {
        sigma2_ <- cbind(sigma2(restricted)[1, 2],
                         sigma2(restricted))
    }
    freq.hd <- assays(mb.subsamp) %>%
        group_by(batch) %>%
        summarize(N=n(),
                  n=sum(likely_deletion)) %>%
        filter(n/N < 0.05)
    if(nrow(freq.hd) > 0){
        loc.scale <- tibble(theta=hdmean,
                            sigma2=sigma2_[, 1],
                            phat=max(p(sb4)[1], 0.05),
                            batch=seq_len(nrow(theta_)))
        loc.scale <- left_join(freq.hd, loc.scale, by="batch") %>%
            select(-N)
        start.index <- length(grep("augment", id(restricted))) + 1
        imp.hd <- impute(restricted, loc.scale, start.index=start.index)
        if(any(isSimulated(restricted))){
            imp1 <- filter(assays(restricted), is_simulated)
            imp.hd <- bind_rows(imp.hd, imp1)
        }
        obsdat <- assays(mb.subsamp) %>%
            mutate(is_simulated=FALSE)
        simdat <- bind_rows(obsdat, imp.hd) %>%
            arrange(batch)
    } else {
        imp1 <-filter(assays(restricted), is_simulated)
        simdat <- bind_rows(assays(mb.subsamp),
                            imp1) %>%
            arrange(batch)
    }
    simdat
}

batchLevels <- function(mb){
  labels <- assays(mb)$batch_labels
  levs <- unique(labels)
  levels <- unique(as.integer(factor(labels, levels=levs)))
  sort(levels [ !is.na(levels) ])
}

.mcmcWithHomDel <- function(simdat, mod_2.3){
    mp <- mcmcParams(mod_2.3)
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
    tmp <- tryCatch(posteriorSimulation(mod_1.3),
                    warning=function(w) NULL)
    bad_start <- FALSE
    if(is.null(tmp)){
        bad_start <- TRUE
    }
    if(!is.null(tmp)){
        if(is.nan(log_lik(tmp)))
            bad_start <- TRUE
    }
    if(bad_start){
        adat <- assays(mod_1.3)
        ##fdat <- filter(adat, oned > -1, !is_simulated)
        fdat <- filter(adat, !likely_deletion, !is_simulated) 
        mod_1.3 <- warmup(fdat, model)
    } else mod_1.3 <- tmp
    if(is.null(mod_1.3)) return(NULL)
    internal.count <- flags(mod_1.3)$.internal.counter
    any_dropped <- TRUE
    if(internal.count < 100){
        ##iter(mod_1.3) <- 1000
        mcmcParams(mod_1.3) <- mp
        mod_1.3 <- posteriorSimulation(mod_1.3)
    }
    mod_1.3
}


meanFull_sdRestricted <- function(mb, restricted){
    hdmean <- filter(assays(mb), likely_deletion) %>%
        pull(oned) %>%
        mean(na.rm=TRUE)
    restricted.sd <- sdRestrictedModel(restricted)
    mn_sd <- list(hdmean, restricted.sd)
    mn_sd
}


mcmcWithHomDel <- function(mb, sb,
                           restricted_model,
                           THR,
                           model="MB3"){
    mb.observed <- mb[ !isSimulated(mb) ]        
    if(any(isSimulated(restricted_model))){
        ##
        ## Bind the simulated observations from the restricted model
        ## to the observed data from the full model
        ##
        ##
        mb1 <- filter(assays(restricted_model),
                      is_simulated) %>%
            bind_rows(assays(mb.observed)) %>%
            arrange(batch) %>%
            MultiBatchList(data=.) %>%
            "[["(model)
    } else {
        ## If there were no simulated observtions in the restricted model,
        ## mb1 is just the observed data in the full model
        mb1 <- MultiBatchList(data=assays(mb.observed))[[model]]
    }
    ## If at least one batch has fewer than 5% subjects with
    ## homozygous deletion, augment the data for homozygous deletions
    mn_sd <- meanFull_sdRestricted(mb, restricted_model)
    simdat <- .augment_homozygous(mb1, mn_sd,
                                  phat=max(p(sb)[1], 0.05))  
    full <- .mcmcWithHomDel(simdat, restricted_model)
    full
}


mcmcHomDelOnly <- function(simdat,
                           restricted_model,
                           sb,
                           model){
  mp <- mcmcParams(sb)
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
    singlebatchvar <- sigma2(sb)[, 1]
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
  mcmcParams(mod_1.3) <- mp
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
    rowMaxs
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

hemizygous_cutoff <- function(dat){
  ## is there a peak less than -0.2
  dens <- density(dat$oned[ dat$oned < -0.2] )
  firstderivative <- diff(dens$y)
  changesign <- which(diff(sign(firstderivative)) != 0)
  if(length(changesign) == 0){
    return(-Inf)
  }
  ## choose the one with the maximal y
  index <- changesign[which.max(dens$y[-1][changesign])]
  if(FALSE){
    plot(dens)
    abline(v=dens$x[index])
  }
  peak <- dens$x[index]
  xy <- tibble(x=dens$x, y=dens$y) %>%
    filter(x > peak)
  lowest_point_after_peak <- xy$x[xy$y==min(xy$y)]
  lowest_point_after_peak
}

duplication_cutoff <- function(dat, min_cutoff=0.1){
  ## is there a peak greater than 0.1
  dens <- density(dat$oned[ dat$oned > min_cutoff] , adjust=1/2)
  firstderivative <- diff(dens$y)
  changesign <- which(diff(sign(firstderivative)) != 0)
  if(length(changesign) == 0){
    return(+Inf)
  }
  ## choose the one with the maximal y
  changesign <- changesign[-1]
  index <- changesign[which.max(dens$y[-1][changesign])]
  if(FALSE){
    plot(dens)
    abline(v=dens$x[index])
  }
  peak <- dens$x[index]
  xy <- tibble(x=dens$x, y=dens$y) %>%
    filter(x < peak, x > min_cutoff)
  lowest_point_before_peak <- xy$x[xy$y==min(xy$y)]
  lowest_point_before_peak
}

#' @export
median_summary <- function(se, provisional_batch, assay_index=1, THR){
    medians <- colMedians(assays(se)[[assay_index]], na.rm=TRUE)
    dat <- tibble(id=colnames(se),
                  oned=medians,
                  provisional_batch=provisional_batch,
                  batch=1,
                  batch_labels="1") %>%
        mutate(likely_deletion = oned < THR)
    ## if no homozygous deletions, check for hemizygous deletions
    ##  if(!any(dat$likely_deletion)){
    ##    THR <- hemizygous_cutoff(dat)
    ##    dat$likely_deletion <- dat$oned < THR
    ##  }
    dat
}

resampleFromRareProvisionalBatches <- function(dat, dat.nohd){
  prov.batch <- unique(dat.nohd$provisional_batch)
  all.prov.batch <- unique(dat$provisional_batch)
  notsampled <- all.prov.batch[ !all.prov.batch %in% prov.batch ]
  ## just include all samples from these batches
  ix <- which(dat$provisional_batch %in% notsampled)
  dat2 <- dat[ix, ]
  dat.nohd <- bind_rows(dat.nohd, dat2) %>%
    filter(!likely_deletion)
  dat.nohd
}

#' @export
down_sample <- function(dat, S){
    dat.nohd <- filter(dat, !likely_deletion)
    ##
    ## Group chemistry plates, excluding homozygous deletions
    ##
    S <- min(nrow(dat.nohd), S)
    ix <- sample(seq_len(nrow(dat.nohd)), S, replace=TRUE)
    dat.nohd <- dat.nohd[ix, ]
    ##
    ## ensure all provisional batches were sampled
    ##
    dat.nohd <- resampleFromRareProvisionalBatches(dat, dat.nohd)
    dat.subsampled <- dat.nohd %>%
        bind_rows(filter(dat, likely_deletion)) %>%
        mutate(is_simulated=FALSE)
    dat.subsampled
}

down_sample2 <- function(dat, S, min_size=100){
    ##
    ## Goals:
    ##  -- we do not want to downsample likely deletions; these tend to be rare
    ##  -- we do not want to downsample batches for which there is little data
    ##  -- the number of observations in each batch should roughly reflect the overall proportions
    ##
    ## Based on S and the number per group,
    ## computed the expected number of observations after downsampling
    expected <- group_by(dat, batch) %>%
        tally() %>%
        mutate(p=n/sum(n),
               expected=p*S)  %>%
        select(batch, expected)
    dat2 <- left_join(dat, expected, by="batch")
    holdout1 <- filter(dat2, likely_deletion)
    holdout2 <- filter(dat2, expected < min_size, !likely_deletion)
    if(nrow(holdout2) > 0){
        ## if min_size is greater than number of observations in group,
        ## the result will be silently truncated to the group size
        holdout2.resampled <- holdout2 %>%
            group_by(batch) %>%
            slice_sample(n=min_size) %>%
            ungroup()

    } else holdout2.resampled <- holdout2
    H <- nrow(holdout2.resampled)
    holdouts <- bind_rows(holdout1, holdout2)
    downsamp <- filter(dat2, !id %in% holdouts$id) %>%
        slice_sample(n=S-H)
    downsamp2 <- holdout2.resampled %>%
        bind_rows(downsamp) %>%
        bind_rows(holdout1) %>%
        arrange(batch) %>%
        mutate(batch_labels=as.character(batch)) %>%
        select(-expected)
    downsamp2
}


##kolmogorov_batches <- function(dat, KS_cutoff, THR){
##    mb <- MultiBatch("MB3", data=dat)
##    mb.ks <- mb %>%
##        findSurrogates(KS_cutoff)
##}

kolmogorov_batches <- function(dat, KS_cutoff){
    ##mb <- MultiBatch("MB3", data=dat)
    findSurrogates(dat, KS_cutoff)
}

add_batchinfo <- function(dat, mb){
    batches <- assays(mb) %>%
        group_by(provisional_batch) %>%
        summarize(batch=unique(batch))
    pr.batch <- assays(mb)$provisional_batch
    stopifnot(all(pr.batch %in% dat$provisional_batch))
    dropbatch <- select(dat, -batch)
    dat <- left_join(dropbatch, batches, by="provisional_batch")
    dat
}

add2small_batches <- function(dat, mb, min_size=50){
    ##message("Check batches")
    ##
    batchfreq <- assays(mb) %>%
        group_by(batch) %>%
        summarize(n=n())
    if(all(batchfreq$n >= min_size)) return(mb)
    batchlabels <- group_by(assays(mb), batch) %>%
        summarize(n=n(),
                  batch_labels=unique(batch_labels))
    batchfreq <- filter(batchfreq, n < min_size)
    adat <- assays(mb) %>%
        filter(!batch %in% batchfreq$batch)
    bdat <- filter(dat, batch %in% batchfreq$batch) %>%
        left_join(batchlabels, by="batch") %>%
        mutate(is_simulated=FALSE) %>%
        select(colnames(adat))
    adat2 <- bind_rows(adat, bdat)
    mb <- MultiBatch("MB2", data=adat2)
    mb
}

add_deletion_stats <- function(dat, mb, THR){
    hdmean <- median(dat$oned[dat$likely_deletion])
    ##expect_equal(hdmean, -3.887, tolerance=0.001)
    ##
    if(is.na(hdmean)) THR <- hdmean <- -1
    assays(mb)$homozygousdel_mean <- hdmean
    summaries(mb)$deletion_cutoff <- deletion_midpoint(dat)
    mb
}


#' @export
summarize_region <- function(se,
                             provisional_batch,
                             THR=-1,
                             assay_index=1,
                             KS_cutoff=0.001,
                             S=1000,
                             min_size=50){
    dat <- median_summary(se,
                          provisional_batch,
                          assay_index=assay_index,
                          THR=THR)
    ##dat2 <- down_sample(dat, S)
    dat2 <- kolmogorov_batches(dat, KS_cutoff)
    ##dat <- add_batchinfo(dat, mb)
    ##mb <- add2small_batches(dat, mb)
    ##mb <- add_deletion_stats(dat, mb, THR)
    dat3 <- down_sample2(dat2, S, min_size)
    dat3
}

## select high confidence samples for each component
## - assume component corresponds to same copy number state across batches
select_highconfidence_samples <- function(model, snpdat){
    if(iter(model) == 0) stop("No saved MCMC simulations. Must specify iter > 0")
    model2 <- dropSimulated(model)
    pz <- probz(model2) %>%
        "/"(rowSums(.))
    maxpz <- tibble(prob=rowMaxs(pz),
                    batch=batch(model2),
                    z=map_z(model2))
    cutoffs <- group_by(maxpz, z) %>%
        summarize(cutoff=max(quantile(prob, 0.75), 0.8))
        ##mutate(cutoff=ifelse(cutoff > 0.9, 0.9, cutoff))
    maxpz <- left_join(maxpz, cutoffs, by="z")
    highconf <- group_by(maxpz, z, batch)  %>%
        summarize(total=n(),
                  n=sum(prob >= cutoff))
    if(!all(highconf$n >= 5)){
        mapping(model) <- rep("?", k(model))
        return(model)
    }
    is_highconf <- maxpz$prob >= maxpz$cutoff
    ##keep <- !isSimulated(model) & is_highconf
    gmodel <- model2[ is_highconf ]
    gmodel
}

select_highconfidence_samples2 <- function(model, snpdat){
    if(iter(model) == 0) stop("No saved MCMC simulations. Must specify iter > 0")
    model2 <- dropSimulated(model)
    pz <- probz(model2) %>%
        "/"(rowSums(.))
    maxpz <- tibble(prob=rowMaxs(pz),
                    batch=batch(model2),
                    z=map_z(model2))
##    
##    cutoffs <- group_by(maxpz, z) %>%
##        summarize(cutoff=max(quantile(prob, 0.75), 0.8))
##        ##mutate(cutoff=ifelse(cutoff > 0.9, 0.9, cutoff))
##    maxpz <- left_join(maxpz, cutoffs, by="z")
    highconf <- group_by(maxpz, z, batch)  %>%
        summarize(total=n(),
                  n=sum(prob >= 0.75))
    if(!all(highconf$n >= 1)){
        mapping(model) <- rep("?", k(model))
        return(model)
    }
    is_highconf <- maxpz$prob >= 0.75
    ##keep <- !isSimulated(model) & is_highconf
    gmodel <- model2[ is_highconf ]
    gmodel
}


#' @export
genotype_model <- function(model, snpdat){
  gmodel <- select_highconfidence_samples2(model, snpdat)
  if(all(mapping(gmodel) == "?")) {
    mapping(gmodel) <- rep("2", k(gmodel))
    return(gmodel)
  }
  keep <- !duplicated(id(gmodel))
  gmodel <- gmodel[ keep ]
  snpdat2 <- snpdat[, id(gmodel) ]
  ##identical(colnames(snpdat), id(final.model))
  clist <- CnList(gmodel)
  stats <- baf_loglik(clist, snpdat2)
  mapping(gmodel) <- strsplit(stats$cn.model[1], ",")[[1]]
  mapping(model) <- mapping(gmodel)
  summaries(model)$baf_loglik <- stats
  model
}

join_baf_oned <- function(model, snpdat){
  ## Plot BAFs for samples with high confidence
  ##model2 <- select_highconfidence_samples(model, snpdat)
  xlimit <- range(oned(model))
  if(diff(xlimit) < 4){
    xlimit <- c(-3, 1)
  }
  maxpz <- probz(model) %>%
    "/"(rowSums(.)) %>%
    rowMaxs
  bafdat <- assays(snpdat)[["baf"]] %>%
      as_tibble() %>%
      mutate(rsid=rownames(snpdat)) %>%
      gather("id", "BAF", -rsid)
  ##
  ## tibble of copy number probabilities
  ##
  gmodel <- model
  cndat <- tibble(id=id(gmodel),
                  batch=batch(gmodel),
                  oned=oned(gmodel),
                  pz=maxpz,
                  z=map_z(gmodel)) %>%
    mutate(cn=mapping(gmodel)[z]) %>%
    mutate(cn=factor(cn))
  bafdat <- left_join(bafdat, cndat, by="id") %>%
    filter(!is.na(cn))
  ## Want cutoff that allows plotting of all states but that removes
  ## low quality SNPs
  model2 <- select_highconfidence_samples(model, snpdat)
  bafdat2 <- filter(bafdat, id %in% id(model2))
  bafdat2
}

#' @export
mixture_plot <- function(model, snpdat, xlimit=c(-4, 1), bins=100){
  bafdat <- join_baf_oned(model, snpdat)
  figs <- list_mixture_plots(model, bafdat, xlimit=xlimit, bins=bins)
  figs
}

component_labels <- function(model){
  cnlabels <- paste0(seq_len(k(model)),
                     "%->%",
                     mapping(model))
  labs <- as.expression(parse(text=cnlabels))
  labs
}

list_mixture_plots <- function(model, bafdat,
                               xlimit=c(-4, 1), bins=100){
  rg <- range(theta(model)[, 1])
  if(min(rg) < xlimit[1]){
    s <- mean(sigma(model)[, 1])
    xlimit[1] <- min(rg) - 2*s
  }
  A <- ggMixture(model, bins=bins) +
    xlab(expression(paste("Median ", log[2], " R ratio"))) +
    ylab("Density\n") +
    xlim(xlimit)
  ## predictive densities excluding simulated data
  A2 <- ggMixture(model[ !isSimulated(model) ], bins=bins) +
    xlab(expression(paste("Median ", log[2], " R ratio"))) +
    ylab("Density\n") +
    xlim(xlimit)
  ##
  ## PLot BAFs
  ##
  labs <- component_labels(model)
  xlab <- expression(paste("\n",
                           "Mixture component"%->%"Copy number"))
  legtitle <- "Mixture\ncomponent\nprobability"
  B <- ggplot(bafdat, aes(factor(z), BAF)) +
    geom_hline(yintercept=c(0, 1/3, 0.5, 2/3, 1), color="gray95") +
    geom_jitter(aes(color=pz), width=0.1, size=0.3) +
    scale_y_continuous(expand=c(0, 0.05)) +
    scale_x_discrete(breaks=seq_len(k(model)),
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

#' @export
mixture_layout <- function(figure_list, augmented=TRUE,
                           newpage=TRUE){
  if(augmented) {
    A <- figure_list[["augmented"]]
  } else A <- figure_list[["observed"]]
  B <- figure_list[["baf"]]
  if(newpage) grid.newpage()
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
  popViewport(2)
}

fit_restricted <- function(mb, sb, model="MBP2",
                           use_restricted=FALSE){
    ##fdat <- filter(assays(mb), oned > THR)  %>%
    fdat <- filter(assays(mb), !likely_deletion) %>%
        mutate(is_simulated=FALSE)
    warm <- warmup(fdat, model)
    mcmcParams(warm) <- mcmcParams(mb)
    restricted <- posteriorSimulation(warm)
    ##
    ## why check sb and not restricted?
    ##
    ok <- ok_hemizygous(sb)
    if(!ok) {
        restricted <- augment_rarecomponent(restricted,
                                            sb,
                                            use_restricted=use_restricted)
        restricted <- posteriorSimulation(restricted)
    }
    restricted
}

fit_restricted2 <- function(mb, model="MBP2", ...){
    ##    fdat <- filter(assays(mb), !likely_deletion) %>%
    ##        mutate(is_simulated=FALSE)
    warm <- warmup(assays(mb), model, ...)
    mcmcParams(warm) <- mcmcParams(mb)
    restricted <- posteriorSimulation(warm)
    restricted
}

explore_multibatch <- function(sb, model="MBP2", ...){
    mb <- revertToMultiBatch(sb)
    mbr <- assays(mb) %>%
        filter(!likely_deletion) %>%
        MultiBatch(data=.)
    mcmcParams(mbr) <- mcmcParams(mb)
    ok <- ok_hemizygous(sb)
    if(!ok) {
        mbr <- augment_rarecomponent(mbr, sb)
    }
    restricted <- fit_restricted2(mbr, model=model, ...)
    full <- mcmcWithHomDel(mb, sb, restricted)
    ok <- ok_model(full, restricted)
    if(!ok){
        model <- gsub("P", "", modelName(full))
        full <- mcmcHomDelOnly(assays(full), restricted, sb, model)
    }
    full
}

##explore_multibatch2 <- function(sb, model="MBP2"){
##    mb <- revertToMultiBatch(sb)
##    restricted <- fit_restricted2(mb, sb, model=model)
##    
##    ## augment_rareduplication?
##    ##message("Fitting full model")
##    full <- mcmcWithHomDel(mb, sb, restricted, THR)
##    ok <- ok_model(full, restricted)
##    if(!ok){
##        model <- gsub("P", "", modelName(full))
##        full <- mcmcHomDelOnly(assays(full), restricted, sb, model)
##    }
##    full
##}

few_hemizygous <- function(model){
  pz <- probz(model) %>%
    "/"(rowSums(.)) %>%
    rowMaxs
  tmp <- tibble(batch=batch(model), z=map_z(model), pz=pz) %>%
    group_by(batch) %>%
    summarize(N=n(),
              n=sum(z==2 & pz > 0.9))
  any(tmp$n < 10)
}

explore_restricted <- function(mb, sb, THR=-1, model="MB3"){
  if(few_hemizygous(mb))
    return(mb)
  mbp <- fit_restricted(mb, sb)
  restricted <- mcmcHomDelOnly(assays(mb),
                               restricted_model=mbp,
                               sb,
                               model=model)
  restricted
}

as_tibble.density <- function(x, ...){
  dens <- x
  result <- tibble(index=seq_along(dens$x),
                   x=dens$x,
                   y=dens$y)
}

compute_density <- function(mb, thr){
    batches <- batchLevels(mb)
    densities <- vector("list", length(batches))
    for(i in seq_along(batches)){
        ##
        ## Coarse
        ##
        index <- which(batch(mb)==i & oned(mb) > thr)
        if(length(index) > 2){
            tmp <- density(oned(mb)[index])
            tmp <- as_tibble(tmp)
        } else tmp <- NULL
        ##
        ## Fine
        ##
        tmp2 <- density(oned(mb)[batch(mb)==i & oned(mb) < thr]) %>%
            as_tibble()
        densities[[i]] <- list(coarse=as_tibble(tmp),
                               fine=as_tibble(tmp2))
    }
    densities
}

compute_modes <- function(densities){
  primary_mode <- sapply(densities, function(dens.list){
    dens <- dens.list[["fine"]]
    dens$x[which.max(dens$y)]
  })
  primary_mode
}

sample_homdel_from_density <- function(densities, modes, N=5){
  isnull_dens <- sapply(densities, function(x) nrow(x[[1]])==0)
  homdel_dens <- densities[!isnull_dens]
  ## just sample one if multiple available
  ix <- sample(seq_along(homdel_dens), 1)
  homdel_dens <- homdel_dens[[ix]][["coarse"]]
  imputed <- list()
  for(i in seq_along(densities)){
    if(!isnull_dens[[i]]) next()
    y <- sample(homdel_dens$x, N, prob=homdel_dens$y)
    adjust <- modes[i] - modes[ix]
    y <- y - adjust
    imputed[[i]] <- tibble(x=y, batch=i)
  }
  imputed <- do.call(rbind, imputed)
  return(imputed)
}

sample_hemdel_from_density <- function(dens, N=5){
  peaks <- hemizygous_peaks(dens)
  fine <- dens[["fine"]] %>%
    filter(x < max(peaks$maxx))
  x <- sample(fine$x, N, prob=fine$y)
  x
}

sample_from_distribution <- function(dat, modes, THR, N){
  del_dat <- filter(dat, oned < THR)
  mn <- mean(del_dat$oned)
  s <-  sd(del_dat$oned)*2
  if(missing(s)) s <- 0.2
  B <- length(unique(dat$batch))
  total <- N * B
  mn <- rep(mn, B)
  delta <- modes - median(modes)
  mn <- rep(mn - delta, each=N)
  x <- rnorm(length(mn), mn, s)
  tibble(x=x, batch=rep(B, each=N))
}

impute_needed <- function(mb, THR){
  dat <- assays(mb)
  if(!any(dat$oned < THR)) return(FALSE)
  tab <- filter(dat, oned < THR) %>%
    group_by(batch) %>%
    summarize(n=n())
  if(all(tab$n > 5)) return(FALSE)
  TRUE
}

impute_homozygous <- function(mb, modes, densities, THR){
  need2impute <- impute_needed(mb, THR)
  dat <- assays(mb)
  if(!need2impute) return(dat)
  isnull_dens <- sapply(densities, function(x) nrow(x[[1]])==0)
  if(any(!isnull_dens)) sample_from_density <- TRUE
  if(sample_from_density){
    x <- sample_homdel_from_density(densities, modes, N=5)
  } else {
    x <- sample_from_distribution(dat, THR)
  }
  nsim <- sum(isSimulated(mb))
  impdat <- tibble(id=paste0("augment_", seq_len(nrow(x)) + nsim),
                   oned=x$x,
                   provisional_batch=NA,
                   likely_deletion=TRUE,
                   batch=x$batch,
                   batch_labels=NA,
                   is_simulated=TRUE)
  dat2 <- bind_rows(dat, imputed) %>%
    arrange(batch)
  dat2
}

.peak <- function(dens, h){
  exceeds_h <- dens$y >= h
  regions <- cumsum(c(0, diff(exceeds_h) != 0))
  dens$region <- regions
  dens$h <- h
  peaks <- filter(dens, y > h) %>%
    group_by(region) %>%
    summarize(minx=min(x),
              maxx=max(x),
              miny=unique(h),
              maxy=max(y)) %>%
    mutate(h=h)
}

find_peaks <- function(dens, max_number=5, min_number=3){
  miny <- min(dens$y[dens$y > 0]) ## smallest y greater than 0
  maxy <- max(dens$y)
  heights <- seq(miny, maxy, length.out=100)
  peak_list <- vector("list", length(heights))
  for(i in seq_along(heights)){
    peak_list[[i]] <- .peak(dens, heights[i])
  }
  npeaks <- elementNROWS(peak_list)
  peak_list[ npeaks <= max_number & npeaks >= min_number]
}

hemizygous_peaks <- function(dens){
  fine <- dens[["fine"]]
  if(FALSE){
    ggplot(fine, aes(x, y, ymin=0, ymax=y)) +
      geom_ribbon()
  }
  peaks <- find_peaks(fine, 50, 2)
  npeaks <- elementNROWS(peaks)
  peaks <- peaks[npeaks == max(npeaks)]
  ix <- sapply(peaks, function(x) which.max(x$maxy))
  peaks <- peaks[ix > 1]
  peaks <- peaks[[length(peaks)]]
  peaks$index <- seq_len(nrow(peaks))
  diploid_index <- which(peaks$maxy==max(peaks$maxy))
  hemizygous_index <- peaks$index < diploid_index
  peaks[hemizygous_index, ]
}

homozygous_peaks <- function(dens){
  coarse <- dens[["coarse"]]
  if(FALSE){
    ggplot(fine, aes(x, y, ymin=0, ymax=y)) +
      geom_ribbon()
  }
  peaks <- find_peaks(coarse, 50, 1)
  npeaks <- elementNROWS(peaks)
  peaks <- peaks[npeaks == max(npeaks)]
  peaks <- peaks[[1]]
}

duplication_peaks <- function(dens){
  fine <- dens[["fine"]]
  if(FALSE){
    ggplot(fine, aes(x, y, ymin=0, ymax=y)) +
      geom_ribbon()
  }
  peaks <- find_peaks(fine, 50, 2)
  npeaks <- elementNROWS(peaks)
  peaks <- peaks[npeaks == max(npeaks)]
  peaks <- peaks[[1]]
  peaks$index <- seq_len(nrow(peaks))
  diploid_index <- which(peaks$maxy==max(peaks$maxy))
  duplication_index <- peaks$index > diploid_index
  peaks[duplication_index, ]
}

handle_missing <- function(fun, dens){
  peaks <- tryCatch(fun(dens), error=function(e) NULL,
                    warning=function(w) NULL)
  if(is.null(peaks)){
    tab <-  tibble(min=NA, max=NA)
    return(tab)
  }
  tab <- tibble(min=min(peaks$minx),
                max=max(peaks$maxx))
  tab
}

homozygous_cutoff2 <- function(dens){
  handle_missing(homozygous_peaks, dens)
}

hemizygous_cutoff2 <- function(dens){
  handle_missing(hemizygous_peaks, dens)
}

duplication_cutoff2 <- function(dens){
  handle_missing(duplication_peaks, dens)
}

impute_hemizygous <- function(mb, modes, densities, N){
  ##hemizygous_cutoff2(densities[[1]])
  cutoffs <- round(sapply(densities, hemizygous_cutoff2),
                   3)
  cuts <- tibble(batch=seq_along(densities),
                 hemizygous_cutoff=cutoffs)
  dat2 <- assays(mb) %>%
    left_join(cuts, by="batch")
  tab <- dat2 %>%
    group_by(batch) %>%
    summarize(n=sum(oned < hemizygous_cutoff))
  ## if TRUE, no need to impute
  if(all(tab$n > N)) return(assays(mb))
  batches <- batchLevels(mb)
  imputed_list <- vector("list", sum(tab$n <= N))
  for(i in seq_along(batches)){
    freq <- tab$n[tab$batch==i]
    if(freq > N) next()
    dens <- densities[[i]]
    if(!is.null(dens)){
      x <- sample_hemdel_from_density(dens, N=N)
    } else{
      browser()
      ##imputed <- sample_hemdel_from_distribution(dat, THR)
    }
    imputed <- tibble(x=x, batch=i)
    imputed_list[[i]] <- imputed
  }
  imputed <- do.call(rbind, imputed_list)
  index <- sum(isSimulated(mb)) + 1
  impdat <- tibble(id=paste0("augment_", index),
                   oned=imputed$x,
                   provisional_batch=NA,
                   likely_deletion=TRUE,
                   batch=imputed$batch,
                   batch_labels=NA,
                   is_simulated=TRUE)
  dat <- bind_rows(assays(mb), imputed) %>%
    arrange(batch)
  dat
}

summarize_blocks <- function(blocks){
  blocks2 <- blocks %>%
    group_by(type) %>%
    summarize(min_mean=mean(min, na.rm=TRUE),
              max_mean=mean(max, na.rm=TRUE))
  blocks3 <- left_join(blocks, blocks2, by="type") %>%
    mutate(min=ifelse(is.na(min) | !is.finite(min), min_mean, min),
           max=ifelse(is.na(max) | !is.finite(max), max_mean, max)) %>%
    select(-c(min_mean, max_mean)) %>%
    filter(!is.na(min))
  blocks4 <- blocks3 %>%
    mutate(mn=(min+max)/2)
  ##blocks4 <- blocks3 %>%
  blocks3
}

make_blocks <- function(densities){
  homdel_cutoffs <- lapply(densities, homozygous_cutoff2) %>%
    do.call(rbind, .) %>%
    mutate(batch=paste("Batch", seq_along(densities)),
           type="homozygous deletion")
  hemdel_cutoffs <- lapply(densities, hemizygous_cutoff2) %>%
    do.call(rbind, .) %>%
    mutate(batch=paste("Batch", seq_along(densities)),
           type="hemizygous deletion")
  dup_cutoffs <- lapply(densities, duplication_cutoff2) %>%
    do.call(rbind, .) %>%
    mutate(batch=paste("Batch", seq_along(densities))) %>%
    mutate(type="duplication")
  cuts <- bind_rows(homdel_cutoffs, hemdel_cutoffs) %>%
    bind_rows(dup_cutoffs) %>%
    mutate(ymin=0,
           ymax=Inf)
  blocks <- summarize_blocks(cuts)
  blocks
}

equivalent_variance <- function(model){
  s <- colMeans(sigma(chains(model)))
  th <- colMeans(theta(chains(model)))
  L <- length(s)
  if(L >= 3 && th[1] < -1){
    s_other <- mean(s[-c(1, L)])
    fc <- s[L] / s_other
    if(fc >= 1.5) return(FALSE) else return(TRUE)
  }
  NA
}

evaluate_sb3 <- function(mb, mp, ...){
    sb3 <- warmup(assays(mb), "SBP3", "SB3", ...)
    mcmcParams(sb3) <- mp
    sb3 <- posteriorSimulation(sb3)
    sb3
}

#' @export
homdel_model <- function(mb, mp, augment=TRUE, ...){
    if(augment){
        assays(mb) <- augment_homozygous(mb)
    } 
    sb3 <- evaluate_sb3(mb, mp, ...)
    if(stop_early(sb3) || numBatch(mb) == 1) return(sb3)
    final <- explore_multibatch(sb3, ...)
    final
}

hemdel_model <- function(mb.subsamp, mp, ...){
    sb <- warmup(assays(mb.subsamp), "SBP2", "SB2", ...)
    mcmcParams(sb) <- mp
    sb <- posteriorSimulation(sb)
    finished <- stop_early(sb)
    if(finished) return(sb)
    mb <- warmup(assays(mb.subsamp), "MBP2", "MB1")
    mcmcParams(mb) <- mp
    mb <- posteriorSimulation(mb)
    mb
}

hd4comp <- function(mod_2.4, simdat2, mb.subsamp, mp){
    model <- incrementK(mod_2.4) %>%
        gsub("P", "", .)
    ##THR <- summaries(mb.subsamp)$deletion_cutoff
    THR <- deletion_midpoint(mb.subsamp)
    mod_1.4 <- MultiBatchList(data=simdat2)[[ model ]]
    hdmean <- homozygousdel_mean(mb.subsamp, THR)
    hdvar <- homozygousdel_var(mb.subsamp, THR)
    theta(mod_1.4) <- cbind(hdmean, theta(mod_2.4))
    V <- matrix(sigma2(mod_2.4)[, 1],
                nrow(sigma2(mod_2.4)), 3,
                byrow=FALSE)
    sigma2(mod_1.4) <- cbind(hdvar, V)
    mcmcParams(mod_1.4) <- mp
    ##mod_1.4 <- mcmc_homozygous(mod_1.4)
    mod_1.4 <- posteriorSimulation(mod_1.4)
    mod_1.4
}

hd3comp <- function(restricted, simdat, mb.subsamp, mp){
  model <- incrementK(restricted) %>%
    gsub("P", "", .)
  ##THR <- summaries(mb.subsamp)$deletion_cutoff
  THR <- deletion_midpoint(mb.subsamp)
  mod_1.3 <- MultiBatchList(data=simdat)[[ model ]]
  hdmean <- homozygousdel_mean(mb.subsamp, THR)
  hdvar <- homozygousdel_var(mb.subsamp, THR)
  theta(mod_1.3) <- cbind(hdmean, theta(restricted))
  V <- matrix(sigma2(restricted)[, 1],
              nrow(sigma2(restricted)), 3,
              byrow=FALSE)
  sigma2(mod_1.3) <- cbind(hdvar, V)
  mcmcParams(mod_1.3) <- mp
  mod_1.3 <- mcmc_homozygous(mod_1.3)
  mod_1.3
}

homdeldup_model <- function(mb, mp, augment=TRUE, ...){
    if(augment){
        assays(mb) <- augment_homozygous(mb)
    }
    sb <- warmup(assays(mb), "SBP4", "SB4", ...)
    mcmcParams(sb) <- mp
    sb <- posteriorSimulation(sb)
    if(substr(modelName(sb), 3, 3) == "P"){
        ##
        ## if 4th component appears diploid in pooled variance model, don't proceed
        ##
        appears_diploid <- not_duplication(sb)
        if(appears_diploid) return(sb)
    }
    finished <- stop_early(sb, 0.98, 0.98)
    if(finished) return(sb)
    ##
    ## 4th component variance is much too big
    ##
    mb.subsamp <- mb
    fdat <- filter(assays(mb.subsamp), !likely_deletion)
    mb <- warmup(fdat, "MBP3", ...)
    mcmcParams(mb) <- mp
    message("Fitting restricted model")
    mod_2.4 <- restricted_homhemdup(mb, mb.subsamp, mp, ...)
    message("Data augmentation for homozygous deletions")
    simdat2 <- augment_rarehomdel(mod_2.4, sb, mb.subsamp)
    mod_1.4 <- hd4comp(mod_2.4, simdat2, mb.subsamp, mp)
    mod_1.4
}

setMethod("bic", "MultiBatchP", function(object){
  object2 <- object
  object2 <- dropSimulated(object2)
  ## number of free parameters to estimate (counting just top level of model)
  ## tau^2_k, mu_k, nu.0, sigma2.0, pi
  K <- length(tau2(object2)) +
    length(mu(object2)) +
    length(nu.0(object2)) +
    length(sigma2.0(object2)) +
    length(p(object2)[1, ])
  n <- length(oned(object2))
  ll <- .compute_loglik(object2)
  bicstat <- -2*(ll + logPrior(object2)) + K*(log(n) - log(2*pi))
  bicstat
})

setMethod("bic", "MultiBatch", function(object){
  object <- dropSimulated(object)
  K <- length(tau2(object)) +
    length(mu(object)) +
    length(nu.0(object)) +
    length(sigma2.0(object)) +
    length(p(object)[1, ])
  n <- length(oned(object))
  ll <- .compute_loglik(object)
  bicstat <- -2*(ll + logPrior(object)) + K*(log(n) - log(2*pi))
  bicstat
})

setMethod(".compute_loglik", "MultiBatch", function(object){
  model <- as(object, "MultiBatchModel")
  .compute_loglik(model)
})

setMethod(".compute_loglik", "MultiBatchP", function(object){
  model <- as(object, "MultiBatchPooled")
  .compute_loglik(model)
})

validModel <- function(model){
  model <- model[ !isSimulated(model) ]
  !any(is.na(theta(model))) &&
    !any(is.na(sigma(model))) &&
     (ncol(theta(model)) == k(model))
}

model_checks <- function(models){
  var.check <- sapply(models, equivalent_variance)
  dip.check <- sapply(models, same_diploid_component)
  distinct <- sapply(models, distinct_components)
  model_names <- sapply(models, modelName)
  tib <- tibble(model=model_names,
                equiv_variance=var.check,
                same_diploid_comp=dip.check,
                distinct_comp=distinct)
}

deletion_models <- function(mb, snp_se, mp, THR){
##  if(missing(THR)){
##    if("deletion_cutoff" %in% names(summaries(mb))){
##      THR <- summaries(mb)$deletion_cutoff
##    } else stop("THR not provided")
    ##  } else summaries(mb)$deletion_cutoff <- THR
  THR <- deletion_midpoint(mb)
  if(!any(oned(mb) < THR)) stop("No observations below deletion cutoff")
  model3 <- homdel_model(mb, mp)
  model4 <- homdeldup_model(mb, mp)
  if(!is.null(snp_se)){
    model3 <- genotype_model(model3, snp_se)
    model4 <- genotype_model(model4, snp_se)
  }
  ##compare bic without data augmentation
  model.list <- list(model3, model4)
  model.list
}

hemideletion_models <- function(mb.subsamp, snp_se, mp, 
                                augment=TRUE, ...){
    ##assays(mb.subsamp)$deletion_cutoff <- THR
    mb1 <- hemdel_model(mb.subsamp, mp, ...)
    mb2 <- hemdeldup_model2(mb.subsamp, mp, ...)
    if(is.null(snp_se)){
        model.list <- list(mb1, mb2)
        model.list
        return(model.list)
    }
    if(nrow(snp_se) > 0){
        mb1 <- genotype_model(mb1, snp_se)
        if(!is.null(mb2))
            mb2 <- genotype_model(mb2, snp_se)
    }
    model.list <- list(mb1, mb2)
    model.list
}

posthoc_checks <- function(model.list){
  checks <- model_checks(model.list)
  bics <- sapply(model.list, bic)
  checks$bic <- bics
  checks
}

hemdeldup_model <- function(mb.subsamp, mp, THR=-0.25){
  THR <- summaries(mb.subsamp)$deletion_cutoff <- THR
  simdat <- augment_homozygous(mb.subsamp)
  sb <- warmup(assays(mb.subsamp),
               "SBP3",
               "SB3")
  mcmcParams(sb) <- mp
  sb <- posteriorSimulation(sb)
  finished <- stop_early(sb, 0.99, 0.99)
  if(is.na(finished)) finished <- FALSE
  if(finished) return(sb)
  ##THR <- summaries(mb.subsamp)$deletion_cutoff
  mb <- explore_multibatch(sb, simdat)
  return(mb)
}

hemdeldup_model2 <- function(mb.subsamp, mp, ...){
    sb <- warmup(assays(mb.subsamp),
                 "SBP3",
                 "SB3", ...)
    mcmcParams(sb) <- mp
    sb <- tryCatch(posteriorSimulation(sb),
                   warning=function(w) NULL)
    if(is.null(sb)) return(NULL)
    finished <- stop_early(sb, 0.98, 0.98)
    if(is.na(finished)) finished <- FALSE
    if(finished) return(sb)
    ##    while(sum(oned(mb.subsamp) < THR) < 5){
    ##        THR <- THR + 0.05
    ##    }
    sb.meds <- colMedians(theta(chains(sb)))
    ##thr <- deletion_midpoint(mb.subsamp)
    thr <- theta(sb3)[, 3] - 2*sigma(sb3)[1, 1]    
    densities <- compute_density(mb.subsamp, thr)
    diploid_modes <- compute_modes(densities)
    dist <- c(sb.meds[2] - sb.meds[1],
              sb.meds[3] - sb.meds[2])
    hemdel <- diploid_modes - dist[1]
    dup <- diploid_modes + dist[2]
    ## standard deviations will be inflated in SB model
    if(modelName(sb)=="SBP3"){
        s <- rep(sigma(sb)[1,1]/2, 2)
    } else s <- sigma(sb)[1, c(1, 3)]/2
    B <- numBatch(mb.subsamp)
    model_names <- rep(c("hemidel", "dup"), each=B)
    means <- c(hemdel, dup)
    sds <- rep(s, each=B)
    tab <- tibble(component=model_names,
                  mean=means,
                  sd=sds)
    x <- vector("list", nrow(tab))
    for(i in seq_len(nrow(tab))){
        x[[i]] <- rnorm(10, tab$mean[i], tab$sd[i])
    }
    x <- unlist(x)
    likely_deletion <- c(rep(TRUE, B*10), rep(FALSE, B*10))
    sdat <- tibble(id=paste0("augment_", seq_along(x)),
                   oned=x,
                   provisional_batch=NA,
                   likely_deletion=likely_deletion,
                   is_simulated=TRUE,
                   batch=rep(rep(seq_len(B), each=10), 2),
                   homozygousdel_mean=NA)
                   ##likely_hd=NA)
    tmp <- bind_rows(assays(mb.subsamp), sdat) %>%
      arrange(batch)
    mb <- MultiBatchList(data=tmp)[["MBP3"]]
    mb <- warmup(tmp, "MBP3", "MB3", ...)
    mcmcParams(mb) <- mp
    mb <- posteriorSimulation(mb)
    return(mb)
}

restricted_homhemdup <- function(mb, mb.subsamp, mp, ...){
    mod_2.4 <- suppressWarnings(posteriorSimulation(mb))
    is_flagged <- mod_2.4@flags$.internal.counter > 40
    if(!is_flagged) return(mod_2.4)
    sb3 <- warmup(assays(mb), "SBP3", ...)
    mcmcParams(sb3) <- mp
    sb3 <- posteriorSimulation(sb3)
    ## simulate from the pooled model for each batch
    full.dat <- assays(mb.subsamp)
    ##    ggMixture(mod_2.4) +
    ##        geom_vline(xintercept=theta(sb3)[, 3] - 2*sigma(sb3)[1, 1])
    thr <- theta(sb3)[, 3] - 2*sigma(sb3)[1, 1]
    simdat <- augment_rareduplication(sb3,
                                      mod_2.4,
                                      full_data=full.dat,
                                      THR=thr)
    mod_2.4.2 <- MultiBatchList(data=simdat)[["MBP3"]]
    simdat2 <- augment_rarehemdel(sb3,
                                  mod_2.4.2,
                                  full_data=full.dat)
    filtered.dat <- filter(simdat2, !likely_deletion)
    mb <- warmup(filtered.dat, "MBP3", ...)
    mcmcParams(mb) <- mp
    mod_2.4 <- posteriorSimulation(mb)
    mod_2.4
}

dropSimulated <- function(model) model[!isSimulated(model)]

distinct_components <- function(model){
  if(numBatch(model) == 1) return(TRUE)
  p <- probz(model)
  nearly_equivalent <- rowSums(p > 0.4 & p < 0.6) > 0
  if(sum(nearly_equivalent) == 0) return(TRUE)
  colnames(p) <- paste0("comp", seq_len(ncol(p)))
  b <- batch(model)
  tib <- as_tibble(p) %>%
    mutate(nearly_equivalent=nearly_equivalent,
           batch=b) %>%
    filter(nearly_equivalent)
  batches <- seq_len(numBatch(model))
  nbatch <- tibble(batch=batches, N=as.numeric(table(b)))
  equiv <- tib %>%
    left_join(nbatch, by="batch") %>%
    group_by(batch) %>%
    summarize(number_equivalent=n(),
              N=unique(N)) %>%
    mutate(prop_equivalent=number_equivalent/N) %>%
    filter(prop_equivalent > 1/3)
  ifelse(nrow(equiv) > 0, FALSE, TRUE)
}

same_diploid_component<- function(model){
  if(numBatch(model) == 1 | k(model) == 1) return(TRUE)
  pmix <- p(model)
  ranks <- apply(pmix, 1, order, decreasing=TRUE)
  diploid.ranks <- ranks[1, ]
  ifelse(length(unique(diploid.ranks)) > 1, FALSE, TRUE)
}

not_duplication <- function(model){
  pmix <- p(model)
  if(nrow(pmix) > 1){
    sds <- colSds(pmix)
    if(any(sds > 0.1)){
      ## suggests overfitting
      return(TRUE)
    }
  }
  k_ <- ncol(pmix)
  ## most samples belong to 4th component
  ranks <- apply(pmix, 1, order, decreasing=TRUE)
  diploid.ranks <- ranks[1, ]
  appears_diploid <- any(diploid.ranks == k_ & pmix[, k_] > 0.6)
  notdup <- appears_diploid
  notdup
}

##singlebatch_no_duplication <- function(model){
##  mname <- modelName(model)
##  is_pooled <- substr(mname, 3, 3) == "P"
##  if(is_pooled){
##    appears_diploid <- not_duplication(model)
##  }else{
##    ## variances not pooled
##    ## check if 3rd or 4th component has large variance
##    ## if so, we need to explore multibatch models
##    pred <- predictive(chains(model))
##    maxsd <- max(sd(pred[, 3]), sd(pred[, 4]))
##    
##  }
##  notdup
##}

batchLabels <- function(object){
  assays(object)$batch_labels
}

anyMissingBatchLabels <- function(object) any(is.na(batchLabels(object)))

duplication_models <- function(mb.subsamp, snpdat, mp, THR=-0.25){
  ## duplication model
  sb <- warmup(assays(mb.subsamp),
               "SBP2",
               "SB2")
  mcmcParams(sb) <- mp
  sb <- tryCatch(posteriorSimulation(sb),
                 warning=function(w) NULL)
  if(!is.null(sb)){
    appears_diploid <- not_duplication(sb)
    if(!is.null(snpdat)){
      sb <- genotype_model(sb, snpdat)
    }
    ## probability > 0.98 for 99% or more of participants
    finished <- stop_early(sb, 0.98, 0.99)
    if(finished){
      return(list(sb))
    }
  }
  ##
  ## Try MultiBatch
  ##
  mb <- warmup(assays(mb.subsamp), "MBP2")
  mcmcParams(mb) <- mp
  mb <- tryCatch(posteriorSimulation(mb), warning=function(w) NULL)
  if(!is.null(mb) && !is.null(snpdat)){
    mb <- genotype_model(mb, snpdat)
  }
  list(sb, mb)
}

select_models <- function(mb){
  minlogr <- min(oned(mb), na.rm=TRUE)
  if(minlogr < -1){
    model <- deletion_models
    return(model)
  }
  number <- sum(oned(mb) >= -1  & oned(mb) < -0.25)
  if(number > 1 && minlogr >= -1  && minlogr < -0.25){
    model <- hemideletion_models
    return(model)
  }
  duplication_models
}

#' @export
use_cutoff <- function(mb){
  minlogr <- min(oned(mb), na.rm=TRUE)
  if(minlogr < -1){
    cutoff <- -1
  }
  if(minlogr >= -1  & minlogr < -0.25){
    cutoff <- -0.25
  }
  if(minlogr >= -0.25){
    cutoff <- 0
  }
  cutoff
}


choose_model <- function(model.list, mb){
  is_null <- sapply(model.list, is.null)
  model.list <- model.list[ !is_null ]
  if(length(model.list) == 1) return(model.list[[1]])
  if(length(model.list) == 2) {
    posthoc <- posthoc_checks(model.list)
    appears_diploid <- not_duplication(model.list[[2]])
    if(appears_diploid){
      model <- model.list[[1]]
    } else {
      ix <- which.min(posthoc$bic)
      model <- model.list[[ix]]
    }
    return(model)
  }
  if(length(model.list) == 0){
    ## return single component model
    mb <- MultiBatchList(data=assays(mb))[["MB1"]]
    mcmcParams(mb) <- mp
    model <- posteriorSimulation(mb)
  }
  model
}

preliminary_checks <- function(mb, grange){
  if(any(is.na(assays(mb)$batch_labels))) {
    message("Missing data in `assays(mb)$batch_labels`")
    return(FALSE)
  }
  if(is.null(grange)) return(TRUE)
  if(length(grange) > 1){
    warning("Multiple elements for `grange` object. Only using first")
    return(TRUE)
  }
  TRUE
}

#' @export
cnv_models <- function(mb,
                       grange,
                       snp_se,
                       mp=McmcParams(iter=400, burnin=500),
                       THR){
  ok <- preliminary_checks(mb, grange)
  stopifnot(ok)
  grange <- grange[1]
  if(!is.null(snp_se))
    snpdat <- subsetByOverlaps(snp_se, grange)
  modelfun <- select_models(mb)
  if(missing(THR))
    THR <- use_cutoff(mb)
  model.list <- modelfun(mb, snpdat, mp, THR)
  model <- choose_model(model.list, mb)
  model
}

#' @export
upsample <- function(model, se, provisional_batch){
    dat <- getData2(se[1, ], provisional_batch, model)
    pred <- predictiveDist(model)
    dat2 <- predictiveProb(pred, dat) %>%
        mutate(copynumber=mapping(model)[inferred_component])
    if(nrow(dat2) == 0) return(dat2)
    ix <- which(colnames(dat2) %in% as.character(0:4))
    ##
    ## multiple components can map to the same copy number state
    ## -- add probabilities belonging to same component
    select <- dplyr::select
    tmp <- dat2[, ix] %>%
        mutate(id=dat2$id) %>%
        gather("state", "prob", -c(id)) %>%
        mutate(component_index=as.numeric(state) + 1) %>%
        mutate(copynumber=mapping(model)[component_index]) %>%
        group_by(id, copynumber) %>%
        summarize(prob=sum(prob)) %>%
        select(c(id, prob, copynumber)) %>%
        spread(copynumber, prob)
    dat2 <- dat2[, -ix] %>%
      left_join(tmp, by="id")
    dat3 <- tidy_cntable(dat2)
    dat3
}

tidy_cntable <- function(dat){
    tmp <- tibble(id=dat$id,
                  `0`=rep(0, nrow(dat)),
                  `1`=rep(0, nrow(dat)),
                  `2`=rep(0, nrow(dat)),
                  `3`=rep(0, nrow(dat)),
                  `4`=rep(0, nrow(dat)))
    dat2 <- left_join(dat, tmp, by="id") 
    dropcols <- paste0(0:4, ".y")
    dropcols <- dropcols[ dropcols %in% colnames(dat2)]
    dat2 <- select(dat2, -dropcols)
    renamecols <- paste0(0:4, ".x")
    renamecols <- renamecols[renamecols %in% colnames(dat2)]
    renameto <- gsub("\\.x", "", renamecols)
    colnames(dat2)[match(renamecols, colnames(dat2))] <- renameto
    colnames(dat2)[match(0:4, colnames(dat2))] <- paste0("cn_", 0:4)
    keep <- c("id", "batch", "copynumber", paste0("cn_", 0:4))
    dat3 <- select(dat2, keep) %>%
        mutate(copynumber=as.integer(copynumber))
    dat3
}
