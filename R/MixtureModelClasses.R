#' @include AllGenerics.R
NULL

#' An object for running MCMC simulations.
#'
#' BatchModel and MarginalModel both inherit from this class.
#' @slot data a tibble with one-dimensional summaries, id, batches, and is_sampled indicator
#' @slot parameters list of parameters
#' @slot chains object of class McmcChains
#' @slot current_values current value of each chain
#' @slot summaries list of empirical data and model summaries
#' @slot flags list of model flags
#' @export
setClass("MultiBatch", representation(data="tbl_df",
                                      parameters="list",
                                      chains="McmcChains",
                                      current_values="list",
                                      summaries="list",
                                      flags="list"))


model_spec <- function(model, data) {
  models <- c("SB", "SBP", "MB", "MBP")
  K <- 1:5
  avail.models <- lapply(models, paste0, K) %>%
    unlist
  if(missing(model)) model <- "MB3" else model <- model[1]
  if(!model %in% avail.models) stop("Model not recognized")
  number_batches <- max(1, length(unique(data$batch)))
  k <- substr(model, nchar(models), nchar(models)) %>%
    as.integer
  is_SB <- substr(model, 1, 2) == "SB"
  number_batches <- ifelse(is_SB, 1, number_batches)
  number_obs <- sum(data$is_sampled)
  tab <- tibble(model=model,
                k=k,
                number_batches=number_batches,
                number_obs=number_obs)
  tab
}

listChains1 <- function(model_specs, parameters){
  S <- iter(parameters[["mp"]])
  num.chains <- nStarts(parameters[["mp"]])
  B <- model_specs$number_batches
  K <- model_specs$k
  mc.list <- vector("list", num.chains)
  for(i in seq_len(num.chains)){
    mc.list[[i]] <- initialize_mcmc(K[i], S, B[i])
  }
  mc.list
}

##
## Constructor
##
MultiBatch <- function(model="MB3",
                       data=modelData(),
                       specs=modelSpecs(model, data),
                       parameters=modelParameters(),
                       chains=listChains1(specs, parameters),
                       current_values=modelValues(),
                       summaries=modelSummaries(),
                       flags=modelFlags()){
  new("MultiBatch",
      data=data,
      parameters=parameters,
      chains=chains,
      current_values=current_values,
      summaries=summaries,
      flags=flags)
}


setMethod("show", "MultiBatch", function(object){
  ##callNextMethod()
  cls <- class(object)
  cat(paste0("An object of class ", cls), "\n")
  cat("     n. obs              :", nrow(assays(object)), "\n")
  cat("     n. batches          :", nBatch(object), "\n")
  cat("     k                   :", k(object), "\n")
  cat("     nobs/batch          :", table(batch(object)), "\n")
  cat("     log lik (s)         :", round(log_lik(object), 1), "\n")
  cat("     log prior (s)       :", round(logPrior(object), 1), "\n")
  cat("     log marginal lik (s):", round(marginal_lik(object), 1), "\n")
})

setClass("MultiBatchPooled2", contains="MultiBatch")

setClass("GenotypeModel", contains="MultiBatch",
         representation(mapping="character"))

setClass("GenotypePooled", contains="MultiBatch",
         representation(mapping="character"))

modelParameters <- function(hp=Hyperparameters(),
                            mp=McmcParams()){
  list(hp=hp, mp=mp)
}

extractParameters <- function(old){
  list(hp=hyperParams(old),
       mp=mcmcParams(old))
}

modelValues <- function(theta=matrix(),
                        sigma2=matrix(),
                        nu.0=numeric(),
                        sigma2.0=numeric(),
                        p=numeric(),
                        mu=numeric(),
                        tau2=numeric(),
                        u=numeric(),
                        z=integer(),
                        logprior=numeric(),
                        loglik=numeric(),
                        probz=matrix()){
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

modelSummaries <- function(data.mean=matrix(),
                           data.prec=matrix(),
                           zfreq=integer(),
                           marginal_lik=numeric(),
                           modes=list()){
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

modelFlags <- function(.internal.constraint=numeric(),
                       .internal.counter=numeric(),
                       label_switch=logical()){
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
                      batch=integer(),
                      is_sampled=logical()){
  tibble(id=id,
         oned=oned,
         batch=batch,
         is_sampled=is_sampled)
}

extractData <- function(old){
  tibble(id=seq_along(y(old)),
         oned=y(old),
         batch=batch(old),
         is_sampled=TRUE)
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
  MultiBatch(data=data,
             parameters=params,
             chains=chains(from),
             current_values=values,
             summaries=summaries,
             flags=flags)
})


##
## Accessors
##
setMethod("values", "MultiBatch", function(x, ...){
  x@current_values
})

setMethod("theta", "MultiBatch", function(object){
  th <- values(object)[["theta"]]
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

setReplaceMethod("chains", c("MultiBatch", "McmcChains"), function(object, value){
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
                   nStarts(mcmcParams(object)) <- as.integer(value)
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

setReplaceMethod("batch", c("MultiBatch", "numeric"), function(object, value){
  assays(object)[["batch"]] <- as.integer(value)
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

##MB2 <- function(data=tibble(),
##                hp=HyperparametersMultiBatch(),
##                mp=McmcParams(iter=1000, thin=10,
##                              burnin=1000, nStarts=4)){
##  if(nrow(dat) == 0){
##    return(.empty_mb(hp, mp))
##  }
##  iter <- 0
##  validZ <- FALSE
##  mp.tmp <- McmcParams(iter=0, burnin=burnin(mp), thin=1, nStarts=1)
##  while(!validZ){
##    ##
##    ## Burnin with MB model
##    ##
##    mb <- .MB(dat, hp, mp.tmp, batches)
##    mb <- runBurnin(mb)
##    tabz1 <- table(batch(mb), z(mb))
##    tabz2 <- table(z(mb))
##    validZ <- length(tabz2) == k(hp) && all(tabz1 > 1)
##    iter <- iter + 1
##    if(iter == 50) {
##      message("Trouble initializing a valid model. The number of components is likely too large")
##      return(NULL)
##    }
##  }
##  mb2 <- sortComponentLabels(mb)
##  mcmcParams(mb2) <- mp
##  chains(mb2) <- McmcChains(mb2)
##  ## coerce back to MultiBatchModel2
##  mb2
##}
