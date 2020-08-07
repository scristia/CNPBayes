#' @include AllGenerics.R
NULL

setValidity("MixtureModel", function(object){
  msg <- TRUE
  if(ncol(p(object)) != k(object)){
    msg <- "Mixture probability vector must be the same length as k"
    return(msg)
  }
  if(k(object)!=k(hyperParams(object))){
    msg <- "disagreement of k in hyperparams and model"
    return(msg)
  }
  if(length(y(object))!=length(u(object))){
    msg <- "u-vector must be same length as data"
    return(msg)
  }
  if(iter(object) != iter(chains(object))){
    msg <- "number of iterations not the same between chains and model"
    return(msg)
  }
  th.len <- prod(dim(theta(object)))
  pr.len <- length(object@predictive)
  if(th.len != pr.len){
    msg <- "predictive slot in current values should have length K x B"
    return(msg)
  }
  zz <- z(object)
  if(length(zz) == 0) return(msg)
  if(any(is.na(zz))){
    msg <- "NAs in z slot"
    return(msg)
  }
  ## z is an index -- it must be greater than 0
  if(min(zz) < 1){
    msg <- "The minimum z-value must be 1"
    return(msg)
  }
  msg
})

setMethod("hyperParams", "MixtureModel", function(object) object@hyperparams)

setReplaceMethod("hyperParams", c("MixtureModel", "Hyperparameters"),
  function(object, value) {
    object@hyperparams <- value
    object
})

setReplaceMethod("batch", "MixtureModel", function(object, value){
  object@batch <- value
  object
})

observed <- function(object) object@data

sigma <- function(object) sqrt(sigma2(object))

tau <- function(object) sqrt(tau2(object))


setMethod("nu.0", "MixtureModel", function(object) object@nu.0)


setMethod("sigma2.0", "MixtureModel", function(object) object@sigma2.0)

sigma.0 <- function(object) sqrt(sigma2.0(object))


setMethod("p", "MixtureModel", function(object){
  object@pi
})

nComp <- function(object) length(p(object))

setMethod("dataMean", "MixtureModel", function(object) object@data.mean)
setMethod("dataPrec", "MixtureModel", function(object) object@data.prec)
setMethod("dataSd", "MixtureModel", function(object) sqrt(1/dataPrec(object)))


setMethod("k", "MixtureModel", function(object) object@k)


setReplaceMethod("k", "MixtureModel",
                 function(object, value) {
                   browser()
        k <- as.integer(value)
        hypp <- hyperParams(object)
        hypp@k <- k
        hypp@alpha <- rep(1, k)
        hyperParams(object) <- hypp
        object@k <- k
        object@pi <- matrix(1, nrow(theta(object)), k)
        object@probz <- matrix(0, length(y(object)), k)
        object <- startingValues(object)
        object
    }
)

setReplaceMethod("z", "MixtureModel", function(object, value){
  ##object@z <- factor(value, levels=seq_len(k(object)))
  object@z <- value
  object
})

setReplaceMethod("theta", "MixtureModel", function(object, value){
  object@theta <- value
  object
})

setReplaceMethod("sigma2", "MixtureModel", function(object, value){
  object@sigma2 <- value
  object
})

setReplaceMethod("p", "MixtureModel", function(object, value){
  object@pi <- value
  object
})

setReplaceMethod("nu.0", "MixtureModel", function(object, value){
  object@nu.0 <- value
  object
})

setReplaceMethod("sigma2.0", "MixtureModel", function(object, value){
  object@sigma2.0 <- value
  object
})

setMethod("chains", "MixtureModel", function(object) object@mcmc.chains)

setReplaceMethod("chains", "MixtureModel", function(object, value){
  object@mcmc.chains <- value
  object
})

setMethod("show", "MixtureModel", function(object){
  cat("An object of class '", class(object), "'\n")
  cat("  data: \n")
  cat("     n           :", length(y(object)), "\n")
  cat("     k           :", nComp(object), "\n")
  cat("     table(z)    :", paste(tablez(object), collapse=", "), "\n")
  cat("     mix prob (s):", paste(round(p(object), 2), collapse=", "), "\n")
  sds <- showSigmas(object)
  mns <- showMeans(object)
  cat("     sigma (s)   :", sds, "\n")
  cat("     theta (s)   :", mns, "\n")
  cat("     sigma2.0 (s):", round(sigma2.0(object), 2), "\n")
  cat("     nu.0 (s)    :", nu.0(object), "\n")
  cat("     logprior(s):", round(logPrior(object), 2), "\n")
  cat("     loglik (s)  :", round(log_lik(object), 2), "\n")
  cat("     log marginal lik (s)  :", round(log_lik(object), 2), "\n")
})

setMethod("show", "SingleBatchModel", function(object){
  object <- as(object, "MultiBatchModel")
  show(object)
})

setMethod("alpha", "MixtureModel", function(object) alpha(hyperParams(object)))

setMethod("y", "MixtureModel", function(object) object@data)

setReplaceMethod("y", "MixtureModel", function(object, value){
  object@data <- value
  object
})

setMethod("oned", "MixtureModel", function(object) object@data)


setMethod("batch", "MixtureModel", function(object) object@batch)

setMethod("z", "MixtureModel", function(object) object@z)


setMethod("computePrior", "SingleBatchModel", function(object){
  compute_logprior(object)
})


.computeProbZ <- function(object){
  pZ <- probz(object)
  zz <- as.integer(z(object))
  for(j in seq_len(k(object))){
    pZ[, j] <- pZ[, j] + as.integer(zz==j)
  }
  pZ
}


setMethod("numberObs", "MixtureModel", function(model) length(y(model)))


setReplaceMethod("probz", "MixtureModel", function(object, value){
  object@probz <- value
  object
})

setMethod("probz", "MixtureModel", function(object) {
  if(iter(object) > 0){
    pz <- object@probz/iter(object)
    return(pz)
  }
  object@probz
})


selectByLogLik <- function(model.list){
  lp <- sapply(model.list, log_lik)
  isfin <- is.finite(lp)
  isnonzero <- lp != 0
  is.candidate <- isfin & isnonzero
  if(!any(is.candidate)){
    stop("Bad starting values.")
  }
  lp <- lp[is.candidate]
  model.list <- model.list[is.candidate]
  select <- which.max(lp)
  if(length(select) == 0 ) stop("No model selected")
  model <- model.list[[select]]
  model
}

psParams <- function(warnings=TRUE,
                     returnNULLonWarnings=FALSE){
  list(warnings=warnings)
}

setMethod("isOrdered", "MixtureModel", function(object){
  identical(order(theta(object)), seq_along(theta(object)))
})

setMethod("isOrdered", "MultiBatchModel", function(object){
  .ordered_thetas_multibatch(object)
})

setReplaceMethod("dataMean", "MixtureModel", function(object, value){
  object@data.mean <- value
  object
})

setReplaceMethod("dataPrec", "MixtureModel", function(object, value){
  object@data.prec <- value
  object
})

setMethod("mu.0", "MixtureModel", function(object) mu.0(hyperParams(object)))
setMethod("mu.0", "Hyperparameters", function(object) object@mu.0)

setMethod("tau2.0", "MixtureModel", function(object) tau2.0(hyperParams(object)))
setMethod("tau2.0", "Hyperparameters", function(object) object@tau2.0)

tau.0 <- function(object) sqrt(tau2.0(object))


setMethod("eta.0", "MixtureModel", function(object) eta.0(hyperParams(object)))

setMethod("eta.0", "Hyperparameters", function(object) object@eta.0)

setMethod("m2.0", "MixtureModel", function(object) m2.0(hyperParams(object)))


setMethod("m2.0", "Hyperparameters", function(object) object@m2.0)

setReplaceMethod("eta.0", "MixtureModel", function(object, value){
  eta.0(hyperParams(object)) <- value
  object
})

setReplaceMethod("m2.0", "MixtureModel", function(object, value){
  m2.0(hyperParams(object)) <- value
  object
})

setReplaceMethod("eta.0", "Hyperparameters", function(object, value){
  object@eta.0 <- value
  object
})

setReplaceMethod("m2.0", "Hyperparameters", function(object, value){
  object@m2.0 <- value
  object
})




##
## the batch names tend to be much too long
##
makeUnique <- function(x, nchar=8){
  ub <- unique(x)
  ##names(ub) <- ub
  maxchar <- pmin(nchar(ub), nchar)
  abbrv <- setNames(make.unique(substr(ub, 1, maxchar)), ub)
  as.character(abbrv[x])
}

#' Calculate the maximum a posteriori estimate of latent variable assignment.
#'
#' @examples
#'      head(map_z(SingleBatchModelExample))
#' @param object an object of class MixtureModel.
#' @return map estimate of latent variable assignment for each observation
#' @export
map_z <- function(object) {
  max.col(probz(object))
}

setMethod("thetac", "MixtureModel", function(object) theta(chains(object)))

setMethod("thetaMean", "MixtureModel", function(object) colMeans(thetac(object)))

setMethod("sigmaMean", "MixtureModel", function(object) colMeans(sigmac(object)))

logLikc <- function(object) log_lik(chains(object))


sigmac <- function(object) sigma(chains(object))

pic <- function(object) p(chains(object))

setMethod("pMean", "MixtureModel", function(object){
  colMeans(pic(object))
})

muc <- function(object) mu(chains(object))

muMean <- function(object) {
  x <- muc(object)
  if(is(object, "MultiBatchModel")){
    return(colMeans(x))
  }
  mean(x)
}

tauc <- function(object) sqrt(tau2(chains(object)))

tauMean <- function(object){
  x <- tauc(object)
  if(is(object, "MultiBatchModel")){
    return(colMeans(x))
  }
  mean(x)
}


setMethod("modes", "MixtureModel", function(object) object@modes)

setReplaceMethod("modes", "MixtureModel", function(object, value) {
  object@modes <- value
  object
})

setMethod("log_lik", "MixtureModel", function(object){
  object@loglik
})

setReplaceMethod("log_lik", "MixtureModel", function(object, value){
  object@loglik <- value
  object
})

argMax <- function(object){
  ll <- log_lik(chains(object))
  if(length(ll) == 1) return(1)
  ##lp <- logPrior(chains(object))
  ##p <- ll+lp
  p <- ll
  p <- p[is.finite(p)]
  if(length(p) == 0){
    return(1)
  }
  maxp <- max(p)
  which(p == maxp)[1]
}

setMethod("isSB", "SingleBatchModel", function(object) TRUE)

setMethod("isSB", "MultiBatchModel", function(object) FALSE)

startAtTrueValues <- function(model, truth){
  theta(model) <- theta(truth)
  sigma2(model) <- sigma2(truth)
  p(model) <- p(truth)
  z(model) <- z(truth)
  mu(model) <- mu(truth)
  tau2(model) <- tau2(truth)
  dataMean(model) <- computeMeans(model)
  dataPrec(model) <- 1/computeVars(model)
  zFreq(model) <- as.integer(table(z(model)))
  log_lik(model) <- computeLoglik(model)
  model
}

restartAtChainIndex <- function(model, index){
  ch <- chains(model)
  nb <- nrow(theta(model))
  theta(model) <- matrix(theta(ch)[index, ], nrow=nb)
  sigma2(model) <- matrix(sigma2(ch)[index, ], nrow=nb)
  p(model) <- matrix(p(ch)[index, ], nrow=nb)
  mu(model) <- mu(ch)[index]
  tau2(model) <- tau2(ch)[index]
  sigma2.0(model) <- sigma2.0(ch)[index]
  nu.0(model) <- nu.0(ch)[index]
  zFreq(model) <- as.integer(table(z(model)))
  dataMean(model) <- computeMeans(model)
  dataPrec(model) <- 1/computeVars(model)
  model
}

setMethod("zFreq", "MixtureModel", function(object) object@zfreq)

setReplaceMethod("zFreq", "MixtureModel", function(object, value){
  object@zfreq <- value
  object
})


setMethod("mcmcParams", "MixtureModel", function(object) object@mcmc.params )


setMethod("iter", "MixtureModel", function(object) iter(mcmcParams(object)))


setMethod("nStarts", "MixtureModel", function(object) nStarts(mcmcParams(object)))


setReplaceMethod("nStarts", "MixtureModel", function(object, value){
  mcmcParams(object)@nstarts <- as.integer(value)
  object
})


setMethod("thin", "MixtureModel", function(object) thin(mcmcParams(object)))


setMethod("burnin", "MixtureModel", function(object) burnin(mcmcParams(object)))

setMethod("numBatch", "MixtureModel", function(object) numBatch(chains(object)))


setReplaceMethod("burnin", "MixtureModel", function(object, value){
  burnin(mcmcParams(object)) <- value
  object
})

setReplaceMethod("iter", "MixtureModel", function(object, value){
  mp <- mcmcParams(object)
  iter(mp) <- as.integer(value)
  mcmcParams(object) <- mp
  iter(chains(object)) <- as.integer(value)
  object
})

setMethod("logPrior", "MixtureModel", function(object) object@logprior)

setReplaceMethod("logPrior", "MixtureModel", function(object, value) {
  object@logprior <- value
  object
})

setMethod("paramUpdates", "MixtureModel", function(x) paramUpdates(mcmcParams(x)))
setReplaceMethod("paramUpdates", "MixtureModel", function(x, value){
  paramUpdates(mcmcParams(x)) <- value
  x
})

nu0c <- function(object) nu.0(chains(object))
sigma20c <- function(object) sigma2.0(chains(object))

setReplaceMethod("mcmcParams", "MixtureModel", function(object, value){
  S <- iter(value)
  B <- numBatch(object)
  K <- k(object)
  if(iter(object) != S){
    if(S > iter(object)){
      object@mcmc.params <- value
      ## create a new chain
      mcmc_chains <- initialize_mcmc(K, S, B)
    } else {
      object@mcmc.params <- value
      index <- seq_len(S)
      mcmc_chains <- chains(object)[index, ]
      mcmc_chains@iter <- S
      mcmc_chains@B <- B
      mcmc_chains@k <- K
    }
    chains(object) <- mcmc_chains
    return(object)
  }
  ## if we've got to this point, it must be safe to update mcmc.params
  ## (i.e., size of chains is not effected)
  object@mcmc.params <- value
  object
})


setReplaceMethod("mcmcParams", "list",
                 function(object, value){
                   for(i in seq_along(object)){
                     mcmcParams(object[[i]]) <- value
                   }
                   object
                 })

setMethod("zChain", "MixtureModel", function(object) chains(object)@z)

mapModel <- function(model){
  model2 <- restartAtChainIndex(model, argMax(model))
  model2
}

mapCnProbability <- function(model){
  ## Cardin et al. : calculate probabilistic copy number assignments
  ## using Bayes Rule applied at the MAP estimates of the cluster
  ## mean, variance, and class proportion parameters
  map_model <- mapModel(model)
  p <- updateMultinomialProb(map_model)
  return(p)
}


setReplaceMethod("thin", c("MixtureModel", "numeric"), function(object, value){
  mp <- mcmcParams(object)
  mp@thin <- as.integer(value)
  mcmcParams(object) <- mp
  object
})


setMethod("mcmcParams", "list", function(object){
  mcmcParams(object[[1]])
})


setMethod("label_switch", "MixtureModel", function(object) object@label_switch)

setReplaceMethod("label_switch", c("MixtureModel", "logical"),
                 function(object, value){
                   object@label_switch <- value
                   object
                 })

setMethod("marginal_lik", "MixtureModel", function(object){
  object@marginal_lik
})

setReplaceMethod("marginal_lik", c("MixtureModel", "numeric"),
                 function(object, value){
                   object@marginal_lik <- value
                   object
                 })


dst <- dlocScale_t


upSample2 <- function(orig.data,
                      model, ## the downsampled model
                      up_sample=TRUE){
  model2 <- useModes(model)
  ## if we do not upSample, we should be able to recover the original probabilities
  ##orig.data$ix <- seq_len(nrow(orig.data))
  ##orig.data <- orig.data[order(orig.data$batch_index), ]
  if(up_sample){
    y(model2) <- orig.data$medians
    if(length(unique_batch(batch(model))) > 1) {
      batch(model2) <- orig.data$batch_index
    } else batch(model2) <- rep(1L, nrow(orig.data))
    model2@u <- rchisq(length(y(model2)), df=dfr(model))
  }
  thetas <- theta(model2)
  sigmas <- sigma(model2)
  if(!is.matrix(thetas)) {
    thetas <- matrix(thetas, nrow=1)
    sigmas <- matrix(sigmas, nrow=1)
  }
  p.comp <- p(model2)
  df <- dfr(model2)
  pooled <- class(model) %in% c("SingleBatchPooled", "MultiBatchPooled")
  K <- seq_len(k(model2))
  ##B <- unique(orig.data$batch_index)
  B <- unique_batch(batch(model2))
  pz <- matrix(NA, length(y(model2)), max(K))
  for(b in B){
    j <- which(batch(model2) == b)
    m <- thetas[b, ]
    if(pooled){
      ss <- sigmas[b]
    } else {
      ss <- sigmas[b, ]
    }
    yy <- y(model2)[j]
    pb <- p.comp[b, ]
    for(k in K){
      if(pooled){
        temp <- pb[k] * dst(yy, df=df, mu=m[k], sigma=ss)
      } else {
        temp <- pb[k] * dst(yy, df=df, mu=m[k], sigma=ss[k])
      }
      pz[j, k] <- temp
    }
  }
  pz2 <- pz/rowSums(pz)
  ## the probz slot expects a frequency
  freq <- pz2 * (iter(model) - 1)
  freq2 <- matrix(as.integer(freq), nrow=nrow(freq), ncol=ncol(freq))
  probz(model2) <- freq2
  ## update z's
  if(class(model2)=="MultiBatchPooled"){
    z(model2) <- z_multibatch_pvar(model2)
  }
  if(class(model2)=="MultiBatchModel"){
    z(model2) <- update_z(model2)
  }
  model2
}

setMethod("u", "MixtureModel", function(object) object@u )


setMethod("dfr", "MixtureModel", function(object) object@hyperparams@dfr )


setReplaceMethod("dfr", c("MixtureModel", "numeric"),
                 function(object, value) {
                   object@hyperparams@dfr <- value
                   object@u <- rchisq(length(y(object)), value)
                   object
                 })

## i can be logical or numeric


#' @aliases "[",MixtureModel,logical,ANY,ANY
#' @rdname subsetting-methods
setMethod("[", c("MixtureModel", "logical"),
          function(x, i, j, ..., drop=TRUE){
            x@data <- oned(x)[i]
            x@z <- z(x)[i]
            x@probz <- x@probz[i, , drop=FALSE]
            x@zfreq <- as.integer(table(z(x)))
            x@batch <- batch(x)[i]
            x@batchElements <- as.integer(table(batch(x)))
            return(x)
})

#' @aliases "[",MixtureModel,numeric,ANY,ANY
#' @rdname subsetting-methods
setMethod("[", c("MixtureModel", "numeric"),
          function(x, i, j, ..., drop=TRUE){
            x@data <- oned(x)[i]
            x@z <- z(x)[i]
            x@probz <- x@probz[i, , drop=FALSE]
            x@zfreq <- as.integer(table(z(x)))
            x@batch <- batch(x)[i]
            x@batchElements <- as.integer(table(batch(x)))
            return(x)
          })

setMethod("useModes", "MixtureModel", function(object){
  m2 <- object
  theta(m2) <- modes(object)[["theta"]]
  sigma2(m2) <- modes(object)[["sigma2"]]
  tau2(m2) <- modes(object)[["tau2"]]
  nu.0(m2) <- modes(object)[["nu0"]]
  sigma2.0(m2) <- modes(object)[["sigma2.0"]]
  p(m2) <- modes(object)[["mixprob"]]
  zFreq(m2) <- as.integer(modes(object)[["zfreq"]])
  log_lik(m2) <- modes(object)[["loglik"]]
  logPrior(m2) <- modes(object)[["logprior"]]
  ##
  ## update z using the modal values from above
  ##
  z(m2) <- updateZ(m2)
  m2
})

setMethod("compute_marginal_lik", "MixtureModel", function(object, params){
  ##
  ## evaluate marginal likelihood. Relax default conditions
  ##
  if(missing(params)){
    params <- mlParams(root=1/2,
                       reject.threshold=exp(-100),
                       prop.threshold=0.5,
                       prop.effective.size=0)
  }
  ml <- tryCatch(marginalLikelihood(object, params),
                 warning=function(w) NULL)
  if(!is.null(ml)){
    marginal_lik(object) <- ml
    message("     marginal likelihood: ", round(marginal_lik(object), 2))
  } else {
    ##warning("Unable to compute marginal likelihood")
    message("Unable to compute marginal likelihood")
  }
  object
})
