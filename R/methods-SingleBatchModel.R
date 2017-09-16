#' @include methods-MixtureModel.R
NULL

#' Constructor for list of single-batch models
#'
#' An object of class MarginalModel is constructed for each k, creating a list of
#' MarginalModels.
#'
#' @param data numeric vector of average log R ratios
#' @param k numeric vector indicating the number of mixture components for each model
#' @param mcmc.params an object of class \code{McmcParams}
#' @param ... additional arguments passed to \code{Hyperparameters}
#' @seealso \code{\link{MarginalModel}} \code{\link{BatchModelList}}
#' @return a list. Each element of the list is a \code{BatchModel}
#' @examples
#' mlist <- MarginalModelList(data=y(MarginalModelExample), k=1:4)
#' mcmcParams(mlist) <- McmcParams(iter=1, burnin=1, nStarts=0)
#' mlist2 <- posteriorSimulation(mlist)
#' @export
SingleBatchModelList <- function(data=numeric(), k=numeric(),
                              mcmc.params=McmcParams(),
                              ...){
  model.list <- vector("list", length(k))
  for(i in seq_along(k)){
    hypp <- Hyperparameters(k=k[i], ...)
    model.list[[i]] <- SingleBatchModel(data=data, k=k[i], mcmc.params=mcmc.params,
                                     hypp=hypp)
  }
  model.list
}


#' Create an object for running single batch MCMC simulations.
#' @examples
#'      model <- SingleBatchModel(data=rnorm(10), k=1)
#' @param data the data for the simulation.
#' @param k An integer value specifying the number of latent classes.
#' @param hypp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mcmc.params An object of class 'McmcParams'
#' @return An object of class 'SingleBatchModel'
#' @export
SingleBatchModel <- function(data=numeric(), k=3, hypp, mcmc.params){
  if(missing(hypp)) hypp <- Hyperparameters(k=3)
  batch <- rep(1L, length(data))
  if(missing(k)){
    k <- k(hypp)
  } else{
    k(hypp) <- k
  }
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  if(missing(hypp)) hypp <- HyperparametersSingleBatch(k=k)
  nbatch <- setNames(as.integer(table(batch)), levels(batch))
  zz <- sample(seq_len(k), length(data), replace=TRUE)
  zfreq <- as.integer(table(zz))
  object <- new("SingleBatchModel",
                k=as.integer(k),
                hyperparams=hypp,
                theta=numeric(k),
                sigma2=numeric(k),
                mu=numeric(1),
                tau2=numeric(1),
                nu.0=1L,
                sigma2.0=1L,
                pi=rep(1/k, k),
                data=data,
                data.mean=numeric(k),
                data.prec=numeric(k),
                z=zz,
                zfreq=zfreq,
                probz=matrix(0, length(data), k),
                logprior=numeric(1),
                loglik=numeric(1),
                mcmc.chains=McmcChains(),
                batch=batch,
                batchElements=nbatch,
                modes=list(),
                mcmc.params=mcmc.params,
                label_switch=FALSE,
                marginal_lik=as.numeric(NA),
                .internal.constraint=5e-4,
                .internal.counter=0L)
  object <- startingValues(object)
}


SingleBatchModel2 <- function(dat=numeric(), hp=Hyperparameters(),
                           mp=McmcParams(iter=1000, burnin=1000,
                                         thin=10, nStarts=4)){
  K <- k(hp)
  ##mu <- rnorm(1, mu.0(hp), sqrt(tau2.0(hp)))
  mu <- rnorm(1, median(dat), sd(dat)) 
  tau2 <- 1/rgamma(1, 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
  p <- rdirichlet(1, alpha(hp))[1, ]
  theta <- sort(rnorm(k(hp), mu, sqrt(tau2)))
  nu.0 <- 3.5
  sigma2.0 <- 0.25
  sigma2 <- 1/rgamma(k(hp), 0.5 * nu.0, 0.5 * nu.0 * sigma2.0)
  object <- new("SingleBatchModel",
                k=as.integer(K),
                hyperparams=hp,
                theta=theta,
                sigma2=sigma2,
                mu=mu,
                tau2=tau2,
                nu.0=nu.0,
                sigma2.0=sigma2.0,
                pi=p,
                data=dat,
                data.mean=numeric(K),
                data.prec=numeric(K),
                z=integer(length(dat)),
                zfreq=integer(K),
                probz=matrix(0, length(dat), K),
                logprior=numeric(1),
                loglik=numeric(1),
                mcmc.chains=McmcChains(),
                batch=rep(1L, length(dat)),
                batchElements=1L,
                modes=list(),
                mcmc.params=mp,
                label_switch=FALSE,
                marginal_lik=as.numeric(NA),
                .internal.constraint=5e-4,
                .internal.counter=0L)
  chains(object) <- McmcChains(object)
  object
}


getK <- function(object){
  hypp <- hyperParams(object)
  getK(hypp)
}

#' @rdname mu-method
#' @aliases mu,SingleBatchModel-method
#' @export
setMethod("mu", "SingleBatchModel", function(object) object@mu)



#' @rdname tau2-method
#' @aliases tau2,SingleBatchModel-method
setMethod("tau2", "SingleBatchModel", function(object) object@tau2)

setMethod("show", "SingleBatchModel", function(object) callNextMethod())

setMethod("computeMeans", "SingleBatchModel", function(object){
  compute_means(object)
})

setMethod("computeVars", "SingleBatchModel", function(object){
  compute_vars(object)
})



setReplaceMethod("tau2", "SingleBatchModel", function(object, value){
  object@tau2 <- value
  object
})


setReplaceMethod("mu", "SingleBatchModel", function(object, value){
  object@mu <- value
  object
})

#' @rdname bic-method
#' @aliases bic,SingleBatchModel-method
setMethod("bic", "SingleBatchModel", function(object){
  object <- useModes(object)
  ## K: number of free parameters to be estimated
  ##   - component-specific parameters:  theta, sigma2   (3 x k(model))
  ##   - mixing probabilities:  k-1
  ##   - length-one parameters: mu, tau2, sigma2.0, nu.0             +4
  K <- 2*k(object) + (k(object)-1) + 4
  n <- length(y(object))
  -2*(log_lik(object) + logPrior(object)) + K*(log(n) - log(2*pi))
})

#' @rdname theta-method
#' @aliases theta,SingleBatchModel-method
setMethod("theta", "SingleBatchModel", function(object) object@theta)

#' @rdname sigma2-method
#' @aliases sigma2,SingleBatchModel-method
setMethod("sigma2", "SingleBatchModel", function(object) object@sigma2)

newMarginalModel <- function(object){
  mp <- mcmcParams(object)
  object2 <- SingleBatchModel(y(object), k=k(object), mcmc.params=mp,
                           hypp=hyperParams(object))
  theta(object2) <- theta(object)
  sigma2(object2) <- sigma2(object)
  p(object2) <- p(object)
  z(object2) <- z(object)
  nu.0(object2) <- nu.0(object)
  mu(object2) <- mu(object)
  tau2(object2) <- tau2(object)
  zFreq(object2) <- zFreq(object)
  probz(object2) <- probz(object)
  sigma2.0(object2) <- sigma2.0(object)
  dataMean(object2) <- dataMean(object)
  dataPrec(object2) <- dataPrec(object)
  log_lik(object2) <- log_lik(object)
  logPrior(object2) <- logPrior(object)
  modes(object2) <- modes(object)
  object2
}

newSingleBatchModel <- function(object){
  mp <- mcmcParams(object)
  object2 <- SingleBatchModel(y(object), k=k(object), mcmc.params=mp,
                           hypp=hyperParams(object))
  theta(object2) <- theta(object)
  sigma2(object2) <- sigma2(object)
  p(object2) <- p(object)
  z(object2) <- z(object)
  nu.0(object2) <- nu.0(object)
  mu(object2) <- mu(object)
  tau2(object2) <- tau2(object)
  zFreq(object2) <- zFreq(object)
  probz(object2) <- probz(object)
  sigma2.0(object2) <- sigma2.0(object)
  dataMean(object2) <- dataMean(object)
  dataPrec(object2) <- dataPrec(object)
  log_lik(object2) <- log_lik(object)
  logPrior(object2) <- logPrior(object)
  modes(object2) <- modes(object)
  object2
}

setMethod("relabel", "SingleBatchModel", function(object, zindex){
  object <- newSingleBatchModel(object)
  if(identical(zindex, seq_len(k(object)))) return(object)
  ##
  ## Permute the labels for the components
  ##
  zz <- factor(z(object), levels=zindex)
  zz <- as.integer(zz)
  z(object) <- zz
  zFreq(object) <- as.integer(table(zz))
  dataMean(object) <- dataMean(object)[zindex]
  dataPrec(object) <- dataPrec(object)[zindex]
  object
})

.computeModesMarginal <- function(object){
  i <- argMax(object)
  mc <- chains(object)
  thetamax <- theta(mc)[i, ]
  sigma2max <- sigma2(mc)[i,]
  pmax <- p(mc)[i, ]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mu(mc)[i],
                tau2=tau2(mc)[i],
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i],
                zfreq=zFreq(mc)[i, ],
                loglik=log_lik(mc)[i],
                logprior=logPrior(mc)[i])
  modes
}

setMethod("computeModes", "SingleBatchModel", function(object){
  .computeModesMarginal(object)
})

setMethod("showMeans", "SingleBatchModel", function(object){
  paste(round(theta(object), 2), collapse=", ")
})


setMethod("showSigmas", "SingleBatchModel", function(object){
  paste(round(sqrt(sigma2(object)), 2), collapse=", ")
})

setMethod("tablez", "SingleBatchModel", function(object) table(z(object)))

modelOtherModes <- function(model, maxperm=5){
  kperm <- permnK(k(model), maxperm)
  model.list <- vector("list", nrow(kperm))
  for(i in seq_along(model.list)){
    model.list[[i]] <- relabel(model, kperm[i, ])
  }
  model.list
}

#' Compute the log bayes factor between models.
#'
#' Models of varying component sizes are compared. The log bayes factor is 
#' calculated comparing each set of two models by marginal likelihood, 
#' as computed by \code{marginalLikelihood}.
#' @param x the result of a call to \code{computeMarginalLik}.
#' @return Log Bayes factor comparing the two models with highest likelihood.
#' @export
logBayesFactor <- function(x) {
    k <- length(x)
    mat <- matrix(0, nrow=k, ncol=k, dimnames=list(names(x), names(x)))
    for (i in seq_len(length(x))) {
        x_i <- x[c(i:k, 0:(i-1))]
        diff <- x_i[1] - x_i
        mat[i, names(diff)] <- diff
    }

    return(mat)
}

setMethod("updateMultinomialProb", "SingleBatchModel", function(object){
  update_multinomialPr(object)
})

setMethod("computeLoglik", "SingleBatchModel", function(object){
  loglik(object)
})

reorderSingleBatchChains <- function(model){
  K <- k(model)
  if(K < 2) return(model)
  ##pstar=.blockUpdates(model, mcmcParams(mlist2[[2]]))
  ch <- chains(model)
  theta.ch <- thetac(model)
  ix <- apply(theta.ch, 1, function(x) paste0(order(x), collapse=","))
  tab.ix <- table(ix)
  ord <- names(tab.ix)[which.max(tab.ix)]
  ord <- as.integer(unlist(strsplit(ord, ",")))
  if(identical(ord, seq_len(K))){
    return(model)
  }
  theta.ch <- theta.ch[, ord]
  sigma.ch <- sigmac(model)
  sigma.ch <- sigma.ch[, ord]
  p.ch <- pic(model)
  p.ch <- p.ch[, ord]
  zfreq.ch <- zFreq(ch)

  ch@theta <- theta.ch
  ch@sigma2 <- sigma.ch^2
  ch@pi <- p.ch
  ch@zfreq <- zfreq.ch[, ord]
  chains(model) <- ch
  model
}
