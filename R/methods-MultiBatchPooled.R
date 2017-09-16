#' @include methods-MultiBatchModel.R
NULL

MultiBatchPooled <- function(dat=numeric(),
                             hp=HyperparametersMultiBatch(),
                             mp=McmcParams(iter=1000, burnin=1000,
                                           thin=10, nStarts=4),
                             batches=integer()){
  obj <- MultiBatchModel2(dat, hp, mp, batches)
  ## average variances across components
  obj@sigma2 <- rowMeans(sigma2(obj))
  ch <- chains(obj)
  nb <- length(unique(batches))
  ch@sigma2 <- matrix(NA, iter(obj), nb)
  obj@mcmc.chains <- ch
  as(obj, "MultiBatchPooled")
}

#' @rdname sigma2-method
#' @aliases sigma2,MultiBatchModel-method
setMethod("sigma2", "MultiBatchPooled", function(object) {
  s2 <- object@sigma2
  ##s2 <- matrix(s2, nBatch(object), k(object))
  names(s2) <- uniqueBatch(object)
  s2
})

setReplaceMethod("sigma2", "MultiBatchPooled", function(object, value){
  names(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})


setMethod("sigmaMean", "MultiBatchPooled", function(object) {
  mns <- colMeans(sigmac(object))
  ##mns <- matrix(mns, nBatch(object), k(object))
  names(mns) <- uniqueBatch(object)
  mns
})

.modesMultiBatchPooled <- function(object){
  i <- argMax(object)
  mc <- chains(object)
  B <- nBatch(object)
  K <- k(object)
  thetamax <- matrix(theta(mc)[i, ], B, K)
  sigma2max <- sigma2(mc)[i, ]
  pmax <- p(mc)[i, ]
  mumax <- mu(mc)[i, ]
  tau2max <- tau2(mc)[i,]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mumax,
                tau2=tau2max,
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i],
                zfreq=zFreq(mc)[i, ],
                loglik=log_lik(mc)[i],
                logprior=logPrior(mc)[i])
  modes
}

setMethod("computeModes", "MultiBatchPooled", function(object){
  .modesMultiBatchPooled(object)
})

setMethod("computeLoglik", "MultiBatchModel", function(object){
  loglik_multibatch_pvar(object)
})
