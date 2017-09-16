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
setMethod("sigma2", "MultiBatchModel", function(object) {
  s2 <- object@sigma2
  ##s2 <- matrix(s2, nBatch(object), k(object))
  names(s2) <- uniqueBatch(object)
  s2
})
