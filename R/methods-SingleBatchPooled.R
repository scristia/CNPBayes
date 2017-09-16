#' @include methods-SingleBatchModel.R
NULL

SingleBatchPooledVar <- function(data=numeric(), k=2, hypp, mcmc.params){
  obj <- SingleBatchModel(data, k, hypp, mcmc.params)
  obj@sigma2 <- mean(sigma2(obj))
  ch <- chains(obj)
  ch@sigma2 <- matrix(NA, iter(obj), 1)
  obj@mcmc.chains <- ch
  as(obj, "SingleBatchPooledVar")
}

## ad-hoc constructor
SingleBatchPooledVar2 <- function(data=numeric(), hp=Hyperparameters(),
                                  mp=McmcParams(iter=1000, burnin=1000,
                                                thin=10, nStarts=4)){
  obj <- SingleBatchModel2(data, hp, mp)
  obj@sigma2 <- mean(sigma2(obj))
  ch <- chains(obj)
  ch@sigma2 <- matrix(NA, iter(obj), 1)
  obj@mcmc.chains <- ch
  as(obj, "SingleBatchPooledVar2")
}

setValidity("SingleBatchPooledVar", function(object){
  s2 <- sigma2(object)
  if(length(s2) != 1){
    return("sigma2 slot should be length-one numeric vector")
  }
  TRUE
})
