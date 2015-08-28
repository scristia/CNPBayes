#' Compute the hierarchical likelihood of a converged model.
#' @param model An object of class \code{BatchModel}, or a list of 
#'        \code{BatchModel}'s.
#' @param niter The number of iterations for the reduced Gibb's sampler.
#' @export
batchLikelihood <- function(model, niter=1000L){
  if (!is.list(model)) {
    mlist <- list(model)
  } else {
    mlist <- model
  }

  sapply(mlist, 
         function(model, niter) {
           mp <- McmcParams(iter=niter)
           logLik <- modes(model)[["loglik"]] ## includes 2nd stage
           model2 <- useModes(model)
           stage2.loglik <- stageTwoLogLikBatch(model2)
           logPrior <- modes(model)[["logprior"]]
           pstar <- blockUpdatesBatch(model, mp)
           m.y <- logLik + stage2.loglik + logPrior - sum(pstar) + 
                  log(factorial(k(model)))
           m.y
         }, niter=niter)
}
