##failSmallPstar <- function(ptheta.star, params=mlParams()){
##  ptheta.star <- ptheta.star[!is.nan(ptheta.star)]
##  reject.threshold <- params$reject.threshold
##  prop.threshold <- params$prop.threshold
##  small.theta.red <- mean(ptheta.star < reject.threshold)
##  ## might small values inflate the marginal likelihood
##  small.theta.red >= prop.threshold
##}
##
##.message_theta <- function(model, params, logprobs){
##  reject.threshold <- params$reject.threshold
##  prop.threshold <- params$prop.threshold
##  ignore.small.pstar <- params$ignore.small.pstar
##  warnings <- params$warnings
##  if(paramUpdates(model)[["theta"]]==0) {
##    ignore.small.pstar <- TRUE
##  }
##  if(!ignore.small.pstar){
##    failed <- failSmallPstar(logprobs$theta, params)
##    if (failed) {
##      msg <- paste("The model for k=", k(model), " may be overfit.",
##                   "This can lead to an incorrect marginal likelihood")
##      if(warnings) warning(msg)
##      return(FALSE)
##    }
##  }
##  TRUE
##}

.blockUpdatesBatch <- function(model, params){
  model.reduced <- model
  logprobs <- tibble(theta=marginal_theta_batch(model.reduced))
  ##continue <- .message_theta(model, params, logprobs)
  ##if(!continue) return(matrix(NA))
  logprobs$sigma <- reduced_sigma_batch(model.reduced)
  stopifnot(identical(modes(model.reduced), modes(model)))
  ##logprobs$pi <- reduced_pi_batch(model.reduced)
  logprobs$pi <- reduced_pi_batch2(model.reduced)
  ## Block updates for stage 2 parameters
  ##
  logprobs$mu <- reduced_mu_batch(model.reduced)
  stopifnot(identical(modes(model.reduced), modes(model)))
  ##
  ## with theta and mu fixed, nothing stochastic -- no need to do reduced sampler
  ##logprobs$tau2 <- reduced_tau_batch(model.reduced)
  logprobs$tau2 <- log_prob_tau2(model.reduced)
  identical(modes(model.reduced), modes(model))
  logprobs$nu0 <- reduced_nu0_batch(model.reduced)
  identical(modes(model.reduced), modes(model))
  ## nothing stochastic at this point in the reduced gibbs sampler
  logprobs$s20 <- log_prob_s20(model.reduced) ;
  probs <- exp(logprobs)
  probs
}

.blockUpdatesMultiBatchPooled <- function(model, params=mlParams()){
  model.reduced <- model
  logprobs <- tibble(theta=marginal_theta_pooled(model.reduced))
  ##continue <- .message_theta(model, params, logprobs)
  logprobs$sigma2 <- reduced_sigma_pooled(model.reduced)
  logprobs$pi <- reduced_pi_pooled2(model.reduced)
  ##
  ## Block updates for stage 2 parameters
  ##
  logprobs$mu <- reduced_mu_pooled(model.reduced)
  ## nothing stochastic here -- return single value from density
  logprobs$tau2 <- log_prob_tau2(model.reduced)
  logprobs$nu0 <- reduced_nu0_pooled(model.reduced)
  logprobs$s20 <- log_prob_s20p(model.reduced) ;
  probs <- exp(logprobs)
  probs
}

blockUpdates <- function(reduced_gibbs, root) {
  pstar <- apply(reduced_gibbs, 2, function(x) log(mean(x^(root), na.rm=TRUE)))
}

#' Parameters for evaluating marginal likelihood
#'
#' @param root length-one numeric vector. We exponentiate \code{p(theta* | ...)}
#'   by the value of \code{root}. Values less than one reduce the influence of
#'   extreme observations.
#' @param reject.threshold length-one numeric vector between 0 and 1.
#'   Probabilities in the reduced Gibbs model for the thetas that are less than
#'   this threshold are flagged.
#' @param prop.threshold length-one numeric vector between 0 and 1. If more than
#'   \code{prop.threshold} are flagged, the marginal likelihood is not
#'   evaluated.
#' @param prop.effective.size Logical. If the effective size / total iterations
#'   is less than \code{prop.effective.size}, the marginal likelihood is not
#'   evaluated (unless \code{ignore.effective.size} is \code{TRUE}).
#' @param ignore.effective.size Logical. By default, if the effective size of
#'   any theta chain is less than 0.02, the marginal likelihood is not
#'   calculated. If this parameter is set to TRUE, the effective size is
#'   ignored. Occasionally, the effective size is misleading. See details.
#' @param ignore.small.pstar Logical. Flags from the \code{reject.threshold}
#'   parameter are ignored and the marginal likelihood is calculated.
#'
#' @param warnings Logical. If FALSE, warnings are not issued. This is FALSE by
#'   default for the marginalLikelihood-list method, and TRUE otherwise.
#'
#' @details
#'
#'
#'  For mixture models, a low effective size of one or more theta chains can
#'  occur for the following reasons:
#'
#' A. the model has not yet converged
#'
#' B. the model is overfit and there is lots of mixing (label swapping )between
#'   some of the chains
#'
#' C. the model is not overfit but there is a lot of mixing of the thetas
#'
#' For both (A) and (B) it is desirable to return NAs. While (C) can also occur,
#' it can be easily diagnosed by visual inspection of the chains. To the extent
#' that (C) occurs, the correction factor may not be needed.
#'
#' @examples
#' mlParams()
#'
#'
#' @return a list of parameters to be passed to \code{marginalLikelihood}.
#' @seealso \code{\link[coda]{effectiveSize}} \code{\link{marginalLikelihood}}
#' @export
mlParams <- function(root=1/10,
                     reject.threshold=exp(-10),
                     prop.threshold=0.5,
                     prop.effective.size=0.05,
                     ignore.effective.size=FALSE,
                     ignore.small.pstar=FALSE,
                     warnings=TRUE){
  list(root=root,
       reject.threshold=reject.threshold,
       prop.threshold=prop.threshold,
       prop.effective.size=prop.effective.size,
       ignore.effective.size=ignore.effective.size,
       ignore.small.pstar=ignore.small.pstar,
       warnings=warnings)
}

## used for debugging
## computeML <- function(model, params=mlParams()){
##   ## calculate p(x|theta)
##   logLik <- modes(model)[["loglik"]] ## includes 2nd stage
##   stage2.loglik <- stageTwoLogLik(model)
## 
##   ## calculate log p(theta)
##   logPrior <- modes(model)[["logprior"]]
## 
##   ## calculate log p(theta|x)
##   ##mp <- McmcParams(iter=niter)
##   mp <- McmcParams(iter=iter(model))
##   red_gibbs <- .blockUpdates(model, params)
##   pstar <- blockUpdates(red_gibbs, root=1)
## 
##   ## calculate p(x|model)
##   setNames(c(logLik, stage2.loglik, logPrior,
##              sum(pstar), log(factorial(k(model)))),
##            c("loglik", "stage2lik", "logprior",
##              "total.pstar", "correction"))
## }

thetaEffectiveSize <- function(model){
  thetas <- as.data.frame(thetac(model))
  eff_size <- sapply(thetas, effectiveSize)
  eff_size
}

##
##  For mixture models, a low effective size of one or more theta chains can
##  occur for the following reasons:
##
## A. the model has not yet converged
##
## B. the model is overfit and there is lots of mixing (label swapping )between
##   some of the chains
##
## C. the model is not overfit but there is a lot of mixing of the thetas
##
## For both (A) and (B) it is desirable to return NAs. While (C) can also occur,
## it can be easily diagnosed by visual inspection of the chains. To the extent
## that (C) occurs, the correction factor may not be needed
##
##
failEffectiveSize <- function(model, params=mlParams()){
  ##if(params$ignore.effective.size) return(FALSE)
  min.size <- min(thetaEffectiveSize(model))
  THR <- params$prop.effective.size
  (min.size / iter(model) < THR) && paramUpdates(model)[["theta"]] > 0
}

effectiveSizeWarning <- function(model){
  paste("The effective size of one or more theta chains is less than 2%.\n",
        "See ?coda::effectiveSize")
}

.ml_batchmodel <- function(model, params=mlParams()){
  niter <- iter(model)
  root <- params$root
  reject.threshold <- params$reject.threshold
  prop.threshold <- params$prop.threshold
  ## calculate p(x|theta)
  logLik <- modes(model)[["loglik"]] ## includes 2nd stage
  model2 <- useModes(model)
  ## calculate log p(theta)
  logPrior <- modes(model)[["logprior"]]
  mp <- McmcParams(iter=niter)
  red_gibbs <- .blockUpdatesBatch(model2, params)
  pstar <- blockUpdates(red_gibbs, root)
  correction.factor <- log(factorial(k(model)))
  ## calculate p(x|model)
  m.y <- logLik + logPrior - sum(pstar) +
    correction.factor
  if(length(unique(batch(model))) == 1){
    names(m.y) <- paste0("SB", k(model))
  } else {
    names(m.y) <- paste0("MB", k(model))
  }
  m.y
}

.ml_multibatch_pooled <- function(model, params=mlParams()){
  ## calculate effective size of thetas and check against threshold
##  warnings <- params$warnings
##  if (failEffectiveSize(model, params)) {
##    if(warnings) warning(effectiveSizeWarning(model))
##    naresult <- setNames(NA, paste0("MB", k(model)))
##    if(!params$ignore.effective.size)
##      return(naresult)
##  }
  ## get parameters from list params
  niter <- iter(model)
  root <- params$root
  reject.threshold <- params$reject.threshold
  prop.threshold <- params$prop.threshold
  logLik <- modes(model)[["loglik"]] ## includes 2nd stage
  model2 <- useModes(model)
  logPrior <- modes(model)[["logprior"]]
  mp <- McmcParams(iter=niter)
  red_gibbs <- .blockUpdatesMultiBatchPooled(model2, params)
  pstar <- blockUpdates(red_gibbs, root)
##  if(failEffectiveSize(model, params)){
##    ## this can fail because the model is mixing beteen components
##    ## and the correction factor is not needed
##    correction.factor <- 0
  ##  } else
  correction.factor <- log(factorial(k(model)))
  ## calculate p(x|model)
  m.y <- logLik + logPrior - sum(pstar) +
    correction.factor
  if(length(unique(batch(model))) == 1){
    names(m.y) <- paste0("SBP", k(model))
  } else {
    names(m.y) <- paste0("MBP", k(model))
  }
  m.y
}

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,MultiBatchModel-method marginalLikelihood,MultiBatchModel,ANY-method
setMethod("marginalLikelihood", "MultiBatchModel",
          function(model, params=mlParams()){
            .ml_batchmodel(model, params)
          })

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,MultiBatchPooled-method marginalLikelihood,MultiBatchModel,ANY-method
setMethod("marginalLikelihood", "MultiBatchPooled",
          function(model, params=mlParams()){
            .ml_multibatch_pooled(model, params)
          })

.ml_list <- function(model.list, params=mlParams(warnings=FALSE)){
  ml <- sapply(model.list, marginalLikelihood, params=params)
  return(ml)
}

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,list-method marginalLikelihood,list,ANY-method
setMethod("marginalLikelihood", "list",
          function(model, params=mlParams(warnings=FALSE)){
            .ml_list(model, params)
          })
