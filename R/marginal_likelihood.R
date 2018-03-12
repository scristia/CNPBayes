failSmallPstar <- function(ptheta.star, params=mlParams()){
  ptheta.star <- ptheta.star[!is.nan(ptheta.star)]
  reject.threshold <- params$reject.threshold
  prop.threshold <- params$prop.threshold
  small.theta.red <- mean(ptheta.star < reject.threshold)
  ## might small values inflate the marginal likelihood
  small.theta.red >= prop.threshold
}

.blockUpdates <- function(model, params){
  ignore.small.pstar <- params$ignore.small.pstar
  warnings <- params$warnings
  model.reduced <- model
  ## Let theta1 = theta
  ## and theta2 = [sigma2, tau2, mu, nu.0]  (all parameters except theta)
  ##
  ## p(theta* | y) = sum_g^G p(theta* | y, theta2^(g), z^(g))
  ptheta.star <- marginal_theta(model)
  ##
  ## ptheta.star
  ##
  if(paramUpdates(model)[["theta"]]==0) {
    ignore.small.pstar <- TRUE
  }
  if(!ignore.small.pstar){
    failed <- failSmallPstar(ptheta.star, params)
    if (failed) {
      msg <- paste("Small values for p(theta* | ...) in the reduced Gibbs.")
      if(warnings) warning(msg)
      return(NA)
    }
  }
  ##
  ## p(sigma2* | y, theta1*) = sum_g=1^G p(sigma2* | y, theta1*, z^(g), theta2^(g))
  ##
  model.psigma2 <- reduced_sigma(model.reduced)
  psigma.star <- p_sigma_reduced(model.psigma2)

  model.pistar <- reduced_pi(model.reduced)
  p.pi.star <- p_pmix_reduced(model.pistar)

  ##
  ## Block updates for stage 2 parameters
  ##
  model.mustar <- reduced_mu(model.reduced)
  p.mustar <- p_mu_reduced(model.mustar)

  model.taustar <- reduced_tau(model.reduced)
  p.taustar <- p_tau_reduced(model.mustar)

  model.nu0star <- reduced_nu0(model.reduced)
  p.nu0star <- p_nu0_reduced(model.nu0star)

  model.s20star <- reduced_s20(model.reduced)
  p.s20star <- p_s20_reduced(model.s20star)

  reduced_gibbs <- cbind(ptheta.star, psigma.star, p.mustar, p.pi.star,
                         p.taustar, p.nu0star, p.s20star)

  colnames(reduced_gibbs) <- c("theta", "sigma", "pi", "mu",
                               "tau", "nu0", "s20")
  reduced_gibbs
}

.blockUpdatesPooledVar <- function(model, params){
  reject.threshold <- params$reject.threshold
  prop.threshold <- params$prop.threshold
  ignore.small.pstar <- params$ignore.small.pstar
  warnings <- params$warnings
  model.reduced <- model
  ##
  ## Block updates for stage 1 parameters
  ##
  ptheta.star <- full_theta_pooled(model)
  if(paramUpdates(model)[["theta"]]==0) {
    ignore.small.pstar <- TRUE
  }
  if(!ignore.small.pstar){
    failed <- failSmallPstar(ptheta.star, params)
    if (failed) {
      msg <- paste("The model for k=", k(model), " may be overfit.",
                   "This can lead to an incorrect marginal likelihood")
      if(warnings) warning(msg)
      return(NA)
    }
  }
  model.psigma2 <- reduced_sigma_pooled(model.reduced)
  identical(modes(model.psigma2), modes(model))
  psigma.star <- p_sigma_reduced_pooled(model.psigma2)

  model.pistar <- reduced_pi_pooled(model.reduced)
  identical(modes(model.pistar), modes(model))
  p.pi.star <- p_pmix_reduced_pooled(model.pistar)

  ##
  ## Block updates for stage 2 parameters
  ##
  model.mustar <- reduced_mu_pooled(model.reduced)
  stopifnot(identical(modes(model.mustar), modes(model)))
  p.mustar <- p_mu_reduced_pooled(model.mustar)

  model.taustar <- reduced_tau_pooled(model.reduced)
  identical(modes(model.taustar), modes(model))
  p.taustar <- p_tau_reduced_pooled(model.mustar)

  model.nu0star <- reduced_nu0_pooled(model.reduced)
  identical(modes(model.nu0star), modes(model))
  p.nu0star <- p_nu0_reduced_pooled(model.nu0star)

  model.s20star <- reduced_s20_pooled(model.reduced)
  p.s20star <- p_s20_reduced_pooled(model.s20star)

  reduced_gibbs <- cbind(ptheta.star, psigma.star, p.mustar, p.pi.star,
                         p.taustar, p.nu0star, p.s20star)
  colnames(reduced_gibbs) <- c("theta", "sigma", "pi", "mu",
                               "tau", "nu0", "s20")
  reduced_gibbs
}

reduced <- function(model, params=mlParams()){
  .blockUpdatesBatch(model, params)
}

.blockUpdatesBatch <- function(model, params){
  reject.threshold <- params$reject.threshold
  prop.threshold <- params$prop.threshold
  ignore.small.pstar <- params$ignore.small.pstar
  warnings <- params$warnings
  model.reduced <- model
  logprobs <- tibble(theta=marginal_theta_batch(model))
  if(paramUpdates(model)[["theta"]]==0) {
    ignore.small.pstar <- TRUE
  }
  ##tmp <- cbind(p.theta=ptheta.star, theta=thetac(model)) %>%
  ##as.tibble
  if(!ignore.small.pstar){
    failed <- failSmallPstar(log_ptheta, params)
    if (failed) {
      msg <- paste("The model for k=", k(model), " may be overfit.",
                   "This can lead to an incorrect marginal likelihood")
      if(warnings) warning(msg)
      return(matrix(NA))
    }
  }
  ##model.psigma2 <- reduced_sigma_batch(model.reduced)
  log_psigma <- reduced_sigma_batch(model.reduced)
  stopifnot(identical(modes(model.reduced), modes(model)))

  ##psigma.star <- p_sigma_reduced_batch(model.psigma2)
  log_pmix <- reduced_pi_batch(model.reduced)
  ##identical(modes(model.pistar), modes(model))
  ##p.pi.star <- p_pmix_reduced_batch(model.pistar)

  ##
  ## Block updates for stage 2 parameters
  ##
  log_pmu <- reduced_mu_batch(model.reduced)
  ##model.mustar <- reduced_mu_batch(model.reduced)
  stopifnot(identical(modes(model.reduced), modes(model)))
  ##p.mustar <- p_mu_reduced_batch(model.mustar)

  log_ptau2 <- reduced_tau_batch(model.reduced)
  identical(modes(model.reduced), modes(model))
  ##p.taustar <- p_tau_reduced_batch(model.mustar)

  log_pnu0 <- reduced_nu0_batch(model.reduced)
  identical(modes(model.reduced), modes(model))
  ##p.nu0star <- p_nu0_reduced_batch(model.nu0star)

  log_ps20 <- reduced_s20_batch(model.reduced)
  reduced_gibbs <- tibble(theta=log_ptheta,
                          sigma=log_psigma,
                          mu=log_pmu,
                          pmix=log_pmix,
                          tau2=log_ptau2,
                          nu0=log_pnu0,
                          s20=log_ps20)
  reduced_gibbs
}

.blockUpdatesMultiBatchPooled <- function(model, params=mlParams()){
  reject.threshold <- params$reject.threshold
  prop.threshold <- params$prop.threshold
  ignore.small.pstar <- params$ignore.small.pstar
  warnings <- params$warnings

  model.reduced <- model
  ##ptheta.star <- marginal_theta_batch(model)
  ptheta.star <- theta_multibatch_pvar_red(model)
  if(paramUpdates(model)[["theta"]]==0) {
    ignore.small.pstar <- TRUE
  }
  if(!ignore.small.pstar){
    failed <- failSmallPstar(ptheta.star, params)
    if (failed) {
      msg <- paste("The model for k=", k(model), " may be overfit.",
                   "This can lead to an incorrect marginal likelihood")
      if(warnings) warning(msg)
      return(matrix(NA))
    }
  }
  model.psigma2 <- sigma_multibatch_pvar_red(model.reduced)
  psigma.star <- psigma_multibatch_pvar_red(model.psigma2)

  model.pistar <- pi_multibatch_pvar_red(model.reduced)
  p.pi.star <- p_pmix_reduced_batch(model.pistar)

  ##
  ## Block updates for stage 2 parameters
  ##
  model.mustar <- mu_multibatch_pvar_red(model.reduced)
  p.mustar <- p_mu_reduced_batch(model.mustar)

  model.taustar <- tau_multibatch_pvar_red(model.reduced)
  p.taustar <- p_tau_reduced_batch(model.mustar)

  model.nu0star <- nu0_multibatch_pvar_red(model.reduced)
  p.nu0star <- pnu0_multibatch_pvar_red(model.nu0star)

  model.s20star <- s20_multibatch_pvar_red(model.reduced)
  p.s20star <- ps20_multibatch_pvar_red(model.s20star)

  reduced_gibbs <- cbind(ptheta.star, psigma.star, p.mustar, p.pi.star,
                         p.taustar, p.nu0star, p.s20star)

  colnames(reduced_gibbs) <- c("theta", "sigma", "pi", "mu",
                               "tau", "nu0", "s20")

  reduced_gibbs
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

.ml_singlebatch <- function(model, params=mlParams()){
  reject.threshold <- params$reject.threshold
  prop.threshold <- params$prop.threshold
  warnings <- params$warnings
  if (failEffectiveSize(model, params)) {
    msg <- effectiveSizeWarning(model)
    if(warnings) warning(msg)
    naresult <- setNames(NA, paste0("SB", k(model)))
    if(!params$ignore.effective.size)
      return(naresult)
  }
  ## get parameters from list params
  ##niter <- params$niter
  niter <- iter(model)
  root <- params$root

  ## calculate p(x|theta)
  logLik <- modes(model)[["loglik"]] ## includes 2nd stage
  model2 <- useModes(model)
  ##stage2.loglik <- stageTwoLogLik(model2)

  ## calculate log p(theta)
  logPrior <- modes(model)[["logprior"]]

  ## calculate log p(theta|x)
  mp <- McmcParams(iter=niter)
  red_gibbs <- .blockUpdates(model2, params)
  if(length(red_gibbs) == 1){
    ##if(is.na(red_gibbs[1])){
    names(red_gibbs) <- paste0("SB", k(model))
    return(red_gibbs)
  }
  pstar <- blockUpdates(red_gibbs, root)
  if(failEffectiveSize(model, params)){
    ## this can fail because the model is mixing beteen components
    ## and the correction factor is not needed
    correction.factor <- 0
  } else correction.factor <- log(factorial(k(model)))

  ## calculate log p(x|model)
  m.y <- logLik + logPrior - sum(pstar) +
    correction.factor
  names(m.y) <- paste0("SB", k(model))
  m.y
}



## used for debugging
computeML <- function(model, params=mlParams()){
  ## calculate p(x|theta)
  logLik <- modes(model)[["loglik"]] ## includes 2nd stage
  stage2.loglik <- stageTwoLogLik(model)

  ## calculate log p(theta)
  logPrior <- modes(model)[["logprior"]]

  ## calculate log p(theta|x)
  ##mp <- McmcParams(iter=niter)
  mp <- McmcParams(iter=iter(model))
  red_gibbs <- .blockUpdates(model, params)
  pstar <- blockUpdates(red_gibbs, root=1)

  ## calculate p(x|model)
  setNames(c(logLik, stage2.loglik, logPrior,
             sum(pstar), log(factorial(k(model)))),
           c("loglik", "stage2lik", "logprior",
             "total.pstar", "correction"))
}

thetaEffectiveSize <- function(model){
  thetas <- as.data.frame(thetac(model))
  eff_size <- sapply(thetas, effectiveSize)
  ##  nr <- nrow(thetas)
  ##  chunks <- sort(rep(1:5, length.out=nr))
  ##  thetas.list <- split(thetas, chunks)
  ##  ## protects against high autocorrelation induced by label switching
  ##  eff_size <- t(sapply(thetas.list, effectiveSize))
  ##  eff_size_total <- colSums(eff_size)
  ##  ##eff_size_theta <- min(effectiveSize(theta(chains(model))))
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

.ml_pooled <- function(model, params=mlParams()){
  warnings <- params$warnings
  if (failEffectiveSize(model, params)) {
    if(warnings) warning(effectiveSizeWarning(model))
    naresult <- setNames(NA, paste0("SBP", k(model)))
    if(!params$ignore.effective.size){
      return(naresult)
    }
  }
  # get parameters from list params
  niter <- iter(model)
  root <- params$root
  reject.threshold <- params$reject.threshold
  prop.threshold <- params$prop.threshold

  # calculate p(x|theta)
  logLik <- modes(model)[["loglik"]] ## includes 2nd stage
  model2 <- useModes(model)

  # calculate log p(theta)
  logPrior <- modes(model)[["logprior"]]

  # calculate log p(theta|x)
  mp <- McmcParams(iter=niter)
  red_gibbs <- .blockUpdatesPooledVar(model2, params)
  pstar <- blockUpdates(red_gibbs, root)
  if(failEffectiveSize(model, params)){
    ## this can fail because the model is mixing beteen components
    ## and the correction factor is not needed
    correction.factor <- 0
  } else correction.factor <- log(factorial(k(model)))
  ## calculate p(x|model)
  m.y <- logLik + logPrior - sum(pstar) +
    correction.factor
  names(m.y) <- paste0("SBP", k(model))
  m.y
}



.ml_batchmodel <- function(model, params=mlParams()){
  ## calculate effective size of thetas and check against threshold
  warnings <- params$warnings
  if (failEffectiveSize(model, params)) {
    if(warnings) warning(effectiveSizeWarning(model))
    naresult <- setNames(NA, paste0("MB", k(model)))
    if(!params$ignore.effective.size)
      return(naresult)
  }
  ## get parameters from list params
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
  if(failEffectiveSize(model, params)){
    ## this can fail because the model is mixing beteen components
    ## and the correction factor is not needed
    correction.factor <- 0
  } else correction.factor <- log(factorial(k(model)))
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
  warnings <- params$warnings
  if (failEffectiveSize(model, params)) {
    if(warnings) warning(effectiveSizeWarning(model))
    naresult <- setNames(NA, paste0("MB", k(model)))
    if(!params$ignore.effective.size)
      return(naresult)
  }
  ## get parameters from list params
  niter <- iter(model)
  root <- params$root
  reject.threshold <- params$reject.threshold
  prop.threshold <- params$prop.threshold

  ## calculate p(x|theta)
  logLik <- modes(model)[["loglik"]] ## includes 2nd stage
  model2 <- useModes(model)
  ##stage2.loglik <- stageTwoLogLik_pooled(model2)

  ## calculate log p(theta)
  logPrior <- modes(model)[["logprior"]]

  mp <- McmcParams(iter=niter)
  red_gibbs <- .blockUpdatesMultiBatchPooled(model2, params)
  pstar <- blockUpdates(red_gibbs, root)
  if(failEffectiveSize(model, params)){
    ## this can fail because the model is mixing beteen components
    ## and the correction factor is not needed
    correction.factor <- 0
  } else correction.factor <- log(factorial(k(model)))
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
#' @aliases marginalLikelihood,SingleBatchModel-method marginalLikelihood,SingleBatchModel,ANY-method
setMethod("marginalLikelihood", "SingleBatchModel",
          function(model, params=mlParams()) {
            .ml_singlebatch(model, params)
          })

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,SingleBatchPooled-method marginalLikelihood,SingleBatchPooled,ANY-method
setMethod("marginalLikelihood", "SingleBatchPooled",
          function(model, params=mlParams()){
            ## calculate effective size of thetas and check against threshold
            .ml_pooled(model, params)
          })

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


## A copy of marginal_theta_batch for debugging
r_theta_multibatch <- function(model){
  model2 <- useModes(model)
  thetastar <- theta(model2)
  tau2c <- tau2(chains(model))
  muc <- mu(chains(model))
  Z <- z(chains(model))
  model.tmp <- model
  sigma2c <- sigma2(chains(model))
  B <- nrow(theta(model))
  K <- k(model)
  df <- dfr(model)
  S <- nrow(theta(chains(model)))
  p_theta <- rep(NA, S)
  ## we need to save the u-chains
  for(s in seq_len(S)){
    tauc <- sqrt(tau2c[s, ])
    zz <- Z[s, ]
    z(model.tmp) <- zz
    ## make a B x k matrix of the counts
    n_hb <- tableBatchZ(model);
    data_mean <- compute_heavy_sums_batch(model);
    sumu <- compute_u_sums_batch(model) ;
    mu <- muc[s, ]
    tau2 <- tau2c[s, ]
    tau2_tilde = 1.0 / tau2;
    sigma2 <- matrix(sigma2c[s, ], B, K)
    ##invs2 = 1.0 / sigma2(s, Rcpp::_);    // this is a vector of length B*K
    invs2 <- 1.0/sigma2
    ##sigma2_tilde = Rcpp::as<Rcpp::NumericVector>(toMatrix(invs2, B, K));
    sigma2_tilde <- matrix(invs2, B, K)
    prod <- 1
    heavyn <- 0
    heavy_mean <- 0
    for (b in seq_len(B)) {
      for(k in seq_len(K)){
        heavyn <- n_hb[b, k] * sumu[b, k] / df;
        heavy_mean = data_mean[b, k] / df;
        post_prec <- 1.0/tau2[k] + heavyn*1.0/sigma2[b, k] ;
        tau_n <- sqrt(1.0/post_prec) ;
        w1 = (1.0/tau2[k])/post_prec ;
        w2 = (n_hb[b, k] * 1.0/sigma2[b, k])/post_prec ;
        mu_n = w1*mu[k] + w2*heavy_mean ;
        tmp = dnorm(thetastar[b, k], mu_n, tau_n);
        prod <- prod * tmp;
      }
    }
    p_theta[s] = prod;
  }
}

## R version of marginal_theta
r_marginal_theta <- function(model){
  model2 <- useModes(model)
  thetastar <- theta(model2)
  S <- nrow(theta(chains(model)))
  Z <- z(chains(model))
  ##U <- u(chains(model))
  K <- k(model)
  N <- length(y(model))
  df <- dfr(model)
  tau2c <- tau2(chains(model))
  sigma2 <- sigma2(chains(model))
  muc <- mu(chains(model))
  p_theta <- rep(NA, S)
  for(s in seq_len(S)){
    zz <- Z[s, ]
    ##zz = Z(s, _) ;
    ##uu = U(s, _) ;
    ##uu <-  rchisq(N, df) ;
    ##uu <- U[s, ]
    z(model) <- zz
    model@u <- uu
    ##model.slot("z") = zz ;
    ##model.slot("u") = uu ;
    counts <- tableZ(K, zz) ;
    ##NumericVector data_mean =  compute_heavy_means(xmod) ;
    data_mean <- compute_heavy_means(model)/df
    ##data_mean =  data_mean/df ;
    ##NumericVector sumu = compute_u_sums(xmod) ;
    ##sumu <- compute_u_sums(model)/df
    ##sumu = sumu/df ;
    ##nn = as<NumericVector>(counts) * sumu ;
    nn <- counts
    ##nn <- counts*sumu
    nn <- sumu
    tau2_tilde = 1/tau2c[s]
    ##CHECK: length-k vector?
    s2 <- sigma2[s, ]
    ##sigma2_tilde = 1.0/sigma2[s, ] ;
    ##CHECK: length-k vector?
    sigma2_tilde <- 1.0/s2
    ##tmp = dnorm(thetastar, muc[s], tauc[s]) ;
    ##double prod = 1.0;
    prod <- 1.0
    ##for(int k = 0; k < K; ++k) {
    for(k in seq_len(K)){
      post_prec <- tau2_tilde + sigma2_tilde[k] * nn[k];
      tau_n = sqrt(1/post_prec);
      w1 = tau2_tilde/post_prec;
      w2 = nn[k]*sigma2_tilde[k]/post_prec;
      mu_n = w1*muc[s] + w2*data_mean[k];
      tmp = dnorm(thetastar[k], mu_n, tau_n) ;
      prod = prod * tmp[1] ;
    }
    p_theta[s] = prod ;
  }
}

##p_theta_batch2 <- function(model, thetastar){
##  prod <- 1.0
##  B <- length(unique(batch(model)))
##  K <- k(model)
##  n_hb <- tableBatchZ(model)
##  data_mean = compute_heavy_sums_batch(model);
##  sumu = compute_u_sums_batch(model) ;
##  df <- dfr(model)
##  ##invs2 <- 1.0 / sigma2(s, Rcpp::_);    // this is a vector of length B*K
##  invs2 <- 1.0/sigma2(model)
##  ##sigma2_tilde <- toMatrix(invs2, B, K);
##  sigma2_tilde <- invs2
##  for (b in seq_len(B)){
##    for(k in seq_len(K)){
##      heavyn <- n_hb[b, k] * sumu[b, k] / df;
##      heavy_mean <- data_mean[b, k] / df;
##      post_prec <- 1.0/tau2[k] + heavyn*1.0 * sigma2_tilde[b, k] ;
##      tau_n <- sqrt(1.0/post_prec) ;
##      w1 <- (1.0/tau2[k])/post_prec ;
##      w2 <- (n_hb[b, k] * sigma2_tilde[b, k])/post_prec ;
##      mu_n <- w1*mu[k] + w2*heavy_mean ;
##      theta <- thetastar[b, k];
##      tmp <- dnorm(theta, mu_n, tau_n);
##      prod <- prod*tmp[0];
##    }
##  }
##  prod
##}
##
##r_marginal_theta_batch <- function(model){
##  thetastar <- modes(model)[["theta"]];
##  xmod <- model
##  S <- iter(model)
##  p_theta <- rep(NA, S)
##  for(i in seq_len(S)){
##    ##model.slot("z") = update_z_batch(xmod) ;
##    z(model) <- update_z_batch(xmod) ;
##    ##model.slot("zfreq") = tableZ(K, model.slot("z")) ;
##    zFreq(model) <- tableZ(K, z(model)) 
##    ##model.slot("data.mean") = compute_means_batch(xmod) ;
##    dataMean(model) <- compute_means_batch(xmod) ;
##    ##model.slot("data.prec") = compute_prec_batch(xmod) ;
##    ##dataPrec(model) <- compute_prec_batch(xmod) ;
##    ##model.slot("theta") = update_theta_batch(xmod) ;
##    theta(model) <- update_theta_batch(xmod) ;
##    model.slot("sigma2") = update_sigma2_batch(xmod) ;
##    model.slot("mu") = update_mu_batch(xmod) ;
##    model.slot("tau2") = update_tau2_batch(xmod) ;
##    model.slot("sigma2.0") = update_sigma20_batch(xmod) ;
##    model.slot("nu.0") = update_nu0_batch(xmod) ;
##    model.slot("pi") = update_p_batch(xmod) ;
##    model.slot("u") = Rcpp::rchisq(N, df) ;
##    p_theta[i] <- p_theta_batch(model, thetastar)
##  }
##  p_theta
##
##  for(s in seq_len(S)){
##    tauc <- sqrt(tau2c[s, ])
##    ##tauc = sqrt(tau2c(s, Rcpp::_));
##    ##// extract z from the chain
##    ##zz = Z(s, Rcpp::_);
##    zz <- Z[s, ]
##    ##// add z to the current slot in the model
##    model@z <- zz
##    ## make a B x k matrix of the counts
##    n_hb = tableBatchZ(model);
##    data_mean <- compute_means_batch(model)
##    ##data_mean = compute_heavy_sums_batch(model);
##    ##sumu = compute_u_sums_batch(model) ;
##    mu = muc(s, Rcpp::_);
##    tau2 = tau2c(s, Rcpp::_);
##    tau2_tilde = 1.0 / tau2;
##    invs2 = 1.0 / sigma2(s, Rcpp::_);    // this is a vector of length B*K
##    sigma2_tilde = Rcpp::as<Rcpp::NumericVector>(toMatrix(invs2, B, K));
##    double prod = 1.0;
##    double heavyn = 0.0;
##    double heavy_mean = 0.0;
##    for (int b = 0; b < B; ++b) {
##      for(int k = 0; k < K; ++k){
##        heavyn = n_hb(b, k) * sumu(b, k) / df;
##        heavy_mean = data_mean(b, k) / df;
##        post_prec = 1.0/tau2[k] + heavyn*1.0 * sigma2_tilde(b, k) ;
##        if (post_prec == R_PosInf) {
##          throw std::runtime_error("Bad simulation. Run again with different start.");
##        }
##        tau_n = sqrt(1.0/post_prec) ;
##        w1 = (1.0/tau2[k])/post_prec ;
##        w2 = (n_hb(b, k) * sigma2_tilde(b, k))/post_prec ;
##        // mu_n = w1*mu[k] + w2*data_mean(b, k) ;
##        mu_n = w1*mu[k] + w2*heavy_mean ;
##        //theta_new(b, k) = as<double>(rnorm(1, mu_n, tau_n)) ;
##        theta[0] = thetastar(b, k);
##        tmp = dnorm(theta, mu_n, tau_n);
##        prod *= tmp[0];
##      }
##    }
##    p_theta[s] = prod;
##  }
##  return p_theta;
##}
