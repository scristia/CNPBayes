blockUpdatesBatch <- function(model, mp){
  ##
  ## Block updates for stage 1 parameters
  ##
  pstar <- setNames(rep(NA, 7),
                    c("theta",
                      "sigma",
                      "pi",
                      "mu",
                      "tau",
                      "nu0",
                      "s20"))

  ptheta.star <- marginal_theta_batch(model)
  pstar["theta"] <- log(mean(ptheta.star))

  model.reduced <- model
  mcmcParams(model.reduced, force=TRUE) <- mp

  model.psigma2 <- reduced_sigma_batch(model.reduced)
  identical(modes(model.psigma2), modes(model))
  psigma.star <- p_sigma_reduced_batch(model.psigma2)
  pstar["sigma"] <- log(mean(psigma.star))


  model.pistar <- reduced_pi_batch(model.reduced)
  identical(modes(model.pistar), modes(model))
  p.pi.star <- p_pmix_reduced_batch(model.pistar)
  pstar["pi"] <- log(mean(p.pi.star))
  ##
  ## Block updates for stage 2 parameters
  ##
  model.mustar <- reduced_mu_batch(model.reduced)
  stopifnot(identical(modes(model.mustar), modes(model)))
  p.mustar <- p_mu_reduced_batch(model.mustar)
  pstar["mu"] <- log(mean(p.mustar))

  model.taustar <- reduced_tau_batch(model.reduced)
  identical(modes(model.taustar), modes(model))
  p.taustar <- p_tau_reduced_batch(model.mustar)
  pstar["tau"] <- log(p.taustar)

  model.nu0star <- reduced_nu0_batch(model.reduced)
  identical(modes(model.nu0star), modes(model))
  p.nu0star <- p_nu0_reduced_batch(model.nu0star)
  pstar["nu0"] <- log(mean(p.nu0star))

  model.s20star <- reduced_s20_batch(model.reduced)
  p.s20star <- p_s20_reduced_batch(model.s20star)
  pstar["s20"] <- log(p.s20star)

  pstar
}

#' Compute the hierarchical likelihood of a converged model.
#' @param model An object of class \code{BatchModel}, or a list of 
#'        \code{BatchModel}'s.
#' @param niter The number of iterations for the reduced Gibb's sampler.
#' @export
batchLikelihood <- function(model, niter=1000L){
  if (!is.list(model)) {
    mlist <- list(model)
  }

  sapply(mlist, 
         function(model, niter) {
           mp <- McmcParams(iter=niter)
           logLik <- modes(model)[["loglik"]] ## includes 2nd stage
           model2 <- useModes(model)
           stage2.loglik <- stageTwoLogLik(model2)
           logPrior <- modes(model)[["logprior"]]
           pstar <- blockUpdatesBatch(model, mp)
           m.y <- logLik + stage2.loglik + logPrior - sum(pstar) + 
                  log(factorial(k(model)))
           m.y
         }, niter=niter)
}
