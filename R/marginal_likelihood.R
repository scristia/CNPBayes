completeLogLik <- function(model){
  model <- useModes(model)
  loglik <- modes(model)[["loglik"]]
  n0 <- nu.0(model)
  sigma2.0 <- sigma2.0(model)
  stage2.loglik <- sum(log(dnorm(theta(model), mu(model), tau(model)) *
                               dgamma(1/sigma2(model), 1/2*n0, 1/2*n0*sigma2.0)))
  result <- loglik + stage2.loglik
  result
}

blockUpdates <- function(model, mp){
  ##
  ## Block updates for stage 1 parameters
  ##
  pstar <- setNames(rep(NA, 7),
                    "theta",
                    "sigma",
                    "pi",
                    "mu",
                    "tau",
                    "nu0",
                    "s20")

  ptheta.star <- .Call("marginal_theta", model)
  pstar["theta"] <- log(mean(ptheta.star))

  model.reduced <- model
  mcmcParams(model.reduced, force=TRUE) <- mp

  model.psigma2 <- .Call("reduced_thetafixed", model.reduced)
  identical(modes(model.psigma2), modes(model))
  psigma.star <- .Call("p_sigma2_zpermuted", model.psigma2)
  pstar["sigma"] <- log(mean(psigma.star))


  model.pistar <- .Call("reduced_theta_sigma_fixed", model.reduced)
  identical(modes(model.pistar), modes(model))
  p.pi.star <- .Call("p_pmix_reduced", model.pistar)
  pstar["pi"] <- log(mean(p.pi.star))
  ##
  ## Block updates for stage 2 parameters
  ##
  model.mustar <- .Call("reduced_mu", model.reduced)
  identical(modes(model.mustar), modes(model))
  pstar["mu"] <- .Call("p_mu_reduced", model.mustar)

  model.taustar <- .Call("reduced_tau", model.reduced)
  identical(modes(model.taustar), modes(model))
  p.taustar <- .Call("p_tau_reduced", model.mustar)
  pstar["tau"] <- log(p.taustar)

  model.nu0star <- .Call("reduced_nu0", model.reduced)
  identical(modes(model.nu0star), modes(model))
  p.nu0star <- .Call("p_nu0_reduced", model.nu0star)
  pstar["nu0"] <- log(mean(p.nu0star))

  model.s20star <- .Call("reduced_s20", model.reduced)
  p.s20star <- .Call("p_s20_reduced", model.s20star)
  pstar["s20"] <- log(p.s20star)

  pstar
}

marginalLikelihood <- function(model, mp){
  logLik <- completeLogLik(model) ## includes 2nd stage
  logPrior <- modes(model)[["logprior"]]
  pstar <- blockUpdates(model, mp)
  m.y <- logLik + logPrior - sum(pstar)
  m.y
}
