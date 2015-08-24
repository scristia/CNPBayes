blockUpdates <- function(model, mp){
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
  psigma.star <- p_sigma_reduced(model.psigma2)
  pstar["sigma"] <- log(mean(psigma.star))


  model.pistar <- reduced_pi(model.reduced)
  identical(modes(model.pistar), modes(model))
  p.pi.star <- p_pmix_reduced(model.pistar)
  pstar["pi"] <- log(mean(p.pi.star))
  ##
  ## Block updates for stage 2 parameters
  ##
  model.mustar <- reduced_mu(model.reduced)
  stopifnot(identical(modes(model.mustar), modes(model)))
  p.mustar <- p_mu_reduced(model.mustar)
  pstar["mu"] <- log(mean(p.mustar))

  model.taustar <- reduced_tau(model.reduced)
  identical(modes(model.taustar), modes(model))
  p.taustar <- p_tau_reduced(model.mustar)
  pstar["tau"] <- log(p.taustar)

  model.nu0star <- reduced_nu0(model.reduced)
  identical(modes(model.nu0star), modes(model))
  p.nu0star <- p_nu0_reduced(model.nu0star)
  pstar["nu0"] <- log(mean(p.nu0star))

  model.s20star <- reduced_s20(model.reduced)
  p.s20star <- p_s20_reduced(model.s20star)
  pstar["s20"] <- log(p.s20star)

  pstar
}
