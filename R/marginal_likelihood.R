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

  ptheta.star <- marginal_theta(model)
  pstar["theta"] <- log(mean(ptheta.star))

  model.reduced <- model
  mcmcParams(model.reduced, force=TRUE) <- mp

  model.psigma2 <- reduced_sigma(model.reduced)
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


blockUpdatesPooledVar <- function(model, mp){
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
  ptheta.star <- full_theta_pooled(model)
  pstar["theta"] <- log(mean(ptheta.star))
  ##
  model.reduced <- model
  mcmcParams(model.reduced, force=TRUE) <- mp
  ##
  model.psigma2 <- reduced_sigma_pooled(model.reduced)
  identical(modes(model.psigma2), modes(model))
  psigma.star <- p_sigma_reduced_pooled(model.psigma2)
  pstar["sigma"] <- log(mean(psigma.star))
  ##
  model.pistar <- reduced_pi_pooled(model.reduced)
  identical(modes(model.pistar), modes(model))
  p.pi.star <- p_pmix_reduced_pooled(model.pistar)
  pstar["pi"] <- log(mean(p.pi.star))
  ##
  ## Block updates for stage 2 parameters
  ##
  model.mustar <- reduced_mu_pooled(model.reduced)
  stopifnot(identical(modes(model.mustar), modes(model)))
  p.mustar <- p_mu_reduced_pooled(model.mustar)
  pstar["mu"] <- log(mean(p.mustar))
  ##
  model.taustar <- reduced_tau_pooled(model.reduced)
  identical(modes(model.taustar), modes(model))
  p.taustar <- p_tau_reduced_pooled(model.mustar)
  pstar["tau"] <- log(p.taustar)
  ##
  model.nu0star <- reduced_nu0_pooled(model.reduced)
  identical(modes(model.nu0star), modes(model))
  p.nu0star <- p_nu0_reduced_pooled(model.nu0star)
  pstar["nu0"] <- log(mean(p.nu0star))
  ##
  model.s20star <- reduced_s20_pooled(model.reduced)
  p.s20star <- p_s20_reduced_pooled(model.s20star)
  pstar["s20"] <- log(p.s20star)
  pstar
}

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

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,MarginalModel,integer-method
setMethod("marginalLikelihood", c("MarginalModel", "integer"),
    function(model, niter) {
        mp <- McmcParams(iter=niter)
        logLik <- modes(model)[["loglik"]] ## includes 2nd stage
        model2 <- useModes(model)
        stage2.loglik <- stageTwoLogLik(model2)
        logPrior <- modes(model)[["logprior"]]
        pstar <- blockUpdates(model, mp)
        m.y <- logLik + stage2.loglik + logPrior - sum(pstar) +
               log(factorial(k(model)))
        m.y
    }
          )

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,SingleBatchPooledVar,integer-method
setMethod("marginalLikelihood", c("SingleBatchPooledVar", "integer"),
          function(model, niter) {
            mp <- McmcParams(iter=niter)
            logLik <- modes(model)[["loglik"]] ## includes 2nd stage
            model2 <- useModes(model)
            stage2.loglik <- stageTwoLogLik_pooled(model2)
            logPrior <- modes(model)[["logprior"]]
            pstar <- blockUpdatesPooledVar(model2, mp)
            m.y <- logLik + stage2.loglik + logPrior - sum(pstar) +
                log(factorial(k(model)))
            m.y
          })

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,MarginalModel-method
setMethod("marginalLikelihood", "MarginalModel",
    function(model) {
        marginalLikelihood(model, 1000L)
    }
)

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,BatchModel,integer-method
setMethod("marginalLikelihood", c("BatchModel", "integer"),
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
    }
)

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,BatchModel-method
setMethod("marginalLikelihood", "BatchModel",
    function(model) {
        marginalLikelihood(model, 1000L)
    }
)

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,list,integer-method
setMethod("marginalLikelihood", c("list", "integer"),
    function(model, niter) {
        sapply(model, marginalLikelihood, niter=niter)
    }
)

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,list-method
setMethod("marginalLikelihood", "list",
    function(model, niter) {
        sapply(model, marginalLikelihood, niter=1000L)
    }
)
