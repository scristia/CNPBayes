.blockUpdates <- function(model, mp, reject.threshold, prop.threshold) {
    model.reduced <- model
    mcmcParams(model.reduced, force=TRUE) <- mp
  
    ptheta.star <- marginal_theta(model.reduced)
    small.theta.red <- mean(ptheta.star < reject.threshold)

    if (small.theta.red >= reject.threshold) {
        warning("The model for k=", k(model), " may be overfit.",
                " This can lead to an incorrect marginal likelihood")
        return(matrix(NA))
    }
  
    model.psigma2 <- reduced_sigma(model.reduced)
    identical(modes(model.psigma2), modes(model))
    psigma.star <- p_sigma_reduced(model.psigma2)
  
    model.pistar <- reduced_pi(model.reduced)
    identical(modes(model.pistar), modes(model))
    p.pi.star <- p_pmix_reduced(model.pistar)
  
    ##
    ## Block updates for stage 2 parameters
    ##
    model.mustar <- reduced_mu(model.reduced)
    stopifnot(identical(modes(model.mustar), modes(model)))
    p.mustar <- p_mu_reduced(model.mustar)
  
    model.taustar <- reduced_tau(model.reduced)
    identical(modes(model.taustar), modes(model))
    p.taustar <- p_tau_reduced(model.mustar)
  
    model.nu0star <- reduced_nu0(model.reduced)
    identical(modes(model.nu0star), modes(model))
    p.nu0star <- p_nu0_reduced(model.nu0star)
  
    model.s20star <- reduced_s20(model.reduced)
    p.s20star <- p_s20_reduced(model.s20star)
  
    reduced_gibbs <- cbind(ptheta.star, psigma.star, p.mustar, p.pi.star,
                           p.taustar, p.nu0star, p.s20star)
  
    colnames(reduced_gibbs) <- c("theta", "sigma", "pi", "mu",
                                 "tau", "nu0", "s20")
  
    reduced_gibbs
}

.blockUpdatesPooledVar <- function(model, mp, reject.threshold, 
                                   prop.threshold) {
    model.reduced <- model
    mcmcParams(model.reduced, force=TRUE) <- mp

    ##
    ## Block updates for stage 1 parameters
    ##
    ptheta.star <- full_theta_pooled(model)
    small.theta.red <- mean(ptheta.star < reject.threshold)

    if (small.theta.red >= reject.threshold) {
        warning("The model for k=", k(model), " may be overfit.",
                " This can lead to an incorrect marginal likelihood")
        return(matrix(NA))
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

.blockUpdatesBatch <- function(model, mp, reject.threshold, prop.threshold) {
    model.reduced <- model
    mcmcParams(model.reduced, force=TRUE) <- mp
  
    ptheta.star <- marginal_theta_batch(model)
    small.theta.red <- mean(ptheta.star < reject.threshold)

    if (small.theta.red >= reject.threshold) {
        warning("The model for k=", k(model), " may be overfit.",
                " This can lead to an incorrect marginal likelihood")
        return(matrix(NA))
    }
  
    model.psigma2 <- reduced_sigma_batch(model.reduced)
    identical(modes(model.psigma2), modes(model))
    psigma.star <- p_sigma_reduced_batch(model.psigma2)
  
    model.pistar <- reduced_pi_batch(model.reduced)
    identical(modes(model.pistar), modes(model))
    p.pi.star <- p_pmix_reduced_batch(model.pistar)
  
    ##
    ## Block updates for stage 2 parameters
    ##
    model.mustar <- reduced_mu_batch(model.reduced)
    stopifnot(identical(modes(model.mustar), modes(model)))
    p.mustar <- p_mu_reduced_batch(model.mustar)
  
    model.taustar <- reduced_tau_batch(model.reduced)
    identical(modes(model.taustar), modes(model))
    p.taustar <- p_tau_reduced_batch(model.mustar)
  
    model.nu0star <- reduced_nu0_batch(model.reduced)
    identical(modes(model.nu0star), modes(model))
    p.nu0star <- p_nu0_reduced_batch(model.nu0star)
  
    model.s20star <- reduced_s20_batch(model.reduced)
    p.s20star <- p_s20_reduced_batch(model.s20star)

    reduced_gibbs <- cbind(ptheta.star, psigma.star, p.mustar, p.pi.star,
                           p.taustar, p.nu0star, p.s20star)
  
    colnames(reduced_gibbs) <- c("theta", "sigma", "pi", "mu",
                                 "tau", "nu0", "s20")
  
    reduced_gibbs
}

blockUpdates <- function(reduced_gibbs, root) {
    pstar <- apply(reduced_gibbs, 2, function(x) log(mean(x^(root)))),
}

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,MarginalModel-method marginalLikelihood,MarginalModel,ANY-method
setMethod("marginalLikelihood", "MarginalModel",
    function(model, params=list(niter=1000L,
                                root=(1/10),
                                reject.threshold=1e-50,
                                prop.threshold=0.5)) {
        # calculate effective size of thetas and check against threshold
        eff_size_theta <- min(effectiveSize(theta(chains(model))))
        if (eff_size_theta / iter(model) < 0.05) {
            warning("The model for k=", k(model), " may be overfit.",
                    " This can lead to an incorrect marginal likelihood")
            return(NA)
        }

        # get parameters from list params
        niter <- params$niter
        root <- params$niter
        reject.threshold <- params$reject.threshold
        prop.threshold <- params$prop.threshold

        # calculate p(x|theta)
        logLik <- modes(model)[["loglik"]] ## includes 2nd stage
        model2 <- useModes(model)
        stage2.loglik <- stageTwoLogLik(model2)

        # calculate log p(theta)
        logPrior <- modes(model)[["logprior"]]

        # calculate log p(theta|x)
        mp <- McmcParams(iter=niter)
        red_gibbs <- .blockUpdates(model, mp, reject.threshold,
                                   prop.threshold)
        pstar <- blockUpdates(red_gibbs, root)

        # calculate p(x|model)
        m.y <- logLik + stage2.loglik + logPrior - sum(pstar) +
               log(factorial(k(model)))
        m.y
})

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,SingleBatchPooledVar-method marginalLikelihood,SingleBatchPooledVar,ANY-method
setMethod("marginalLikelihood", "SingleBatchPooledVar",
    function(model, params=list(niter=1000L,
                                root=(1/10),
                                reject.threshold=1e-50,
                                prop.threshold=0.5)) {
        # calculate effective size of thetas and check against threshold
        eff_size_theta <- min(effectiveSize(theta(chains(model))))
        if (eff_size_theta / iter(model) < 0.05) {
            warning("The model for k=", k(model), " may be overfit.",
                    " This can lead to an incorrect marginal likelihood")
            return(NA)
        }

        # get parameters from list params
        niter <- params$niter
        root <- params$niter
        reject.threshold <- params$reject.threshold
        prop.threshold <- params$prop.threshold

        # calculate p(x|theta)
        logLik <- modes(model)[["loglik"]] ## includes 2nd stage
        model2 <- useModes(model)
        stage2.loglik <- stageTwoLogLik_pooled(model2)

        # calculate log p(theta)
        logPrior <- modes(model)[["logprior"]]

        # calculate log p(theta|x)
        mp <- McmcParams(iter=niter)
        red_gibbs <- .blockUpdatesPooledVar(model2, mp, reject.threshold,
                                            prop.threshold)
        pstar <- blockUpdates(red_gibbs, root)

        # calculate p(x|model)
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
        eff_size_theta <- min(effectiveSize(theta(chains(model))))
        if (eff_size_theta / iter(model) < 0.05) {
            warning("The model for k=", k(model), " may be overfit.",
                    " This can lead to an incorrect marginal likelihood")
            return(NA)
        }
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
        marg.list <- sapply(model, marginalLikelihood, niter=niter)
        names <- sapply(model, function(x) paste0(class(x), k(x)))
        names <- gsub("MarginalModel", "SB", names)
        names <- gsub("BatchModel", "MB", names)
        names(marg.list) <- names
        return(marg.list)
    }
)

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,list-method
setMethod("marginalLikelihood", "list",
    function(model, niter) {
        marg.list <- sapply(model, marginalLikelihood, niter=1000L)
        names <- sapply(model, function(x) paste0(class(x), k(x)))
        names <- gsub("MarginalModel", "SB", names)
        names <- gsub("BatchModel", "MB", names)
        names(marg.list) <- names
        return(marg.list)
    }
)
