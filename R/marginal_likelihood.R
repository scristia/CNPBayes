.blockUpdates <- function(model, mp, reject.threshold, prop.threshold) {
    model.reduced <- model
    mcmcParams(model.reduced, force=TRUE) <- mp
  
    ptheta.star <- marginal_theta(model)[1:iter(mp)]
    small.theta.red <- mean(ptheta.star < reject.threshold)
    if(paramUpdates(model)[["theta"]]==0) small.theta.red <- 0
    
    if (small.theta.red >= prop.threshold) {
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
    ptheta.star <- full_theta_pooled(model)[1:iter(mp)]
    ## Doesn't make sense if theta is fixed
    small.theta.red <- mean(ptheta.star < reject.threshold)
    if(paramUpdates(model)[["theta"]]==0) small.theta.red <- 0

    if (small.theta.red >= prop.threshold) {
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
  
    ptheta.star <- marginal_theta_batch(model)[1:iter(mp)]
    small.theta.red <- mean(ptheta.star < reject.threshold)
    if(paramUpdates(model)[["theta"]]==0) small.theta.red <- 0
    
    if (small.theta.red >= prop.threshold) {
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
    pstar <- apply(reduced_gibbs, 2, function(x) log(mean(x^(root))))
}

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,MarginalModel-method marginalLikelihood,MarginalModel,ANY-method
setMethod("marginalLikelihood", "MarginalModel",
    function(model, params=list(niter=1000L,
                                root=(1/10),
                                reject.threshold=1e-50,
                                prop.threshold=0.5)) {
        # check if length of chains is shorter than gibbs niter
        if (iter(model) < params$niter) {
            stop("The number of posterior samples must be greater than or ",
                 "equal to the number of reduced gibbs samples")
        }

        if (isOverfit(model, params)) {
            warning("The model for k=", k(model), " may be overfit.",
                    " This can lead to an incorrect marginal likelihood")
            return(NA)
        }

        # get parameters from list params
        niter <- params$niter
        root <- params$root
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
        red_gibbs <- .blockUpdates(model2, mp, reject.threshold,
                                   prop.threshold)
        pstar <- blockUpdates(red_gibbs, root)

        # calculate p(x|model)
        m.y <- logLik + stage2.loglik + logPrior - sum(pstar) +
               log(factorial(k(model)))
        m.y
    })

isOverfit <- function(model, params){
  eff_size_theta <- min(effectiveSize(theta(chains(model))))
  (eff_size_theta / iter(model) < 0.05) && paramUpdates(model)[["theta"]] > 0
}

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,SingleBatchPooledVar-method marginalLikelihood,SingleBatchPooledVar,ANY-method
setMethod("marginalLikelihood", "SingleBatchPooledVar",
    function(model, params=list(niter=1000L,
                                root=(1/10),
                                reject.threshold=1e-50,
                                prop.threshold=0.5)) {
        # check if length of chains is shorter than gibbs niter
        if (iter(model) < params$niter) {
            stop("The number of posterior samples must be greater than or ",
                 "equal to the number of reduced gibbs samples")
        }

        # calculate effective size of thetas and check against threshold
        if (isOverfit(model, params)) {
            warning("The model for k=", k(model), " may be overfit.",
                    " This can lead to an incorrect marginal likelihood")
            return(NA)
        }

        # get parameters from list params
        niter <- params$niter
        root <- params$root
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
#' @aliases marginalLikelihood,BatchModel-method marginalLikelihood,BatchModel,ANY-method
setMethod("marginalLikelihood", "BatchModel",
    function(model, params=list(niter=1000L,
                                root=(1/10),
                                reject.threshold=1e-50,
                                prop.threshold=0.5)) {
        # check if length of chains is shorter than gibbs niter
        if (iter(model) < params$niter) {
            stop("The number of posterior samples must be greater than or ",
                 "equal to the number of reduced gibbs samples")
        }

        # calculate effective size of thetas and check against threshold
        if (isOverfit(model, params)) {
            warning("The model for k=", k(model), " may be overfit.",
                    " This can lead to an incorrect marginal likelihood")
            return(NA)
        }

        # get parameters from list params
        niter <- params$niter
        root <- params$root
        reject.threshold <- params$reject.threshold
        prop.threshold <- params$prop.threshold

        # calculate p(x|theta)
        logLik <- modes(model)[["loglik"]] ## includes 2nd stage
        model2 <- useModes(model)
        stage2.loglik <- stageTwoLogLikBatch(model2)

        # calculate log p(theta)
        logPrior <- modes(model)[["logprior"]]

        mp <- McmcParams(iter=niter)
        red_gibbs <- .blockUpdatesBatch(model2, mp, reject.threshold,
                                        prop.threshold)
        pstar <- blockUpdates(red_gibbs, root)

        # calculate p(x|model)
        m.y <- logLik + stage2.loglik + logPrior - sum(pstar) +
               log(factorial(k(model)))
        m.y
    }
)

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,list-method marginalLikelihood,list,ANY-method
setMethod("marginalLikelihood", "list",
    function(model, params=list(niter=1000L,
                                root=(1/10),
                                reject.threshold=1e-50,
                                prop.threshold=0.5)) {
        marg.list <- sapply(model, marginalLikelihood, params=params)
        names <- sapply(model, function(x) paste0(class(x), k(x)))
        names <- gsub("MarginalModel", "SB", names)
        names <- gsub("BatchModel", "MB", names)
        names(marg.list) <- names
        return(marg.list)
    }
)
