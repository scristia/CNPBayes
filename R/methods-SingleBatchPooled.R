#' @include methods-SingleBatchModel.R
NULL

setMethod("computeLoglik", "SingleBatchPooled", function(object){
  loglik_pooled(object)
})

## ad-hoc constructor
##.SingleBatchPooled <- function(dat=numeric(),
##                              hp=Hyperparameters(),
##                              mp=McmcParams(iter=1000, burnin=1000,
##                                            thin=10, nStarts=4)){
##  model <- .SingleBatchPooled2(dat=dat, hp=hp, mp=mp)
##  model
##}


.SingleBatchPooled <- function(dat=numeric(),
                               hp=Hyperparameters(),
                               mp=McmcParams(iter=1000, burnin=1000,
                                             thin=10, nStarts=4)){
  if(length(dat) == 0){
    sb <- .empty_singlebatch_model(hp, mp)
    sbp <- as(sb, "SingleBatchPooled")
    return(sbp)
  }
  K <- k(hp)
  mu <- rnorm(1, median(dat), sd(dat))
  tau2 <- 1/rgamma(1, 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
  p <- rdirichlet(1, alpha(hp))[1, ]
  theta <- sort(rnorm(k(hp), mu, sqrt(tau2)))
  ##nu.0 <- 3.5
  ##sigma2.0 <- 0.25
  ##sigma2 <- 1/rgamma(k(hp), 0.5 * nu.0, 0.5 * nu.0 * sigma2.0)
  ##sigma2 <- abs(rnorm(1, 0.1, 0.1))^2
  nu.0 <- rgeom(1, betas(hp))
  sigma2.0 <- rgamma(1, 0.5*a(hp), 0.5*a(hp)*b(hp))
  sigma2 <- 1/rgamma(1, 0.5 * nu.0, 0.5 * nu.0 * sigma2.0)
  u <- rchisq(length(dat), hp@dfr)
  object <- new("SingleBatchPooled",
                k=as.integer(K),
                hyperparams=hp,
                theta=theta,
                sigma2=sigma2,
                mu=mu,
                tau2=tau2,
                nu.0=nu.0,
                sigma2.0=sigma2.0,
                pi=p,
                data=dat,
                u=dat,
                data.mean=numeric(K),
                data.prec=numeric(K),
                z=integer(length(dat)),
                zfreq=integer(K),
                probz=matrix(0, length(dat), K),
                logprior=numeric(1),
                loglik=numeric(1),
                mcmc.chains=McmcChains(),
                batch=rep(1L, length(dat)),
                batchElements=1L,
                modes=list(),
                mcmc.params=mp,
                label_switch=FALSE,
                marginal_lik=as.numeric(NA),
                .internal.constraint=5e-4,
                .internal.counter=0L)
  chains(object) <- McmcChains(object)
  object
}

.computeModesSBP <- function(object){
  i <- argMax(object)
  if(length(i) == 0) i <- iter(object)
  mc <- chains(object)
  thetamax <- theta(mc)[i, ]
  sigma2max <- sigma2(mc)[i, ]
  pmax <- p(mc)[i, ]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mu(mc)[i],
                tau2=tau2(mc)[i],
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i],
                zfreq=zFreq(mc)[i, ],
                loglik=log_lik(mc)[i],
                logprior=logPrior(mc)[i])
  modes
}

setMethod("computeModes", "SingleBatchPooled", function(object){
  .computeModesSBP(object)
})

## We initialize the z's to be all zeros, then z is updated as the first step of
## the mcmc. However, if the first update of z results in some components not
## being observed, then z will not be updated and will stay at zero. This
## creates NaNs for the thetas and several other parameters. To try to
## circumvent this issue, we have a while loop that simulates a model and runs 5
## iterations burnin. If the parameter values are valid, we stop. If not, we
## simulate a new model.
#' @export
#' @rdname SingleBatchModel2
SingleBatchPooled <- function(dat=numeric(),
                              hp=Hyperparameters(),
                              mp=McmcParams(iter=1000, burnin=1000,
                                            thin=10, nStarts=4)){
  if(length(dat) == 0){
    return(.SingleBatchPooled(dat, hp, mp))
  }
  iter <- 0
  validZ <- FALSE
  ##
  ## Burnin with the more flexible SingleBatch model to obtain starting values,
  ## then convert to SingleBatchPooled class
  ##
  mp.tmp <- McmcParams(iter=0, burnin=burnin(mp), thin=1, nStarts=1)
  while(!validZ){
    ##
    ## Burnin with SB model
    ##
    sb <- .SingleBatchModel2(dat, hp, mp.tmp)
    sb <- runBurnin(sb)
    tabz <- table(z(sb))
    if(length(tabz) == k(hp)) validZ <- TRUE
    iter <- iter + 1
    if(iter > 50) stop("Trouble initializing valid model. Try increasing the burnin")
  }
  sbp <- as(sb, "SingleBatchPooled")
  mcmcParams(sbp) <- mp
  sbp <- sortComponentLabels(sbp)
  log_lik(sbp) <- loglik_pooled(sbp)
  sbp
}

#' @export
SBP <- SingleBatchPooled

setValidity("SingleBatchPooled", function(object){
  s2 <- sigma2(object)
  if(length(s2) != 1){
    return("sigma2 slot should be length-one numeric vector")
  }
  TRUE
})

setMethod("updateZ", "SingleBatchPooled", function(object){
  z_pooled(object)
})

combine_singlebatch_pooled <- function(model.list, batches){
  ch.list <- map(model.list, chains)
  . <- NULL
  fun <- function(ch) ch@pi
  prob <- map(ch.list, fun) %>% do.call(rbind, .)
  th <- map(ch.list, theta) %>% do.call(rbind, .)
  s2 <- map(ch.list, sigma2) %>% do.call(rbind, .)
  ll <- map(ch.list, log_lik) %>% unlist
  ##pp <- map(ch.list, prob) %>% do.call(rbind, .)
  n0 <- map(ch.list, nu.0) %>% unlist
  s2.0 <- map(ch.list, sigma2.0) %>% unlist
  logp <- map(ch.list, logPrior) %>% unlist
  zz <- map(ch.list, z) %>% do.call(rbind, .)
  .mu <- map(ch.list, mu) %>% unlist
  .tau2 <- map(ch.list, tau2) %>% unlist
  zfreq <- map(ch.list, zFreq) %>% do.call(rbind, .)
  mc <- new("McmcChains",
            theta=th,
            sigma2=s2,
            pi=prob,
            mu=.mu,
            tau2=.tau2,
            nu.0=n0,
            sigma2.0=s2.0,
            zfreq=zfreq,
            logprior=logp,
            loglik=ll,
            z=zz)
  hp <- hyperParams(model.list[[1]])
  mp <- mcmcParams(model.list[[1]])
  iter(mp) <- nrow(th)
  pm.th <- colMeans(th)
  pm.s2 <- colMeans(s2)
  pm.p <- colMeans(prob)
  pm.n0 <- median(n0)
  pm.mu <- mean(.mu)
  pm.tau2 <- mean(.tau2)
  pm.s20 <- mean(s2.0)
  pm.ll <- mean(ll)
  pz <- map(model.list, probz) %>% Reduce("+", .)
  pz <- pz/length(model.list)
  ## the accessor will divide by number of iterations - 1
  pz <- pz * (iter(mp) - 1)
  zz <- max.col(pz)
  yy <- y(model.list[[1]])
  y_mns <- as.numeric(tapply(yy, zz, mean))
  y_prec <- as.numeric(1/tapply(yy, zz, var))
  zfreq <- as.integer(table(zz))
  any_label_swap <- any(map_lgl(model.list, label_switch))
  ## use mean marginal likelihood in combined model,
  ## or NA if marginal likelihood has not been estimated
  ml <- map_dbl(model.list, marginal_lik)
  if(all(is.na(ml))) {
    ml <- as.numeric(NA)
  } else ml <- mean(ml, na.rm=TRUE)
  model <- new(class(model.list[[1]]),
               k=k(hp),
               hyperparams=hp,
               theta=pm.th,
               sigma2=pm.s2,
               mu=pm.mu,
               tau2=pm.tau2,
               nu.0=pm.n0,
               sigma2.0=pm.s20,
               pi=pm.p,
               data=y(model.list[[1]]),
               u=u(model.list[[1]]),
               data.mean=y_mns,
               data.prec=y_prec,
               z=zz,
               zfreq=zfreq,
               probz=pz,
               logprior=numeric(1),
               loglik=numeric(1),
               mcmc.chains=mc,
               batch=rep(1L, length(yy)),
               batchElements=1L,
               modes=list(),
               mcmc.params=mp,
               label_switch=any_label_swap,
               marginal_lik=ml,
               .internal.constraint=5e-4,
               .internal.counter=0L)
  modes(model) <- computeModes(model)
  log_lik(model) <- computeLoglik(model)
  logPrior(model) <- computePrior(model)
  model
}

finiteLoglik <- function(model){
  is.finite(log_lik(model))
}


gibbs_singlebatch_pooled <- function(hp, mp, dat, max_burnin=32000){
  nchains <- nStarts(mp)
  nStarts(mp) <- 1L ## because posteriorsimulation uses nStarts in a different way
  if(iter(mp) < 500){
    stop("Require at least 500 Monte Carlo simulations")
  }
  while(burnin(mp) < max_burnin && thin(mp) < 100){
    message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
    mod.list <- replicate(nchains, SingleBatchPooled(dat=dat,
                                                     hp=hp,
                                                     mp=mp))
    mod.list <- suppressWarnings(map(mod.list, .posteriorSimulation2))
    label_swapping <- map_lgl(mod.list, label_switch)
    finite_loglik <- map_lgl(mod.list, function(m) is.finite(log_lik(m)))
    nswap <- sum(label_swapping | !finite_loglik)
    if(nswap > 0){
      index <- label_swapping | !finite_loglik
      mp@thin <- as.integer(thin(mp) * 2)
      message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
      mod.list2 <- replicate(nswap,
                             SingleBatchPooled(dat=dat,
                                               mp=mp,
                                               hp=hp))
      mod.list2 <- suppressWarnings(map(mod.list2, .posteriorSimulation2))
      mod.list[ index ] <- mod.list2
      label_swapping <- map_lgl(mod.list, label_switch)
      finite_loglik <- map_lgl(mod.list, function(m) is.finite(log_lik(m)))
      if(any(label_swapping | !finite_loglik)){
        message("  Label switching detected")
        mlist <- mcmcList(mod.list)
        neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
        if(is.null(neff)) neff <- 0
        r <- tryCatch(gelman_rubin(mlist, hp), error=function(e) NULL)
        if(is.null(r)) r <- list(mpsrf=10)
        break()
      }
    }
    mod.list <- mod.list[ selectModels(mod.list) ]
    mlist <- mcmcList(mod.list)
    neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
    if(is.null(neff)) neff <- 0
    r <- tryCatch(gelman_rubin(mlist, hp), error=function(e) NULL)
    if(is.null(r)) r <- list(mpsrf=10)
    message("     r: ", round(r$mpsrf, 2))
    message("     eff size (minimum): ", round(min(neff), 1))
    message("     eff size (median): ", round(median(neff), 1))
    if(all(neff > 500) && r$mpsrf < 1.2) break()
    burnin(mp) <- as.integer(burnin(mp) * 2)
    mp@thin <- as.integer(thin(mp) * 2)
  }
  model <- combine_singlebatch_pooled(mod.list)
  meets_conditions <- all(neff > 500) && r$mpsrf < 2 && !label_switch(model)
  if(meets_conditions){
    model <- compute_marginal_lik(model)
  }
  model
}


gibbsSingleBatchPooled <- function(hp,
                                   mp,
                                   k_range=c(1, 4),
                                   dat,
                                   max_burnin=32000,
                                   reduce_size=TRUE){
  K <- seq(k_range[1], k_range[2])
  hp.list <- map(K, updateK, hp)
  model.list <- map(hp.list,
                    gibbs_singlebatch_pooled,
                    mp=mp,
                    dat=dat,
                    max_burnin=max_burnin)
  names(model.list) <- paste0("SBP", map_dbl(model.list, k))
  ## sort by marginal likelihood
  ##
  ## if(reduce_size) TODO:  remove z chain, keep y in one object
  ##
  ix <- order(map_dbl(model.list, marginal_lik), decreasing=TRUE)
  models <- model.list[ix]
}

reorderPooledVar <- function(model){
  thetas <- theta(model)
  K <- k(model)
  ix <- order(thetas)
  if(identical(ix, seq_len(K))) return(model)
  thetas <- thetas[ix]
  zs <- as.integer(factor(z(model), levels=ix))
  ps <- p(model)[ix]
  theta(model) <- thetas
  p(model) <- ps
  z(model) <- zs
  dataPrec(model) <- 1/computeVars(model)
  dataMean(model) <- computeMeans(model)
  model
}

setMethod("sortComponentLabels", "SingleBatchPooled", function(model){
  reorderPooledVar(model)
})
