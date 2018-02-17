#' @include methods-MultiBatchModel.R
NULL

MultiBatchPooled <- function(dat=numeric(),
                             hp=HyperparametersMultiBatch(),
                             mp=McmcParams(iter=1000, burnin=1000,
                                           thin=10, nStarts=4),
                             batches=integer()){
  if(length(dat) == 0){
    mb <- MultiBatchModel2(dat, hp, mp, batches)
    mbp <- as(mb, "MultiBatchPooled")
    return(mbp)
  }
  iter <- 0
  validZ <- FALSE
  mp.tmp <- McmcParams(iter=0, burnin=burnin(mp), thin=1, nStarts=1)
  while(!validZ){
    ##
    ## Burnin with MB model
    ##
    mb <- MultiBatchModel2(dat, hp, mp.tmp, batches)
    mb <- runBurnin(mb)
    tabz <- table(z(mb))
    if(length(tabz) == k(hp)) validZ <- TRUE
    iter <- iter + 1
    if(iter > 50) stop("Trouble initializing valid model. Try increasing the burnin")
  }
  ## average variances across components
  mbp <- as(mb, "MultiBatchPooled")
  mbp <- sortComponentLabels(mbp)
  mcmcParams(mbp) <- mp
  log_lik(mbp) <- loglik_multibatch_pvar(mbp)
  mbp
}

#' @export
MBP <- MultiBatchPooled

#' @rdname sigma2-method
#' @aliases sigma2,MultiBatchPooled-method
setMethod("sigma2", "MultiBatchPooled", function(object) {
  s2 <- object@sigma2
  ##s2 <- matrix(s2, nBatch(object), k(object))
  names(s2) <- uniqueBatch(object)
  s2
})

#' @rdname sigma2-method
#' @aliases sigma2,MultiBatchCopyNumberPooled-method
setMethod("sigma2", "MultiBatchCopyNumberPooled", function(object) {
  s2 <- object@sigma2
  names(s2) <- uniqueBatch(object)
  s2
})

setReplaceMethod("sigma2", "MultiBatchPooled", function(object, value){
  names(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})


setReplaceMethod("sigma2", "MultiBatchCopyNumberPooled", function(object, value){
  names(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})


setMethod("sigmaMean", "MultiBatchPooled", function(object) {
  mns <- colMeans(sigmac(object))
  ##mns <- matrix(mns, nBatch(object), k(object))
  names(mns) <- uniqueBatch(object)
  mns
})

.modesMultiBatchPooled <- function(object){
  i <- argMax(object)
  mc <- chains(object)
  B <- nBatch(object)
  K <- k(object)
  thetamax <- matrix(theta(mc)[i, ], B, K)
  sigma2max <- sigma2(mc)[i, ]
  pmax <- p(mc)[i, ]
  mumax <- mu(mc)[i, ]
  tau2max <- tau2(mc)[i,]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mumax,
                tau2=tau2max,
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i],
                zfreq=zFreq(mc)[i, ],
                loglik=log_lik(mc)[i],
                logprior=logPrior(mc)[i])
  modes
}

setMethod("computeModes", "MultiBatchPooled", function(object){
  .modesMultiBatchPooled(object)
})

setMethod("computeLoglik", "MultiBatchPooled", function(object){
  loglik_multibatch_pvar(object)
})

setMethod("updateZ", "MultiBatchPooled", function(object){
  z_multibatch_pvar(object)
})

setMethod("updateZ", "MultiBatchCopyNumberPooled", function(object){
  z_multibatch_pvar(object)
})

combine_multibatch_pooled <- function(model.list, batches){
  ch.list <- map(model.list, chains)
  . <- NULL
  th <- map(ch.list, theta) %>% do.call(rbind, .)
  s2 <- map(ch.list, sigma2) %>% do.call(rbind, .)
  ll <- map(ch.list, log_lik) %>% unlist
  pp <- map(ch.list, p) %>% do.call(rbind, .)
  n0 <- map(ch.list, nu.0) %>% unlist
  s2.0 <- map(ch.list, sigma2.0) %>% unlist
  logp <- map(ch.list, logPrior) %>% unlist
  zz <- map(ch.list, z) %>% do.call(rbind, .)
  uu <- map(ch.list, u) %>% do.call(rbind, .)
  .mu <- map(ch.list, mu) %>% do.call(rbind, .)
  .tau2 <- map(ch.list, tau2) %>% do.call(rbind, .)
  zfreq <- map(ch.list, zFreq) %>% do.call(rbind, .)
  mc <- new("McmcChains",
            theta=th,
            sigma2=s2,
            pi=pp,
            mu=.mu,
            tau2=.tau2,
            nu.0=n0,
            sigma2.0=s2.0,
            zfreq=zfreq,
            logprior=logp,
            loglik=ll,
            z=zz,
            u=uu)
  hp <- hyperParams(model.list[[1]])
  mp <- mcmcParams(model.list[[1]])
  iter(mp) <- nrow(th)
  B <- length(unique(batches))
  K <- k(model.list[[1]])
  pm.th <- matrix(colMeans(th), B, K)
  pm.s2 <- colMeans(s2)
  pm.p <- colMeans(pp)
  pm.n0 <- median(n0)
  pm.mu <- colMeans(.mu)
  pm.tau2 <- colMeans(.tau2)
  pm.s20 <- mean(s2.0)
  pz <- map(model.list, probz) %>% Reduce("+", .)
  pz <- pz/length(model.list)
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
  nbatch <- as.integer(table(batch(model.list[[1]])))
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
               batch=batch(model.list[[1]]),
               batchElements=nbatch,
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


gibbs_multibatch_pooled <- function(hp, mp, dat, max_burnin=32000, batches){
  nchains <- nStarts(mp)
  nStarts(mp) <- 1L ## because posteriorsimulation uses nStarts in a different way
  if(iter(mp) < 500){
    warning("Require at least 500 Monte Carlo simulations")
    MIN_EFF <- ceiling(iter(mp) * 0.5)
  } else MIN_EFF <- 500
  while(burnin(mp) < max_burnin && thin(mp) < 100){
    message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
    mod.list <- replicate(nchains, MultiBatchPooled(dat=dat,
                                                    hp=hp,
                                                    mp=mp,
                                                    batches=batches))
    mod.list <- suppressWarnings(map(mod.list, .posteriorSimulation2))
    label_swapping <- map_lgl(mod.list, label_switch)
    nswap <- sum(label_swapping)
    if(nswap > 0){
      mp@thin <- as.integer(thin(mp) * 2)
      if(thin(mp) > 100){
        mlist <- mcmcList(mod.list)
        neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
        if(is.null(neff)) neff <- 0
        r <- tryCatch(gelman_rubin(mlist, hp), error=function(e) NULL)
        if(is.null(r)) r <- list(mpsrf=10)
        break()
      }
      message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
      mod.list2 <- replicate(nswap,
                             MultiBatchPooled(dat=dat,
                                              mp=mp,
                                              hp=hp,
                                              batches=batches))
      mod.list2 <- suppressWarnings(map(mod.list2, .posteriorSimulation2))
      mod.list[ label_swapping ] <- mod.list2
      label_swapping <- map_lgl(mod.list, label_switch)
      if(any(label_swapping)){
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
    if(all(neff > MIN_EFF) && r$mpsrf < 1.2) break()
    burnin(mp) <- as.integer(burnin(mp) * 2)
    mp@thin <- as.integer(thin(mp) * 2)
  }
  model <- combine_multibatch_pooled(mod.list, batches)
  meets_conditions <- all(neff > MIN_EFF) && r$mpsrf < 2 && !label_switch(model)
  if(meets_conditions){
    model <- compute_marginal_lik(model)
  }
  model
}

gibbsMultiBatchPooled <- function(hp,
                                  mp,
                                  k_range=c(1, 4),
                                  dat,
                                  batches,
                                  max_burnin=32000,
                                  reduce_size=TRUE){
  K <- seq(k_range[1], k_range[2])
  hp.list <- map(K, updateK, hp)
  model.list <- map(hp.list,
                    gibbs_multibatch_pooled,
                    mp=mp,
                    dat=dat,
                    batches=batches,
                    max_burnin=max_burnin)
  names(model.list) <- paste0("MBP", map_dbl(model.list, k))
  ## sort by marginal likelihood
  ##
  ## if(reduce_size) TODO:  remove z chain, keep y in one object
  ##
  ix <- order(map_dbl(model.list, marginal_lik), decreasing=TRUE)
  models <- model.list[ix]
}

gibbsPooled <- function(hp.list,
                        mp,
                        dat,
                        batches,
                        k_range=c(1, 4),
                        max_burnin=32000,
                        top=3){
  message("Fitting multi-batch models K=", min(k_range), " to K=", max(k_range))
  mb.models <- gibbsMultiBatchPooled(hp.list[["multi_batch"]],
                                     k_range=k_range,
                                     mp=mp,
                                     dat=dat,
                                     batches=batches,
                                     max_burnin=max_burnin)
  message("Fitting single-batch models K=", min(k_range), " to K=", max(k_range))
  sb.models <- gibbs_K(hp.list[["single_batch"]],
                       k_range=k_range,
                       mp=mp,
                       dat=dat,
                       max_burnin=max_burnin)
  models <- c(mb.models, sb.models)
  ml <- map_dbl(models, marginal_lik)
  ix <- head(order(ml, decreasing=TRUE), top)
  models <- models[ix]
  models
}

reorderMultiBatchPooled <- function(model){
  is_ordered <- .ordered_thetas_multibatch(model)
  if(is_ordered) return(model)
  ## thetas are not all ordered
  thetas <- theta(model)
  ##s2s <- sigma2(model)
  K <- k(model)
  ix <- order(thetas[1, ])
  B <- nBatch(model)
  zlist <- split(z(model), batch(model))
  for(i in seq_len(B)){
    ix.next <- order(thetas[i, ])
    thetas[i, ] <- thetas[i, ix.next]
    zlist[[i]] <- as.integer(factor(zlist[[i]], levels=ix.next))
  }
  zs <- unlist(zlist)
  ps <- p(model)[ix]
  ## sigmas are batch-specific not component-specific, and therefore do not need to be reordered
  mu(model) <- mu(model)[ix]
  tau2(model) <- tau2(model)[ix]
  ##sigma2(model) <- s2s
  theta(model) <- thetas
  p(model) <- ps
  z(model) <- zs
  dataMean(model) <- computeMeans(model)
  dataPrec(model) <- computePrec(model)
  log_lik(model) <- computeLoglik(model)
  model
}

setMethod("sortComponentLabels", "MultiBatchPooled", function(model){
  reorderMultiBatchPooled(model)
})


#' @aliases sigma,MultiBatchCopyNumberPooled-method
#' @rdname sigma2-method
setMethod("sigma", "MultiBatchCopyNumberPooled", function(object){
  s2 <- object@sigma2
  names(s2) <- uniqueBatch(object)
  sqrt(s2)
})
