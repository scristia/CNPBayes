set_param_names <- function(x, nm){
  K <- seq_len(ncol(x))
  set_colnames(x, paste0(nm, K))
}

mcmcList <- function(model.list){
  ch.list <- map(model.list, chains)
  theta.list <- map(ch.list, theta) %>%
    map(set_param_names, "theta")
  sigma.list <- map(ch.list, sigma) %>%
    map(set_param_names, "sigma")
  p.list <- map(ch.list, p) %>%
    map(set_param_names, "p")
  nu0.list <- map(ch.list, nu.0) %>%
    map(as.matrix) %>%
    map(set_param_names, "nu.0")
  s20.list <- map(ch.list, sigma2.0) %>%
    map(as.matrix) %>%
    map(set_param_names, "sigma2.0")
  mu.list <- map(ch.list, mu) %>%
    map(as.matrix) %>%
    map(set_param_names, "mu")
  tau2.list <- map(ch.list, tau2) %>%
    map(as.matrix) %>%
    map(set_param_names, "tau2")
  loglik <- map(ch.list, log_lik) %>%
    map(as.matrix) %>%
    map(set_param_names, "log_lik")
  half <- floor(nrow(theta.list[[1]])/2)
  first_half <- function(x, half){
    x[seq_len(half), , drop=FALSE]
  }
  last_half <- function(x, half){
    i <- (half + 1):(half*2)
    x <- x[i, , drop=FALSE]
    x
  }
  theta.list <- c(map(theta.list, first_half, half),
                  map(theta.list, last_half, half))
  sigma.list <- c(map(sigma.list, first_half, half),
                  map(sigma.list, last_half, half))
  p.list <- c(map(p.list, first_half, half),
              map(p.list, last_half, half))
  nu0.list <- c(map(nu0.list, first_half, half),
                map(nu0.list, last_half, half))
  s20.list <- c(map(s20.list, first_half, half),
                map(s20.list, last_half, half))
  mu.list <- c(map(mu.list, first_half, half),
               map(mu.list, last_half, half))
  tau2.list <- c(map(tau2.list, first_half, half),
                map(tau2.list, last_half, half))
  vars.list <- vector("list", length(p.list))
  for(i in seq_along(vars.list)){
    vars.list[[i]] <- cbind(theta.list[[i]],
                            sigma.list[[i]],
                            p.list[[i]],
                            nu0.list[[i]],
                            s20.list[[i]],
                            mu.list[[i]],
                            tau2.list[[i]])
  }
  vars.list <- map(vars.list, mcmc)
  ##loglik2 <- mcmc(do.call(rbind, loglik))
  mlist <- mcmc.list(vars.list)
  mlist
}

diagnostics <- function(model.list){
  mlist <- mcmcList(model.list)
  neff <- effectiveSize(mlist)
  r <- gelman.diag(mlist[, -9])  
  list(neff=neff, r=r)
}

combineModels <- function(model.list){
  ch.list <- map(model.list, chains)
  ##pz <- map(model.list, probz) %>% Reduce("+", .)
  ##pz <- pz/length(model.list)
  ##z <- map(ch.list, "z") %>%  Reduce("+", .)
  ##z <- max.col(pz)
  th <- map(ch.list, theta) %>% do.call(rbind, .)
  s2 <- map(ch.list, sigma2) %>% do.call(rbind, .)
  ll <- map(ch.list, log_lik) %>% unlist
  pp <- map(ch.list, p) %>% do.call(rbind, .)
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
            pi=pp,
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
  pm.p <- colMeans(pp)
  pm.n0 <- median(n0)
  pm.mu <- mean(.mu)
  pm.tau2 <- mean(.tau2)
  pm.s20 <- mean(s2.0)
  pz <- map(model.list, probz) %>% Reduce("+", .)
  pz <- pz/length(model.list)
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
  model <- new("MarginalModel",
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

selectModels <- function(model.list){
  ## cluster models in two groups by mean of the log likelihood
  ## discard models that cluster in a group of low log likelihoods
  ch.list <- map(model.list, chains)
  ll <- map(ch.list, log_lik) %>%
    map_dbl(mean)
  cl <- kmeans(ll, centers=2)$cluster
  mean.ll <- sort(map_dbl(split(ll, cl), mean))
  keep <- cl == names(mean.ll)[2]
  if(sum(keep) < 3){
    ## keep all models
    keep <- rep(TRUE, length(model.list))
  }
  keep
}

gelman_rubin <- function(mcmc_list, hp){
  r <- tryCatch(gelman.diag(mcmc_list), error=function(e) NULL)
  if(is.null(r)){
    ## happens because p is not positive definite
    pcolumn <- match(paste0("p", k(hp)), colnames(mcmc_list[[1]]))
    f <- function(x, pcolumn){
      x[, -pcolumn]
    }
    mcmc_list <- map(mcmc_list, f, pcolumn)
    r <- gelman.diag(mcmc_list)
  }
  r
}

gibbs <- function(mp, hp, dat, max_burnin=32000){
  nchains <- nStarts(mp)
  nStarts(mp) <- 1L ## because posteriorsimulation uses nStarts in a different way
  while(burnin(mp) < max_burnin){
    message("burnin: ", burnin(mp), ", thin: ", thin(mp))
    mod.list <- replicate(nchains,
                          MarginalModel2(data=dat,
                                         k=k(hp),
                                         mcmc.params=mp,
                                         hypp=hp))
    mod.list <- map(mod.list, posteriorSimulation)
    label_swapping <- map_lgl(mod.list, label_switch)
    if(any(label_swapping)){
      ## try one more time
      mod.list <- replicate(nchains,
                            MarginalModel2(data=dat,
                                           k=k(hp),
                                           mcmc.params=mp,
                                           hypp=hp))
      mod.list <- map(mod.list, posteriorSimulation)
      label_swapping <- map_lgl(mod.list, label_switch)
      message("Label switching detected")
      mlist <- mcmcList(mod.list)
      neff <- effectiveSize(mlist)
      r <- gelman_rubin(mlist, hp)
      break()
    }
    mod.list <- mod.list[ selectModels(mod.list) ]
    mlist <- mcmcList(mod.list)
    neff <- effectiveSize(mlist)
    r <- gelman_rubin(mlist, hp)
    message("   r: ", round(r$mpsrf, 2))
    message("   eff size (minimum): ", round(min(neff), 1))
    message("   eff size (median): ", round(median(neff), 1))
    if(all(neff > 500)){
      if(r$mpsrf < 2) break()
      ## neff is > 500 but mpsrf is big -- more burnin; new starts
    } else {
      ## neff is small-- autocorrelation or label switching
      mp@thin <- thin(mp) * 2
    }
    burnin(mp) <- burnin(mp) * 2
  }
  model <- combineModels(mod.list)
  meets_conditions <- all(neff > 500) && r$mpsrf < 2 && !label_switch(model)
  if(meets_conditions){
    ##
    ## evaluate marginal likelihood
    ##
    marginal_lik(model) <- marginalLikelihood(model)
    message("   marginal likelihood: ", round(marginal_lik(model), 2))
  }
  model
}

gibbs_multipleK <- function(mp, hp, k_range=c(1, 4), dat, max_burnin=32000){
  K <- seq(k_range[1], k_range[2])
  updateK <- function(ncomp, h) {
    k(h) <- ncomp
    h
  }
  hp.list <- map(K, updateK, hp)
  model.list <- vector("list", length(K))
  for(i in seq_along(model.list)){
    model.list[[i]] <- gibbs(mp, hp.list[[i]], dat, max_burnin)
  }
  model.list
}
