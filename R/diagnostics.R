set_param_names <- function(x, nm){
  K <- seq_len(ncol(x))
  set_colnames(x, paste0(nm, K))
}

mcmcList <- function(model.list){
  if(!is(model.list, "list")){
    model.list <- list(model.list)
  }
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
  if(k(model.list[[1]]) == 1){
    ## there is no p chain
    dropP <- function(x) x[, -match("p1", colnames(x))]
    mlist <- map(mlist, dropP)
  }
  mlist
}


diagnostics <- function(model.list){
  mlist <- mcmcList(model.list)
  neff <- effectiveSize(mlist)
  r <- gelman_rubin(mlist, hyperParams(model.list[[1]]))
  list(neff=neff, r=r)
}

combine_batch <- function(model.list, batches){
  ch.list <- map(model.list, chains)
  th <- map(ch.list, theta) %>% do.call(rbind, .)
  s2 <- map(ch.list, sigma2) %>% do.call(rbind, .)
  ll <- map(ch.list, log_lik) %>% unlist
  pp <- map(ch.list, p) %>% do.call(rbind, .)
  n0 <- map(ch.list, nu.0) %>% unlist
  s2.0 <- map(ch.list, sigma2.0) %>% unlist
  logp <- map(ch.list, logPrior) %>% unlist
  zz <- map(ch.list, z) %>% do.call(rbind, .)
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
            z=zz)
  hp <- hyperParams(model.list[[1]])
  mp <- mcmcParams(model.list[[1]])
  iter(mp) <- nrow(th)
  B <- length(unique(batches))
  K <- k(model.list[[1]])
  pm.th <- matrix(colMeans(th), B, K)
  pm.s2 <- matrix(colMeans(s2), B, K)
  pm.p <- colMeans(pp)
  pm.n0 <- median(n0)
  pm.mu <- colMeans(.mu)
  pm.tau2 <- colMeans(.tau2)
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
  cl <- tryCatch(kmeans(ll, centers=2)$cluster, error=function(e) NULL)
  if(is.null(cl)){
    return(rep(TRUE, length(model.list)))
  }
  mean.ll <- sort(map_dbl(split(ll, cl), mean))
  keep <- cl == names(mean.ll)[2]
  if(sum(keep) < 3){
    ## keep all models
    keep <- rep(TRUE, length(model.list))
  }
  keep
}

gelman_rubin <- function(mcmc_list, hp){
  anyNA <- function(x){
    any(is.na(x))
  }
  any_nas <- map_lgl(mcmc_list, anyNA)
  mcmc_list <- mcmc_list[ !any_nas ]
  if(length(mcmc_list) < 2 ) stop("Need at least two MCMC chains")
  r <- tryCatch(gelman.diag(mcmc_list, autoburnin=FALSE), error=function(e) NULL)
  if(is.null(r)){
    ## gelman rubin can fail if p is not positive definite
    pcolumn <- match(paste0("p", k(hp)), colnames(mcmc_list[[1]]))
    f <- function(x, pcolumn){
      x[, -pcolumn]
    }
    mcmc_list <- map(mcmc_list, f, pcolumn) %>%
      as.mcmc.list
    r <- gelman.diag(mcmc_list, autoburnin=FALSE)
  }
  r
}

constructor <- function(nm, hp, mp, batches){
  x <- switch(nm,
              SingleBatchModel=SingleBatchModel2(dat=dat, hp=hp, mp=mp),
              MultiBatchModel=MultiBatchModel(dat=dat, hp=hp, mp=mp, batches=batches))
  x
}

gibbs <- function(hp, mp, dat, max_burnin=32000){
  nchains <- nStarts(mp)
  nStarts(mp) <- 1L ## because posteriorsimulation uses nStarts in a different way
  if(iter(mp) < 500){
    stop("Require at least 500 Monte Carlo simulations")
  }
  while(burnin(mp) < max_burnin){
    message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
    mod.list <- replicate(nchains, SingleBatchModel2(dat=dat,
                                                  hp=hp,
                                                  mp=mp))
    mod.list <- suppressWarnings(map(mod.list, posteriorSimulation))
    label_swapping <- map_lgl(mod.list, label_switch)
    nswap <- sum(label_swapping)
    if(nswap > 0){
      ## try one more time
      mp@thin <- as.integer(thin(mp) * 2)
      mod.list2 <- replicate(nswap,
                             SingleBatchModel2(dat=dat,
                                              mp=mp,
                                              hp=hp))
      mod.list2 <- suppressWarnings(map(mod.list2, posteriorSimulation))
      mod.list[ label_swapping ] <- mod.list2
      label_swapping <- map_lgl(mod.list, label_switch)
      if(any(label_swapping)){
        message("  Label switching detected")
        mlist <- mcmcList(mod.list)
        neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
        if(is.null(neff)) neff <- 0
        r <- gelman_rubin(mlist, hp)
        break()
      }
    }
    mod.list <- mod.list[ selectModels(mod.list) ]
    mlist <- mcmcList(mod.list)
    neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
    if(is.null(neff))  neff <- 0
    r <- gelman_rubin(mlist, hp)
    message("     r: ", round(r$mpsrf, 2))
    message("     eff size (minimum): ", round(min(neff), 1))
    message("     eff size (median): ", round(median(neff), 1))
    if(all(neff > 500) && r$mpsrf < 1.2) break()
    burnin(mp) <- as.integer(burnin(mp) * 2)
    mp@thin <- as.integer(thin(mp) * 2)
  }
  model <- combineModels(mod.list)
  meets_conditions <- all(neff > 500) && r$mpsrf < 2 && !label_switch(model)
  if(meets_conditions){
    model <- compute_marginal_lik(model)
  }
  model
}

compute_marginal_lik <- function(model){
  ##
  ## evaluate marginal likelihood. Relax default conditions
  ##
  params <- mlParams(root=1/2,
                     reject.threshold=exp(-1),
                     prop.threshold=0.5,
                     prop.effective.size=0)
  ml <- tryCatch(marginalLikelihood(model, params), warning=function(w) NULL)
  if(!is.null(ml)){
    marginal_lik(model) <- ml
    message("     marginal likelihood: ", round(marginal_lik(model), 2))
  } else {
    warning("Unable to compute marginal likelihood")
  }
  model
}

gibbs_batch <- function(hp, mp, dat, max_burnin=32000, batches){
  nchains <- nStarts(mp)
  nStarts(mp) <- 1L ## because posteriorsimulation uses nStarts in a different way
  if(iter(mp) < 500){
    stop("Require at least 500 Monte Carlo simulations")
  }
  while(burnin(mp) < max_burnin && thin(mp) < 100){
    message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
    mod.list <- replicate(nchains, MultiBatchModel2(dat=dat,
                                                    hp=hp,
                                                    mp=mp,
                                                    batches=batches))
    mod.list <- suppressWarnings(map(mod.list, posteriorSimulation))
    label_swapping <- map_lgl(mod.list, label_switch)
    nswap <- sum(label_swapping)
    if(nswap > 0){
      mp@thin <- as.integer(thin(mp) * 2)
      message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
      mod.list2 <- replicate(nswap,
                            MultiBatchModel2(dat=dat,
                                            mp=mp,
                                            hp=hp,
                                            batches=batches))
      mod.list2 <- suppressWarnings(map(mod.list2, posteriorSimulation))
      mod.list[ label_swapping ] <- mod.list2
      label_swapping <- map_lgl(mod.list, label_switch)
      if(any(label_swapping)){
        message("  Label switching detected")
        mlist <- mcmcList(mod.list)
        neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
        if(is.null(neff)) neff <- 0
        r <- gelman_rubin(mlist, hp)
        break()
      }
    }
    mod.list <- mod.list[ selectModels(mod.list) ]
    mlist <- mcmcList(mod.list)
    neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
    if(is.null(neff)) neff <- 0
    r <- gelman_rubin(mlist, hp)
    message("     r: ", round(r$mpsrf, 2))
    message("     eff size (minimum): ", round(min(neff), 1))
    message("     eff size (median): ", round(median(neff), 1))
    if(all(neff > 500) && r$mpsrf < 1.2) break()
    burnin(mp) <- as.integer(burnin(mp) * 2)
    mp@thin <- as.integer(thin(mp) * 2)
  }
  model <- combine_batch(mod.list, batches)
  meets_conditions <- all(neff > 500) && r$mpsrf < 2 && !label_switch(model)
  if(meets_conditions){
    model <- compute_marginal_lik(model)
  }
  model
}

updateK <- function(ncomp, h) {
  k(h) <- ncomp
  h
}

gibbs_K <- function(hp,
                    mp,
                    k_range=c(1, 4),
                    dat,
                    max_burnin=32000,
                    reduce_size=TRUE){
  K <- seq(k_range[1], k_range[2])
  hp.list <- map(K, updateK, hp)
  model.list <- map(hp.list,
                    gibbs,
                    mp=mp,
                    dat=dat,
                    max_burnin=max_burnin)
  names(model.list) <- paste0("SB", map_dbl(model.list, k))
  ## sort by marginal likelihood
  ix <- order(map_dbl(model.list, marginal_lik), decreasing=TRUE)
  ##
  ## if(reduce_size) TODO:  remove z chain, keep y in one object
  ##
  model.list[ix]
}

gibbs_batch_K <- function(hp,
                          mp,
                          k_range=c(1, 4),
                          dat,
                          batches,
                          max_burnin=32000,
                          reduce_size=TRUE){
  K <- seq(k_range[1], k_range[2])
  hp.list <- map(K, updateK, hp)
  model.list <- map(hp.list,
                    gibbs_batch,
                    mp=mp,
                    dat=dat,
                    batches=batches,
                    max_burnin=max_burnin)
  names(model.list) <- paste0("MB", map_dbl(model.list, k))
  ## sort by marginal likelihood
  ##
  ## if(reduce_size) TODO:  remove z chain, keep y in one object
  ##
  ix <- order(map_dbl(model.list, marginal_lik), decreasing=TRUE)
  models <- model.list[ix]
}

reload <- function(){
  svpacks(); load_all("CNPBayes"); papp()
}

reload2 <- function(){
  svpacks()
  load_all("CNPBayes")
  document("CNPBayes")
  load_all("CNPBayes")
  papp()
}

#' Evaluate both single-batch and multi-batch models with the specified range for the number of components, returning the top models sorted by marginal likelihood
#' 
#' @param hp.list a list of hyperparameters. See example.
#' @param mp a \code{McmcParams} object
#' @param dat numeric vector of CNP summary statistics (e.g., median log R ratios)
#' @param batches an integer vector of the same length as the data providing an index for the batch
#' @param k_range a length-two integer vector providing the minimum and maximum number of components
#' @param max_burnin a length-one integer vector indicating the maximum number of burnin iterations
#' @param top the number of models to return after ordering by the marginal likelihood
#' @return a list of models
#' @examples
#'
#'  set.seed(100)
#'  nbatch <- 3
#'  k <- 3
#'  means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
#'      0, 0.05), nbatch, k, byrow = FALSE)
#'  sds <- matrix(0.15, nbatch, k)
#'  sds[, 1] <- 0.3
#'  N <- 1000
#'  truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
#'                                                length.out = N),
#'                             p = c(1/10, 1/5, 1 - 0.1 - 0.2),
#'                             theta = means,
#'                             sds = sds)
#'  hp <- HyperparametersMultiBatch(k=3,
#'                             mu=-0.75,
#'                             tau2.0=0.4,
#'                             eta.0=32,
#'                             m2.0=0.5)
#'  hp.sb <- Hyperparameters(tau2.0=0.4,
#'                           mu.0=-0.75,
#'                           eta.0=32,
#'                           m2.0=0.5)
#'  hp.list <- list(single_batch=hp.sb,
#'                  multi_batch=hp)
#'  mp <- McmcParams(iter = 1000,
#'                   burnin = 1000,
#'                   nStarts = 4,
#'                   thin=10)
#' \dontrun{
#'    models <- gibbs_all(hp.list=hp.list, dat=y(truth),
#'                        batches=batch(truth),
#'                        mp=mp,
#'                        top=3)
#' }
gibbs_all <- function(hp.list,
                      mp,
                      dat,
                      batches,
                      k_range=c(1, 4),
                      max_burnin=32000,
                      top=3){
  message("Fitting multi-batch models K=", min(k_range), " to K=", max(k_range))
  mb.models <- gibbs_batch_K(hp.list[["multi_batch"]],
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
