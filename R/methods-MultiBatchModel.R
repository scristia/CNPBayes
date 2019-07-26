.empty_batch_model <- function(hp, mp){
  K <- k(hp)
  B <- 0L
  S <- iter(mp)
  ch <- initialize_mcmc(K, S, B)
  obj <- new("MultiBatchModel",
             k=as.integer(K),
             hyperparams=hp,
             theta=matrix(NA, 0, K),
             sigma2=matrix(NA, 0, K),
             mu=numeric(K),
             tau2=numeric(K),
             nu.0=numeric(1),
             sigma2.0=numeric(1),
             ##pi=numeric(K),
             pi=matrix(NA, 0, K),
             predictive=numeric(K*B),
             zstar=integer(K*B),
             data=numeric(0),
             data.mean=matrix(NA, B, K),
             data.prec=matrix(NA, B, K),
             z=integer(0),
             zfreq=integer(K),
             probz=matrix(0, S, K),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=ch,
             mcmc.params=mp,
             batch=integer(0),
             batchElements=integer(0),
             label_switch=FALSE,
             marginal_lik=as.numeric(NA),
             .internal.constraint=5e-4,
             .internal.counter=0L)
  chains(obj) <- McmcChains(obj)
  obj
}

.ordered_thetas_multibatch <- function(model){
  thetas <- theta(model)
  checkOrder <- function(theta) identical(order(theta), seq_along(theta))
  is_ordered <- apply(thetas, 1, checkOrder)
  all(is_ordered)
}

reorderMultiBatch <- function(model){
  is_ordered <- .ordered_thetas_multibatch(model)
  if(is_ordered) return(model)
  ## thetas are not all ordered
  thetas <- theta(model)
  s2s <- sigma2(model)
  K <- k(model)
  ix <- order(thetas[1, ])
  B <- nBatch(model)
  . <- NULL
  tab <- tibble(z_orig=z(model),
                z=z(model),
                batch=batch(model)) %>%
    mutate(index=seq_len(nrow(.)))
  z_relabel <- NULL
  for(i in seq_len(B)){
    ix.next <- order(thetas[i, ])
    thetas[i, ] <- thetas[i, ix.next]
    s2s[i, ] <- s2s[i, ix]
    index <- which(tab$batch == i)
    tab2 <- tab[index, ] %>%
      mutate(z_relabel=factor(z, levels=ix.next)) %>%
      mutate(z_relabel=as.integer(z_relabel))
    tab$z[index] <- tab2$z_relabel
  }
  ps <- p(model)[, ix, drop=FALSE]
  mu(model) <- mu(model)[ix]
  tau2(model) <- tau2(model)[ix]
  sigma2(model) <- s2s
  theta(model) <- thetas
  p(model) <- ps
  z(model) <- tab$z
  log_lik(model) <- computeLoglik(model)
  model
}

setMethod("sortComponentLabels", "MultiBatchModel", function(model){
  reorderMultiBatch(model)
})

.MB <- function(dat=numeric(),
                hp=HyperparametersMultiBatch(),
                mp=McmcParams(iter=1000, thin=10,
                              burnin=1000, nStarts=4),
                batches=integer()){
  ## If the data is not ordered by batch,
  ## its a little harder to sort component labels
  ##dat2 <- tibble(y=dat, batch=batches)
  ub <- unique(batches)
  nbatch <- setNames(as.integer(table(batches)), ub)
  B <- length(ub)
  N <- length(dat)
  ## move to setValidity
  if(length(dat) != length(batches)) {
    stop("batch vector must be the same length as data")
  }
  K <- k(hp)
  ## mu_k is the average across batches of the thetas for component k
  ## tau_k is the sd of the batch means for component k
  mu <- sort(rnorm(k(hp), mu.0(hp), sqrt(tau2.0(hp))))
  tau2 <- 1/rgamma(k(hp), 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
  p <- rdirichlet(1, alpha(hp))[1, ]
  p <- matrix(p, B, K, byrow=TRUE)
  sim_theta <- function(mu, tau, B) sort(rnorm(B, mu, tau))
  . <- NULL
  thetas <- map2(mu, sqrt(tau2), sim_theta, B) %>%
    do.call(cbind, .) %>%
    apply(., 1, sort) %>%
    t
  if(K == 1) thetas <- t(thetas)
  nu.0 <- 3.5
  sigma2.0 <- 0.25
  sigma2s <- 1/rgamma(k(hp) * B, 0.5 * nu.0, 0.5 * nu.0 * sigma2.0) %>%
    matrix(B, k(hp))
  u <- rchisq(length(dat), hp@dfr)
  S <- iter(mp)
  ch <- initialize_mcmc(K, S, B)
  obj <- new("MultiBatchModel",
             k=as.integer(K),
             hyperparams=hp,
             theta=thetas,
             sigma2=sigma2s,
             mu=mu,
             tau2=tau2,
             nu.0=nu.0,
             sigma2.0=sigma2.0,
             pi=p,
             predictive=numeric(K*B),
             zstar=integer(K*B),
             data=dat,
             u=u,
             data.mean=matrix(NA, B, K),
             data.prec=matrix(NA, B, K),
             z=sample(seq_len(K), N, replace=TRUE),
             zfreq=integer(K),
             probz=matrix(0, N, K),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=ch,
             mcmc.params=mp,
             batch=batches,
             batchElements=nbatch,
             label_switch=FALSE,
             marginal_lik=as.numeric(NA),
             .internal.constraint=5e-4,
             .internal.counter=0L)
  obj
}

MultiBatchModel2 <- function(dat=numeric(),
                             hp=HyperparametersMultiBatch(),
                             mp=McmcParams(iter=1000, thin=10,
                                           burnin=1000, nStarts=4),
                             batches=integer()){
  if(length(dat) == 0){
    return(.empty_batch_model(hp, mp))
  }
  iter <- 0
  validZ <- FALSE
  mp.tmp <- McmcParams(iter=0, burnin=burnin(mp), thin=1, nStarts=1)
  while(!validZ){
    ##
    ## Burnin with MB model
    ##
    mb <- .MB(dat, hp, mp.tmp, batches)
    mb <- runBurnin(mb)
    tabz1 <- table(batch(mb), z(mb))
    tabz2 <- table(z(mb))
    validZ <- length(tabz2) == k(hp) && all(tabz1 > 1)
    iter <- iter + 1
    if(iter == 100) {
      message("Trouble initializing a valid model. The number of components is likely too large")
      return(NULL)
    }
  }
  mb2 <- sortComponentLabels(mb)
  mcmcParams(mb2) <- mp
  chains(mb2) <- McmcChains(mb2)
  probz(mb2)[,] <- 0
  mb2
}

MB <- MultiBatchModel2

ensureAllComponentsObserved <- function(obj){
  zz <- table(batch(obj), z(obj))
  K <- seq_len(k(obj))
  if(any(zz<=1)){
    index <- which(rowSums(zz<=1) > 0)
    for(i in seq_along(index)){
      j <- index[i]
      zup <- z(obj)[batch(obj) == j]
      zfact <- factor(zup, levels=K)
      minz <- as.integer(names(table(zfact))[table(zfact) <= 1])
      ##missingk <- K[!K %in% unique(zup)]
      maxk <- names(table(zfact))[which.max(table(zfact))]
      nreplace <- length(minz)*2
      zup[sample(which(zup == maxk), nreplace)] <- as.integer(minz)
      obj@z[batch(obj) == j] <- as.integer(zup)
    }
  }
  obj
}


setMethod("[", "MultiBatchModel", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    y(x) <- y(x)[i]
    z(x) <- z(x)[i]
    batch(x) <- batch(x)[i]
  }
  x
})

setMethod("bic", "MultiBatchModel", function(object){
  object <- useModes(object)
  ## K: number of free parameters to be estimated
  ##   - component and batch-specific parameters:  theta, sigma2  ( k(model) * nBatch(model))
  ##   - mixing probabilities: (k-1)*nBatch
  ##   - component-specific parameters: mu, tau2                 2 x k(model)
  ##   - length-one parameters: sigma2.0, nu.0                   +2
  K <- 2*k(object)*nBatch(object) + (k(object)-1) + 2*k(object) + 2
  n <- length(y(object))
  ll <- compute_loglik(object)
  bicstat <- -2*(ll + logPrior(object)) + K*(log(n) - log(2*pi))
  bicstat
})

setMethod("bic", "MultiBatchPooled", function(object){
  object <- useModes(object)
  ## K: number of free parameters to be estimated
  ##   - component and batch-specific parameters:  theta  ( k(model) * nBatch(model))
  ##   - batch specific parameter: sigma2
  ##   - mixing probabilities: (k-1)*nBatch
  ##   - component-specific parameters: mu, tau2                 2 x k(model)
  ##   - length-one parameters: sigma2.0, nu.0                   +2
  K <- k(object)*nBatch(object) + nBatch(object) + (k(object)-1) + 2*k(object) + 2
  n <- length(y(object))
  ll <- loglik_multibatch_pvar(object)
  bicstat <- -2*(ll + logPrior(object)) + K*(log(n) - log(2*pi))
  bicstat
})

setMethod(".compute_loglik", "MultiBatchModel", function(object){
  compute_loglik(object)
})

setMethod(".compute_loglik", "MultiBatchPooled", function(object){
  loglik_multibatch_pvar(object)
})


setMethod("collapseBatch", "MultiBatchModel", function(object){
  collapseBatch(y(object), as.character(batch(object)))
})


batchLik <- function(x, p, mean, sd)  p*dnorm(x, mean, sd)


setMethod("computeMeans", "MultiBatchModel", function(object){
  compute_means(object)
})


setMethod("computePrec", "MultiBatchModel", function(object){
  compute_prec(object)
})


setMethod("computePrior", "MultiBatchModel", function(object){
  compute_logprior(object)
})

.computeModesBatch <- function(object){
  mc <- chains(object)
  if(iter(mc) == 0) return(list())
  i <- argMax(object)
  B <- nBatch(object)
  K <- k(object)
  thetamax <- matrix(theta(mc)[i, ], B, K)
  sigma2max <- matrix(sigma2(mc)[i, ], B, K)
  pmax <- matrix(p(mc)[i, ], B, K)
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

setMethod("computeModes", "MultiBatchModel", function(object){
  .computeModesBatch(object)
})

componentVariances <- function(y, z)  v <- sapply(split(y, z), var)


setMethod("computeVars", "MultiBatchModel", function(object){
  compute_vars(object)
})


setMethod("mu", "MixtureModel", function(object) object@mu)


setReplaceMethod("mu", "MultiBatchModel", function(object, value){
  object@mu <- value
  object
})

nBatch <- function(object) length(uniqueBatch(object))

batchElements <- function(object) object@batchElements

setReplaceMethod("p", "MixtureModel", function(object, value){
  object@pi <- value
  object
})

setMethod("pMean", "MultiBatchModel", function(object) {
  mns <- colMeans(pic(object))
  mns
})

setMethod("showMeans", "MixtureModel", function(object){
  thetas <- round(theta(object), 2)
  mns <- c("\n", paste0(t(cbind(thetas, "\n")), collapse="\t"))
  mns <- paste0("\t", mns[2])
  mns <- paste0("\n", mns[1])
  mns
})

setMethod("showSigmas", "MixtureModel", function(object){
  sigmas <- round(sqrt(sigma2(object)), 2)
  sigmas <- c("\n", paste0(t(cbind(sigmas, "\n")), collapse="\t"))
  sigmas <- paste0("\t", sigmas[2])
  sigmas <- paste0("\n", sigmas[1])
  sigmas
})


setReplaceMethod("sigma2", "MixtureModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})

setReplaceMethod("sigma", "MixtureModel", function(object, value){
  sigma2(object) <- value^2
  object
})

setMethod("sigma2", "MixtureModel", function(object) {
  s2 <- object@sigma2
  ##s2 <- matrix(s2, nBatch(object), k(object))
  ##rownames(s2) <- uniqueBatch(object)
  s2
})

setMethod("sigma", "MixtureModel", function(object) {
  s2 <- sigma2(object)
  ##s2 <- matrix(s2, nBatch(object), k(object))
  sqrt(s2)
})

setMethod("tablez", "MultiBatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})

setMethod("sigmaMean", "MultiBatchModel", function(object) {
  mns <- colMeans(sigmac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})


setMethod("tau2", "MixtureModel", function(object) object@tau2)

setReplaceMethod("tau2", "MultiBatchModel", function(object, value){
  object@tau2 <- value
  object
})

setMethod("theta", "MixtureModel", function(object) {
  b <- object@theta
  ##b <- matrix(b, nBatch(object), k(object))
  ##rownames(b) <- uniqueBatch(object)
  b
})


setReplaceMethod("theta", "MultiBatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@theta <- value
  object
})

setMethod("thetaMean", "MultiBatchModel", function(object) {
  mns <- colMeans(thetac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

setMethod("show", "MultiBatchModel", function(object){
  ##callNextMethod()
  cls <- class(object)
  cat(paste0("An object of class ", cls), "\n")
  cat("     n. obs              :", length(y(object)), "\n")
  cat("     n. batches          :", nBatch(object), "\n")
  cat("     k                   :", k(object), "\n")
  cat("     nobs/batch          :", table(batch(object)), "\n")
  cat("     log lik (s)         :", round(log_lik(object), 1), "\n")
  cat("     log prior (s)       :", round(logPrior(object), 1), "\n")
  cat("     log marginal lik (s):", round(marginal_lik(object), 1), "\n")
})

setMethod("tablez", "MixtureModel", function(object){
  tab <- table(batch(object), z(object))
  tab <- tab[uniqueBatch(object), , drop=FALSE]
  tab
})

uniqueBatch <- function(object) unique(batch(object))

setMethod("updateMultinomialProb", "MultiBatchModel", function(object){
  update_multinomialPr(object)
})


setMethod("computeLoglik", "MultiBatchModel", function(object){
  compute_loglik(object)
})

setMethod("updateZ", "MultiBatchModel", function(object){
  update_z(object)
})

setMethod("updateObject", "MultiBatchModel", function(object){
  chains(object) <- updateObject(chains(object))
  object <- callNextMethod(object)
  object
})
