#' @include methods-MultiBatchModel.R
NULL

setValidity("MultiBatchPooled", function(object){
  msg <- TRUE
  if(length(p(object)) != k(object)){
    msg <- "Mixture probability vector must be the same length as k"
    return(msg)
  }
  if(k(object)!=k(hyperParams(object))){
    msg <- "disagreement of k in hyperparams and model"
    return(msg)
  }
  if(length(y(object))!=length(u(object))){
    msg <- "u-vector must be same length as data"
    return(msg)
  }
  if(iter(object) != iter(chains(object))){
    msg <- "number of iterations not the same between chains and model"
    return(msg)
  }
  th.len <- prod(dim(theta(object)))
  pr.len <- length(object@predictive)
  if(th.len != pr.len){
    msg <- "predictive slot in current values should have length K x B"
    return(msg)
  }
  msg
})

reorderMultiBatchPooled <- function(model){
  is_ordered <- .ordered_thetas_multibatch(model)
  if(is_ordered) return(model)
  thetas <- theta(model)
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
    index <- which(tab$batch == i)
    tab2 <- tab[index, ] %>%
      mutate(z_relabel=factor(z, levels=ix.next)) %>%
      mutate(z_relabel=as.integer(z_relabel))
    tab$z[index] <- tab2$z_relabel
  }
  ps <- p(model)[ix]
  mu(model) <- mu(model)[ix]
  tau2(model) <- tau2(model)[ix]
  theta(model) <- thetas
  p(model) <- ps
  z(model) <- tab$z
  log_lik(model) <- computeLoglik(model)
  model
}

setMethod("sortComponentLabels", "MultiBatchPooled", function(model){
  reorderMultiBatchPooled(model)
})

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
  mb <- MB(dat, hp, mp, batches)
  ## average variances across components
  mbp <- as(mb, "MultiBatchPooled")
  mbp <- sortComponentLabels(mbp)
  mcmcParams(mbp) <- mp
  log_lik(mbp) <- loglik_multibatch_pvar(mbp)
  mbp
}

#' Constructor for MultiBatchPooled model
#'
#' Initializes a MultiBatchPooled model, a container for storing data, parameters, and MCMC output for mixture models with batch- and component-specific means. The variance is assumed to be the same for all components, but allowed to differ by batch.
#'
#' @param dat the data for the simulation.
#' @param batches an integer-vector of the different batches
#' @param hp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mp An object of class 'McmcParams'
#' @return An object of class `MultiBatchPooled`
#' @export
#' @examples
#'   model <- MBP(rnorm(10), batch=rep(1L, 10))
#'   set.seed(100)
#'   nbatch <- 3
#'   k <- 3
#'   means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
#'       0, 0.05), nbatch, k, byrow = FALSE)
#'   sds <- matrix(0.15, nbatch, k)
#'   sds[, 1] <- 0.3
#'   N <- 1000
#'   truth <- simulateBatchData(N = 2500,
#'                              batch = rep(letters[1:3],
#'                              length.out = 2500),
#'                              theta = means, sds = sds,
#'                              p = c(1/5, 1/3, 1 - 1/3 - 1/5))
#'   MBP(dat=y(truth), batches=batch(truth),
#'       hp=hpList(k=3)[["MB"]])
MBP <- MultiBatchPooled

## MBP

#' @rdname sigma2-method
#' @aliases sigma2,MultiBatchPooled-method
setMethod("sigma2", "MultiBatchPooled", function(object) {
  s2 <- object@sigma2
  ##s2 <- matrix(s2, nBatch(object), k(object))
  names(s2) <- uniqueBatch(object)
  s2
})

setMethod("sigma", "MultiBatchPooled", function(object) {
  s2 <- sigma2(object)
  sqrt(s2)
})

#' @rdname sigma2-method
#' @aliases sigma2<-,MultiBatchCopyNumberPooled-method
setReplaceMethod("sigma2", "MultiBatchPooled", function(object, value){
  names(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})

setReplaceMethod("sigma", "MultiBatchPooled", function(object, value) {
  sigma2(object) <- value^2
  object
})

## MultiBatch Copy number

#' @rdname sigma2-method
#' @aliases sigma2,MultiBatchCopyNumberPooled-method
setMethod("sigma2", "MultiBatchCopyNumberPooled", function(object) {
  s2 <- object@sigma2
  names(s2) <- uniqueBatch(object)
  s2
})

#' @rdname sigma2-method
#' @aliases sigma,MultiBatchCopyNumberPooled-method
setMethod("sigma", "MultiBatchCopyNumberPooled", function(object) {
  sqrt(sigma2(object))
  s2
})

#' @rdname sigma2-method
#' @aliases sigma2,MultiBatchCopyNumberPooled-method
setReplaceMethod("sigma2", "MultiBatchCopyNumberPooled", function(object, value){
  names(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})

#' @rdname sigma2-method
#' @aliases sigma,MultiBatchCopyNumberPooled-method
setReplaceMethod("sigma", "MultiBatchCopyNumberPooled", function(object, value){
  sigma2(object) <- value^2
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
  .mu <- map(ch.list, mu) %>% do.call(rbind, .)
  .tau2 <- map(ch.list, tau2) %>% do.call(rbind, .)
  zfreq <- map(ch.list, zFreq) %>% do.call(rbind, .)
  pred <- map(ch.list, predictive) %>% do.call(rbind, .)
  zz <- map(ch.list, zstar) %>% do.call(rbind, .)
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
            predictive=pred,
            zstar=zz,
            iter=nrow(pred),
            k=ncol(.mu),
            B=length(unique(batches)))
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
               predictive=predictive(mc)[nrow(th), ],
               zstar=zstar(mc)[nrow(th), ],
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




#' @aliases sigma,MultiBatchCopyNumberPooled-method
#' @rdname sigma2-method
setMethod("sigma_", "MultiBatchCopyNumberPooled", function(object){
  s2 <- object@sigma2
  names(s2) <- uniqueBatch(object)
  sqrt(s2)
})

##setMethod("updateObject", "MultiBatchPooled", function(object){
##  chains(object) <- updateObject(chains(object),
##                                 k=k(object),
##                                 iter=iter(object),
##                                 B=nBatch(object))
##  object
##})
