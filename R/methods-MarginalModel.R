#' Create an object for running marginal MCMC simulations.
#' @param data the data for the simulation.
#' @param k An integer value specifying the number of latent classes.
#' @param hypp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mcmc.params An object of class 'McmcParams'
#' @export
MarginalModel <- function(data=numeric(), k=2, hypp, mcmc.params){
  batch <- rep(1L, length(data))
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  if(missing(hypp)) hypp <- HyperparametersMarginal(k=k)
  nbatch <- setNames(as.integer(table(batch)), levels(batch))
  zz <- sample(seq_len(k), length(data), replace=TRUE)
  zfreq <- as.integer(table(zz))
  object <- new("MarginalModel",
                k=as.integer(k),
                hyperparams=hypp,
                theta=numeric(k),
                sigma2=numeric(k),
                mu=numeric(1),
                tau2=numeric(1),
                nu.0=1L,
                sigma2.0=1L,
                pi=rep(1/k, k),
                data=data,
                data.mean=numeric(k),
                data.prec=numeric(k),
                z=zz,
                zfreq=zfreq,
                probz=matrix(0, length(data), k),
                logprior=numeric(1),
                loglik=numeric(1),
                mcmc.chains=McmcChains(),
                batch=batch,
                batchElements=nbatch,
                hwe=numeric(),
                modes=list(),
                mcmc.params=mcmc.params)
  object <- startingValues(object)
}

getK <- function(object){
  hypp <- hyperParams(object)
  .Call("getK", hypp)
}

#' @rdname mu-method
#' @aliases mu,MarginalModel-method
#' @export
setMethod("mu", "MarginalModel", function(object) object@mu)

#' @rdname tau2-method
#' @aliases tau2,MarginalModel-method
setMethod("tau2", "MarginalModel", function(object) object@tau2)

## compute p(theta) * p(y | theta)
setMethod("computePotential", "MarginalModel", function(object){
  ll <- computeLogLikxPrior(object)
  ll.phi <- .loglikPhiMarginal(object)
  ll+ll.phi
})

setMethod("initializeSigma2.0", "MarginalModel", function(object){
  hypp <- hyperParams(object)
  sum(alpha(hypp)*sigma2(object))/sum(alpha(hypp))
})


## just choose a big number
setMethod("initializeTau2", "MarginalModel", function(object)  1000)

setMethod("posteriorMultinomial", "MarginalModel", function(object){
  .Call("update_multinomialPr", object)
})

setMethod("show", "MarginalModel", function(object) callNextMethod())

setMethod("computeMeans", "MarginalModel", function(object){
  .Call("compute_means", object)
})

setMethod("computeVars", "MarginalModel", function(object){
  .Call("compute_vars", object)
})

setMethod("simulateY", "MarginalModel", function(object){
  zz <- z(object)
  yy <- rnorm(length(zz), mean=theta(object)[zz], sd=sigma(object)[zz])
})

setMethod("updateThetaCpp", "MarginalModel", function(object, constrain) {
  .Call("update_theta", object, constrain=constrain)
})

setMethod("updateTheta", "MarginalModel", function(object) {
  .Call("update_theta", object)
})

setMethod("updateSigma2", "MarginalModel", function(object) {
  .Call("update_sigma2", object)
})

setMethod("updateSigma2.0", "MarginalModel", function(object){
  .Call("update_sigma2_0", object)
})

setReplaceMethod("tau2", "MarginalModel", function(object, value){
  object@tau2 <- value
  object
})

setReplaceMethod("mu", "MarginalModel", function(object, value){
  object@mu <- value
  object
})

##
## For the marginal model, mu has already been initialized in the hyperparameters
##
setMethod("initializeMu", "MarginalModel", function(object)   mu(object))

#' @rdname bic-method
#' @aliases bic,MarginalModel-method
setMethod("bic", "MarginalModel", function(object){
  object <- useModes(object)
  ## K: number of free parameters to be estimated
  ##   - component-specific parameters:  theta, sigma2   (3 x k(model))
  ##   - mixing probabilities:  k-1
  ##   - length-one parameters: mu, tau2, sigma2.0, nu.0             +4
  K <- 2*k(object) + (k(object)-1) + 4
  n <- length(y(object))
  -2*(log_lik(object) + logPrior(object)) + K*(log(n) - log(2*pi))
})

setMethod("initializeTheta", "MarginalModel", function(object){
  initializeTheta(k(object))
})

#' @rdname theta-method
#' @aliases theta,MarginalModel-method
setMethod("theta", "MarginalModel", function(object) object@theta)

#' @rdname sigma2-method
#' @aliases sigma2,MarginalModel-method
setMethod("sigma2", "MarginalModel", function(object) object@sigma2)

setMethod("reorderComponents", "MarginalModel", function(object){
  ##
  ## First, update the model so that the components are ordered by theta
  ##
  ix <- order(theta(object))
  theta(object) <- sort(theta(object))
  sigma2(object) <- sigma2(object)[ix]
  p(object) <- p(object)[ix]
  zz <- z(object)
  ##
  ## Set the labels of the latent variable such that 1=first
  ## components, 2= second component, ...
  ##
  zz <- factor(as.integer(factor(zz, levels=ix)), levels=seq_len(k(object)))
  z(object) <- zz
  dataMean(object) <- dataMean(object)[ix]
  dataPrec(object) <- dataPrec(object)[ix]
  object
})

newMarginalModel <- function(object){
  mp <- mcmcParams(object)
  object2 <- MarginalModel(y(object), k=k(object), mcmc.params=mp,
                           hypp=hyperParams(object))
  theta(object2) <- theta(object)
  sigma2(object2) <- sigma2(object)
  p(object2) <- p(object)
  z(object2) <- z(object)
  nu.0(object2) <- nu.0(object)
  mu(object2) <- mu(object)
  tau2(object2) <- tau2(object)
  zFreq(object2) <- zFreq(object)
  probz(object2) <- probz(object)
  sigma2.0(object2) <- sigma2.0(object)
  dataMean(object2) <- dataMean(object)
  dataPrec(object2) <- dataPrec(object)
  log_lik(object2) <- log_lik(object)
  logPrior(object2) <- logPrior(object)
  modes(object2) <- modes(object)
  object2
}

newBatchModel <- function(object){
  mp <- mcmcParams(object)
  object2 <- BatchModel(y(object), batch=batch(object),
                        k=k(object), mcmc.params=mp,
                        hypp=hyperParams(object))
  theta(object2) <- theta(object)
  sigma2(object2) <- sigma2(object)
  p(object2) <- p(object)
  z(object2) <- z(object)
  nu.0(object2) <- nu.0(object)
  mu(object2) <- mu(object)
  tau2(object2) <- tau2(object)
  zFreq(object2) <- zFreq(object)
  probz(object2) <- probz(object)
  sigma2.0(object2) <- sigma2.0(object)
  dataMean(object2) <- dataMean(object)
  dataPrec(object2) <- dataPrec(object)
  log_lik(object2) <- log_lik(object)
  logPrior(object2) <- logPrior(object)
  modes(object2) <- modes(object)
  object2
}

setMethod("relabel", "MarginalModel", function(object, zindex){
  object <- newMarginalModel(object)
  if(identical(zindex, seq_len(k(object)))) return(object)
  ##
  ## Permute the labels for the components
  ##
  zz <- factor(z(object), levels=zindex)
  zz <- as.integer(zz)
  z(object) <- zz
  zFreq(object) <- as.integer(table(zz))
  dataMean(object) <- dataMean(object)[zindex]
  dataPrec(object) <- dataPrec(object)[zindex]
  object
})

setMethod("relabel", "BatchModel", function(object, zindex){
  object <- newBatchModel(object)
  if(identical(zindex, seq_len(k(object)))) return(object)
  ##
  ## Permute only the latent variables
  ##
  zz <- factor(z(object), levels=zindex)
  zz <- as.integer(zz)
  z(object) <- zz
  zFreq(object) <- as.integer(table(zz))
  dataMean(object) <- dataMean(object)[, zindex, drop=FALSE]
  dataPrec(object) <- dataPrec(object)[, zindex, drop=FALSE]
  object
})

setMethod("sort", "MarginalModel", function(x, decreasing=FALSE, ...){
  mc <- chains(x)
  pot <- logpotential(mc)
  index <- which.max(pot)
  thetas <- theta(mc)[index, ]
  if(identical(thetas, sort(thetas))){
    ## nothing to do
    return(x)
  }
  cn <- order(thetas)
  theta(mc) <- theta(mc)[, cn]
  theta(x) <- theta(x)[cn]

  sigma2(mc) <- sigma2(mc)[, cn]
  sigma2(x) <- sigma2(x)[cn]

  p(mc) <- p(mc)[, cn]
  p(x) <- p(x)[cn]

  mu(x) <- mu(x)[cn]
  tau2(x) <- tau2(x)[cn]

  probz(x) <- probz(x)[, cn]

  zz <- as.integer(z(x))
  z(x) <- factor(as.integer(factor(zz, levels=cn)), levels=sort(unique(zz)))
  dataMean(x) <- dataMean(x)[cn]
  dataPrec(x) <- dataPrec(x)[cn]
  chains(x) <- mc
  x
})

.computeModesMarginal <- function(object){
  i <- argMax(object)
  mc <- chains(object)
  thetamax <- theta(mc)[i, ]
  sigma2max <- sigma2(mc)[i,]
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

setMethod("computeModes", "MarginalModel", function(object){
  .computeModesMarginal(object)
})

sumSquares <- function(object){
  B <- batch(object)
  thetas <- theta(object)
  yy <- y(object)
  K <- k(object)
  zz <- z(object)
  ss <- matrix(NA, nBatch(object), k(object))
  rownames(ss) <- uniqueBatch(object)
  batch.index <- split(seq_len(length(yy)), B)
  zz <- z(object)
  for(b in uniqueBatch(object)){
    k <- batch.index[[b]]
    y <- yy[k]
    cn <- zz[k]
    m <- thetas[b, ]
    ## This could be tricky in C.  It works in R because of the factor to an integer:
    ##  as.integer(factor(c(1, 3), levels=c("1", "2", "3"))) evaluates to 1,3
    m <- m[as.integer(cn)]
    squares <- (y - m)^2
    ss[b, ] <- sapply(split(squares, factor(cn, levels=seq_len(K))), sum)
  }
  ss
}

setMethod("showMeans", "MarginalModel", function(object){
  paste(round(theta(object), 2), collapse=", ")
})

setMethod("showSigmas", "MarginalModel", function(object){
  paste(round(sqrt(sigma2(object)), 2), collapse=", ")
})

setMethod("tablez", "MarginalModel", function(object) table(z(object)))

permuteZ <- function(object, ix){
  zz <- zChain(object)
  zz2 <- zz
  for(i in seq_along(ix)){
    zz2[zz == ix[i]] <- i
  }
  zz2
}


setMethod("permuteModes", "MarginalModel", function(object, ix){
  modes(object)[["theta"]] <- modes(object)[["theta"]][ix]
  modes(object)[["sigma2"]] <- modes(object)[["sigma2"]][ix]
  modes(object)[["mixprob"]] <- modes(object)[["mixprob"]][ix]
  object
})

permuteZandModes <- function(kmod, ix){
  zChain(kmod) <- permuteZ(kmod, ix)
  kmod <- permuteModes(kmod, ix)
  kmod
}

setMethod("pTheta", "MarginalModel", function(object){
  exp(.Call("marginal_theta", object))
})

setMethod("pTheta", "BatchModel", function(object){
  .Call("p_theta_batch", object)
})

setMethod("pTheta_Zfixed", "MarginalModel", function(object){
  exp(.Call("p_theta_zpermuted", object))  ## use permuted
})

setMethod("pTheta_Zfixed", "BatchModel", function(object){
  .Call("p_theta_zfixed_batch", object)  ## use permuted
})

setMethod("reducedGibbsZThetaFixed", "MarginalModel", function(object){
  .Call("permutedz_reduced1", object)
})

setMethod("reducedGibbsZThetaFixed", "BatchModel", function(object){
  .Call("reduced_z_theta_fixed", object)
})

setMethod("pSigma2", "MarginalModel", function(object) {
  exp(.Call("p_sigma2_zpermuted", object))
})

setMethod("pSigma2", "BatchModel", function(object) {
  .Call("p_sigma2_batch", object)
})

## same for marginal and batch models

setMethod("pMixProb", "MixtureModel", function(object) {
  exp(.Call("p_pmix_reduced", object))
})

setMethod("reducedGibbsThetaFixed", "MarginalModel", function(object){
  .Call("simulate_z_reduced1", object)
})

setMethod("reducedGibbsThetaFixed", "BatchModel", function(object){
  .Call("simulate_z_reduced1_batch", object)
})

setMethod("reducedGibbsThetaSigmaFixed", "MarginalModel", function(object){
  .Call("simulate_z_reduced2", object)
})

setMethod("reducedGibbsThetaSigmaFixed", "BatchModel", function(object){
  .Call("simulate_z_reduced2_batch", object)
})


.pthetastar <- function(kmod, kmodZ1, kmodZ2, T2, ix){
  K <- k(kmod)
  if(identical(ix, seq_len(K))){
    ##
    ## Non-permuted modes
    ## estimate p_theta from existing chain
    ##
    p_theta.z <- pTheta(kmod)
  } else {
    ## one of the K! - 1 other modes
    iter(kmod, force=TRUE) <- T2
    iter(kmodZ1, force=TRUE) <- T2
    iter(kmodZ2, force=TRUE) <- T2
    kmod <- permuteZandModes(kmod, ix)
    kmodZ1 <- permuteZandModes(kmodZ1, ix)
    kmodZ2 <- permuteZandModes(kmodZ2, ix)
    kmodZ1 <- reducedGibbsZThetaFixed(kmodZ1)
    p_theta.z <- pTheta_Zfixed(kmod)
  }
  if(FALSE){
    truth <- as.numeric(table(zChain(kmodZ2)[1, ])/2500)
    result <- modes(kmodZ2)[["mixprob"]]
    ## truth and results should be in the same ballpark
    pstar <- result
    zz <- as.numeric(table(zChain(kmodZ2)[1, ]))
    gtools::ddirichlet(pstar, zz)

    Z <- zChain(kmodZ2)
    NN <- length(y(kmodZ2))
    ptab <- t(apply(Z, 1, function(x, NN) table(x)/NN, NN=NN))
    plot.ts(ptab, plot.type="single", col=1:4)
    abline(h=modes(kmodZ2)[["mixprob"]], col=1:4, lwd=2)

    ii <- argMax(kmod)
    checkEquals(pic(kmod)[ii,], modes(kmodZ2)[["mixprob"]])
    checkEquals(thetac(kmod)[ii, ],  modes(kmodZ2)[["theta"]])
    checkEquals(thetac(kmod)[ii, ],  theta(kmodZ2))
    checkEquals(sigmac(kmod)[ii, ],  sigma(kmodZ2))
    Z2 <- table(zChain(kmod)[ii, ])/NN
    checkEquals(as.numeric(Z2), pic(kmod)[ii,], tolerance=0.05)
    checkTrue(log(gtools::ddirichlet(pic(kmod)[ii,], as.numeric(Z2))) > -Inf)

    kmod3 <- reducedGibbsThetaSigmaFixed(kmod)
    Z <- zChain(kmod3)
    ptab3 <- t(apply(Z, 1, function(x, NN) table(x)/NN, NN=NN))
    plot.ts(ptab3, plot.type="single", col=1:4)
    abline(h=modes(kmodZ2)[["mixprob"]], col=1:4, lwd=2)
    p_pmix.z <- pMixProb(kmod3)
  }
  p_sigma2.z <- pSigma2(kmodZ1)
  p_pmix.z <- pMixProb(kmodZ2)
  p_theta <- mean(p_theta.z)
  p_sigma2 <- mean(p_sigma2.z)
  p_pmix <- mean(p_pmix.z)
  pstar <- c(p_theta, p_sigma2, p_pmix)
  pstar
}


setMethod("pThetaStar", "MixtureModel",
  function(kmod, maxperm=5, T2=1/10*iter(kmod)){
    burnin(kmod) <- 0
    nStarts(kmod) <- 1
    permutations <- permnK(k(kmod), Inf)
    if(nrow(permutations) > maxperm){
      index <- sample(2:nrow(permutations), maxperm)
      permutations <- permutations[c(1, index), ]
    }
    NP <- nrow(permutations)
    results <- matrix(NA, NP, 3)
    nms <- c("p(theta*|y)", "p(sigma2*|theta*, y)", "p(p*|theta*, sigma2*, y)")
    colnames(results) <- nms
    ##if(TRUE) results <- results[, 1:2, drop=FALSE]
    kmodZ1 <- reducedGibbsThetaFixed(kmod)
    kmodZ2 <- reducedGibbsThetaSigmaFixed(kmod)
    ##if(k(kmod) == 4) browser()
    for(i in seq_len(nrow(permutations))){
      ##message("entering .pthetastar")
      results[i, ] <- .pthetastar(kmod, kmodZ1, kmodZ2, T2, permutations[i, ])
    }
    results
})

.berkhof <- function(x, model){
  NP <- nrow(x)
  p.I <- prod(x[1, ])  ## Chib's estimate
  if(is.na(log(p.I))) browser()
  if(k(model) == 1 ){
    p.V <- p.I
  } else {
    p.IV <- 1/(NP-1)*sum(rowProds(x[-1, , drop=FALSE]))
    p.V <- 1/NP * p.I + (NP-1)/NP * p.IV
  }
  list(p_theta=x,
       chib=log(p.I),
       berkhof=log(p.V),
       marginal=modes(model)[["loglik"]] + modes(model)[["logprior"]] - log(p.V))
}

berkhofEstimate <- function(model, T2=0.2*iter(model), maxperm=5){
  ##message("entering pThetaStar...")
  x <- pThetaStar(model, T2=T2, maxperm=maxperm)
  ##message("entering .berkhof")
  results <- .berkhof(x, model=model)
  ##message("computing PosteriorSummary")
  PosteriorSummary(p_theta=results$p_theta,
                   chib=results$chib,
                   berkhof=results$berkhof,
                   marginal=results$marginal)
}

# get rid of this function after finalized.
summarizeMarginalEstimates <- function(x){
  chibs <- round(sapply(x, chib), 2)
  berk <- round(sapply(x, berkhof), 2)
  my <- sapply(x, marginal)
  ly <- sapply(my, length)
  if(any(ly == 0)){
    my[[which(ly==0)]] <- NA
    my <- unlist(my)
  }
  my <- round(my, 2)
  xx <- cbind(chibs, berk, my)
  colnames(xx) <- c("chib", "berkhof", "marginal")
  xx
}

summarizeMarginalEstimates <- function(x){
    chibs <- round(chib(x), 2)
    berk <- round(berkhof(x), 2)
    my <- marginal(x)
    xx <- c(chibs, berk, my)
    names(xx) <- c("chib", "berkhof", "marginal")
    xx
}

updateMultipleChains <- function(nchains, modellist, mp){
  if(k==1) nchains <- 1
  if(missing(batch)){
    kmodlist <- replicate(nchains, {
      kmod <- MarginalModel(y, k=k, mcmc.params=mp)
      posteriorSimulation(kmod)
    })
    return(kmodlist)
  }
  ## batch not missing
  kmodlist <- replicate(nchains, {
    kmod <- BatchModel(y, batch=batch, k=k, mcmc.params=mp)
    posteriorSimulation(kmod)
  })
  return(kmodlist)
}

#' Estimate the marginal density of a mixture model with k components
#'
#' This function is used to estimate the marginal density of a mixture
#' model with k components. A list of coverged models is passed by the user
#' to this function. Each model must have already converged and should be of
#' a different component size.
#' @param modlist A list of converged models with different component sizes.
#' @param post.iter number of MCMC iterations after permuting the modes.  See \code{maxperm}.
#' @param maxperm a length-one integer vector.  For a mixture model
#' with K components, there are K! possible modes.  \code{maxperm}
#' indicates the maximum number of permutations to explore.
#' @param method specifies chib or laplace method
#' @export
#' @seealso \code{\link{logBayesFactor}} for computing the bayes
#' factor for two models, \code{\link{orderModels}} for ordering a
#' list of models by decreasing marginal density, and \code{plot} for
#' visualizing the component densities.
computeMarginalLik <- function(modlist,
                                post.iter=200,
                                maxperm=3,
                                method='chib'){
    K <- sapply(modlist, k)
    my <- vector("list", length(K))
    mlist <- vector("list", length(K))
  
    for (i in seq_along(K)) {
        model_lik <- berkhofEstimate(modlist[[i]], T2=post.iter, maxperm=maxperm)
#         log_lik <- CNPBayes:::modalLoglik(model)
        xx <- summarizeMarginalEstimates(model_lik)
        my[[i]] <- xx
#         mlist[[i]] <- modlist[[i]]
    }

    if(is(modlist[[1]], "MarginalModel")){
        nms <- paste0("M", K)
        mlist <- MarginalModelList(modlist, names=nms)
    } else {
        nms <- paste0("B", K)
        mlist <- BatchModelList(modlist, names=nms)
    }

    names(my) <- nms

    results <- list(models=mlist, marginal=my)
    results
}

modelOtherModes <- function(model, maxperm=5){
  kperm <- permnK(k(model), maxperm)
  model.list <- vector("list", nrow(kperm))
  for(i in seq_along(model.list)){
    model.list[[i]] <- relabel(model, kperm[i, ])
  }
  model.list
}

#' Reorder models of varying component sizes.
#'
#' Models are ordered according to marginal likelihood. The marginal likelihood is computed for each chain of each component size model separately. The mean is taken by model, and ordering by this mean marginal is performed. For each model, the difference of marginal likelihoods is calculated for each chain and the range is taken. If the sum of these ranges across models is greater than \code{maxdev}, a NULL is returned and a warning message printed.
#' @param x the result of a call to \code{computeMarginalLik}.
#' @export
orderModels <- function(x, maxdev=5){
  models <- x$models
  ##K <- k(models)
  K <- names(models)
  ##maxdev <- sapply(K, function(x) log(factorial(x))) + 0.5
  marginal.est.list <- x$marginal
  m <- sapply(marginal.est.list, function(x) x["marginal"])
  K <- K[order(m, decreasing=TRUE)]
  ix <- match(K, names(models))
  models <- models[ix]
  return(models)
}

#' Compute the log bayes factor between models.
#'
#' Models of varying component sizes are compared. The log bayes factor is calculated comparing the two models with the highest marginal likelihood, as computed by \code{computeMarginalLik}.
#' @param x the result of a call to \code{computeMarginalLik}.
#' @export
logBayesFactor <- function(x){
  models <- orderModels(x)
  if(length(models) <= 1) {
    return(NA)
  }
  K <- k(models)
  ##nms <- names(orderModels(x))
  nms <- names(models)
  ml <- sapply(x$marginal[nms], function(x) x["marginal"])
  names(ml) <- substr(names(ml), 1, 2)
  bf <- setNames(ml[1]-ml[2], paste0(names(ml[1:2]), collapse="-"))
  bf
}

setMethod("updateMultinomialProb", "MarginalModel", function(object){
  .Call("update_multinomialPr", object)
})

setMethod("updateMultinomialProb", "BatchModel", function(object){
  .Call("update_multinomialPr_batch", object)
})

NestedMarginalModel <- function(model){
  model <- useModes(model)
  if(k(model) == 1) stop("model only has 1 component")
  ##
  ## Transition to K=K-1 model, selecting K-1 modes with the greatest separation
  ##
  th <- setNames(theta(model), paste0("comp", seq_along(theta(model))))
  d <- diff(sort(th))
  drop_component <- names(d)[which.min(d)]
  if(which.min(d) > 1){
    replacement_for_dropped_component <- names(d)[which.min(d)-1]
  } else {
    replacement_for_dropped_component <- names(th)[!names(th) %in% names(d)]
  }
  zlabel_replace <- as.integer(strsplit(replacement_for_dropped_component, "comp")[[1]][2])
  zlabel_dropped <- as.integer(strsplit(drop_component, "comp")[[1]][2])
  zz <- z(model)
  zz[zz == zlabel_dropped] <- zlabel_replace
  ii <- zlabel_dropped
  ## next relabel z such that the maximum value can not be greater
  ## than the number of components
  zz <- as.integer(factor(zz))
  hypp <- hyperParams(model)
  k(hypp) <- k(model)-1L
  kmod <-  new("MarginalModel",
               k=k(model)-1L,
               hyperparams=hypp,
               theta=theta(model)[-ii],
               sigma2=sigma2(model)[-ii],
               mu=mu(model),
               tau2=tau2(model),
               nu.0=nu.0(model),
               sigma2.0=sigma2.0(model),
               pi=p(model)[-ii],
               data=y(model),
               data.mean=dataMean(model)[-ii],
               data.prec=dataPrec(model)[-ii],
               z=zz,
               zfreq=as.integer(table(zz)),
               probz=matrix(0, length(y(model)), k(model)-1),
               mcmc.chains=McmcChains(),
               batch=batch(model),
               batchElements=nBatch(model),
               mcmc.params=mcmcParams(model))
  ## initialize modes
  modes(kmod) <- list(theta=theta(kmod),
                      sigma2=sigma2(kmod),
                      mixprob=p(kmod),
                      mu=mu(kmod),
                      tau2=tau2(kmod),
                      nu0=nu.0(kmod),
                      sigma2.0=sigma2.0(kmod),
                      zfreq=zFreq(kmod),
                      loglik=log_lik(kmod),
                      logprior=logPrior(kmod))
  chains(kmod) <- McmcChains(kmod)
  kmod
}

setMethod("computeLoglik", "BatchModel", function(object){
  .Call("compute_loglik_batch", object)
})

setMethod("computePotential", "BatchModel", function(object){
  computeLogLikxPrior(object)
})

setMethod("computeLoglik", "MarginalModel", function(object){
  .Call("loglik", object)
})

setMethod("computeLogLikxPrior", "MixtureModel", function(object){
  .Call("compute_llxprior", object)
})
