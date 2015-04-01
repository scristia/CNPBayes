setMethod("fullGibbs", "MarginalModel", function(object, mcmcp){
  object=.Call("mcmc_marginal_burnin", object, mcmcp)
  object=.Call("mcmc_marginal", object, mcmcp)
})

setMethod("fullGibbs", "BatchModel", function(object, mcmcp){
  ##.fullGibbs(object, mcmcp)
  object=.Call("mcmc_batch_burnin", object, mcmcp)
  object=.Call("mcmc_batch", object, mcmcp)
})

## TODO: change name to reducedGibbs
.fullGibbs <- function(post, mcmcp){
  up <- paramUpdates(mcmcp)
  if(up["theta"] >= 1L)
    theta(post) <- updateTheta(post)
  if(up["sigma2"] >= 1L)
    sigma2(post) <- updateSigma2(post)
  if(up["mu"] >= 1L)
    mu(post) <- updateMu(post)
  if(up["tau2"] >= 1L)
    tau2(post) <- updateTau2(post)
  if(up["sigma2.0"] >= 1L)
    sigma2.0(post) <- updateSigma2.0(post)
  if(up["nu.0"] >= 1L)
    nu.0(post) <- updateNu.0(post)
  if(up["p"] >= 1L)
    p(post) <- updateMixProbs(post)
  if(up["z"] >= 1L)
    z(post) <- updateZ(post)
  dataMean(post) <- computeMeans(post)
  dataPrec(post) <- 1/computeVars(post)
  logpotential(post) <- computePotential(post)
  post
}

## Sample with full conditionals
## p(sigma2 | y, theta*, z)
## and
## p(z | y, sigma2, theta*)
## reducedGibbsSigma2 <- function(post, move_chain){
##   ##theta(post) <- updateTheta(post)
##   sigma2(post) <- updateSigma2(post)
##   mu(post) <- updateMu(post)
##   tau2(post) <- updateTau2(post)
##   sigma2.0(post) <- updateSigma2.0(post)
##   nu.0(post) <- updateNu.0(post)
##   p(post) <- updateMixProbs(post)
##   z(post) <- updateZ(post)
##   dataMean(post) <- computeMeans(post)
##   dataPrec(post) <- 1/computeVars(post)
##   logpotential(post) <- computePotential(post)
##   post
## }

reducedGibbsZ <- function(post, move_chain){
  theta(post) <- updateTheta(post)
  sigma2(post) <- updateSigma2(post)
  mu(post) <- updateMu(post)
  tau2(post) <- updateTau2(post)
  sigma2.0(post) <- updateSigma2.0(post)
  nu.0(post) <- updateNu.0(post)
  p(post) <- updateMixProbs(post)
  ##z(post) <- updateZ(post)
  dataMean(post) <- computeMeans(post)
  dataPrec(post) <- 1/computeVars(post)
  logpotential(post) <- computePotential(post)
  post
}

reducedGibbsP <- function(post, move_chain){
  ##theta(post) <- updateTheta(post)
  ##sigma2(post) <- updateSigma2(post)
  mu(post) <- updateMu(post)
  tau2(post) <- updateTau2(post)
  sigma2.0(post) <- updateSigma2.0(post)
  nu.0(post) <- updateNu.0(post)
  p(post) <- updateMixProbs(post)
  z(post) <- updateZ(post)
  dataMean(post) <- computeMeans(post)
  dataPrec(post) <- 1/computeVars(post)
  logpotential(post) <- computePotential(post)
  post
}


setMethod("posteriorTheta", "MarginalModel", function(object){
  ##.posteriorTheta(object, mcmcp)
  .Call("marginal_theta", object, mcmcParams(object))
})

setMethod("posteriorTheta", "BatchModel", function(object){
  .posteriorThetaBatch(object)
})

setMethod("posteriorSigma2", "MarginalModel", function(object){
  mcmcp <- mcmcParams(object)
  paramUpdates(mcmcp)["theta"] <- 0L
  .Call("marginal_sigma2", object, mcmcParams(object))
})

setMethod("posteriorSigma2", "BatchModel", function(object){
  paramUpdates(object)["theta"] <- 0L
  .posteriorSigma2(object)
})

## .posteriorTheta <- function(object, mcmcp){
##   ##object <- useModes(object)
##   thetastar <- modes(object)[["theta"]]
##   S <- iter(mcmcp)
##   object = .Call("mcmc_marginal", object, mcmcp)
##   mus <- muc(object)
##   taus <- tauc(object)
##   K <- k(object)
##   p.theta <- matrix(NA, S+1, K)
##   for(k in seq_len(k(object))){
##     p.theta[, k] <- dnorm(thetastar[k], mus, taus)
##   }
##   rowProds(p.theta)
## }

.posteriorThetaBatch <- function(object){
  mcmcp <- mcmcParams(object)
  thetastar <- as.numeric(modes(object)[["theta"]])
  S <- iter(mcmcp)
  object <- .Call("mcmc_batch_burnin", object, mcmcp)
  object <- .Call("mcmc_batch", object, mcmcp)
  mus <- muc(object)
  taus <- tauc(object)
  thetas <- thetac(object)
  K <- k(object)
  B <- nBatch(object)
  p.theta <- matrix(NA, S, B*K)
  j <- 1
  for(k in 1:K){
    for(b in 1:B){
      p.theta[, j] <- dnorm(thetastar[j], mus[, k], taus[, k])
      j <- j+1
    }
  }
  lp.theta <- rowSums(log(p.theta))
  lp.theta
}

mode.numeric <- function(x, na.rm=TRUE){
  if(na.rm){
    notna <- !is.na(x)
    xx <- x[notna]
    dens <- density(xx)
    x_mode <- dens$x[which.max(dens$y)]
  } else {
    dens <- density(x)
    x_mode <- dens$x[which.max(dens$y)]
  }
  return(x_mode)
}

.posteriorSigma2 <- function(object){
  mcmcp <- mcmcParams(object)
  ##paramUpdates(mcmcp)["theta"] <- 0L
  isigma2star <- 1/modes(object)[["sigma2"]]
  object <- fullGibbs(object, mcmcp)
  nuc <- nu.0(mcmcChains(object))
  s20 <- sigma2.0(mcmcChains(object))
  shape.s <- 0.5*nuc
  rate.s <- 0.5*nuc*s20
  K <- k(object)
  post <- matrix(NA, length(nuc), K)
  for(k in seq_len(K)){
    post[, k] <- log(dgamma(isigma2star[k], shape=shape.s, rate=rate.s))
  }
  rowSums(post)
}

posteriorP <- function(object){
  mcmcp <- mcmcParams(object)
  x <- p(object)
  paramUpdates(mcmcp)[c("theta", "sigma2")] <- 0L
  if(is(object, "MarginalModel")){
    object=.Call("mcmc_marginal_burnin", object, mcmcp)
    object=.Call("mcmc_marginal", object, mcmcp)
  } else {
    object=.Call("mcmc_batch_burnin", object, mcmcp)
    object=.Call("mcmc_batch", object, mcmcp)
  }
  ztab <- zFreq(mcmcChains(object))
  ztab <- ztab + alpha(object)
  probs <- apply(ztab, 1, function(alpha, x){
    ddirichlet(x, alpha)
  }, x=x)
  log(probs)
}

.priorProbsBatch <- function(object){
  ## evaluate the prior probabilities at the modal values
  x <- modes(object)[["mu"]]
  p.mu <- prod(dnorm(x, mu.0(object), tau.0(object)))

  thetastar <- modes(object)[["theta"]]
  mustar <- modes(object)[["mu"]]
  tau2star <- modes(object)[["tau2"]]

  mustar <- matrix(mustar, nrow(thetastar), length(mustar), byrow=TRUE)
  tau2star <- matrix(tau2star, nrow(thetastar), length(tau2star), byrow=TRUE)
  p.theta <- prod(dnorm(as.numeric(thetastar), as.numeric(mustar), as.numeric(sqrt(tau2star))))

  sigma2star <- modes(object)[["sigma2"]]
  nu0 <- modes(object)[["nu0"]]
  sigma2.0 <- modes(object)[["sigma2.0"]]
  p.sigma2 <- prod(dgamma(as.numeric(1/sigma2star), shape=0.5*nu0, rate=0.5*nu0*sigma2.0))

  x <- 1/modes(object)[["tau2"]]
  p.tau2 <- prod(dgamma(x, shape=0.5*eta.0(object), rate=0.5*eta.0(object)*m2.0(object)))

  .beta <- hyperParams(object)@beta
  p.nu0 <- dgeom(nu0, .beta)

  x <- modes(object)[["sigma2.0"]]
  hyp <- hyperParams(object)
  p.sigma2.0 <- dgamma(x, shape=a(hyp), rate=b(hyp))

  p.theta*p.sigma2*p.mu*p.tau2*p.nu0*p.sigma2.0
}



setMethod("priorProbs", "MarginalModel", function(object) {
  .priorProbs(object)
})

setMethod("priorProbs", "BatchModel", function(object) {
  .priorProbsBatch(object)
})

.priorProbs <- function(object){
  ## evaluate the prior probabilities at the modal values
  x <- modes(object)[["mu"]]
  p.mu <- dnorm(x, mu.0(object), tau.0(object))

  thetastar <- modes(object)[["theta"]]
  mustar <- modes(object)[["mu"]]
  tau2star <- modes(object)[["tau2"]]
  p.theta <- prod(dnorm(thetastar, mustar, sqrt(tau2star)))

  sigma2star <- modes(object)[["sigma2"]]
  nu0 <- modes(object)[["nu0"]]
  sigma2.0 <- modes(object)[["sigma2.0"]]
  p.sigma2 <- prod(dgamma(1/sigma2star, shape=0.5*nu0, rate=0.5*nu0*sigma2.0))

  x <- 1/tau2star
  p.tau2 <- dgamma(x, shape=0.5*eta.0(object), rate=0.5*eta.0(object)*m2.0(object))

  .beta <- hyperParams(object)@beta
  p.nu0 <- dgeom(nu0, .beta)

  x <- modes(object)[["sigma2.0"]]
  hyp <- hyperParams(object)
  p.sigma2.0 <- dgamma(x, shape=a(hyp), rate=b(hyp))

  p.theta*p.sigma2*p.mu*p.tau2*p.nu0*p.sigma2.0
}

useModes <- function(object){
  m2 <- object
  theta(m2) <- modes(object)[["theta"]]
  sigma2(m2) <- modes(object)[["sigma2"]]
  tau2(m2) <- modes(object)[["tau2"]]
  nu.0(m2) <- modes(object)[["nu0"]]
  sigma2.0(m2) <- modes(object)[["sigma2.0"]]
  p(m2) <- modes(object)[["mixprob"]]
  zFreq(m2) <- as.integer(modes(object)[["zfreq"]])
  logLik(m2) <- modes(object)[["loglik"]]
  logPrior(m2) <- modes(object)[["logprior"]]
  ##
  ## update z using the modal values from above
  ##
  if(is(object, "MarginalModel")){
    z(m2) <- .Call("update_z", m2)
  } else {
    z(m2) <- .Call("update_z_batch", m2)
  }
  m2
}

#' @export
partialGibbs <- function(object){
  log_ptheta <- posteriorTheta(object)
  log_psigma2 <- posteriorSigma2(object)
  log_pmix <- posteriorP(object)
  results <- cbind(log_ptheta, log_psigma2, log_pmix)
  colnames(results) <- c("logtheta", "logsigma2", "logp")
  results
}


#' @export
partialGibbsSummary <- function(model, x){
  log.prior <- modes(model)[["logprior"]]
  log.lik <- modes(model)[["loglik"]]
  results <- setNames(c(log.prior,
                        log.lik,
                        mean(x[, "logtheta"]),
                        mean(x[, "logsigma2"]),
                        mean(x[, "logp"])),
                      c("prior", "loglik", "theta", "sigma2", "pmix"))
  results
}

posteriorPsi <- function(x) sum(x[["theta"]] + x[["sigma2"]] + x[["pmix"]])


#' @export
marginalY <- function(x){
  x[["prior"]] + x[["loglik"]] - posteriorPsi(x)
}

##bayesFactor <- function(m1, m2) exp(marginalY(m1) - marginalY(m2))


pairwiseModels <- function(nmodels){
  sapply(apply(combn(nmodels, 2), 2, list), "[", 1)
}


modelFiles <- function(outdir, cnpids){
  list(batch=file.path(outdir, paste0("mm_batch_", cnpids, ".rds")),
       marginal=file.path(outdir, paste0("mm_marginal_", cnpids, ".rds")))
}

postSummaryFiles <- function(outdir, cnpids){
  model.files <- modelFiles(outdir, cnpids)

  files <- model.files[["batch"]]
  my.batch <- gsub("mm", "margy", files)
  my.batch2 <- gsub("margy", "margy2", my.batch)
  my.batch3 <- gsub("margy", "margy3", my.batch)

  mfiles <- model.files[["marginal"]]
  my.margy <- gsub("mm", "margy", mfiles)
  my.margy2 <- gsub("margy", "margy2", my.margy)
  my.margy3 <- gsub("margy", "margy3", my.margy)

  mfiles <- cbind(my.margy, my.margy2, my.margy3)
  bfiles <- cbind(my.batch, my.batch2, my.batch3)
  list(marginal=mfiles,
       batch=bfiles)
}




readLogLik <- function(files){
  parseLogLik <- function(file){
    tmp <- readRDS(file)
    loglik <- tmp["loglik", ]
  }
  m.loglik <- do.call(rbind, lapply(files[["marginal"]][,1], parseLogLik))
  b.loglik <- do.call(rbind, lapply(files[["batch"]][,1], parseLogLik))
  list(marginal=m.loglik,
       batch=b.loglik)
}


readPThetaPrior <- function(files){
  parsefun <- function(file){
    tmp <- readRDS(file)
    prior <- tmp["prior", ]
  }
  marginal <- do.call(rbind, lapply(files[["marginal"]][,1], parsefun))
  batch <- do.call(rbind, lapply(files[["batch"]][,1], parsefun))
  list(marginal=marginal, batch=batch)
}




listPostTheta <- function(x){
  M1 <- x[, 1, drop=FALSE]
  B1 <- x[, 1, drop=FALSE]
  M2 <- x[, c(2,3)]
  B2 <- x[, c(2,3)+9]
  M3 <- x[, 4:6]
  B3 <- x[, 4:6 + 9]
  M4 <- x[, 7:9]
  B4 <- x[, 7:9 + 9]
  list(M1=M1,
       M2=M2,
       M3=M3,
       M4=M4,
       B1=B1,
       B2=B2,
       B3=B3,
       B4=B4)
}

deltaPTheta <- function(x){
  require(matrixStats)
  maxs <- rowMaxs(x)
  mins <- rowMins(x)
  maxs-mins
}


avgMarginalTheta <- function(files){
  p_theta <- readPostTheta(files)
  D <- lapply(p_theta, deltaPTheta)
  D <- do.call(cbind, D)
  maxD <- rowMaxs(D, na.rm=TRUE)
  ##
  ## If maxD is greater than 5, the estimates of p(theta|y) are not
  ## reliable.  Plug in -Inf to force selection of another model.
  ##
  marginal_theta <- matrix(NA, nrow(D), ncol(D))
  dimnames(marginal_theta) <- list(rownames(se.ea), colnames(D))
  for(j in 1:ncol(D)){
    p <- p_theta[[j]]
    if(ncol(p) == 1){
      marginal_theta[, j] <- p
      next()
    }
    accept <- D[, j] <= 5
    mns <- rowMeans(p)
    mns[!accept] <- -Inf
    marginal_theta[, j] <- mns
  }
  marginal_theta
}

#' @export
computeMarginalProbs <- function(model, mcmcp, maxperm=5){
  mmod <- useModes(model)
  iter(mmod) <- iter(mcmcp)
  burnin(mmod) <- burnin(mcmcp)
  nStarts(mmod) <- 1L
  K <- k(model)
  model.list <- ModelEachMode(mmod, maxperm)
  ##
  ## Run partial Gibbs sampler for each mode
  ##
  results <- matrix(NA, length(model.list), 5)
  for(i in seq_along(model.list)){
    model1 <- model.list[[i]]
    pg <- partialGibbs(model1)
    results[i, ] <- partialGibbsSummary(model1, pg)
  }
  colnames(results) <- c("logprior", "loglik", "logtheta", "logsigma2", "logp")
  marginal.y <- results[, "logprior"] + results[, "loglik"] - results[, "logtheta"] -
      results[, "logsigma2"] - results[, "logp"]
  as.numeric(marginal.y)
}

#' @export
computeMarginalEachK <- function(data, K=1:4, hypp, mcmcp=McmcParams(), MAX.RANGE=5){
  j <- 1
  marginaly <- setNames(rep(NA, length(K)), paste0("M", K))
  if(missing(hypp)) hypp <- Hyperparameters("marginal")
  for(k in K){
    k(hypp) <- k
    if(k == 1){
      mp <- McmcParams(iter=min(iter(mcmcp), 500), burnin=min(burnin(mcmcp), 100))
    } else mp <- mcmcp
    kmod <- MarginalModel(data, k=k, mcmc.params=mp, hypp=hypp)
    kmod <- posteriorSimulation(kmod)
    m.y(kmod) <- computeMarginalProbs(kmod, mp)
    my <- m.y(kmod)
    if( diff(range(my)) < MAX.RANGE ){
      marginaly[j] <- mean(my)
    }
    j <- j+1
  }
  marginaly
}

#' @export
computeMarginalEachK2 <- function(data, batch, K=1:4, mcmcp=McmcParams(),
                                  MAX.RANGE=5,
                                  hypp){
  if(missing(hypp)) hypp <- Hyperparameters("batch")
  j <- 1
  marginaly <- setNames(rep(NA, length(K)), paste0("B", K))
  for(k in K){
    k(hypp) <- k
    my <- .computeMarginal(data, batch, mcmcp, hypp)
    my <- my[!is.nan(my)]
    if( diff(range(my)) < MAX.RANGE ){
      marginaly[j] <- mean(my)
    }
    j <- j+1
  }
  marginaly
}

.computeMarginal <- function(data, batch, mcmcp, hypp){
  if(k(hypp) == 1){
    mp <- McmcParams(iter=min(iter(mcmcp), 500), burnin=min(burnin(mcmcp), 100), nStarts=1)
  } else mp <- mcmcp
  kmod <- BatchModel(data, batch, k=k(hypp), mcmc.params=mp, hypp=hypp)
  ##if(k(kmod) == 3) browser()
  kmod <- posteriorSimulation(kmod)
  m.y(kmod) <- computeMarginalProbs(kmod, mp)
  my <- m.y(kmod)
  my
}


ModelEachMode <- function(model, maxperm=5){
  kperm <- permn(seq_len(k(model)))
  kperm <- kperm[1:min(maxperm, length(kperm))]
  ##
  ##  Reorder the z's
  ##
  model.list <- lapply(kperm, function(zindex, model) relabel(model, zindex), model=model)
  is_marginal <- is(model.list[[1]], "MarginalModel")
  ##
  ## Reorder the modal values from the original model according to z
  ##
  if(length(model.list) == 1) return(model.list)
  for(i in 2:length(model.list)){
    index <- kperm[[i]]
    modal.values <- modes(model.list[[i]])
    if(is_marginal){
      modal.values[["theta"]] <- modal.values[["theta"]][index]
      modal.values[["sigma2"]] <- modal.values[["sigma2"]][index]
    } else {
      modal.values[["theta"]] <- modal.values[["theta"]][, index]
      modal.values[["sigma2"]] <- modal.values[["sigma2"]][, index]
      modal.values[["mu"]] <- modal.values[["mu"]][index]
      modal.values[["tau2"]] <- modal.values[["tau2"]][index]
    }
    modal.values[["p"]] <- modal.values[["p"]][index]
    modes(model.list[[i]]) <- modal.values
  }
  model.list
}

#' @export
computeMarginalPr <- function(model, mcmcp){
  model.list <- ModelEachMode(model)
  ##J <- min(3, length(model.list))
  ##model.list <- model.list[seq_len(J)]
  pg <- lapply(model.list, partialGibbs, mcmcp=mcmcp)
  results <- foreach(model=model.list, x=pg, .combine="cbind") %do% {
    partialGibbsSummary(model, x)
  }
  results <- as.matrix(results, 5, J)
  colnames(results) <- paste0("M", k(model), "_", seq_len(J))
  results
}

#' @export
posteriorRange <- function(x, thr=5){
  msg <- FALSE
  post <- colSums(x[c("theta", "sigma2", "pmix"), , drop=FALSE])
  diff(range(post))
}
