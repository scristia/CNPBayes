fullGibbs <- function(post, move_chain){
  theta(post) <- updateTheta(post)
  sigma2(post) <- updateSigma2(post)
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

## Sample with full conditionals
## p(sigma2 | y, theta*, z)
## and
## p(z | y, sigma2, theta*)
reducedGibbsSigma2 <- function(post, move_chain){
  ##theta(post) <- updateTheta(post)
  sigma2(post) <- updateSigma2(post)
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


setMethod("posteriorTheta", "MarginalModel", function(object, mcmcp){
  .posteriorTheta(object, mcmcp)
})

setMethod("posteriorTheta", "BatchModel", function(object, mcmcp){
  .posteriorThetaBatch(object, mcmcp)
})

setMethod("posteriorSigma2", "MarginalModel", function(object, mcmcp){
  .posteriorSigma2(object, mcmcp)
})

setMethod("posteriorSigma2", "BatchModel", function(object, mcmcp){
  .posteriorSigma2(object, mcmcp)
})

.posteriorTheta <- function(object, mcmcp){
  ##object <- useModes(object)
  thetastar <- modes(object)[["theta"]]
  S <- 2:(savedIterations(mcmcp)+1)
  message("Running an additional ", savedIterations(mcmcp), " simulations from full Gibbs to estimate p(theta*|y)")
  do_thin <- thin(mcmcp) > 1
  T <- seq_len(thin(mcmcp))
  p.theta <- rep(NA, length(S))
  p.theta[1] <- prod(dnorm(thetastar, mu(object), tau(object)))
  for(s in S){
    object <- fullGibbs(object, TRUE)
    object <- moveChain(object, s)
    p.theta[s] <- prod(dnorm(thetastar, mu(object), tau(object)))
    ## update without moving chain
    if(do_thin){
      for(t in T) object <- fullGibbs(object, FALSE)
    }
  }
  p.theta
}

.posteriorThetaBatch <- function(object, mcmcp){
  ##object <- useModes(object)
  thetastar <- as.numeric(modes(object)[["theta"]])
  S <- 2:(savedIterations(mcmcp)+1)
  message("Running an additional ", savedIterations(mcmcp), " simulations from full Gibbs to estimate p(theta*|y)")
  do_thin <- thin(mcmcp) > 1
  T <- seq_len(thin(mcmcp))
  p.theta <- rep(NA, length(S))
  mus <- as.numeric(matrix(mu(object), nBatch(object), k(object), byrow=TRUE))
  taus <- as.numeric(matrix(tau(object), nBatch(object), k(object), byrow=TRUE))
  p.theta[1] <- prod(dnorm(thetastar, mus, taus))
  for(s in S){
    object <- fullGibbs(object, TRUE)
    object <- moveChain(object, s)
    mus <- as.numeric(matrix(mu(object), nBatch(object), k(object), byrow=TRUE))
    taus <- as.numeric(matrix(tau(object), nBatch(object), k(object), byrow=TRUE))
    p.theta[s] <- prod(dnorm(thetastar, mu(object), tau(object)))
    ## update without moving chain
    if(do_thin){
      for(t in T) object <- fullGibbs(object, FALSE)
    }
  }
  p.theta
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

.posteriorSigma2 <- function(object, mcmcp){
  ##object <- useModes(object)
  isigma2star <- as.numeric(1/sigma2(object))
  S <- 2:(savedIterations(mcmcp)+1)
  message("Running an additional ", savedIterations(mcmcp), " simulations from reduced Gibbs to estimate p(sigma2*|y)")
  do_thin <- thin(mcmcp) > 1
  T <- seq_len(thin(mcmcp))
  post <- rep(NA, length(S))
  post[1] <- sum(log(dgamma(isigma2star, shape=0.5*nu.0(object), rate=0.5*nu.0(object)*sigma2.0(object))))
  for(s in S){
    object <- reducedGibbsSigma2(object, TRUE)
    object <- moveChain(object, s)
    shape.s <- 0.5*nu.0(object)
    rate.s <- 0.5*nu.0(object)*sigma2.0(object)
    post[s] <- sum(log(dgamma(isigma2star, shape=shape.s, rate=rate.s)))
    if(do_thin){
      for(t in T) object <- reducedGibbsSigma2(object, FALSE)
    }
  }
  post
}

posteriorP <- function(object, mcmcp){
  ##object <- useModes(object)
  x <- p(object)
  S <- 2:(savedIterations(mcmcp)+1)
  message("Running an additional ", savedIterations(mcmcp), " simulations from reduced Gibbs to estimate p(pi*|y)")
  do_thin <- thin(mcmcp) > 1
  T <- seq_len(thin(mcmcp))
  post <- rep(NA, length(S))
  ##
  ## The ddirichlet function uses a very strict definition requires
  ## sum(p) ==1 instead of allowing a tolerance
  ##
  if(length(x) > 1){
    x[length(x)] <- 1-sum(x[-length(x)])
  }
  alpha.n <- updateAlpha(object)
  post[1] <- ddirichlet(x, alpha.n)
  for(s in S){
    object <- reducedGibbsP(object, TRUE)
    object <- moveChain(object, s)
    alpha.n <- updateAlpha(object)
    post[s] <- ddirichlet(x, alpha.n)
    ##if(post[s] > 10e3) browser()
    if(do_thin){
      for(t in T) object <- reducedGibbsP(object, FALSE)
    }
  }
  ## some values from dirichlet seem to be very extreme
  ##median(post)
  post
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

##computeLoglik2 <- function(object){
##  thetastar <- modes(object)[["theta"]]
##  sigmastar <- sqrt(modes(object)[["sigma2"]])
##  yy <- y(object)
##  pp <- modes(object)[["mixprob"]]
##  ## We should keep the z-values for the MCMC iteration that has the
##  ## highest likelihood.  Here, I use the maximum a posteriori
##  ## estimate as a plug-in, but this is not quite correct
##  zz <- map(object)
##  lik <- pp[zz]*dnorm(yy, thetastar[zz], sigmastar[zz])
##  sum(log(lik))
##}

useModes <- function(object){
  m2 <- object
  theta(m2) <- modes(object)[["theta"]]
  sigma2(m2) <- modes(object)[["sigma2"]]
  z(m2) <- factor(map(object), levels=seq_len(k(object)))
  tau2(m2) <- modes(object)[["tau2"]]
  nu.0(m2) <- modes(object)[["nu0"]]
  sigma2.0(m2) <- modes(object)[["sigma2.0"]]
  p(m2) <- modes(object)[["mixprob"]]
  m2
}


##computeBayesFactor <- function(model1, model2, mcmcp){
##  log.m.y1 <- marginalY(model1, mcmcp)
##  log.m.y2 <- marginalY(model2, mcmcp)
##  exp(log.m.y1-log.m.y2)
##}

#' @export
partialGibbs <- function(model, mcmcp){
  m <- useModes(model)
  mcmcChains(m) <- McmcChains(m, mcmcp)
  ptheta <- posteriorTheta(m, mcmcp)
  modes(m) <- modes(model)
  lpsigma2 <- posteriorSigma2(m, mcmcp)
  modes(m) <- modes(model)
  pmix <- posteriorP(m, mcmcp)
  results <- cbind(ptheta, lpsigma2, pmix)
  colnames(results) <- c("theta", "sigma2", "p")
  results
}


#' @export
partialGibbsSummary <- function(model, x){
  p.priors <- priorProbs(model)
  ll <- logLik(mcmcChains(model))
  ll <- ll[is.finite(ll)]
  log.lik <- mean(ll)
  results <- setNames(c(log(p.priors),
                        log.lik,
                        log(mean(x[, "theta"])),
                        mode.numeric(x[, "sigma2"]),
                        log(mean(x[, "p"]))),
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
computeMarginalPr <- function(model, mcmcp){
  kperm <- permn(seq_len(k(model)))
  model <- useModes(model)
  model.list <- lapply(kperm, function(zindex, model) relabel(model, zindex), model=model)
  J <- min(3, length(model.list))
  model.list <- model.list[seq_len(J)]
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
