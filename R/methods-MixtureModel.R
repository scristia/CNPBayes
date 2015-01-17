setMethod("hyperParams", "MixtureModel", function(object) object@hyperparams)

setReplaceMethod("hyperParams", c("MixtureModel", "Hyperparameters"),
                 function(object, value) {
                   object@hyperparams <- value
                   object
                 })


setMethod("batch", "MixtureModel", function(object) object@batch)

setReplaceMethod("batch", "MixtureModel", function(object, value){
  object@batch <- value
  object
})


observed <- function(object) object@data



#' @export
sigma <- function(object) sqrt(sigma2(object))

#' @export
tau <- function(object) sqrt(tau2(object))

setMethod("nu.0", "MixtureModel", function(object) object@nu.0)

setMethod("sigma2.0", "MixtureModel", function(object) object@sigma2.0)

setMethod("y", "MixtureModel", function(object) object@data)

setMethod("z", "MixtureModel", function(object) object@z)

#' @export
p <- function(object) object@pi

nComp <- function(object) length(p(object))
dataMean <- function(object) object@data.mean
dataPrec <- function(object) object@data.prec
dataSd <- function(object) sqrt(1/dataPrec(object))

#' @export
logpotential <- function(object) object@logpotential

setMethod("k", "MixtureModel", function(object) k(hyperParams(object)))

setReplaceMethod("z", "MixtureModel", function(object, value){
  object@z <- value
  object
})

setReplaceMethod("theta", "MixtureModel", function(object, value){
  object@theta <- value
  object
})

setReplaceMethod("sigma2", "MixtureModel", function(object, value){
  object@sigma2 <- value
  object
})

setReplaceMethod("p", "MixtureModel", function(object, value){
  object@pi <- value
  object
})





setReplaceMethod("nu.0", "MixtureModel", function(object, value){
  object@nu.0 <- value
  object
})

setReplaceMethod("sigma2.0", "MixtureModel", function(object, value){
  object@sigma2.0 <- value
  object
})

setReplaceMethod("logpotential", "MixtureModel", function(object, value){
  object@logpotential <- value
  object
})

setMethod("mcmcChains", "MixtureModel", function(object) object@mcmc.chains)

setMethod("dat", "MixtureModel", function(object) object@data)

setReplaceMethod("dat", "MixtureModel", function(object, value) {
  object@data <- value
  object
})

setReplaceMethod("mcmcChains", "MixtureModel", function(object, value){
  object@mcmc.chains <- value
  object
})

setMethod("showMeans", "MarginalModel", function(object){
  paste(round(theta(object), 2), collapse=", ")
})

setMethod("showSigmas", "MarginalModel", function(object){
  paste(round(sqrt(sigma2(object)), 2), collapse=", ")
})

setMethod("showMeans", "BatchModel", function(object){
  thetas <- round(theta(object), 2)
##  dimnames(thetas) <- list(paste0("batch", uniqueBatch(object)),
##                           paste0("component", seq_len(k(object))))
  ##
  mns <- c("\n", paste0(t(cbind(thetas, "\n")), collapse="\t"))
  mns <- paste0("\t", mns[2])
  mns <- paste0("\n", mns[1])
  mns
})

setMethod("showSigmas", "BatchModel", function(object){
  sigmas <- round(sqrt(sigma2(object)), 2)
##  dimnames(sigmas) <- list(paste0("batch", uniqueBatch(object)),
##                           paste0("component", seq_len(k(object))))
  sigmas <- c("\n", paste0(t(cbind(sigmas, "\n")), collapse="\t"))
  sigmas <- paste0("\t", sigmas[2])
  sigmas <- paste0("\n", sigmas[1])
  sigmas
})



setMethod("show", "MixtureModel", function(object){
  cat("An object of class '", class(object), "'\n")
  cat("  data: \n")
  cat("     n           :", length(y(object)), "\n")
  cat("     k           :", nComp(object), "\n")
  cat("     table(z)    :", paste(tablez(object), collapse=", "), "\n")
  cat("     mix prob (s):", paste(round(p(object), 2), collapse=", "), "\n")
  ##post.sd <- paste(round(sqrt(sigma2(object)), 2), collapse=", ")
  sds <- showSigmas(object)
  mns <- showMeans(object)
  cat("     sigma (s)   :", sds, "\n")
  cat("     theta (s)   :", mns, "\n")
  cat("     sigma2.0 (s):", round(sigma2.0(object), 2), "\n")
  cat("     nu.0 (s)    :", nu.0(object), "\n")
  ##cat("     tau2 (s)    :", round(tau2(object), 2),  "\n")
  ##cat("     mu (s)      :", round(mu(object), 2), "\n")
  cat("     potential(s):", round(logpotential(object), 2), "\n")
})

setMethod("tablez", "MarginalModel", function(object) table(z(object)))

setGeneric("updateMixProbs", function(object) standardGeneric("updateMixProbs"))

setMethod("updateMixProbs", "MixtureModel", function(object){
  alpha.n <- updateAlpha(object)
  as.numeric(rdirichlet(1, alpha.n))
})

setMethod("updateMixProbs", "BatchModel", function(object){
  ##With batch-specific, there's not any borrowing
  ## of strength between batches even in the higher levels of the
  ## model.  Unless we have a prior on the parameters for the
  ## Dirichlet.
  ##
  alpha.n <- updateAlpha(object)
  as.numeric(rdirichlet(1, alpha.n))
})

setMethod("alpha", "MixtureModel", function(object) alpha(hyperParams(object)))

##setMethod("alpha", "MixtureModel", function(object) alpha(hyperParams(object)))

setMethod("updateAlpha", "MixtureModel", function(object){
  alpha(object) + table(z(object))
})

##setMethod("updateAlpha", "BatchModel", function(object){
##  alpha(object) + tablez(object)
##})


updateAll <- function(post, move_chain, check_labels=FALSE){
  ##
  ## sample new value of (theta_h, sigma2_h)
  ##   - note these are independent in the prior
  ##
  theta(post) <- updateTheta(post)
  sigma2(post) <- updateSigma2(post)
  ##
  ## update (mu, tau2)
  ##
  mu(post) <- updateMu(post)
  tau2(post) <- updateTau2(post)
  ##
  ## update pi
  ##
  p(post) <- updateMixProbs(post)
  ##
  ##  update auxillary variable
  ##
  z(post) <- updateZ(post)
  ##
  dataMean(post) <- computeMeans(post)
  dataPrec(post) <- 1/computeVars(post)
  #
  ##
  ##  update (nu.0, sigma2.0)
  ##
  sigma2.0(post) <- updateSigma2.0(post)
  nu.0(post) <- updateNu.0(post)
  logpotential(post) <- computePotential(post)
  if(check_labels){
    post <- .updateLabels(post)
  }
  if(move_chain){
    probz(post) <- .computeProbZ(post)
  }
  post
}

.computeProbZ <- function(object){
  pZ <- probz(object)
  zz <- as.integer(z(object))
  for(j in seq_len(k(object))){
    pZ[, j] <- pZ[, j] + as.integer(zz==j)
  }
  pZ
}

setReplaceMethod("probz", "MixtureModel", function(object, value){
  object@probz <- value
  object
})

setMethod("probz", "MixtureModel", function(object) object@probz)

runBurnin <- function(object, mcmcp){
  if(burnin(mcmcp)==0){
    object <- moveChain(object, 1)
    probz(object) <- .computeProbZ(object)
    return(object)
  }
  message("Running burnin...")
  mcmc_burnin <- McmcParams(iter=burnin(mcmcp), burnin=0)
  mcmcChains(object) <- McmcChains(object, mcmc_burnin)
  object <- moveChain(object, 1)
  for(b in seq_len(burnin(mcmcp))){
    ## only store updates in the Posterior object slots; do not
    ## store anything in the chains
    object <- updateAll(object, move_chain=TRUE)
    object <- moveChain(object, b+1)
  }
  ## For recording prob(z) during the burnin, should we just use the
  ## z from the last iteration of the burnin?
  probz(object) <- probz(object)/rowSums(probz(object))
  object
}

runMcmc <- function(object, mcmcp){
  if(burnin(mcmcp) > 0) message("Burnin finished. Running additional MCMC simulations...")
  S <- 2:(savedIterations(mcmcp)+1)
  do_thin <- thin(mcmcp) > 1
  T <- seq_len(thin(mcmcp))
  for(s in S){
    object <- updateAll(object, TRUE, check_labels=checkLabels(mcmcp))
    object <- moveChain(object, s)
    ## update without moving chain
    if(do_thin){
      for(t in T) object <- updateAll(object, FALSE, check_labels=checkLabels(mcmcp))
    }
  }
  object
}


multipleStarts <- function(object, mcmcp){
  if(k(object)==1) return(object)
  message("Running ", nStarts(mcmcp), " chains")
  ##mp <- mcmcp[1]
  ##if(is(object, "BatchModel")){
  mp <- McmcParams(burnin=0, iter=nStartIter(mcmcp))
  mparams <- ModelParams("marginal", y=y(object), k=k(object),
                         mcmc.params=mp)
  mmod <- replicate(nStarts(mcmcp), initializeModel(mparams))
  mmodels <- suppressMessages(lapply(mmod, posteriorSimulation, mp))
  message("Selecting chain with largest log likelihood")
  lp <- sapply(mmodels, function(x) max(logpotentialc(x)))
  mmodel <- mmodels[[which.max(lp)]]
  if(is(object, "MarginalModel")) return(mmodel)
  params <- ModelParams("batch", y=y(object), k=k(object),
                        batch=batch(object),
                        mcmc.params=mcmcp)
  bmodel <- initializeBatchModel(params, z(mmodel))
  bmodel
}

#' @export
posteriorSimulation <- function(post, mcmcp){
  if(nStarts(mcmcp) > 1){
    post <- multipleStarts(post, mcmcp)
  }
  niter <- iter(mcmcp)
  ##
  ## Burn-in
  ##
  post <- runBurnin(post, mcmcp)
  if(niter==0) return(post)
  mcmcChains(post) <- McmcChains(post, mcmcp)
  ## make the first iteration in the stored chain the last iteration from the mcmc
  post <- moveChain(post, 1)
  ##
  ## if the chain moves quickly from the burnin, it might be related
  ## to the constrain imposed after burnin
  ##
  ##probz(post) <- .computeProbZ(post)
  ##
  ## Post burn-in
  ##
  post <- runMcmc(post, mcmcp)
  probz(post) <- probz(post)/(savedIterations(mcmcp)+1)
  ##if(!check_labels){
    ## compute the modes (otherwise, we have assume the modes are known)
  modes(post) <- computeModes(post)
  ##}
  ##  if(check_labels){
  ##
  ##  }
  ##  if(check_labels){
  ##    post <- switchLabels(post)
  ##    ##
  ##    ## We assume that computeModes provides estimates of the modes for the posterior
  ##    ## - for each subsequent iteration, we order the updates to minimize the distance to the modal estimates
  ##    ##
  ##    post <- runMcmc2(post, mcmcp, switch_labels=TRUE)
  ##  }
  orderTheta(post) <- getThetaOrder(post)
  post
}


setMethod("hist", "MixtureModel", function(x, ...){
  op <- par(las=1)
  L <- length(y(x))
  hist(y(x), breaks=L/10,
       col="gray", border="gray", main="", xlab="y",
       freq=FALSE)
  par(op)
})




.class_constructor <- function(class, ...){
  new(class, ...)
}

homozygousComponent <- function(y){
  ##mean(y < -1.5) > 0.005
  sum(y < -1.5) >= 3
}







setMethod("startingValues", "UnivariateMarginalModel", function(object){
  theta(object) <- median(y(object), na.rm=TRUE)
  sigma2(object) <- var(y(object), na.rm=TRUE)
  object
})



setMethod("computeMeans", "UnivariateMarginalModel", function(object){
  median(y(object), na.rm=TRUE)
})

setMethod("computeVars", "UnivariateMarginalModel", function(object){
  var(y(object), na.rm=TRUE)
})



setReplaceMethod("dataMean", "MixtureModel", function(object, value){
  object@data.mean <- value
  object
})

setReplaceMethod("dataPrec", "MixtureModel", function(object, value){
  object@data.prec <- value
  object
})

setMethod("mu.0", "MixtureModel", function(object) mu.0(hyperParams(object)))
setMethod("mu.0", "Hyperparameters", function(object) object@mu.0)

setMethod("tau2.0", "MixtureModel", function(object) tau2.0(hyperParams(object)))
setMethod("tau2.0", "Hyperparameters", function(object) object@tau2.0)

setMethod("eta.0", "MixtureModel", function(object) eta.0(hyperParams(object)))
setMethod("eta.0", "Hyperparameters", function(object) object@eta.0)

setMethod("m2.0", "MixtureModel", function(object) m2.0(hyperParams(object)))
setMethod("m2.0", "Hyperparameters", function(object) object@m2.0)


##updateWithPosteriorMeans <- function(object){
##  mc <- mcmcChains(object)
##  theta(object) <- colMeans(theta(mc))
##  sigma2(object) <- colMeans(sigma2(mc))
##  p(object) <- colMeans(p(mc))
##  mu(object) <- mean(mu(mc))
##  tau2(object) <- mean(tau2(mc))
##  nu.0(object) <- median(nu.0(mc))
##  sigma2.0(object) <- mean(sigma2.0(object))
##  ##...
##  logpotential(object) <- computePotential(object)
##  object
##}

simulateZ <- function(N, p){
  P <- rdirichlet(N, p)
  cumP <- t(apply(P, 1, cumsum))
  u <- runif(N)
  zz <- rep(NA, N)
  zz[u < cumP[, 1]] <- 1
  k <- 2
  while(k <= ncol(P)){
    zz[u < cumP[, k] & u >= cumP[, k-1]] <- k
    k <- k+1
  }
  zz
}

cumProbs <- function(p, k){
  pcum <- list()
  cols <- 2:(k-1)
  for(j in seq_along(cols)){
    g <- cols[j]
    pcum[[j]] <- rowSums(p[, 1:g, drop=FALSE])
  }
  pcum2 <- cbind(p[, 1], do.call(cbind, pcum), 1)
}

.updateZ <- function(p){
  ## generalize to any k, k >= 1
  cumP <- t(apply(p, 1, function(x) cumsum(x)))
  N <- nrow(p)
  u <- runif(N)
  zz <- rep(NA, N)
  zz[u < cumP[, 1]] <- 1
  k <- 2
  while(k <= ncol(p)){
    zz[u < cumP[, k] & u >= cumP[, k-1]] <- k
    k <- k+1
  }
  if(any(is.na(zz))) stop("missing values in zz")
  return(zz)
}

setMethod("updateZ", "MixtureModel", function(object){
  p <- posteriorMultinomial(object)
  ##zz <- simulateZ(length(y(object)), p)
  zz <- .updateZ(p)
  factor(zz, levels=seq_len(k(object)))
})



setMethod("initializeTheta", "numeric", function(object){
  means <- switch(paste0("k", object),
                  k1=0,
                  k2=c(-0.5, 0),
                  k3=c(-2, -0.5, 0),
                  k4=c(-2, -0.5, 0, 0.5),
                  k5=c(-2, -0.5, 0, 0.5, 1),
                  k6=c(-2, -0.5, -0.2, 0.2, 0.5, 1),
                  k7=c(-2, -0.5, -0.2, 0, 0.2, 0.5, 1),
                  NULL)
  if(is.null(means)) stop("k needs to be 1-7")
  means
})

ksTest <- function(object){
  B <- batch(object)
  yy <- y(object)
  ##ks <- rep(NA, choose(nBatch(object), 2))
  ks <- matrix(NA, choose(nBatch(object), 2), 4)
  i <- 1
  uB <- unique(B)
  for(j in seq_along(uB)){
    for(k in seq_along(uB)){
      if(k <= j) next()
      ##cat(j, k, sep="")
      ##cat(" ")
      b1 <- uB[j]
      b2 <- uB[k]
      stat <- suppressWarnings(ks.test(yy[B==b1], yy[B==b2]))
      ks[i, ] <- c(b1, b2, stat$statistic, stat$p.value)
      i <- i+1
    }
  }
  colnames(ks) <- c("batch1", "batch2", "stat", "pval")
  ks
}

##.collapseBatch <- function(object){
.collapseBatch <- function(yy, B){
  #B <- batch(object)
  uB <- unique(B)
  #yy <- y(object)
  ## One plate can pair with many other plates.
  for(j in seq_along(uB)){
    for(k in seq_along(uB)){
      if(k <= j) next() ## next k
      b1 <- uB[j]
      b2 <- uB[k]
      stat <- suppressWarnings(ks.test(yy[B==b1], yy[B==b2]))
      if(stat$p.value < 0.1) next()
      b <- paste(b1, b2, sep=",")
      B[B %in% b1 | B %in% b2] <- b
      ## once we've defined a new batch, return the new batch to the
      ## calling function
      return(B)
    }
  }
  B
}



##
## the batch names tend to be much too long
##
makeUnique <- function(x){
  ub <- unique(x)
  ##names(ub) <- ub
  abbrv <- setNames(make.unique(substr(ub, 1, 8)), ub)
  as.character(abbrv[x])
}


setMethod("collapseBatch", "BatchModel", function(object){
  N <- choose(nBatch(object), 2)
  cond2 <- TRUE
  while(N > 1 && cond2){
    cat('.')
    B <- batch(object)
    batch(object) <- .collapseBatch(y(object), batch(object))
    cond2 <- !identical(B, batch(object))
    N <- nBatch(object)
  }
  makeUnique(batch(object))
})

setMethod("collapseBatch", "SummarizedExperiment", function(object){
  N <- choose(length(unique(object$plate)), 2)
  cond2 <- TRUE
  while(N > 1 && cond2){
    cat('.')
    B <- object$plate
    object$plate <- .collapseBatch(copyNumber(object)[1,], object$plate)
    cond2 <- !identical(B, object$plate)
    N <- choose(length(unique(object$plate)), 2)
  }
  makeUnique(object$plate)
})


HardyWeinberg <- function(object){
  warning("requires cn inference...")
  tab <- table(map(object))
  ##freq <- setNames(table(zz), c("AA", "AB", "BB"))
  if(length(tab) == 1) return(NA)
  if(length(tab)==2){
    tab <- setNames(c(sort(tab, decreasing=TRUE), 0), c("AA", "AB", "BB"))
  }
  if(length(tab)==3){
    tab <- setNames(tab, c("AA", "AB", "BB"))
  }
  if(length(tab) > 3){
    tab <- setNames(tab[1:3], c("AA", "AB", "BB"))
  }
  HWChisq(tab)
}

setMethod("hwe", "MixtureModel", function(object) object@hwe)

map <- function(object) apply(probz(object), 1, which.max)

setMethod("fitMixtureModels", "numeric", function(object, mcmcp, K=1:5){
  message("Fitting ", length(K), " mixture models")
  fit <- vector("list", length(K))
  for(j in seq_along(K)){
    cat(".")
    mp <- mcmcp[j]
    kk <- K[j]
    params <- ModelParams("marginal", y=object, k=kk, mcmc.params=mp)
    modelk <- initializeModel(params)
    fit[[j]] <- suppressMessages(posteriorSimulation(modelk, mp))
  }
  fit
})

setMethod("thetac", "MixtureModel", function(object) theta(mcmcChains(object)))

setMethod("thetaMean", "MixtureModel", function(object) colMeans(thetac(object)))



setMethod("sigmaMean", "MixtureModel", function(object) colMeans(sigmac(object)))




logpotentialc <- function(object) logpotential(mcmcChains(object))

#' @export
sigmac <- function(object) sigma(mcmcChains(object))

#' @export
pic <- function(object) p(mcmcChains(object))

setMethod("pMean", "MixtureModel", function(object){
  colMeans(pic(object))
})

#' @export
muc <- function(object) mu(mcmcChains(object))

#' @export
muMean <- function(object) colMeans(muc(object))

#' @export
tauc <- function(object) sqrt(tau2(mcmcChains(object)))

#' @export
tauMean <- function(object) colMeans(tauc(object))

setMethod("modes", "MixtureModel", function(object) object@modes)

setReplaceMethod("modes", "MixtureModel", function(object, value) {
  object@modes <- value
  object
})
