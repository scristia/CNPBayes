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

setGeneric("showMeans", function(object) standardGeneric("showMeans"))
setGeneric("showSigmas", function(object) standardGeneric("showSigmas"))

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

tablez <- function(post){
  zz <- map(post)
  ##
  ## This is safer than doing table(map(post))
  ##
  tab <- rep(NA, k(post))
  for(j in seq_len(k(post))){
    tab[j] <- sum(zz == j)
  }
  tab
}


updateAll <- function(post, move_chain, constrainTheta=TRUE){
  ##
  ## sample new value of (theta_h, sigma2_h)
  ##   - note these are independent in the prior
  ##
  theta(post) <- updateTheta(post, constrain=constrainTheta)
  sigma2(post) <- updateSigma2(post)
  ##
  ## update (mu, tau2)
  ##
  mu(post) <- updateMu(post)
  tau2(post) <- updateTau2(post)
  ##
  ## update pi
  ##
  hypp <- hyperParams(post)
  ##alpha.n <- alpha(hypp) + tablez(post)
  alpha.n <- alpha(hypp) + table(z(post))
  mixprobs <- as.numeric(rdirichlet(1, alpha.n))
  p(post) <- mixprobs
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
  if(move_chain){
    pZ <- probz(post)
    zz <- as.integer(z(post))
    for(j in seq_len(k(post))){
      pZ[, j] <- pZ[, j] + as.integer(zz==j)
    }
    probz(post) <- pZ
  }
  post
}

setReplaceMethod("probz", "MixtureModel", function(object, value){
  object@probz <- value
  object
})

setMethod("probz", "MixtureModel", function(object) object@probz)

runBurnin <- function(object, mcmcp){
  if(burnin(mcmcp)==0){
    object <- moveChain(object, 1)
    return(object)
  }
  mcmc_burnin <- McmcParams(iter=burnin(mcmcp), burnin=0)
  mcmcChains(object) <- McmcChains(object, mcmc_burnin)
  object <- moveChain(object, 1)
  for(b in seq_len(burnin(mcmcp))){
    ## only store updates in the Posterior object slots; do not
    ## store anything in the chains
    object <- updateAll(object, move_chain=TRUE, constrainTheta=FALSE)
    object <- moveChain(object, b+1)
  }
  ## Reorder the component indices by the iteration during burnin that
  ## had the largest potential.  This requires storing the chain
  ## during burin.  Reordering is important for the sigma2s (i.e,
  ## sigma2_j = sigma2_j' for all j and j' greater than 1 browser()
  ## object <- sort(object)
  ##mcmcChains(object) <- McmcChains(object, mcmcp)
  object
}

runMcmc <- function(object, mcmcp){
  S <- 2:(savedIterations(mcmcp)+1)
  do_thin <- thin(mcmcp) > 1
  T <- seq_len(thin(mcmcp))
  do_constrain <- constrainTheta(mcmcp)
  for(s in S){
    object <- updateAll(object, TRUE, constrainTheta=do_constrain)
    object <- moveChain(object, s)
    ## update without moving chain
    if(do_thin){
      for(t in T) object <- updateAll(object, FALSE, constrainTheta=do_constrain)
    }
  }
  object
}

#' @export
posteriorSimulation <- function(post, mcmcp){
  ##mcmcChains(post) <- McmcChains(post, mcmcp)
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
  probz(post) <- matrix(0, length(y(post)), k(post))
  ##
  ## Post burn-in
  ##
  post <- runMcmc(post, mcmcp)
  probz(post) <- probz(post)/(savedIterations(mcmcp)+1)
  post <- sort(post)
  post
}


setMethod("hist", "MixtureModel", function(x, ...){
  op <- par(las=1)
  hist(y(x), breaks=1/10*length(y(x)),
       col="gray", border="gray", main="", xlab="y",
       freq=FALSE)
  par(op)
})




.class_constructor <- function(class, ...){
  new(class, ...)
}



constructModel <- function(type, data, k, batch){
  nbatch <- length(unique(batch))
  if(k > 1){
    if(type=="batch"){
      object <- BatchModel(data, k, batch)
    } else{
      object <- MarginalModel(data, k, batch)
    }
  } else {
    if(type=="batch"){
      object <- UnivariateBatchModel(data, 1, batch)
    } else  {
      object <- UnivariateMarginalModel(data, 1 , batch)
    }
  }
  object
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

setMethod("computeVars", "MarginalModel", function(object){
  sapply(split(y(object), z(object)), var, na.rm=TRUE)
})

setReplaceMethod("dataMean", "MixtureModel", function(object, value){
  object@data.mean <- value
  object
})

setReplaceMethod("dataPrec", "MixtureModel", function(object, value){
  object@data.prec <- value
  object
})

#' @export
muList <- function(k){
  switch(paste0("k", k),
         k1=list(-2, -01, -0.5, 0, 0.5, 1),
         k2=list(c(-2, -0.5),
           c(-2, 0),
           c(-1, 0),
           c(-0.5, 0),
           c(0, 0.5),
           c(0.5, 1)),
         k3=list(c(-3, -0.5, 0),
           c(-2, -0.4, 0),
           c(-1, -0.5, 0),
           c(-0.5, 0, 0.5),
           c(0, 0.5, 1),
           c(-3, 0, 0.5)),
         k4=list(c(-3, -0.5, 0, 0.5),
           c(-2, -0.4, 0, 0.5),
           c(-1, -0.5, 0, 0.5),
           c(-0.5, 0, 0.5, 1),
           c(-3, 0, 0.5, 1),
           c(-3, -1, -0.4, 0)),
         k5=replicate(6, sort(rnorm(5, mean=c(-2.5, -0.4, 0, 0.4, 1), sd=0.2)),
           simplify=FALSE),
         k6=replicate(6, sort(rnorm(6, mean=c(-2.5, -0.4, -0.1, 0.10, .4, 1), sd=0.2)),
           simplify=FALSE)
         )
}

initialTau2 <- function(k) rep(0.2, k)


#' @export
##initializeModel <- function(type, yy, k=2, batch, .alpha, .theta, .sigma){
##setGeneric("initializeDataMean", function(object) standardGeneric("initializeDataMean"))
##setMethod("initilalizeDataMean", "MarginalModel", function(object){
##  mu(object)
##})
##
##setMethod("initilalizeDataMean", "MarginalModel", function(object){
##  mu(object)
##})


initializeModel <- function(params, .alpha){ ##.theta, .sigma){
  if(missing(.alpha)) .alpha <- rep(1, k(params))
  object <- constructModel(type(params), data=y(params), k=k(params),
                           batch=batch(params))
  hypp <- hyperParams(object) <- Hyperparameters(type(params), k=k(params))
  ##  if(missing(.sigma)) .sigma <- rep(mad(y(params), na.rm=TRUE), k(params))
  ##  if(missing(.theta)) .theta <- mu(object)
  ##
  ## effects precision of component means
  ##
  object <- startingValues(object)
  ##
  ## simulate from priors
  ##
  tau2(object) <- initializeTau2(object) ## redundant?
  p(object) <- as.numeric(rdirichlet(1, alpha(hypp))) ## rows are platform, columns are components
  ## initalize sigma2.0 as weighted average of component variances
  sigma2.0(object) <- initializeSigma2.0(object)
  nu.0(object) <- max(rgeom(1, betas(hypp)), 1)
  ##dataPrec(object) <- tau2(object)
  dataPrec(object) <- 1/sigma2(object)
  mu(object) <- initializeMu(object)
  tau2(object) <- initializeTau2(object)
  dataMean(object) <- theta(object)
  z(object) <- updateZ(object)
  ##
  ##
  ##
  logpotential(object) <- computePotential(object)
  ##
  ## given starting values, simulate z
  ##
  hypp <- hyperParams(object)
  alpha(hypp) <- .alpha
  hyperParams(object) <- hypp
  ##
  alpha.n <- alpha(hypp) + table(z(object))
  p(object) <- as.numeric(rdirichlet(1, alpha.n))
  pisum <- cumsum(p(object))
  u <- runif(N(params), 0, 1)
  ## Simulate z according to pi
  ##zz <- ifelse(u < pisum[1], 1, ifelse(u < pisum[2], 2, 3))
  ##z(object) <- factor(zz, levels=seq_len(k(hypp)))
  mcmcChains(object) <- McmcChains(object, mcmcParams(params))
  object
}


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

#' @export
simulateBatchData <- function(N=2500, .k=3, .batch, .alpha, means, sds){
  if(missing(.alpha)) .alpha <- rep(1, .k)
  if(missing(.batch)) .batch <- as.character(rep(1, N))
  yy <- rnorm(N)
  params <- ModelParams("batch", y=yy, k=.k, batch=.batch)
  object <- constructModel(type(params), data=y(params),
                           k=k(params), batch=batch(params))
  hypp <- hyperParams(object) <- Hyperparameters(type(params), k=k(params))
  object <- initializeModel(params, .alpha=.alpha)
  theta(object) <- means
  sigma2(object) <- sds^2
  rownames(sigma2(object)) <- rownames(theta(object)) <- unique(.batch)
  ##
  ## Simulate data | Z
  ##
  hypp <- hyperParams(object)
  stopif(is.na(a(hypp)))
  stopif(is.na(b(hypp)))
  stopif(is.na(nu.0(object)))
  dat(object) <- simulateY(object)
  dataMean(object) <- computeMeans(object)
  dataPrec(object) <- 1/computeVars(object)
  logpotential(object) <- computePotential(object)

  zz <- z(object)
  tabz <- table(zz)
  p(object) <- as.numeric(tabz/length(zz))
  object
}


#' @export
simulateData <- function(N=2500, .k=3, .batch, .alpha, means, sds){
  if(missing(.alpha)) .alpha <- rep(1, .k)
  if(missing(.batch)) .batch <- as.character(rep(1, N))
  yy <- rnorm(N)
  params <- ModelParams("marginal", y=yy, k=.k, batch=.batch)
  object <- initializeModel(params, .alpha=.alpha)
  ##
  ## Simulate data | Z
  ##
  theta(object) <- means
  sigma2(object) <- sds^2
  hypp <- hyperParams(object)
  stopif(is.na(a(hypp)))
  stopif(is.na(b(hypp)))
  stopif(is.na(nu.0(object)))

  dat(object) <- simulateY(object)
  dataMean(object) <- computeMeans(object)
  dataPrec(object) <- 1/computeVars(object)
  logpotential(object) <- computePotential(object)

  tabz <- table(z(object))
  p(object) <- as.numeric(tabz/length(z(object)))
  object
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

.updateZ <- function(p, k){
  ## generalize to any k, k >= 1
  pcum <- cumProbs(p, k)
  z2 <- rep(NA, nrow(p))
  u <- runif(nrow(p), 0, 1)
  ## update z for first and last component
  z2[u <= pcum[, 1]] <- 1
  z2[u > pcum[, ncol(pcum)-1]] <- k
  ## update z for components 2, ..., k-1
  for(j in 1:(ncol(p)-2)){
    z2[u > pcum[, j] & u <= pcum[, j+1]] <- j+1
  }
  factor(z2, levels=seq_len(k))
}

setMethod("updateZ", "MixtureModel", function(object){
  p <- posteriorMultinomial(object)
  .updateZ(p, k(object))
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

#' @export
setGeneric("collapseBatch", function(object) standardGeneric("collapseBatch"))

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
  fit <- foreach(j = seq_along(K)) %do% {
    cat(".")
    kk <- K[j]
    params <- ModelParams("marginal", y=object, k=kk, mcmc.params=mcmcp)
    modelk <- initializeModel(params)
    modelk <- posteriorSimulation(modelk, mcmcp)
    modelk
  }
  fit
})

#' @export
thetac <- function(object) theta(mcmcChains(object))

#' @export
sigmac <- function(object) sigma(mcmcChains(object))

#' @export
pic <- function(object) p(mcmcChains(object))

muc <- function(object) mu(mcmcChains(object))
