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

#' @export
nu.0 <- function(object) object@nu.0

#' @export
sigma2.0 <- function(object) object@sigma2.0

#' @export
y <- function(object) object@data

#' @export
z <- function(object) object@z

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


updateAll <- function(post, move_chain, s, constrainTheta=TRUE){
  ##
  ## sample new value of (theta_h, sigma2_h)
  ##   - note these are independent in the prior
  ##
  theta(post) <- updateTheta(post, constrain=constrainTheta)
  ##theta(post) <- updateTheta(post, constrain=s)
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
  post <- moveChain(post, s)
  post
}

setReplaceMethod("probz", "MixtureModel", function(object, value){
  object@probz <- value
  object
})

setMethod("probz", "MixtureModel", function(object) object@probz)

#' @export
posteriorSimulation <- function(post, mcmcp){
  niter <- iter(mcmcp)
  ##
  ## Record initial values
  ##
  post <- moveChain(post, 1)
  post2 <- post
  ##
  ## Burn-in
  ##
  if( burnin(mcmcp) > 0){
    for(b in 1:burnin(mcmcp)){
      ## only store updates in the Posterior object slots; do not
      ## store anything in the chains
      post <- updateAll(post, move_chain=FALSE, b, constrainTheta=FALSE)
    }
    ##
    ## The assumption here is that the ordering of the last update is
    ## sufficient.  This could be more robust by using posterior
    ## means.
    ##
    if(!identical(theta(post), sort(theta(post)))) {
      post <- reorderComponents(post)
    }
    ## make the first iteration in the stored chain the last iteration from the mcmc
    moveChain(post, 1)
  }
  ##browser()
  if(niter==0) return(post)
  ##browser()
  if(FALSE){
    table(z(post2), z(post))
    cbind(theta(post2), theta(post))
  }
  probz(post) <- matrix(0, length(y(post)), k(post))
  ##
  ## Post burn-in
  ##
  ## By starting at 2, the last iteration from the burnin or (if no
  ## burnin) the starting values are the first element in the chains
  S <- 2:savedIterations(mcmcp)
  ##browser()
  for(s in S){
    post <- updateAll(post, TRUE, s, constrainTheta=constrainTheta(mcmcp))
    if(thin(mcmcp) > 1){
      for(t in seq_len(thin(mcmcp))){
        post <- updateAll(post, FALSE, s, constrainTheta=constrainTheta(mcmcp))
      }
    }
  }
  probz(post) <- probz(post)/(savedIterations(mcmcp))
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


initializeModel <- function(params, .alpha){ ##.theta, .sigma){
  if(missing(.alpha)) .alpha <- rep(1, k(params))
  object <- constructModel(type(params), data=y(params), k=k(params), batch=batch(params))
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
  dataPrec(object) <- tau2(object)
  mu(object) <- initializeMu(object)
  tau2(object) <- initializeTau2(object)
  dataMean(object) <- mu(object)
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
      stat <- ks.test(yy[B==b1], yy[B==b2])
      ks[i, ] <- c(b1, b2, stat$statistic, stat$p.value)
      i <- i+1
    }
  }
  colnames(ks) <- c("batch1", "batch2", "stat", "pval")
  ks
}

.collapseBatch <- function(object){
  B <- batch(object)
  uB <- unique(B)
  yy <- y(object)
  ## One plate can pair with many other plates.
  for(j in seq_along(uB)){
    for(k in seq_along(uB)){
      if(k <= j) next() ## next k
      b1 <- uB[j]
      b2 <- uB[k]
      stat <- ks.test(yy[B==b1], yy[B==b2])
      if(stat$p.value < 0.1) next()
      b <- paste(b1, b2, sep=",")
      B[B %in% b1 | B %in% b2] <- b
      ## once we've defined a new batch, the unique batches has
      ## changed and we must exit this function
      return(B)
    }
  }
  B
}

#' @export
collapseBatch <- function(object, mcmcp){
  ##ks <- ksTest(object)
  ##b <- .collapseBatch(model, ks)
  i <- 0
  N <- choose(nBatch(object), 2)
  cond2 <- TRUE
  while(N > 1 && cond2){
    cat('.')
    B <- batch(object)
    batch(object) <- .collapseBatch(object)
    cond2 <- !identical(B, batch(object))
    N <- nBatch(object)
    i <- i+1
  }
  params <- ModelParams("batch", y=y(object),
                        k=k(object),
                        batch=batch(object),
                        mcmc.params=mcmcp)
  object <- initializeModel(params)
  object
}


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
