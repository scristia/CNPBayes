setMethod("hyperParams", "MixtureModel", function(object) object@hyperparams)

setReplaceMethod("hyperParams", c("MixtureModel", "Hyperparameters"),
                 function(object, value) {
                   object@hyperparams <- value
                   object
                 })


setMethod("batch", "MixtureModel", function(object) object@batch)


observed <- function(object) object@data

#' @export
theta <- function(object) object@theta

#' @export
sigma2 <- function(object) object@sigma2

#' @export
sigma <- function(object) sqrt(object@sigma2)

#' @export
mu <- function(object) object@mu

#' @export
tau2 <- function(object) object@tau2

#' @export
tau <- function(object) sqrt(object@tau2)

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

setReplaceMethod("mu", "MixtureModel", function(object, value){
  object@mu <- value
  object
})

setReplaceMethod("tau2", "MixtureModel", function(object, value){
  object@tau2 <- value
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
  cat("     theta (s)   :", "\t", mns, "\n")
  cat("     sigma2.0 (s):", round(sigma2.0(object), 2), "\n")
  cat("     nu.0 (s)    :", nu.0(object), "\n")
  cat("     tau2 (s)    :", round(tau2(object), 2),  "\n")
  cat("     mu (s)      :", round(mu(object), 2), "\n")
  cat("     potential(s):", round(logpotential(object), 2), "\n")
})

tablez <- function(post) table(z(post))


updateAll <- function(post, move_chain, s){
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
  hypp <- hyperParams(post)
  alpha.n <- alpha(hypp) + tablez(post)
  p(post) <- as.numeric(rdirichlet(1, alpha.n))
  ##
  ##  update auxillary variable
  ##
  z(post) <- updateZ(post)
  ##
  ##  update (nu.0, sigma2.0)
  ##
  sigma2.0(post) <- updateSigma2.0(post)
  nu.0(post) <- updateNu.0(post)
  if(move_chain) {
    logpotential(post) <- computePotential(post)
    post <- moveChain(post, s+1)
  }
  post
}


#' @export
posteriorSimulation <- function(post, mcmcp){
  niter <- iter(mcmcp)
  ## Count the initial values as the first iteration
  S <- seq_len(niter-1)
  if(length(S) > 0){
    save_iter <- seq(1, niter-1, thin(mcmcp))
    move_chain <- S %in% save_iter
  }
  ##
  ## Record initial values
  ##
  post <- moveChain(post, 1)
  ##
  ## Burn-in
  ##
  if(burnin(mcmcp) > 0){
    for(b in 1:burnin(mcmcp)){
      ## only store updates in the Posterior object slots; do not
      ## store anything in the chains
      post <- updateAll(post, move_chain=FALSE, b)
    }
  }
  if(k(post)==1) return(post)
  ##
  ## Post burn-in
  ##
  for(s in S){
    post <- updateAll(post, move_chain[s], s)
  }
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



constructModel <- function(data, k, batch){
  nbatch <- length(unique(batch))
  if(k > 1){
    if(nbatch > 1){
      object <- BatchModel(data, k, batch)
    } else{
      object <- MarginalModel(data, k, batch)
    }
  } else {
    if(nbatch > 1){
      stop("constructModel for this class not specified")
    } else  {
      object <- UnivariateMarginalModel(data, k , batch)
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
initializeModel <- function(yy, k=2, batch){
  if(missing(batch)) batch <- rep(1, length(yy))
  object <- constructModel(data=yy, k=k, batch=batch)
  hypp <- hyperParams(object) <- Hyperparameters(k=k)
  ##
  ## effects precision of component means
  ##
  object <- startingValues(object)
  ##
  ## simulate from priors
  ##
  tau2(object) <- initializeTau2(object)
  p(object) <- as.numeric(rdirichlet(1, alpha(hypp))) ## rows are platform, columns are components
  ## initalize sigma2.0 as weighted average of component variances
  sigma2.0(object) <- initializeSigma2.0(object)
  nu.0(object) <- rgeom(1, beta(hypp))
  z(object) <- updateZ(object)
  mu(object) <- updateMu(object)

  ##
  ## given starting values, simulate z
  ##
  dataMean(object) <- computeMeans(object)
  dataPrec(object) <- 1/computeVars(object)
  logpotential(object) <- computePotential(object)
  object
}

updateWithPosteriorMeans <- function(object){
  mc <- mcmcChains(object)
  theta(object) <- colMeans(theta(mc))
  sigma2(object) <- colMeans(sigma2(mc))
  p(object) <- colMeans(p(mc))
  mu(object) <- mean(mu(mc))
  tau2(object) <- mean(tau2(mc))
  nu.0(object) <- median(nu.0(mc))
  sigma2.0(object) <- mean(sigma2.0(object))
  ##...
  logpotential(object) <- computePotential(object)
  object
}

#' @export
setMethod("BIC", "MarginalModel", function(object, ...){
  if(k(object) > 1){
    object <- updateWithPosteriorMeans(object)
  }
  K <- 3*2*k(object)
  -2*logpotential(object) + K*(log(length(y(object))) - log(2*pi))
})



#' @export
simulateData <- function(N=2500, .k=3, .theta, .sigma, .batch, .alpha){
  if(missing(.alpha)) .alpha <- rep(1, .k)
  if(missing(.batch)) .batch <- rep(1, N)
  yy <- rnorm(N)
  object <- initializeModel(yy, k=.k, batch=.batch)
  hypp <- hyperParams(object)
  alpha(hypp) <- .alpha
  hyperParams(object) <- hypp
  ##
  theta(object) <- .theta
  mu(object) <- initializeMu(object)
  tau2(object) <- initializeTau2(object)
  sigma2(object) <- .sigma^2
  ##p(post) <- as.numeric(rdirichlet(1, alpha.n))
  p(object) <- as.numeric(rdirichlet(1, alpha(hypp)))
  pisum <- cumsum(p(object))
  u <- runif(N, 0, 1)
  ## Simulate z according to pi
  zz <- ifelse(u < pisum[1], 1, ifelse(u < pisum[2], 2, 3))
  z(object) <- factor(zz, levels=seq_len(k(hypp)))
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

.updateZ <- function(p, k){
  ## generalize to any k, k >= 1
  pcum <- list()
  cols <- 2:(k-1)
  u <- runif(nrow(p), 0, 1)
  for(j in seq_along(cols)){
    g <- cols[j]
    pcum[[j]] <- rowSums(p[, 1:g, drop=FALSE])
  }
  pcum2 <- cbind(p[, 1], do.call(cbind, pcum), 1)
  z2 <- rep(NA, nrow(p))
  z2[u <= pcum2[, 1]] <- 1
  z2[u > pcum2[, ncol(pcum2)-1]] <- k
  for(j in 1:(ncol(p)-2)){
    z2[u > pcum2[, j] & u <= pcum2[, j+1]] <- j+1
  }
  factor(z2, levels=seq_len(k))
}

setMethod("updateZ", "MixtureModel", function(object){
  p <- posteriorMultinomial(object)
  .updateZ(p, k(object))
})
