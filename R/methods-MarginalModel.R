MarginalModel <- function(data, k, batch){
  new("MarginalModel",
      hyperparams=HyperparametersMarginal(k=k),
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
      z=factor(numeric(k)),
      probz=matrix(0, length(data), k),
      logpotential=numeric(1),
      mcmc.chains=McmcChains(),
      batch=batch,
      hwe=numeric(),
      modes=list())
}

UnivariateMarginalModel <- function(data, k=1, batch){
  new("UnivariateMarginalModel",
      hyperparams=HyperparametersMarginal(k=k),
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
      z=factor(numeric(k)),
      probz=matrix(1, length(data), 1),
      logpotential=numeric(1),
      mcmc.chains=McmcChains(),
      batch=batch,
      hwe=numeric(),
      modes=list())
}

#' @export
setMethod("mu", "MarginalModel", function(object) object@mu)

#' @export
setMethod("tau2", "MarginalModel", function(object) object@tau2)

setMethod("updateZ", "UnivariateMarginalModel", function(object){
  factor(rep(1, length(y(object))))
})


setMethod("computePotential", "MarginalModel", function(object){
  hypp <- hyperParams(object)
  K <- k(hypp)
  yy <- y(object)
  zz <- z(object)
  thetas <- theta(object)
  sigma2s <- sigma2(object)
  mus <- mu(object)
  tau2s <- tau2(object)
  pp <- p(object)
  ##p.mu <- dnorm(mus, mu.0(hypp), sqrt(tau2.0(hypp)))
  ##p.mu <- 1
  p.sigma2.0 <- dgamma(sigma2.0(object), a(hypp), b(hypp))
  p.nu.0 <- dgeom(nu.0(object), betas(hypp))
  p.theta <- dnorm(thetas, mus, sqrt(tau2s))
  p.sigma2 <- dgamma(1/sigma2s, 1/2*nu.0(object), 1/2*nu.0(object)*sigma2.0(object))
  pot <- list()
  for(i in seq_len(K)){
    pot[[i]] <- pp[i]*dnorm(yy, thetas[i], sqrt(sigma2s)[i])
  }
  pot <- do.call("cbind", pot)
  pot <- rowSums(pot)
  total_pot <- sum(log(pot)) + sum(log(p.theta)) + sum(log(p.sigma2)) + log(p.nu.0) + log(p.sigma2.0)
  total_pot
})

##
## For the marginal model, mu and tau2 are hyper-parameters.  There is
## no update.
##
setMethod("updateMu", "MarginalModel", function(object){
  ##tau2.0, tau2, eta.0, mu.0, thetas, k, nn){
  ##hypp <- hyperParams(object)
  ##  mu.k <- .updateMu(tau2.0(hypp), tau2(object), k(object), z(object),
  ##                    theta(object), mu.0(hypp))
  ##  mu.k
  .updateMu(object)
})

.updateMu <- function(object){
  hypp <- hyperParams(object)
  tau2.0.tilde <- 1/tau2.0(hypp)
  tau2.tilde <- 1/tau2(object)
  tau2.k.tilde <- tau2.0.tilde + k(object)*tau2.tilde
  ##browser()
  nn <- tablez(object)
  theta.bar <- sum(nn*theta(object))/sum(nn)
  mu.k <- tau2.0.tilde/(tau2.0.tilde + k(object)*tau2.tilde)*mu.0(hypp) +
    k(object)*tau2.tilde/(tau2.0.tilde + k(object)*tau2.tilde)*theta.bar
  if(any(is.na(mu.k))) stop("NAs in mu update")
  if(length(mu.k) > 1) stop("length > 1")
  mu.k
}

setMethod("startingValues", "MarginalModel", function(object){
  hypp <- hyperParams(object)
##  tmp.file <- tempfile()
##  sink(tmp.file)
##  mmfit <- normalmixEM(y(object), arbvar = FALSE,
##                       epsilon = 1e-03, k=k(hypp), maxit=2000)
##  sink()
##  unlink(tmp.file)
##  mus <- mmfit$mu
##  vars <- (mmfit$sigma[order(mmfit$mu)])^2
##  if(FALSE){
##    if(k(object) >= 3 && any(y(object) < - 1.5)) {
##      ## set first component to be near the expected mean for homozygous deletions
##      mus[1] <- median(y(object)[y(object) < -1.5])
##      ## set first component to be near the expected mean for hemizygous deletions
##      mus[2] <- median(y(object)[y(object) > -1.5 & y(object) < -0.15])
##      ##
##      ## the variance may be overestimated
##      vars[1] <- var(y(object)[y(object) < -1.5])
##      vars[2:k(object)] <- var(y(object)[y(object) > -1.5 & y(object) < -0.15])
##      ## initialize using the posterior
##    }
  ##  }
  mu(object) <- rnorm(1, mu.0(hypp), tau.0(hypp))
  tau2.tilde <- rgamma(1, shape=1/2*nu.0(object), rate=1/2*nu.0(object)*m2.0(hypp))
  tau2(object) <- 1/tau2.tilde
  theta(object) <- rnorm(k(object), mu(object), tau(object))
  sigma2(object) <- rep(tau2(object), k(object))
  object
})



setMethod("initializeSigma2.0", "MarginalModel", function(object){
  hypp <- hyperParams(object)
  sum(alpha(hypp)*sigma2(object))/sum(alpha(hypp))
})


## just choose a big number
setMethod("initializeTau2", "MarginalModel", function(object)  1000)


##setMethod("initializeMu", "MarginalModel", function(object){
##  means <- switch(paste0("k", k(object)),
##                  k1=0,
##                  k2=c(-0.5, 0),
##                  k3=c(-2, -0.5, 0),
##                  k4=c(-2, -0.5, 0, 0.5),
##                  k5=c(-2, -0.5, 0, 0.5, 1),
##                  k6=c(-2, -0.5, -0.2, 0.2, 0.5, 1),
##                  k7=c(-2, -0.5, -0.2, 0, 0.2, 0.5, 1),
##                  NULL)
##  if(is.null(means)) stop("k needs to be 1-7")
##  means
##})


setMethod("posteriorMultinomial", "MarginalModel", function(object){
  .multinomial_probs <- .posteriorMultinomial(y(object),
                                              theta(object),
                                              sqrt(sigma2(object)),
                                              p(object))

})

.posteriorMultinomial <- function(y, theta, sd, pi){
  K <- seq_len(length(pi))
  numerator <- vector("list", length(K))
  for(j in K){
    numerator[[j]] <- pi[j]*dnorm(y, theta[j], sd[j])
  }
  numerator <- do.call(cbind, numerator)
  denominator <- rowSums(numerator)
  numerator/denominator
}

setMethod("show", "MarginalModel", function(object) callNextMethod())

setMethod("computeMeans", "MarginalModel", function(object){
  means <- sapply(split(y(object), z(object)), mean, na.rm=TRUE)
  if(any(is.nan(means))) {
    means[is.nan(means)] <- rnorm(sum(is.nan(means)), mu(object), tau(object))
  }
  means
})

setMethod("computeVars", "MarginalModel", function(object){
  vars <- sapply(split(y(object), z(object)), var, na.rm=TRUE)
  if(any(is.nan(vars))){
    vars[is.nan(vars)] <- 1/rgamma(sum(is.nan(vars)), shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  }
  vars
})



#' @export
setMethod("plot", "MarginalModel", function(x, y, use.current=FALSE, ...){
  hist(x)
  pi <- p(x)
  xx <- seq(min(observed(x)), max(observed(x)),  length.out=10e3)
  mc <- mcmcChains(x)
  if(!use.current){
    thetas <- colMeans(theta(mc))
    sds <- colMeans(sqrt(sigma2(mc)))
  } else {
    ## use current value
    thetas <- theta(x)
    sds <- sqrt(sigma2(x))
  }
  cols <- brewer.pal(max(length(pi), 3),  "Set1")
  for(j in seq_along(pi)){
    lines(xx, pi[j]*dnorm(xx,
                          mean=thetas[j],
                          sd=sds[j]),
          col=cols[j], lwd=2)
  }
})

setMethod("simulateY", "MarginalModel", function(object){
  zz <- z(object)
  yy <- rnorm(length(zz), mean=theta(object)[zz], sd=sigma(object)[zz])
})

setMethod("moveChain", "MarginalModel", function(object, s){
  mcmc <- mcmcChains(object)
  theta(mcmc)[s, ] <- as.numeric(theta(object))
  sigma2(mcmc)[s, ] <- as.numeric(sigma2(object))
  p(mcmc)[s, ] <- p(object)
  ##
  mu(mcmc)[s] <- mu(object)
  tau2(mcmc)[s] <- tau2(object)
  ##
  nu.0(mcmc)[s] <- nu.0(object)
  sigma2.0(mcmc)[s] <- sigma2.0(object)
  logpotential(mcmc)[s] <- logpotential(object)
  logLik(mcmc)[s] <- computeLoglik(object)
  mcmcChains(object) <- mcmc
  object
})

setMethod("updateThetaCpp", "MarginalModel", function(object, constrain) {
  .Call("update_theta", object, constrain=constrain)
})
#
#.updateThetaCpp <- function(obj, constrain=TRUE) {
#    theta <- .Call("update", obj, constrain)
#    theta
#}


setMethod("updateTheta", "MarginalModel", function(object) .updateTheta(object))

.updateTheta <- function(object){
  theta.last <- theta(object)
  tau2.tilde <- 1/tau2(object)
  sigma2.tilde <- 1/sigma2(object)
  n.h <- tablez(object)
  n.h <- pmax(n.h, 1)
  tau2.n.tilde <- tau2.tilde + n.h*sigma2.tilde
  tau2.n <- 1/tau2.n.tilde
  denom <- tau2.tilde + n.h*sigma2.tilde
  w1 <- tau2.tilde/denom
  w2 <- n.h*sigma2.tilde/denom
  mu.n <- w1*mu(object) + w2*dataMean(object)
  thetas <- rnorm(k(object), mu.n, sqrt(tau2.n))
  if(any(is.na(thetas))) stop("NAs in theta update")
  if(length(thetas) != length(theta.last)) stop("check thetas")
  thetas
}

setMethod("updateSigma2Cpp", "MarginalModel", function(object) {
  .Call("update_sigma2", object)
})

setMethod("updateSigma2", "MarginalModel", function(object) {
  .updateSigma2(object)
  ##.updateSigma2_2(object)
})


##.updateSigma2 <- function(data.list, thetas, nu.0, sigma2.0, n.h){
.updateSigma2 <- function(object){
  data.list <- split(y(object), z(object))
  thetas <- theta(object)
  nu.0 <- nu.0(object)
  sigma2.0 <- sigma2.0(object)
  n.h <- tablez(object)
  n.h <- pmax(n.h, 1)

  nu.n <- nu.0+n.h

  zz <- z(object)
  m <- thetas[as.integer(zz)]
  squares <- (y(object) - m)^2
  ss <- sapply(split(squares, zz), sum)

  ## weighted average of sums of squares
  sigma2.nh <- 1/nu.n*(nu.0*sigma2.0 + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- rgamma(k(object), shape=shape, rate=rate)
  sigma2.h <- 1/sigma2.h.tilde
  stopif(any(is.nan(sigma2.h)))
  sigma2.h
}

## variance for componets k > 1 assumed to be the same
.updateSigma2_2 <- function(object){
  nz <- nonZeroCopynumber(object)
  if(length(unique(nz)) ==1){
    ## guard against zeros
    nz[which.min(y(object))] <- 0
  }
  n.h <- table(nz)
  nu.n <- nu.0(object)+n.h

  thetas <- theta(object)


  zz <- z(object)
  m <- thetas[as.integer(zz)]
  squares <- (y(object) - m)^2
  ss <- sapply(split(squares, nz), sum)

  ## weighted average of sums of squares
  sigma2.nh <- 1/nu.n*(nu.0(object)*sigma2.0(object) + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- rgamma(2, shape=shape, rate=rate)
  ##tmp <- rgamma(1000, shape=1/2*nu.n[1], rate=1/2*nu.n[1]*sigma2.nh[1])
  sigma2.h <- 1/sigma2.h.tilde
  ##stopif(any(is.nan(sigma2.h)))
  s2 <- c(sigma2.h[1], rep(sigma2.h[2], k(object)-1))
  s2
}


setMethod("updateSigma2.0", "MarginalModel", function(object){
  hypp <- hyperParams(object)
  .updateSigma2.0(a=a(hypp), b=b(hypp), nu.0=nu.0(object),
                  sigma2.h=sigma2(object), k=k(hypp))
})

.updateSigma2.0 <- function(a, b, nu.0, sigma2.h, k){
  a.k <- a+1/2*k*nu.0
  b.k <- b+1/2*sum(1/sigma2.h)
  sigma2.0 <- rgamma(1, shape=a.k, rate=b.k)
  stopif(is.nan(sigma2.0))
  sigma2.0
}

setMethod("updateNu.0", "MarginalModel", function(object){
  hypp <- hyperParams(object)
  .updateNu.0(beta=betas(hypp), sigma2.0=sigma2.0(object), sigma2.h=sigma2(object),
              nu.0=nu.0(object), k=k(object))
})

.updateNu.0 <- function(NUMAX=100, beta, sigma2.0, sigma2.h, nu.0, k){
  x <- seq_len(NUMAX)
  lpnu0 <- k * (0.5 * x * log(sigma2.0 * x/2)-lgamma(x/2)) +
    (x/2 - 1) * sum(log(1/sigma2.h)) +
      -x * (beta + 0.5 * sigma2.0 * sum(1/sigma2.h))
  prob <- exp(lpnu0 - max(lpnu0))
  nu0 <- sample(x, 1, prob=prob)
  nu0 <- max(1, nu0)
  nu0
}

##
## For the marginal model, mu and tau2 are hyper-parameters.  There is
## no update.
##
setMethod("updateTau2", "MarginalModel", function(object){
  ##hypp <- hyperParams(object)
  ##.updateTau2(eta.0(hypp), m2.0(hypp), theta(object), mu(object), k(object))
  ##tau2(object)
  .updateTau2(object)
})

.updateTau2 <- function(object){
  hypp <- hyperParams(object)
  ##eta.0, m2.0, theta, mu, k){
  eta.k <- eta.0(hypp)+k(object)
  s2.k <- sum((theta(object)-mu(object))^2)
  m2.k <- 1/eta.k*(eta.0(hypp) * m2.0(hypp) + s2.k)
  tau2 <- 1/rgamma(1, shape=1/2 * eta.k, rate=1/2 * eta.k * m2.k)
  if(is.nan(tau2) || !is.finite(tau2)){
    tau2 <- 1/rgamma(1, shape=1/2*eta.0(hypp), rate=1/2*eta.0(hypp)*m2.0(hypp))
  }
  if(length(tau2) > 1) stop("tau2 should have length 1")
  tau2
}


setReplaceMethod("tau2", "MarginalModel", function(object, value){
  ##hyperParams(object)@tau2 <- value
  object@tau2 <- value
  object
})

setReplaceMethod("mu", "MarginalModel", function(object, value){
  ##hyperParams(object)@mu <- value
  object@mu <- value
  object
})

##
## For the marginal model, mu has already been initialized in the hyperparameters
##
setMethod("initializeMu", "MarginalModel", function(object)   mu(object))

#' @export
setMethod("bic", "MarginalModel", function(object, ...){
  if(k(object) > 1){
    object <- updateWithPosteriorMeans(object)
  }
  ## K: number of free parameters to be estimated
  ##   - component-specific parameters:  theta, sigma2   (3 x k(model))
  ##   - mixing probabilities:  k-1
  ##   - length-one parameters: mu, tau2, sigma2.0, nu.0             +4
  K <- 2*k(object) + (k(object)-1) + 4
  n <- length(y(object))
  -2*logpotential(object) + K*(log(n) - log(2*pi))
})

setMethod("initializeTheta", "MarginalModel", function(object){
  initializeTheta(k(object))
})

setMethod("theta", "MarginalModel", function(object) object@theta)

setMethod("sigma2", "MarginalModel", function(object) object@sigma2)

setMethod("reorderComponents", "MarginalModel", function(object){
  ix <- order(theta(object))
  theta(object) <- sort(theta(object))
  sigma2(object) <- sigma2(object)[ix]
  p(object) <- p(object)[ix]
  zz <- z(object)
  ## trick to reset the augmentation var
  zz <- factor(as.integer(factor(zz, levels=ix)), levels=seq_len(k(object)))
  z(object) <- zz
  dataMean(object) <- dataMean(object)[ix]
  dataPrec(object) <- dataPrec(object)[ix]
  object
})

setMethod("updateWithPosteriorMeans", "MarginalModel", function(object){
  mc <- mcmcChains(object)
  theta(object) <- colMeans(theta(mc))
  sigma2(object) <- colMeans(sigma2(mc))
  p(object) <- colMeans(p(mc))
  nu.0(object) <- median(nu.0(mc))
  mu(object) <- mean(mu(object))
  tau2(object) <- mean(tau2(object))
  sigma2.0(object) <- mean(sigma2.0(object))
  logpotential(object) <- computePotential(object)
  z(object) <- factor(map(object), levels=seq_len(k(object)))
  object
})

setMethod("sort", "MarginalModel", function(x, decreasing=FALSE, ...){
  mc <- mcmcChains(x)
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
  mcmcChains(x) <- mc
  x
})

.computeModesMarginal <- function(object){
  mc <- mcmcChains(object)
  th <- theta(mc)
  pot <- logpotential(mc)
  i <- which.max(pot)
  thetamax <- theta(mc)[i, ]
  sigma2max <- sigma2(mc)[i,]
  pmax <- p(mc)[i, ]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax)
  modes
}

setMethod("computeModes", "MarginalModel", function(object){
  .computeModesMarginal(object)
})



.computeDistanceMarginal <- function(th, s2, P, modes, param.sds){
  thetad <- .absoluteDistance(th, modes[["theta"]])
  sigma2d <- .absoluteDistance(s2, modes[["sigma2"]])
  pid <- .absoluteDistance(P, modes[["mixprob"]])
  ##
  ## sum the distances for each parameter matrix and standardize the
  ## total distance by the standard deviation of the modal parameter
  ## estimates
  ##
  thd <- rowSums(thetad)/param.sds[["theta"]]
  s2d <- rowSums(sigma2d)/param.sds[["sigma2"]]
  pid <- rowSums(pid)/param.sds[["mixprob"]]
  if(param.sds[["sigma2"]] > 0){
    tot <- thd+s2d+pid
  } else tot <- thd+pid
  tot
}

.absoluteDistance <- function(x, y) abs(t(t(x)-y))

.updateLabels <- function(object){
  modal.params <- modes(object)
  param.sds <- sapply(modal.params, sd)
  th <- theta(object)
  s2 <- sigma2(object)
  P <- p(object)
  ix <- permutations(k(object), k(object))##
  nc <- nrow(ix)
  D <- rep(NA, nc)
  ##
  ## Compute distance to the modes for the current ordering (1,2,3)
  ## at each iteration of the chain.
  ##
  ## subtrace a vector from each row of chain matrix
  D[1] <- .computeDistanceMarginal(th, s2, P, modal.params, param.sds)
  if(FALSE) plot.ts(D[,1], col="gray")
  for(j in 2:nc){
    J <- ix[j, ]
    browser()
    D[j] <- .computeDistanceMarginal(th[J], s2[J], P[J], modal.params, param.sds)
  }
  reordering <- ix[which.min(D)]
  theta(object) <- th[reordering]
  sigma2(object) <- s2[reordering]
  p(object) <- P[reordering]
  z(object) <- factor(z(object), levels=reordering)
  dataMean(object) <- dataMean(object)[reordering]
  dataPrec(object) <- dataPrec(object)[reordering]
  object
}

setGeneric("updateLabels", function(object) standardGeneric("updateLabels"))

setMethod("updateLabels", "MarginalModel", function(object){
  .updateLabels(object)
})


setMethod("computeDistance", "MarginalModel", function(object){
  modal.params <- modes(object)
  param.sds <- sapply(modal.params, sd)
  mc <- mcmcChains(object)
  th <- theta(mc)
  s2 <- sigma2(mc)
  P <- p(mc)
  ix <- permutations(k(object), k(object))##
  nr <- nrow(th)
  nc <- nrow(ix)
  D <- matrix(NA, nr, nc)
  ##
  ## Compute distance to the modes for the current ordering (1,2,3)
  ## at each iteration of the chain.
  ##
  ## subtrace a vector from each row of chain matrix
  D[, 1] <- .computeDistanceMarginal(th, s2, P, modal.params, param.sds)
  if(FALSE) plot.ts(D[,1], col="gray")
  for(j in 2:nrow(ix)){
    J <- ix[j, ]
    D[, j] <- .computeDistanceMarginal(th[, J], s2[, J], P[, J], modal.params, param.sds)
  }
  D
})

setMethod("switchLabels", "MarginalModel", function(object){
  D <- computeDistance(object)
  ordering_index <- apply(D, 1, which.min)
  if(all(ordering_index == 1)) return(object)
  warning("Label switching occurred. Posterior probabilities for z may be incorrect")
  mc <- mcmcChains(object)
  perms <- permutations(k(object), k(object))
  tab <- as.integer(names(table(ordering_index)))
  tab <- tab[tab!=1]
  for(i in seq_along(tab)){
    mcmc.index <- which(ordering_index == tab[i])
    j <- perms[tab[i], ]
    theta(mc)[mcmc.index, ] <- theta(mc)[mcmc.index, j]
    sigma2(mc)[mcmc.index, ] <- sigma2(mc)[mcmc.index, j]
    p(mc)[mcmc.index, ] <- p(mc)[mcmc.index, j]
  }
  mcmcChains(object) <- mc
  object
})

##logLik(mcmc)[s] <- computeLoglik(object)
##  mcmcChains(object) <- mcmc
##  object
##})

##
## TODO: pass arguments to .updateThetaBatch to make it clear what
## parameters the theta update depends on
##
setMethod("updateTheta", "BatchModel", function(object) {
  .updateThetaBatch(object)
})

.updateThetaBatch <- function(object){
  ##  if(constrain==2) browser()
  ##  if(constrain!=2) constrain <- TRUE
  tau2.tilde <- 1/tau2(object)
  sigma2.tilde <- 1/sigma2(object)
  K <- k(object)
  ## should be a vector of length K
  tau2.n.tilde <- rep(NA, K)
  n.hp <- tablez(object)
  ##
  ## Guard against zero-components
  ##
  n.hp <- pmax(n.hp, 1)
  ##
  ## mu and tau2 are not batch-specific
  tau2.tilde <- matrix(tau2.tilde, nBatch(object), k(object), byrow=TRUE)
  ##
  ## Make mu the same dimension to make the arithmetic obvious
  ##
  mus <- matrix(mu(object), nBatch(object), k(object), byrow=TRUE)
  ##
  tau2.n.tilde <- tau2.tilde + n.hp * sigma2.tilde
  tau2.n <- 1/tau2.n.tilde
  tau.n <- sqrt(tau2.n)
  ##
  denom <- tau2.tilde + n.hp*sigma2.tilde
  w1 <- tau2.tilde/denom
  w2 <- n.hp*sigma2.tilde/denom
  ##
  ## when a component has 0 observations, mu.n should just be w1*mu
  ##
  ybar <- dataMean(object)
  if(any(is.nan(ybar))){
    ybar[is.nan(ybar)] <- 0
  }
  mu.n <- w1*mus + w2*ybar
  rownames(tau.n) <- rownames(mu.n) <- uniqueBatch(object)
  ##
  ##  If there are very few observations, we will be sampling from a
  ##  normal distribution with mean equal to the marginal mean. This
  ##  will result in thetas that do not satisfy the order constraints
  ##
  ##
  thetas <- matrix(NA, nBatch(object), k(object))
  rownames(thetas) <- uniqueBatch(object)
  for(b in uniqueBatch(object)){
    thetas[b, ] <- rnorm(K, mu.n[b, ], tau.n[b, ])
  }
  if(any(is.na(thetas))) stop("NAs in thetas")
  return(thetas)
}

truncNormLower <- function(theta.last, epsilon){
  a <- theta.last[1] + epsilon
}
truncNormUpper <- function(theta.last, epsilon, iequalsk){
  if(!iequalsk){
    b <- theta.last - epsilon
  } else b <- 5
  b
}

.update_equal_s2 <- function(object){
  sigma2.current <- sigma2(object)
  n.hp <- table(batch(object))[uniqueBatch(object)]
  ##
  ## guard against zeros
  ##
  nu.n <- nu.0(object) + n.hp
  thetas <- theta(object)
  ##
  ## assume variance for all components is the same
  ##
  ss <- sumSquares2(object)
  ##
  ## weighted average of sums of squares
  ##
  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- setNames(rep(NA, nBatch(object)), uniqueBatch(object))
  for(b in uniqueBatch(object)){
    sigma2.h.tilde[b] <- rgamma(1, shape=shape[b], rate=rate[b])
  }
  sigma2.h <- 1/sigma2.h.tilde
  v <- var(y(object), na.rm=TRUE) + 0.05
  if(any(sigma2.h > v)) {
    sigma2.current <- sigma2.current[, 1]
    tmp <- tryCatch(sigma2.h[sigma2.h > v] <- sigma2.current[sigma2.h > v], error=function(e) NULL)
    if(is.null(tmp)) browser()
  }
  ##
  ## return matrix of original dimension
  ##
  s2 <- matrix(sigma2.h, nBatch(object), k(object))
  rownames(s2) <- uniqueBatch(object)
  s2
}

##
## assumes all variance components are the same
## - (the UnivariateBatchModel should also be able to use this update)
##setMethod("updateSigma2", "BatchModelNoHom", function(object){
##  .update_equal_s2(object)
##})
##
####
#### Allows the first component to have a different variance, but
#### currently no constraint on the relationship between the first
#### component variance and the other components
####
##setMethod("updateSigma2", "BatchModelPlusHom", function(object){
##  .update_nonzero_s2(object)
##})
##
##.update_nonzero_s2 <- function(object){
##  sigma2.current <- sigma2(object)
##  nz <- nonZeroCopynumber(object)
##  if(length(unique(nz)) ==1){
##    ## guard against zeros
##    if(all(nz > 0)){
##      nz[which.min(y(object))] <- 0
##    } else nz[which.max(y(object))] <- 1
##  }
##  n.hp <- table(batch(object), nz)
##  n.hp <- n.hp[uniqueBatch(object), ]
##  ##
##  ## guard against zeros
##  ##
##  n.hp <- pmax(n.hp, 1)
##  nu.n <- nu.0(object) + n.hp
##  thetas <- theta(object)
##  ##
##  ## assume variance for copy number 1-k is the same
##  ##
##  ss <- sumSquares(object)
##  ##
##  ## weighted average of sums of squares
##  ##
##  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
##  shape <- 1/2*nu.n
##  rate <- shape*sigma2.nh
##  sigma2.h.tilde <- matrix(NA, nBatch(object), 2)
##  rownames(sigma2.h.tilde) <- rownames(thetas)
##  for(b in uniqueBatch(object)){
##    sigma2.h.tilde[b, ] <- rgamma(2, shape=shape[b, ], rate=rate[b, ])
##  }
##  sigma2.h <- 1/sigma2.h.tilde
##  v <- var(y(object), na.rm=TRUE) + 0.05
##  if(any(sigma2.h > v)) {
##    sigma2.current <- sigma2.current[, 1:2]
##    tmp <- tryCatch(sigma2.h[sigma2.h > v] <- sigma2.current[sigma2.h > v], error=function(e) NULL)
##    if(is.null(tmp)) browser()
##  }
##  ##
##  ## return matrix of original dimension
##  ##
##  s2 <- cbind(sigma2.h[, 1], matrix(sigma2.h[, 2], nBatch(object), k(object)-1))
##  s2
##}

sumSquares <- function(object){
  ss <- matrix(NA, nBatch(object), k(object))
  rownames(ss) <- uniqueBatch(object)
  B <- batch(object)
  thetas <- theta(object)
  yy <- y(object)
  zz <- z(object)
  for(b in uniqueBatch(object)){
    y <- yy[B==b]
    cn <- zz[B==b]
    ##nonzero <- nz[B==b]
    m <- thetas[b, ]
    ## This could be tricky in C.  It works in R because of the factor to an integer:
    ##  as.integer(factor(c(1, 3), levels=c("1", "2", "3"))) evaluates to 1,3
    m <- m[as.integer(cn)]
    squares <- (y - m)^2
    ss[b, ] <- sapply(split(squares, cn), sum)
  }
  ss
}

.update_sigma2 <- function(object){
  sigma2.current <- sigma2(object)
  n.hp <- tablez(object)
  ##
  ## guard against zeros
  ##
  n.hp <- pmax(n.hp, 1)
  nu.n <- nu.0(object) + n.hp
  ss <- sumSquares(object)
  ##
  ## Zeros in sums of squares occurs for batches with no observations
  ##
  ## should handle this by polymorphism
  if(k(object) == 1) ss <- ss[, 1, drop=FALSE]
  ##
  ## weighted average of sums of squares
  ##
  sigma2.nh <- 1/nu.n*(nu.0(object) * sigma2.0(object) + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- matrix(NA, nBatch(object), k(object))
  rownames(sigma2.h.tilde) <- uniqueBatch(object)
  for(b in uniqueBatch(object)){
    sigma2.h.tilde[b, ] <- rgamma(k(object), shape=shape[b, ], rate=rate[b, ])
  }
  sigma2.h <- 1/sigma2.h.tilde
  stopif(any(is.nan(sigma2.h)))
  sigma2.h
}

## special case when there is only one component
setMethod("updateSigma2", "UnivariateBatchModel", function(object){
  .update_sigma2(object)
})

setMethod("updateSigma2", "BatchModel", function(object){
  .update_sigma2(object)
})

nonZeroCopynumber <- function(object) as.integer(as.integer(z(object)) > 1)



## This is a more parsimonious model.  There are only 2 variance
## estimates for each batch: the variance of the first component and
## the variance of components k>1
##.updateSigma2Batch_2 <- function(object){
##
##}

sumSquares2 <- function(object){
  ss <- setNames(rep(NA, nBatch(object)), uniqueBatch(object))
  B <- batch(object)
  thetas <- theta(object)
  ##nz <- nonZeroCopynumber(object)
  for(b in uniqueBatch(object)){
    yy <- y(object)[B==b]
    zz <- z(object)[B==b]
    ##nonzero <- nz[B==b]
    m <- thetas[b, ]
    m <- m[as.integer(zz)]
    ss[b] <- sum((yy - m)^2)
  }
  ss
}

##
## If K is 2, assume there is no homozygous deletion component
##  Constrain the variances to be the same so that one component will
##  not have heavier tails and capture the outliers
##.updateSigma2Batch_samevar <- function(object){
##
##}

setMethod("updateMu", "BatchModel", function(object){
  .updateMuBatch(object)
})

##.updateMu <- function(tau2.0, tau2, k, z, theta, mu.0){
.updateMuBatch <- function(object){
  hypp <- hyperParams(object)
  tau2.0.tilde <- 1/tau2.0(hypp)
  tau2.tilde <- 1/tau2(object)
  P <- nBatch(object)
  tau2.P.tilde <- tau2.0.tilde + P*tau2.tilde
  n.h <- tablez(object)
  ## guard against components with zero observations
  n.h <- pmax(n.h, 1)
  thetas <- theta(object)
  ##
  ## weights for within-component average of thetas
  w1 <- tau2.0.tilde/(tau2.0.tilde + P*tau2.tilde)
  w2 <- P*tau2.tilde/(tau2.0.tilde + P*tau2.tilde)
  ##
  ## average thetas, giving more weight to batches with more
  ## observations
  ##
  theta.bar <- colSums(n.h*thetas)/colSums(n.h)
  ## when the prior is zero, mu.P is shrunk towards zero
  mu.P <- w1*mu.0(hypp) + w2*theta.bar
  mu.P
}

setMethod("updateTau2", "BatchModel", function(object){
  .updateTau2Batch(object)
})

##.updateTau2Batch <- function(eta.0, m2.0, theta, mu, k){
.updateTau2Batch <- function(object){
  hypp <- hyperParams(object)
  P <- nBatch(object)
  eta.P <- eta.0(hypp)+P
  mus <- mu(object)
  mus <- matrix(mus, P, k(object), byrow=TRUE)
  thetas <- theta(object)
  s2.P <- colSums((thetas-mus)^2)
  m2.P <- 1/eta.P * (eta.0(hypp) * m2.0(hypp) + s2.P)
  tau2 <- 1/rgamma(k(object), shape=1/2 * eta.P, rate=1/2 * eta.P * m2.P)
  tau2
}

setMethod("updateSigma2.0", "BatchModel", function(object){
  .updateSigma2.0Batch(object)
})

##.updateSigma2.0 <- function(a, b, nu.0, sigma2.h, k){
.updateSigma2.0Batch <- function(object){
  hypp <- hyperParams(object)
  P <- nBatch(object)
  sigma2s <- as.numeric(sigma2(object))
  a.k <- a(hypp) + 1/2 * (k(object) * P) * nu.0(object)
  b.k <- b(hypp) + 1/2 * sum(1/sigma2s)
  sigma2.0 <- rgamma(1, shape=a.k, rate=b.k)
  stopifnot(sigma2.0 > 0)
  stopif(is.nan(sigma2.0))
  sigma2.0
}

setMethod("updateNu.0", "BatchModel", function(object){
  .updateNu.0Batch(object)
})

##.updateNu.0Batch <- function(NUMAX=100, beta, sigma2.0, sigma2.h, nu.0, k){
.updateNu.0Batch <- function(object){
  NUMAX <- 100
  hypp <- hyperParams(object)
  x <- seq_len(NUMAX)
  k <- k(object)
  P <- nBatch(object)
  sigma2s <- as.numeric(sigma2(object))
  lpnu0 <- (k*P) * (0.5 * x * log(sigma2.0(object) * x/2)-lgamma(x/2)) +
    (x/2 - 1) * sum(log(1/sigma2s)) +
      -x * (betas(hypp) + 0.5 * sigma2.0(object) * sum(1/sigma2s))
  prob <- exp(lpnu0 - max(lpnu0))
  nu0 <- sample(x, 1, prob=prob)
  nu0
}

#' @export
setMethod("mu", "BatchModel", function(object) object@mu)

#' @export
setMethod("tau2", "BatchModel", function(object) object@tau2)

setReplaceMethod("tau2", "BatchModel", function(object, value){
  object@tau2 <- value
  object
})

setReplaceMethod("mu", "BatchModel", function(object, value){
  object@mu <- value
  object
})



##setMethod("initializeSigma2", "UnivariateBatchModel", function(object){
##  tmp <- t(replicate(nBatch(object), 1/rgamma(k(object), shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))))
##  matrix(tmp, ncol=1)
##})

##setMethod("initializeTheta", "UnivariateBatchModel", function(object){
##  tmp <- t(replicate(nBatch(object), sort(rnorm(k(object), mu(object), tau2(object)))))
##  tmp <- matrix(tmp, ncol=1)
##  tmp
##})

##componentCapturesTails <- function(object){
##  ix <- order(y(object))
##  if(any(head(z(object)[ix]) %in% tail(z(object)[ix]))){
##    is_outlier <- TRUE
##  } else is_outlier <- FALSE
##  is_outlier
##}

setMethod("bic", "BatchModel", function(object, ...){
  if(k(object) > 1){
    object <- updateWithPosteriorMeans(object)
  }
  ## K: number of free parameters to be estimated
  ##   - component and batch-specific parameters:  theta, sigma2  ( k(model) * nBatch(model))
  ##   - mixing probabilities: (k-1)*nBatch
  ##   - component-specific parameters: mu, tau2                 2 x k(model)
  ##   - length-one parameters: sigma2.0, nu.0                   +2
  K <- 2*k(object)*nBatch(object) + (k(object)-1) + 2*k(object) + 2
  n <- length(y(object))
  bicstat <- -2*logpotential(object) + K*(log(n) - log(2*pi))
  bicstat
})


##setMethod("bic", "BatchModelNoHom", function(object, ...){
##  if(k(object) > 1){
##    object <- updateWithPosteriorMeans(object)
##  }
##  ## K: number of free parameters to be estimated
##  ##   - component and batch-specific parameters:  theta  ( k(model) * nBatch(model))
##  ##   - 1 variance estimate for each batch:  nBatch(model)
##  ##   - component-specific parameters: mu, tau2                 2 x k(model)
##  ##   - length-one parameters: sigma2.0, nu.0                   +2
##  ##   - mixture probs:  +3
##  K <- k(object)*nBatch(object) + nBatch(object) +  2*k(object) + 2 + 3
##  ## Experimental: extra penalty
##  K <- K + 2
##  n <- length(y(object))
##  bicstat <- -2*logpotential(object) + K*(log(n) - log(2*pi))
##  bicstat
##})

setMethod("theta", "BatchModel", function(object) {
  b <- object@theta
  b <- matrix(b, nBatch(object), k(object))
  rownames(b) <- uniqueBatch(object)
  b
})

setReplaceMethod("theta", "BatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@theta <- value
  object
})



setReplaceMethod("sigma2", "BatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})

setMethod("sigma2", "BatchModel", function(object) {
  s2 <- object@sigma2
  s2 <- matrix(s2, nBatch(object), k(object))
  rownames(s2) <- uniqueBatch(object)
  s2
})

setMethod("reorderComponents", "BatchModel", function(object){
  thetas <- theta(object)
  sigma2s <- sigma2(object)
  for(i in seq_len(nrow(thetas))){
    x <- thetas[i, ]
    j <- order(x)
    thetas[i, ] <- x[j]
    sigma2s[i, ] <- sigma2s[i, j]
    if(i == 1){
      p(object) <- p(object)[j]
    }
    dataMean(object)[i, ] <- dataMean(object)[i, j]
    dataPrec(object)[i, ] <- dataPrec(object)[i, j]
  }
  theta(object) <- thetas
  sigma2(object) <- sigma2s
  ##ix <- order(theta(object)[1, ])
  ##theta(object) <- theta(object)[, ix]
  ##sigma2(object) <- sigma2(object)[, ix]
  ##p(object) <- p(object)[ix]
  zz <- z(object)
  zz <- factor(as.integer(factor(zz, levels=j)), levels=seq_len(k(object)))
  z(object) <- zz
  ##  dataMean(object) <- dataMean(object)[, ix]
  ##  dataPrec(object) <- dataPrec(object)[, ix]
  object
})



setAs("BatchModel", "SummarizedExperiment", function(from, to){
  cnmat <- matrix(y(from), 1, length(y(from)))
  cnmat <- oligoClasses::integerMatrix(cnmat, 1000)
  se <- SummarizedExperiment(assays=SimpleList(medr=cnmat),
                             colData=DataFrame(plate=batch(from)))
  se
})

#' @export
setMethod("sort", "BatchModel", function(x, decreasing=FALSE, ...){
  mc <- mcmcChains(x)
  pot <- logpotential(mc)
  index <- which.max(pot)
  if(FALSE){
    modal.params <- list(p=pic(x)[index, ],
                         theta=thetac(x)[index, ],
                         sigma=sigmac(x)[index, ])
    ##
    ## TODO: foreach iteration, order so that the distance to the modal
    ## parameters is minimized
    ##
  }
  thetas <- matrix(theta(mc)[index, ], nBatch(x), k(x))
  thetas <- thetas[1, ]
  if(identical(thetas, sort(thetas))){
    ## nothing to do
    return(x)
  }
  B <- nBatch(x); K <- k(x)
  cn <- order(thetas)

  ## figure out the appropriate indices to sort
  tmp <- matrix(seq_len(B*K), B, K)
  tmp <- tmp[, cn]
  ix <- as.numeric(tmp)

  theta(mc) <- theta(mc)[, ix]
  theta(x) <- theta(x)[, cn]

  sigma2(mc) <- sigma2(mc)[, ix]
  sigma2(x) <- sigma2(x)[, cn]

  p(mc) <- p(mc)[, cn]
  p(x) <- p(x)[cn]

  mu(mc) <- mu(mc)[, cn]
  tau2(mc) <- tau2(mc)[, cn]
  mu(x) <- mu(x)[cn]
  tau2(x) <- tau2(x)[cn]

  probz(x) <- probz(x)[, cn]

  zz <- as.integer(z(x))
  z(x) <- factor(as.integer(factor(zz, levels=cn)), levels=sort(unique(zz)))
  dataMean(x) <- dataMean(x)[, cn]
  dataPrec(x) <- dataPrec(x)[, cn]
  mcmcChains(x) <- mc
  x
})

.computeDistOneBatch <- function(th, s2, P, mus, tau2s, modes, param.sds){
  thetad <- .absoluteDistance(th, modes[["theta"]])
  sigma2d <- .absoluteDistance(s2, modes[["sigma2"]])
  pid <- .absoluteDistance(P, modes[["mixprob"]])
  mud <- .absoluteDistance(mus, modes[["mu"]])
  tau2d <- .absoluteDistance(tau2s, modes[["tau2"]])
  ##
  ## sum the distances for each parameter matrix and standardize the
  ## total distance by the standard deviation of the modal parameter
  ## estimates
  ##
  thd <- rowSums(thetad)/param.sds[["theta"]]
  s2d <- rowSums(sigma2d)/param.sds[["sigma2"]]
  pid <- rowSums(pid)/param.sds[["mixprob"]]
  mud <- rowSums(mud)/param.sds[["mu"]]
  tau2d <- rowSums(tau2d)/param.sds[["tau2"]]
  if(param.sds[["sigma2"]] > 0){
    tot <- thd+s2d+pid+mud+tau2d
  } else tot <- thd+pid+mud+tau2d
  tot
}

## compute distance for a given permutation of columns
.computeDistanceOnePerm <- function(mc, column.permutation, modes){
  param.sds <- sapply(modes, sd)
  .computeDistOneBatch(th=theta(mc)[, column.permutation],
                       s2=sigma2(mc)[, column.permutation],
                       mus=mu(mc)[, column.permutation],
                       tau2s=tau2(mc)[, column.permutation],
                       P=p(mc)[, column.permutation],
                       modes=modes,
                       param.sds=param.sds)
}

setMethod("computeDistance", "BatchModel", function(object){
  modal.params <- modes(object)
  ##param.sds <- sapply(modal.params, function(x) sd(x[1,]))
  mc <- mcmcChains(object)
  th <- theta(mc)
  s2 <- sigma2(mc)
  ## ix is the number of possible orderings for each batch
  ix <- permutations(k(object), k(object))##
  nr <- nrow(th)
  nc <- nrow(ix)
  Dlist <- vector("list", nBatch(bmodel))
  ## iterate over batches
  for(b in seq_len(nBatch(object))){
    mc2 <- mc
    batch.index <- seq(b, nBatch(object)*k(object), by=nBatch(object))
    theta(mc2) <- th[, batch.index]
    sigma2(mc2) <- s2[, batch.index]
    D <- matrix(NA, nr, nc)
    m.params <- list(theta=modal.params[["theta"]][b,],
                     sigma2=modal.params[["sigma2"]][b,],
                     mu=modal.params[["mu"]],
                     tau2=modal.params[["tau2"]],
                     mixprob=modal.params[["mixprob"]])
    ## iterate over all possible permutations
    for(j in 1:nrow(ix)){
      J <- ix[j, ]
      D[, j] <- .computeDistanceOnePerm(mc=mc2, column.permutation=J,
                                        modes=m.params)
    }
    Dlist[[b]] <- D
  }
  Dlist
})

setMethod("switchLabels", "BatchModel", function(object){
  Dlist <- computeDistance(object)
  mc <- mcmcChains(object)
  warn <- FALSE
  for(b in seq_along(Dlist)){
    D <- Dlist[[b]]
    ordering_index <- apply(D, 1, which.min)
    if(all(ordering_index == 1)) next()
    warn <- TRUE
    batch.index <- seq(b, nBatch(object)*k(object), by=nBatch(object))
    mc2 <- mc
    theta(mc2) <- theta(mc)[, batch.index]
    sigma2(mc2) <- sigma2(mc)[, batch.index]
    perms <- permutations(k(object), k(object))
    tab <- as.integer(names(table(ordering_index)))
    tab <- tab[tab!=1]
    for(i in seq_along(tab)){
      mcmc.index <- which(ordering_index == tab[i])
      j <- perms[tab[i], ]
      ## rewrite batch.index in mc from the permuted index in mc2
      theta(mc)[mcmc.index, batch.index] <- theta(mc2)[mcmc.index, j]
      sigma2(mc)[mcmc.index, batch.index] <- sigma2(mc2)[mcmc.index, j]
      p(mc)[mcmc.index, ] <- p(mc)[mcmc.index, j]
      mu(mc)[mcmc.index,] <- mu(mc)[mcmc.index, j]
      tau2(mc)[mcmc.index,] <- tau2(mc)[mcmc.index, j]
    }
    mcmcChains(object) <- mc
  }
  if(warn) warning("Label switching occurred. Posterior probabilities for z may be incorrect")
  object
})


.computeModesBatch <- function(object){
  mc <- mcmcChains(object)
  th <- theta(mc)
  nr <- nrow(th)
  nc <- ncol(th)
  pot <- logpotential(mc)
  i <- which.max(pot)
  nb <- nBatch(object)
  kk <- k(object)
  thetamax <- matrix(theta(mc)[i, ], nb, kk)
  sigma2max <- matrix(sigma2(mc)[i, ], nb, kk)
  pmax <- p(mc)[i, ]
  mumax <- mu(mc)[i, ]
  tau2max <- tau2(mc)[i,]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mumax,
                tau2=tau2max)
  modes
}

setMethod("computeModes", "BatchModel", function(object){
  .computeModesBatch(object)
})

setMethod("tracePlot", "BatchModel", function(object, name, ...){
  ilist <- foreach(j=1:nBatch(object)) %do% seq(j, nBatch(object)*k(object), nBatch(object))
  uB <- uniqueBatch(object)
  if(name=="theta"){
    ##op <- par(mfrow=c(3, 3), las=1)
    foreach(k=1:nBatch(object)) %do% {
      plot.ts(thetac(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k], ...)
    }
    ##par(op)
  }
  if(name=="sigma"){
    ##op <- par(mfrow=c(nBatch(object)/3, 3), las=1)
    foreach(k=1:nBatch(object)) %do% {
      plot.ts(sigmac(object)[, ilist[[k]]], ylab="", xlab="",
              plot.type="single", main=uB[k],...)
    }
    ##par(op)
  }
  if(name=="p"){
    ##op <- par(mfrow=c(1, k(object)), las=1)
    ##foreach(k=1:nBatch(object)) %do% {
    plot.ts(pic(object),  ...)
    ##plot.ts(pic(object), col="gray", ...)
    ##par(op)
  }
  if(name=="mu"){
    ##op <- par(mfrow=c(1, k(object)), las=1)
    plot.ts(muc(object),  ...)
    ##par(op)
  }
  if(name=="tau"){
    ##op <- par(mfrow=c(1, k(object)), las=1)
    plot.ts(tauc(object),  ...)
    ##par(op)
  }
})

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})

setMethod("updateZ", "UnivariateBatchModel", function(object){
  factor(rep(1, length(y(object))))
})

setMethod("sigmaMean", "BatchModel", function(object) {
  mns <- colMeans(sigmac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

setMethod("thetaMean", "BatchModel", function(object) {
  mns <- colMeans(thetac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

setMethod("pMean", "BatchModel", function(object) {
  mns <- colMeans(pic(object))
  mns
})


setMethod("[", "BatchModel", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    y(object) <- y(object)[i]
    z(object) <- z(object)[i]
    batch(objct) <- batch(object)[i]
  }
  object
})

setMethod("computeLoglik", "MarginalModel", function(object){
  mn <- theta(object)
  ss <- sigma(object)
  yy <- y(object)
  zz <- as.integer(z(object))
  pp <- p(object)
  lik <- pp[zz]*dnorm(yy, mn[zz], ss[zz])
  sum(log(lik))
})
