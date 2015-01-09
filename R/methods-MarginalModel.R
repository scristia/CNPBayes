MarginalModel <- function(data, k, batch){
  new("MarginalModel",
      hyperparams=HyperparametersMarginal(k=k),
      theta=numeric(k),
      sigma2=numeric(k),
      ##mu=numeric(k),
      ##tau2=numeric(k),
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
      hwe=numeric())
}

UnivariateMarginalModel <- function(data, k=1, batch){
  new("UnivariateMarginalModel",
      hyperparams=HyperparametersMarginal(k=k),
      theta=numeric(k),
      sigma2=numeric(k),
      ##mu=numeric(1),
      ##tau2=numeric(1),
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
      hwe=numeric())
}

#' @export
setMethod("mu", "MarginalModel", function(object) mu(hyperParams(object)))

#' @export
setMethod("tau2", "MarginalModel", function(object) tau2(hyperParams(object)))

setMethod("updateZ", "UnivariateMarginalModel", function(object){
  factor(rep(1, length(y(object))))
})

setMethod("updateZ", "UnivariateBatchModel", function(object){
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
  hypp <- hyperParams(object)
  ##  mu.k <- .updateMu(tau2.0(hypp), tau2(object), k(object), z(object),
  ##                    theta(object), mu.0(hypp))
  ##  mu.k
  mu(object)
})

.updateMu <- function(tau2.0, tau2, k, z, theta, mu.0){
  tau2.0.tilde <- 1/tau2.0
  tau2.tilde <- 1/tau2
  tau2.k.tilde <- tau2.0.tilde + k*tau2.tilde
  nn <- table(z)
  theta.bar <- sum(nn*theta)/sum(nn)
  mu.k <- tau2.0.tilde/(tau2.0.tilde + k*tau2.tilde)*mu.0 +
    k*tau2.tilde/(tau2.0.tilde + k*tau2.tilde)*theta.bar
  mu.k
}

setMethod("startingValues", "MarginalModel", function(object){
  hypp <- hyperParams(object)
  tmp.file <- tempfile()
  sink(tmp.file)
  mmfit <- normalmixEM(y(object), arbvar = FALSE,
                       epsilon = 1e-03, k=k(hypp), maxit=2000)
  sink()
  unlink(tmp.file)
  mus <- mmfit$mu
  vars <- (mmfit$sigma[order(mmfit$mu)])^2
  if(FALSE){
  if(k(object) >= 3 && any(y(object) < - 1.5)) {
    ## set first component to be near the expected mean for homozygous deletions
    mus[1] <- median(y(object)[y(object) < -1.5])
    ## set first component to be near the expected mean for hemizygous deletions
    mus[2] <- median(y(object)[y(object) > -1.5 & y(object) < -0.15])
    ##
    ## the variance may be overestimated
    vars[1] <- var(y(object)[y(object) < -1.5])
    vars[2:k(object)] <- var(y(object)[y(object) > -1.5 & y(object) < -0.15])

    ## initialize using the posterior
  }
}
  theta(object) <- mus
  sigma2(object) <- vars
  object
})



setMethod("initializeSigma2.0", "MarginalModel", function(object){
  hypp <- hyperParams(object)
  sum(alpha(hypp)*sigma2(object))/sum(alpha(hypp))
})

## 0 + 3sd = 0.4
setMethod("initializeTau2", "MarginalModel", function(object){
  ##hypp <- hyperParams(object)
  ##1/rgamma(1, shape=1/2*eta.0(hypp), rate=1/2*eta.0(hypp)*m2.0(hypp))
  ##rep(0.67, k(object))
  ##rep((0.4/3)^2, k(object))
  rep(0.68^2, k(object))
})

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
  means[is.nan(means)] <- mu(object)[is.nan(means)]
  means
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
  yy <- rnorm(length(zz), mean=theta(object)[zz], sd=sqrt(sigma2(object)[zz]))
})

setMethod("moveChain", "MarginalModel", function(object, s){
  mcmc <- mcmcChains(object)
  theta(mcmc)[s, ] <- as.numeric(theta(object))
  sigma2(mcmc)[s, ] <- as.numeric(sigma2(object))
  p(mcmc)[s, ] <- p(object)
  ##
  ## mu and tau2 are now hyper-parameters
  ##mu(mcmc)[s] <- mu(object)
  ##tau2(mcmc)[s] <- tau2(object)
  ##
  nu.0(mcmc)[s] <- nu.0(object)
  sigma2.0(mcmc)[s] <- sigma2.0(object)
  logpotential(mcmc)[s] <- logpotential(object)
  mcmcChains(object) <- mcmc
  object
})


setMethod("updateTheta", "MarginalModel", function(object, constrain) {
  .updateTheta(mu(object), tau2(object),
               sigma2(object),
               dataMean(object),
               n.h=table(z(object)),
               theta.last=theta(object),
               constrain=constrain)
})

.updateTheta <- function(mu, tau2,
                         sigma2,
                         data.mean,
                         n.h,
                         theta.last,
                         constrain=TRUE){
  tau2.tilde <- 1/tau2
  sigma2.tilde <- 1/sigma2
  tau2.n.tilde <- tau2.tilde + n.h*sigma2.tilde
  ##tau2.n.tilde <- posteriorPrecisionConjugateNormal(tau2.tilde, data.precision)
  tau2.n <- 1/tau2.n.tilde
  denom <- tau2.tilde + n.h*sigma2.tilde
  w1 <- tau2.tilde/denom
  w2 <- n.h*sigma2.tilde/denom
  mu.n <- w1*mu + w2*data.mean
  k <- length(tau2.n)
  thetas <- rnorm(k, mu.n, sqrt(tau2.n))
  stopif(any(is.na(thetas)))
  ##
  ## Do not constrain the thetas during burnin
  ##
  if(!constrain) return(thetas)
  if(identical(thetas, sort(thetas)))  return(thetas)
  ##
  ## sorted thetas are not the same
  ##
  ## - constrain updates according to theta.last
  epsilon <- 0.01
  thetas[1] <- rtruncnorm(1, a=-Inf, b=theta.last[2]-epsilon,
                          mean=theta.last[1], sd=sqrt(tau2.n[1]))
  for(i in 2:k){
    a <- thetas[i-1] + epsilon
    ##if(i==4) stop()
    if(i < k){
      ##b <- theta.last[i+1] - epsilon
      b <- max(mu.n[i] + min(3*sqrt(tau2.n[i]), 0.1), a+0.2)
    } else b <- 5
    thetas[i] <- rtruncnorm(1, a=a, b=b,
                            mean=mu.n[i], sd=sqrt(tau2.n[i]))
  }
  stopif(any(is.na(thetas)))
  thetas
}

setMethod("updateSigma2", "MarginalModel", function(object) {
  ##.updateSigma2(object)
  .updateSigma2_2(object)
})


##.updateSigma2 <- function(data.list, thetas, nu.0, sigma2.0, n.h){
.updateSigma2 <- function(object){
  data.list <- split(y(object), z(object))
  thetas <- theta(object)
  nu.0 <- nu.0(object)
  sigma2.0 <- sigma2.0(object)

  nu.n <- nu.0+n.h
  k <- length(nu.n)
  ss <- sumOfSquares(data.list, thetas)
  ## weighted average of sums of squares
  sigma2.nh <- 1/nu.n*(nu.0*sigma2.0 + ss)
  shape <- 1/2*nu.n
  rate <- shape*sigma2.nh
  sigma2.h.tilde <- rgamma(k, shape=shape, rate=rate)
  ##tmp <- rgamma(1000, shape=1/2*nu.n[1], rate=1/2*nu.n[1]*sigma2.nh[1])
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
  hypp <- hyperParams(object)
  ##.updateTau2(eta.0(hypp), m2.0(hypp), theta(object), mu(object), k(object))
  tau2(object)
})

.updateTau2 <- function(eta.0, m2.0, theta, mu, k){
  eta.k <- eta.0+k
  s2.k <- sum((theta-mu)^2)
  m2.k <- 1/eta.k*(eta.0 * m2.0 + s2.k)
  tau2 <- 1/rgamma(1, shape=1/2 * eta.k, rate=1/2 * eta.k * m2.k)
  stopif(is.nan(tau2))
  tau2
}


setReplaceMethod("tau2", "MarginalModel", function(object, value){
  hyperParams(object)@tau2 <- value
  object
})

setReplaceMethod("mu", "MarginalModel", function(object, value){
  hyperParams(object)@mu <- value
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
  ##   - component-specific parameters:  theta, pi  (2 x k(model))
  ##   - 2 variance estimates (non-zero  components assumed to be the same) +2
  ##   - length-one parameters: sigma2.0, nu.0             +2
  K <- 2*k(object) + 2 + 2
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
