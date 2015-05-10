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
                m.y=numeric(1),
                mcmc.params=mcmc.params)
  object <- startingValues(object)
}

#' @export
UnivariateMarginalModel <- function(data, k=1, batch, hypp){
  if(missing(hypp)) hypp <- HyperparametersMarginal(k=k)
  if(!missing(batch)){
    nbatch <- setNames(as.integer(table(batch)), levels(batch))
  } else nbatch <- length(data)
  zz <- as.integer(factor(numeric(K)))
  zfreq <- as.integer(table(zz))
  new("UnivariateMarginalModel",
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
      probz=matrix(1, length(data), 1),
      logprior=numeric(1),
      loglik=numeric(1),
      mcmc.chains=McmcChains(),
      batch=batch,
      ##uniqueBatch=unique(batch),
      batchElements=nbatch,
      hwe=numeric(),
      modes=list(),
      m.y=numeric(1))
}

#' @export
getK <- function(object){
  hypp <- hyperParams(object)
  .Call("getK", hypp)
}

#' @export
setMethod("mu", "MarginalModel", function(object) object@mu)

#' @export
setMethod("tau2", "MarginalModel", function(object) object@tau2)

setMethod("updateZ", "UnivariateMarginalModel", function(object){
  factor(rep(1, length(y(object))))
})

## compute p(theta) * p(y | theta)
setMethod("computePotential", "MarginalModel", function(object){
  ll <- computeLogLikxPrior(object)
  ll.phi <- .loglikPhiMarginal(object)
  ll+ll.phi
})

##
## For the marginal model, mu and tau2 are hyper-parameters.  There is
## no update.
##
setMethod("updateMu", "MarginalModel", function(object){
  .Call("update_mu", object)
  ##.updateMu(object)
})

.updateMu <- function(object){
  hypp <- hyperParams(object)
  tau2.0.tilde <- 1/tau2.0(hypp)
  tau2.tilde <- 1/tau2(object)
  tau2.k.tilde <- tau2.0.tilde + k(object)*tau2.tilde
  nn <- zFreq(object)
  theta.bar <- sum(nn*theta(object))/sum(nn)
  w1 <- tau2.0.tilde/tau2.k.tilde
  w2 <- k(object)*tau2.tilde/tau2.k.tilde
  mu.k <- w1*mu.0(hypp) + w2*theta.bar
  sd <- sqrt(1/tau2.k.tilde)
  mu_new <- rnorm(1, mu.k, sd)
  mu_new
}





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
##    .multinomial_probs <- .posteriorMultinomial(y(object),
##                                                theta(object),
##                                                sqrt(sigma2(object)),
##                                                p(object))
##
  .Call("update_multinomialPr", object)
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
  .Call("compute_means", object)

})

.computeMeans <- function(object){
  means <- sapply(split(y(object), z(object)), mean, na.rm=TRUE)
  if(any(is.nan(means))) {
    means[is.nan(means)] <- rnorm(sum(is.nan(means)), mu(object), tau(object))
  }
  means
}

setMethod("computeVars", "MarginalModel", function(object){

  .Call("compute_vars", object)
})

.computeVars <- function(object){
  vars <- sapply(split(y(object), z(object)), var, na.rm=TRUE)
  if(any(is.nan(vars))){
    vars[is.nan(vars)] <- 1/rgamma(sum(is.nan(vars)), shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  }
  vars
}



#' @export

setMethod("simulateY", "MarginalModel", function(object){
  zz <- z(object)
  yy <- rnorm(length(zz), mean=theta(object)[zz], sd=sigma(object)[zz])
})

setMethod("moveChain", "MarginalModel", function(object, s){
  mcmc <- mcmcChains(object)
  K <- k(object)
  theta(mcmc)[s, ] <- as.numeric(theta(object))
  sigma2(mcmc)[s, ] <- as.numeric(sigma2(object))
  p(mcmc)[s, ] <- p(object)
  ##
  mu(mcmc)[s] <- mu(object)
  tau2(mcmc)[s] <- tau2(object)
  ##
  nu.0(mcmc)[s] <- nu.0(object)
  sigma2.0(mcmc)[s] <- sigma2.0(object)
  ##logpotential(mcmc)[s] <- logpotential(object)
  logLik(mcmc)[s] <- logLik(object)
  ##zz <- factor(z(object), levels=seq_len(K))
  zFreq(mcmc)[s, ] <- as.integer(table(z(object)))
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


setMethod("updateTheta", "MarginalModel", function(object) {
  ##.updateTheta(object)
  .Call("update_theta", object)
})

.updateTheta <- function(object){
  theta.last <- theta(object)
  tau2.tilde <- 1/tau2(object)
  sigma2.tilde <- 1/sigma2(object)
  n.h <- tablez(object)
  n.h <- pmax(n.h, 1)
  ##tau2.n.tilde <- tau2.tilde + n.h*sigma2.tilde
  post.prec <- tau2.tilde + n.h*sigma2.tilde
  tau2.n <- 1/post.prec
  ##denom <- tau2.tilde + n.h*sigma2.tilde
  w1 <- tau2.tilde/post.prec
  w2 <- n.h*sigma2.tilde/post.prec
  mu.n <- w1*mu(object) + w2*dataMean(object)
  thetas <- rnorm(k(object), mu.n, sqrt(tau2.n))
  if(any(is.na(thetas))) stop("NAs in theta update")
  if(length(thetas) != length(theta.last)) stop("check thetas")
  thetas
}

## setMethod("updateSigma2Cpp", "MarginalModel", function(object) {
##   .Call("update_sigma2", object)
## })

setMethod("updateSigma2", "MarginalModel", function(object) {
  ##.updateSigma2(object)
  .Call("update_sigma2", object)
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
  .Call("update_sigma2_0", object)
##     hypp <- hyperParams(object)
##     .updateSigma2.0(a=a(hypp), b=b(hypp), nu.0=nu.0(object),
##                     sigma2.h=sigma2(object), k=k(hypp))
})

.updateSigma2.0 <- function(a, b, nu.0, sigma2.h, k){
  a.k <- a+1/2*k*nu.0
  b.k <- b+1/2*sum(1/sigma2.h)
  sigma2.0 <- rgamma(1, shape=a.k, rate=b.k)
  stopif(is.nan(sigma2.0))
  sigma2.0
}

setMethod("updateNu.0", "MarginalModel", function(object){
  .Call("update_nu0", object)
##   hypp <- hyperParams(object)
##   .updateNu.0(beta=betas(hypp), sigma2.0=sigma2.0(object), sigma2.h=sigma2(object),
##               nu.0=nu.0(object), k=k(object))
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
  ##  .updateTau2(object)
  .Call("update_tau2", object)
})

.updateTau2 <- function(object){
  hypp <- hyperParams(object)
  ##eta.0, m2.0, theta, mu, k){
  eta.k <- eta.0(hypp)+k(object)
  s2.k <- sum((theta(object)-mu(object))^2)
  m2.k <- 1/eta.k * (eta.0(hypp) * m2.0(hypp) + s2.k)

  tau2 <- 1/rgamma(1, shape=1/2 * eta.k, rate=1/2 * eta.k * m2.k)
##   if(is.nan(tau2) || !is.finite(tau2)){
##     tau2 <- 1/rgamma(1, shape=1/2*eta.0(hypp), rate=1/2*eta.0(hypp)*m2.0(hypp))
##   }
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
  logLik(object2) <- logLik(object)
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
  logLik(object2) <- logLik(object)
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
##  theta(object) <- theta(object)[, zindex, drop=FALSE]
##  sigma2(object) <- sigma2(object)[, zindex, drop=FALSE]
##  mu(object) <- mu(object)[zindex]
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
  i <- argMax(object)
  mc <- mcmcChains(object)
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
                loglik=logLik(mc)[i],
                logprior=logPrior(mc)[i])
  modes
}

setMethod("computeModes", "MarginalModel", function(object){
  .computeModesMarginal(object)
})


##logLik(mcmc)[s] <- computeLoglik(object)
##  mcmcChains(object) <- mcmc
##  object
##})

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


##permuteZ <- function(object){
##  zz <- z(object)
##  i <- permn(seq_len(k(object)))[[2]]
##  zz <- factor(as.integer(factor(zz, levels=i)), levels=seq_len(k(object)))
##  z(object) <- zz
##  object
##}

setMethod("showMeans", "MarginalModel", function(object){
  paste(round(theta(object), 2), collapse=", ")
})

setMethod("showSigmas", "MarginalModel", function(object){
  paste(round(sqrt(sigma2(object)), 2), collapse=", ")
})

setMethod("tablez", "MarginalModel", function(object) table(z(object)))

#' @export
marginalModel1 <- function(object, hyp.list,
                           data, mcmcp.list=mcmcpList(),
                           save.it=FALSE, test=FALSE){
  if(test){
    message("Testing with just a few mcmc iterations")
    mcmcp.list <- mcmcpList(TRUE)
    save.it <- FALSE
  }
  mcmcp <- mcmcp.list[[1]]
##   mp.list <- ModelParamList(hyp.list[[1]],
##                             K=1:4,
##                             data=data,
##                             mcmcp=mcmcp)
  modlist <- foreach(hypp=hyp.list) %do% {
    MarginalModel(data=y(truth), hypp=hypp, mcmc.params=mcmcp)
    ##initializeModel(params=param, hypp=hypp)
  }
  models <- foreach(k=1:4, model=modlist) %do% posteriorSimulation(model, mcmcp[k])
  file.out <- postFiles(object)[1]
  x <- lapply(models, computeMarginalPr, mcmcp=mcmcp.list[[2]])
  if(save.it) saveRDS(x, file=file.out)
  xx <- do.call(cbind, lapply(x, rowMeans))
  post.range <- unlist(lapply(x, posteriorRange))
  marginals <- computeMarginal(xx)
  for(i in seq_along(models)){
    m.y(models[[i]]) <- marginals[i]
  }
  if(save.it) saveRDS(models, file=model(object))
  ## important to return NULL -- otherwise memory will skyrocket
  NULL
}

#' @export
marginalExperiment <- function(object,
                               outdir,
                               ##mcmcp.list=mcmcpList(),
                               mcmcp=McmcParams(iter=1000, burnin=200, nStarts=20),
                               hypp,
                               marginaly=TRUE,
                               test=FALSE){
  M <- getFiles(outdir, rownames(object), "marginal")
  if(missing(hypp)){
    hp.list <- HyperParameterList(K=1:4, HyperparametersMarginal(tau2.0=1000))
  } else{
    hp.list <- HyperParameterList(K=1:4, hypp)
  }
  cn <- copyNumber(object)
  J <- seq_len(nrow(object)); j <- NULL
  x <- foreach(j = J, .packages=c("CNPBayes", "foreach")) %dopar% {
    models <- marginalModel1(M[j], data=cn[j, ],
                             hyp.list=hp.list,
                             mcmcp.list=mcmcp.list, save.it=TRUE,
                             test=test)
  }
  TRUE
}

setMethod("marginal", c("SummarizedExperiment", "missing"),
          function(object, batch, mcmc.params, K=1:4, maxperm=5, ...){
            cn <- copyNumber(object)[1, ]
            notna <- !is.na(cn)
            cn <- cn[ notna ]
            mlist <- computeMarginalEachK(cn, K=K,
                                          mcmcp=mcmc.params,
                                          maxperm=maxperm, ...)
            mlist
          })

setMethod("marginal", "MarginalModelList", ## a ModelList
          function(object, batch, mcmc.params, K=1:4, maxperm=5, ...){
            mlist <- computeMarginalEachK(object, maxperm=maxperm)
          })

setMethod("rowMarginal", c("SummarizedExperiment", "missing"),
          function(object, batch, mcmc.params, model.files, K=1:4, maxperm=5, ...){
            rowM_yGivenK(object=object,  mcmc.params=mcmc.params,
                         K=K, maxperm=maxperm, ...)
          })

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

##setMethod("reducedGibbsZThetaSigmaFixed", "MarginalModel", function(object){
##  .Call("permutedz_reduced2", object)
##})

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

    ptab <- matrix(NA, nrow(zChain(kmodZ2)), 3)
    Z <- zChain(kmodZ2)
    NN <- length(y(kmodZ2))
    ptab <- t(apply(Z, 1, function(x, NN) table(x)/NN, NN=NN))
    plot.ts(ptab, plot.type="single", col=1:3)
    abline(h=modes(kmodZ2)[["mixprob"]], col=1:3)
    i <- argMax(kmod)
    pic(kmod)[i,]
    modes(kmodZ2)[["mixprob"]]
    Z2 <- table(zChain(kmod)[i, ])
  }
  p_sigma2.z <- pSigma2(kmodZ1)
  p_pmix.z <- pMixProb(kmodZ2)
  p_theta <- mean(p_theta.z)
  p_sigma2 <- mean(p_sigma2.z)
  p_pmix <- mean(p_pmix.z)
  pstar <- c(p_theta, p_sigma2, p_pmix)
  if(any(is.na(pstar)) || any(!is.finite(pstar))) browser()
  lpstar <- log(pstar)
  if(any(is.na(lpstar)) || any(!is.finite(lpstar))) browser()
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
            kmodZ1 <- reducedGibbsThetaFixed(kmod)
            kmodZ2 <- reducedGibbsThetaSigmaFixed(kmod)
            for(i in seq_len(nrow(permutations))){
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
  x <- pThetaStar(model, T2=T2, maxperm=maxperm)
  results <- .berkhof(x, model=model)
  PosteriorSummary(p_theta=results$p_theta,
                   chib=results$chib,
                   berkhof=results$berkhof,
                   marginal=results$marginal)
}

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

simulateMultipleChains <- function(nchains, y, batch, k, mp){
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

modalLoglik <- function(x) modes(x)[["loglik"]]

#' @export
computeMarginalLik <- function(y, batch, K=1:4,
                               T=1000, burnin=200,
                               T2=200,
                               maxperm=5,
                               nchains=2){
  my <- vector("list", length(K))
  mlist <- vector("list", length(K))
  mp <- McmcParams(iter=T, nStarts=1, burnin=burnin)
  for(i in seq_along(K)){
    k <- K[i]
    kmodlist <- simulateMultipleChains(nchains, y, batch, k, mp)
    mlik <- lapply(kmodlist, berkhofEstimate, T2=T2, maxperm=maxperm)
    ll <- sapply(kmodlist, modalLoglik)
    mlist[[i]] <- kmodlist[[which.max(ll)]]
    xx <- summarizeMarginalEstimates(mlik)
    my[[i]] <- xx
  }
  if(is(mlist[[1]], "MarginalModel")){
    nms <- paste0("M", K)
    mlist <- MarginalModelList(mlist, names=nms)
  } else {
    nms <- paste0("B", K)
    mlist <- BatchModelList(mlist, names=nms)
  }
  names(my) <- nms
  results <- list(models=mlist, marginal=my)
  results
}


updateMarginalLik <- function(modellist, T=1000, burnin=200,
                              T2=200,
                              maxperm=5,
                              nchains=2){
  my <- vector("list", length(modellist))
  mp <- McmcParams(iter=T, nStarts=1, burnin=burnin)
  for(i in seq_along(modellist)){
    model <- modellist[[i]]
    mcmcParams(model, force=TRUE) <- mp
    kmodlist <- modelOtherModes(model, nchains)
    kmodlist <- lapply(kmodlist, posteriorSimulation)
    mlik <- lapply(kmodlist, berkhofEstimate, T2=T2, maxperm=maxperm)
    ll <- sapply(kmodlist, modalLoglik)
    modellist[[i]] <- kmodlist[[which.max(ll)]]
    xx <- summarizeMarginalEstimates(mlik)
    my[[i]] <- xx
  }
  names(my) <- names(modellist)
  results <- list(models=modellist, marginal=my)
  results
}

## setMethod("computeMarginalLik", c("numeric", "vector"),
##           function(y, batch, K=1:4,
##                    T=1000, burnin=200,
##                    T2=200,
##                    maxperm=5,
##                    nchains=2){
##             my <- vector("list", length(K))
##             mlist <- vector("list", length(K))
##             mp <- McmcParams(iter=T, nStarts=1, burnin=burnin)
##             modalLoglik <- function(x) modes(x)[["loglik"]]
##             for(i in seq_along(K)){
##               k <- K[i]
##               ##kmod <- posteriorSimulation(kmod)
##               if( k == 1 ) ns <- 1 else ns <- nchains
##               kmodlist <- replicate(ns, {
##                 kmod <- MarginalModel(y, k=k, mcmc.params=mp)
##                 posteriorSimulation(kmod)
##               })
##               mlik <- lapply(kmodlist, berkhofEstimate, T2=T2, maxperm=maxperm)
##               ll <- sapply(kmodlist, modalLoglik)
##               mlist[[i]] <- kmodlist[[which.max(ll)]]
##               xx <- summarizeMarginalEstimates(mlik)
##               my[[i]] <- xx
##             }
##             nms <- paste0("M", K)
##             names(my) <- nms
##             results <- list(models=MarginalModelList(model_list=mlist,
##                                 names=nms),
##                             marginal=my)
##             results
##           })

#' @export
orderModels <- function(x){
  x <- .trimNA(x)
  models <- x$models
  K <- k(models)
  ##maxdev <- sapply(K, function(x) log(factorial(x))) + 0.5
  maxdev <- 5
  marginal.est.list <- x$marginal
  m <- lapply(marginal.est.list, function(x) x[, "marginal"])
  my <- sapply(m, max)
  d <- sapply(m, function(x) diff(range(x)))
  if(sum(d < maxdev, na.rm=TRUE) >= 1){
    K <- K[ d < maxdev ]
  } else {
    warning("Marginal likelihood estimates are dissimilar for all models.  More iterations are needed")
    return(NULL)
  }
  K <- K[order(my[K], decreasing=TRUE)]
  models[ K ]
}

.trimNA <- function(object){
  mm <- object$marginal
  mm <- lapply(mm, function(x) x[rowSums(is.na(x)) == 0, , drop=FALSE] )
  object$marginal <- mm
  object
}

#' @export
logBayesFactor <- function(x){
  models <- orderModels(x)
  if(length(models) <= 1) {
    return(NA)
  }
  K <- k(models)
  ml <- sapply(x$marginal[K], function(x) mean(x[, "marginal"]))
  bf <- setNames(ml[1]-ml[2], paste0(names(ml[1:2]), collapse="-"))
  bf
}

setMethod("updateMultinomialProb", "MarginalModel", function(object){
  .Call("update_multinomialPr", object)
})

setMethod("updateMultinomialProb", "BatchModel", function(object){
  .Call("update_multinomialPr_batch", object)
})
