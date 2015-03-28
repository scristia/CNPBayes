#' @export
BatchModel <- function(data=numeric(), k=2L, batch, hypp, mcmc.params){
  if(missing(batch)) batch <- as.integer(factor(rep("a", length(data))))
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  mcmc.chains <- McmcChains()
  bf <- factor(batch)
  batch <- as.integer(bf)
  ub <- unique(batch)
  ix <- order(batch)
  nbatch <- setNames(as.integer(table(batch)), levels(bf))
  B <- length(ub)
  if(B==1 && length(data) > 0){
    if(missing(hypp)) hypp <- HyperparametersMarginal(k=k)
    zz <- as.integer(factor(numeric(k)))
    zfreq <- as.integer(table(zz))
    obj <- MarginalModel(data, k, batch, hypp, mcmc.params)
    return(obj)
  }
  if(k == 1) {
    if(missing(hypp)) hypp <- HyperparametersBatch(k=1)
    obj <- UnivariateBatchModel(data, k, batch, hypp, mcmc.params)
    return(obj)
  }
  ##browser()
  ##B <- length(ub)
  if(missing(hypp)) hypp <- HyperparametersBatch(k=k)
  zz <- integer(length(data))
  zfreq <- as.integer(table(zz))
  obj <- new("BatchModel",
             hyperparams=hypp,
             theta=matrix(NA, B, k),
             sigma2=matrix(NA, B, k),
             mu=numeric(k),
             tau2=numeric(k),
             nu.0=numeric(1),
             sigma2.0=numeric(1),
             pi=numeric(k),
             data=data[ix],
             data.mean=matrix(NA, B, k),
             data.prec=matrix(NA, B, k),
             z=zz,
             zfreq=zfreq,
             probz=matrix(0, length(data), k),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=mcmc.chains,
             mcmc.params=mcmc.params,
             batch=batch[ix],
             batchElements=nbatch,
             hwe=numeric(),
             theta_order=numeric(B*k),
             m.y=numeric(1))
  obj <- startingValues(obj)
  obj
}


##
## Multiple batches, but only 1 component
##
UnivariateBatchModel <- function(data, k=1, batch, hypp, mcmc.params){
  mcmc.chains <- McmcChains()
  batch <- factor(batch)
  ix <- order(batch)
  B <- length(levels(batch))
  nbatch <- elementLengths(split(batch, batch))
  zz <- integer(length(data))
  zfreq <- as.integer(table(zz))
  if(missing(hypp)) hypp <- HyperparametersBatch(k=1)
  obj <- new("UnivariateBatchModel",
             hyperparams=hypp,
             theta=matrix(NA, B, 1),
             sigma2=matrix(NA, B, 1),
             mu=numeric(k),
             tau2=numeric(k),
             nu.0=numeric(1),
             sigma2.0=numeric(1),
             pi=numeric(k),
             ##pi=matrix(NA, B, 1),
             data=data[ix],
             data.mean=matrix(NA, B, 1),
             data.prec=matrix(NA, B, 1),
             z=factor(numeric(length(data))),
             zfreq=zfreq,
             probz=matrix(1, length(data), 1),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=mcmc.chains,
             mcmc.params=mcmc.params,
             batch=batch[ix],
             batchElements=nbatch,
             hwe=numeric(),
             m.y=numeric(1))
  obj <- startingValues(obj)
  obj
}



setMethod("[", "BatchModel", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    y(object) <- y(object)[i]
    z(object) <- z(object)[i]
    batch(objct) <- batch(object)[i]
  }
  object
})


setMethod("batchCorrect", "BatchModel", function(object){
  B <- factor(batch(object))
  yy <- dat(object)
  names(yy) <- seq_along(yy)
  ybatch <- split(yy, B)
  zbatch <- split(z(object), B)
  thetas <- theta(object)
  sigmas <- sqrt(sigma2(object))
  i <- NULL
  ystar <- foreach(i=seq_along(ybatch), .combine='c') %do% {
    yy <- ybatch[[i]]
    zz <- zbatch[[i]]
    s <- sigmas[i, ]
    m <- thetas[i, ]
    (yy - m[zz])/s[zz]
  }
  ystar[names(yy)]
})




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

setMethod("collapseBatch", "BatchModel", function(object){
  collapseBatch(y(object), as.character(batch(object)))
})


batchLik <- function(x, p, mean, sd)  p*dnorm(x, mean, sd)


setMethod("computeMeans", "BatchModel", function(object){
  .Call("compute_means_batch", object)

})

setMethod("computePrec", "BatchModel", function(object){
  .Call("compute_prec_batch", object)
})

.computeMeansBatch <- function(object){
   if(length(y(object))==0) return(numeric())
   B <- batch(object)
   ubatch <- uniqueBatch(object)
   ybatch <- split(y(object), B)
   zbatch <- split(z(object), B)
   ymeans <- foreach(y=ybatch, z=zbatch, .combine='rbind') %do%{
     mns <- sapply(split(y, z), mean, na.rm=TRUE)
     mns
   }
   rownames(ymeans) <- names(ybatch)
   ymeans[ubatch, ]
 }

.computeModesBatch <- function(object){
  i <- argMax(object)
  mc <- mcmcChains(object)
  B <- nBatch(object)
  K <- k(object)
  thetamax <- matrix(theta(mc)[i, ], B, K)
  sigma2max <- matrix(sigma2(mc)[i, ], B, K)
  pmax <- p(mc)[i, ]
  mumax <- mu(mc)[i, ]
  tau2max <- tau2(mc)[i,]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mumax,
                tau2=tau2max,
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i],
                zfreq=zFreq(mc)[i, ],
                loglik=logLik(mc)[i],
                logprior=logPrior(mc)[i])
  modes
}

setMethod("computeModes", "BatchModel", function(object){
  .computeModesBatch(object)
})

componentVariances <- function(y, z)  v <- sapply(split(y, z), var)

setMethod("computeVars", "BatchModel", function(object){
  .Call("compute_vars_batch", object)
})

.computeVarsBatch <- function(object){
   if(length(y(object))==0) return(numeric())
   ubatch <- uniqueBatch(object)
   B <- batch(object)
   ybatch <- split(y(object), B)
   zbatch <- split(z(object), B)
   yvars <- foreach(y=ybatch, z=zbatch, .combine='rbind') %do%{
     componentVariances(y, z)
   }
   rownames(yvars) <- names(ybatch)
   yvars[ubatch, ]
}


setMethod("getThetaOrder", "MarginalModel", function(object) order(thetaMean(object)))
setMethod("getThetaOrder", "BatchModel", function(object){
  .getThetaOrdering(object)
})

.getThetaOrdering <- function(object){
  ##th <- thetaMean(object)
  th <- matrix(colMeans(thetac(object)), nBatch(object), k(object))
  ##ix <- object@theta_order
  ix <- matrix(NA, nBatch(object), k(object))
  index <- matrix(seq_len(nBatch(object)*k(object)), nBatch(object), k(object))
  for(j in 1:nrow(th)){
    ix[j, ] <- order(th[j,])
    index[j, ] <- index[j, ix[j, ]]
  }
  index <- as.numeric(index)
  index
}


setMethod("initializeSigma2.0", "BatchModel", function(object){
  hypp <- hyperParams(object)
  sum(alpha(hypp)*colMeans(sigma2(object)))/sum(alpha(hypp))
})


setMethod("moveChain", "BatchModel", function(object, s){
  mcmc <- mcmcChains(object)
  theta(mcmc)[s, ] <- as.numeric(theta(object))
  sigma2(mcmc)[s, ] <- as.numeric(sigma2(object))
  p(mcmc)[s, ] <- p(object)
  mu(mcmc)[s, ] <- mu(object)
  tau2(mcmc)[s, ] <- tau2(object)
  nu.0(mcmc)[s] <- nu.0(object)
  sigma2.0(mcmc)[s] <- sigma2.0(object)
  ##logpotential(mcmc)[s] <- computeLogLikxPrior(object)
  logPrior(mcmc)[s] <- computePrior(object)
  logLik(mcmc)[s] <- computeLoglik(object)
  mcmcChains(object) <- mcmc
  object
})



## y: length n_b vector  (number of samples in batch b)
## theta: length K vector for batch b
## sd:  length K vector for batch b
## pi: length K vector for batch b.
## Returns:  n_b x K matrix for batch b
.multBatchSpecific <- function(y, theta, sd, pi){
  K <- seq_len(length(pi))
  result <- matrix(NA, length(y), length(theta))
  for(j in K){
    result[, j] <- pi[j]*dnorm(y, theta[j], sd[j])
  }
  mix.probs <- result/rowSums(result)
  mix.probs
}



#' @export
setMethod("mu", "BatchModel", function(object) object@mu)

setReplaceMethod("mu", "BatchModel", function(object, value){
  object@mu <- value
  object
})

nBatch <- function(object) length(uniqueBatch(object))

batchElements <- function(object) object@batchElements

setReplaceMethod("p", "BatchModel", function(object, value){
  object@pi <- value
  object
})

setMethod("pMean", "BatchModel", function(object) {
  mns <- colMeans(pic(object))
  mns
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

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  ##tab
  tab[uniqueBatch(object), , drop=FALSE]
})


setMethod("sigmaMean", "BatchModel", function(object) {
  mns <- colMeans(sigmac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

#' @export
setMethod("tau2", "BatchModel", function(object) object@tau2)

setReplaceMethod("tau2", "BatchModel", function(object, value){
  object@tau2 <- value
  object
})


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



setMethod("thetaMean", "BatchModel", function(object) {
  mns <- colMeans(thetac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})


setMethod("show", "BatchModel", function(object){
  ##callNextMethod()
  cls <- class(object)
  cat(paste0("An object of class ", cls), "\n")
  cat("     n. obs      :", length(y(object)), "\n")
  cat("     n. batches  :", nBatch(object), "\n")
  cat("     k           :", k(object), "\n")
  cat("     nobs/batch  :", table(batch(object)), "\n")
  cat("     loglik (s)  :", round(logLik(object), 1), "\n")
  cat("     logprior (s):", round(logPrior(object), 1), "\n")
})

setMethod("simulateY", "BatchModel", function(object){
##  browser()
  zz <- setNames(z(object), seq_along(z(object)))
  B <- batch(object)
  zbatch <- split(zz, B)
  sigmas <- sqrt(sigma2(object))
  thetas <- theta(object)
  ysim <- foreach(b = uniqueBatch(object), .combine="c") %do% {
    m <- thetas[b, ]
    s <- sigmas[b, ]
    Z <- zbatch[[b]]
    setNames(rnorm(length(Z), mean=m[Z], s=s[Z]), names(Z))
  }
  ysim[names(zz)]
})

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab <- tab[uniqueBatch(object), , drop=FALSE]
  tab
})

#' @export
uniqueBatch <- function(object) unique(batch(object))

#' @export
batchExperiment <- function(object, outdir, mcmcp.list=mcmcpList(),
                            hypp, marginaly=TRUE, test=FALSE){
  B <- getFiles(outdir, rownames(object), "batch")
  ## shouldn't batch files be named with respect to the rownames?
  batch.files <- paste0(dirname(model(B)), "/", rownames(object), "_batch.rds")
  if(missing(hypp)){
    hp.list <- HyperParameterList(K=1:4, HyperparametersBatch(tau2.0=1000))
  } else{
    hp.list <- HyperParameterList(K=1:4, hypp)
  }
  J <- seq_len(nrow(object)); j <- NULL
  x <- foreach(j = J, .packages=c("CNPBayes", "foreach")) %dopar% {
    cn <- copyNumber(object)[j, ]
    notna <- !is.na(cn)
    bt <- saveBatch(object[j, notna], batch.files[j])
    ub <- unique(bt)
    if(length(ub) == 1){
      ## create a NULL experiment
      saveRDS(NULL, model(B)[j])
      return(NULL)
    }
    models <- batchModel1(B[j], data=cn[notna],
                          hyp.list=hp.list,
                          mcmcp.list=mcmcp.list,
                          batch=bt, save.it=TRUE, test=test,
                          marginaly=marginaly)
  }
  TRUE
}


#' @export
batchModel1 <- function(object, hyp.list, data, mcmcp.list,
                        batch,
                        save.it=FALSE, test=FALSE,
                        marginaly=TRUE){
  if(test){
    message("Testing with just a few mcmc iterations")
    mcmcp.list <- mcmcpList(TRUE)
    save.it <- FALSE
  }
  mcmcp <- mcmcp.list[[1]]
  b <- factor(batch)
  ix <- order(b)
  mp.list <- ModelParamList(hyp.list[[1]],
                            K=1:4,
                            data=data[ix],
                            batch=b[ix],
                            mcmcp=mcmcp)
  modlist <- foreach(hypp=hyp.list, param=mp.list) %do% {
    initializeBatchModel(params=param, hypp=hypp)
  }
  models <- foreach(k=1:4, model=modlist) %do% posteriorSimulation(model, mcmcp[k])
  if(marginaly){
    x <- lapply(models, computeMarginalPr, mcmcp=mcmcp.list[[2]])
    if(save.it) saveRDS(x, file=postFiles(object)[1])
  }
  marginaly.file <- postFiles(object)[1]
  if(file.exists(marginaly.file)){
    x <- readRDS(marginaly.file)
    xx <- do.call(cbind, lapply(x, rowMeans))
    post.range <- unlist(lapply(x, posteriorRange))
    marginals <- computeMarginal(xx)
    for(i in seq_along(models)){
      m.y(models[[i]]) <- marginals[i]
    }
  }
  if(save.it) saveRDS(models, file=model(object))
  ## important to return NULL -- otherwise memory will skyrocket
  NULL
}
