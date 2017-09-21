#' @include plot-functions.R
#' @include methods-MixtureModel.R
NULL


### From AllClasses.R

#' An object for running MCMC simulations.
#'
#' Run hierarchical MCMC for batch model.
#' @slot k An integer value specifying the number of latent classes.
#' @slot hyperparams An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @slot theta the means of each component and batch
#' @slot sigma2 the variances of each component and batch
#' @slot nu.0 the shape parameter for sigma2
#' @slot sigma2.0 the rate parameter for sigma2
#' @slot pi mixture probabilities which are assumed to be the same for all batches
#' @slot mu means from batches, averaged across batches
#' @slot tau2 variances from batches,  weighted by precisions
#' @slot data the data for the simulation.
#' @slot data.mean the empirical means of the components
#' @slot data.prec the empirical precisions
#' @slot z latent variables
#' @slot zfreq table of latent variables
#' @slot probz n x k matrix of probabilities
#' @slot logprior log likelihood of prior: log(p(sigma2.0)p(nu.0)p(mu))
#' @slot loglik log likelihood: \eqn{\sum p_k \Phi(\theta_k, \sigma_k)}
#' @slot mcmc.chains an object of class 'McmcChains' to store MCMC samples
#' @slot batch a vector of the different batch numbers
#' @slot batchElements a vector labeling from which batch each observation came from
#' @slot modes the values of parameters from the iteration which maximizes log likelihood and log prior
#' @slot label_switch length-one logical vector indicating whether label-switching occurs (possibly an overfit model)
#' @slot mcmc.params An object of class 'McmcParams'
#' @slot .internal.constraint Constraint on parameters. For internal use only.
setClass("BatchModel", contains="MixtureModel")

#' The 'MarginalModel' class
#'
#' Run marginal MCMC simulation
#' @slot k An integer value specifying the number of latent classes.
#' @slot hyperparams An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @slot theta the means of each component and batch
#' @slot sigma2 the variances of each component and batch
#' @slot nu.0 the shape parameter for sigma2
#' @slot sigma2.0 the rate parameter for sigma2
#' @slot pi mixture probabilities which are assumed to be the same for all batches
#' @slot mu overall mean
#' @slot tau2 overall variance
#' @slot data the data for the simulation.
#' @slot data.mean the empirical means of the components
#' @slot data.prec the empirical precisions
#' @slot z latent variables
#' @slot zfreq table of latent variables
#' @slot probz n x k matrix of probabilities
#' @slot logprior log likelihood of prior: log(p(sigma2.0)p(nu.0)p(mu))
#' @slot loglik log likelihood: \eqn{\sum p_k \Phi(\theta_k, \sigma_k)}
#' @slot mcmc.chains an object of class 'McmcChains' to store MCMC samples
#' @slot batch a vector of the different batch numbers
#' @slot batchElements a vector labeling from which batch each observation came from
#' @slot modes the values of parameters from the iteration which maximizes log likelihood and log prior
#' @slot mcmc.params An object of class 'McmcParams'
#' @slot label_switch length-one logical vector indicating whether label-switching occurs (possibly an overfit model)
#' @slot .internal.constraint Constraint on parameters. For internal use only.
setClass("MarginalModel", contains="MixtureModel")


#' Constructor for list of single-batch models
#'
#' An object of class MarginalModel is constructed for each k, creating a list of
#' MarginalModels.
#'
#' @param data numeric vector of average log R ratios
#' @param k numeric vector indicating the number of mixture components for each model
#' @param mcmc.params an object of class \code{McmcParams}
#' @param ... additional arguments passed to \code{Hyperparameters}
#' @seealso \code{\link{MarginalModel}} \code{\link{BatchModelList}}
#' @return a list. Each element of the list is a \code{BatchModel}
#' @examples
#' mlist <- MarginalModelList(data=y(MarginalModelExample), k=1:4)
#' mcmcParams(mlist) <- McmcParams(iter=1, burnin=1, nStarts=0)
#' mlist2 <- posteriorSimulation(mlist)
#' @export
MarginalModelList <- function(data=numeric(), k=numeric(),
                              mcmc.params=McmcParams(),
                              ...){
  .Deprecated("See SingleBatchModelList")
  model.list <- vector("list", length(k))
  for(i in seq_along(k)){
    hypp <- Hyperparameters(k=k[i], ...)
    model.list[[i]] <- MarginalModel(data=data, k=k[i], mcmc.params=mcmc.params,
                                     hypp=hypp)
  }
  model.list
}

#' Create an object for running marginal MCMC simulations.
#' @examples
#'      model <- MarginalModel(data=rnorm(10), k=1)
#' @param data the data for the simulation.
#' @param k An integer value specifying the number of latent classes.
#' @param hypp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mcmc.params An object of class 'McmcParams'
#' @return An object of class 'MarginalModel'
#' @export
MarginalModel <- function(data=numeric(), k=3, hypp, mcmc.params){
  .Deprecated("See SingleBatchModel")
  if(missing(hypp)) hypp <- Hyperparameters(k=3)
  batch <- rep(1L, length(data))
  if(missing(k)){
    k <- k(hypp)
  } else{
    k(hypp) <- k
  }
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
                modes=list(),
                mcmc.params=mcmc.params,
                label_switch=FALSE,
                marginal_lik=as.numeric(NA),
                .internal.constraint=5e-4,
                .internal.counter=0L)
  object <- startingValues(object)
}

#' @rdname mu-method
#' @aliases mu,MarginalModel-method
#' @export
setMethod("mu", "MarginalModel", function(object){
  .Deprecated("see SingleBatchModel")
  object@mu
})

#' @rdname tau2-method
#' @aliases tau2,MarginalModel-method
setMethod("tau2", "MarginalModel", function(object) object@tau2)

setMethod("show", "MarginalModel", function(object) callNextMethod())

setMethod("computeMeans", "MarginalModel", function(object){
  .Deprecated("See SingleBatchModel")
  compute_means(object)
})

setMethod("computeVars", "MarginalModel", function(object){
  .Deprecated("See SingleBatchModel")
  compute_vars(object)
})

setReplaceMethod("tau2", "MarginalModel", function(object, value){
  .Deprecated("See SingleBatchModel")
  object@tau2 <- value
  object
})

setReplaceMethod("mu", "MarginalModel", function(object, value){
  .Deprecated("See SingleBatchModel")
  object@mu <- value
  object
})


#' @rdname bic-method
#' @aliases bic,MarginalModel-method
setMethod("bic", "MarginalModel", function(object){
  .Deprecated("See SingleBatchModel")
  object <- useModes(object)
  ## K: number of free parameters to be estimated
  ##   - component-specific parameters:  theta, sigma2   (3 x k(model))
  ##   - mixing probabilities:  k-1
  ##   - length-one parameters: mu, tau2, sigma2.0, nu.0             +4
  K <- 2*k(object) + (k(object)-1) + 4
  n <- length(y(object))
  -2*(log_lik(object) + logPrior(object)) + K*(log(n) - log(2*pi))
})

#' @rdname theta-method
#' @aliases theta,MarginalModel-method
setMethod("theta", "MarginalModel", function(object) object@theta)

#' @rdname sigma2-method
#' @aliases sigma2,MarginalModel-method
setMethod("sigma2", "MarginalModel", function(object) object@sigma2)

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
  object
})

setMethod("computeModes", "MarginalModel", function(object){
  .computeModesMarginal(object)
})


setMethod("showMeans", "MarginalModel", function(object){
  paste(round(theta(object), 2), collapse=", ")
})

setMethod("showSigmas", "MarginalModel", function(object){
  paste(round(sqrt(sigma2(object)), 2), collapse=", ")
})

setMethod("tablez", "MarginalModel", function(object) table(z(object)))

setMethod("updateMultinomialProb", "MarginalModel", function(object){
  update_multinomialPr(object)
})

setMethod("updateMultinomialProb", "BatchModel", function(object){
  update_multinomialPr_batch(object)
})

setMethod("computeLoglik", "BatchModel", function(object){
  compute_loglik_batch(object)
})

setMethod("computeLoglik", "MarginalModel", function(object){
  loglik(object)
})

setReplaceMethod("sigma2", "BatchModel", function(object, value){
  rownames(value) <- uniqueBatch(object)
  object@sigma2 <- value
  object
})

setMethod("showSigmas", "BatchModel", function(object){
  sigmas <- round(sqrt(sigma2(object)), 2)
  sigmas <- c("\n", paste0(t(cbind(sigmas, "\n")), collapse="\t"))
  sigmas <- paste0("\t", sigmas[2])
  sigmas <- paste0("\n", sigmas[1])
  sigmas
})

#' @rdname sigma2-method
#' @aliases sigma2,BatchModel-method
setMethod("sigma2", "BatchModel", function(object) {
  s2 <- object@sigma2
  ##s2 <- matrix(s2, nBatch(object), k(object))
  rownames(s2) <- uniqueBatch(object)
  s2
})


setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab[uniqueBatch(object), , drop=FALSE]
})


setMethod("sigmaMean", "BatchModel", function(object) {
  mns <- colMeans(sigmac(object))
  mns <- matrix(mns, nBatch(object), k(object))
  rownames(mns) <- uniqueBatch(object)
  mns
})

#' @rdname tau2-method
#' @aliases tau2,BatchModel-method
setMethod("tau2", "BatchModel", function(object) object@tau2)

setReplaceMethod("tau2", "BatchModel", function(object, value){
  object@tau2 <- value
  object
})

#' @rdname theta-method
#' @aliases theta,BatchModel-method
setMethod("theta", "BatchModel", function(object) {
  b <- object@theta
  ##b <- matrix(b, nBatch(object), k(object))
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
  cat("     loglik (s)  :", round(log_lik(object), 1), "\n")
  cat("     logprior (s):", round(logPrior(object), 1), "\n")
})

setMethod("tablez", "BatchModel", function(object){
  tab <- table(batch(object), z(object))
  tab <- tab[uniqueBatch(object), , drop=FALSE]
  tab
})

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,MarginalModel-method marginalLikelihood,MarginalModel,ANY-method
setMethod("marginalLikelihood", "MarginalModel",
          function(model, params=mlParams()) {
            .ml_singlebatch(model, params)
          })

#' @rdname marginalLikelihood-method
#' @aliases marginalLikelihood,BatchModel-method marginalLikelihood,BatchModel,ANY-method
setMethod("marginalLikelihood", "BatchModel",
          function(model, params=mlParams()){
            .ml_batchmodel(model, params)
          })

setMethod("runBurnin", "MarginalModel", function(object){
  mcmc_marginal_burnin(object, mcmcParams(object))
})

setMethod("runBurnin", "BatchModel", function(object){
  mcmc_batch_burnin(object, mcmcParams(object))
})

setMethod("runMcmc", "MarginalModel", function(object){
  mcmc_marginal(object, mcmcParams(object))
})

setMethod("runMcmc", "BatchModel", function(object){
  mcmc_batch(object, mcmcParams(object))
})

setMethod("computeModes", "BatchModel", function(object){
  .computeModesBatch(object)
})

setMethod("sortComponentLabels", "MarginalModel", function(model){
  reorderSingleBatch(model)  
})

setMethod("isOrdered", "BatchModel", function(object){
  .ordered_thetas_multibatch(object)
})


#' @rdname ggplot-functions
setMethod("ggMultiBatch", "BatchModel", function(model, bins){
  .gg_multibatch(model, bins)
})

setMethod("densitiesCluster", "BatchModel", function(object){
  .Deprecated("See MultiBatchModel")
  dens <- densities(object)
  modes <- dens[["modes"]]
  if(length(modes(object)) > 0){
    object <- useModes(object)
  }
  ix <- order(mu(object))
  thetas <- mu(object)[ix]
  names(thetas) <- ix
  if(length(modes) == 1) {
    km <- setNames(rep(1, length(thetas)), seq_along(thetas))
  } else {
    if(length(thetas) > length(modes)){
      km <- kmeans(thetas, centers=modes)$cluster
      names(km) <- seq_along(thetas)
    } else {
      ## length modes = length thetas
      km <- setNames(seq_along(thetas), seq_along(thetas))
    }
  }
  nclusters <- length(unique(km))
  batch <- dens[["batch"]]
  batch.list <- split(batch, km)
  batch2 <- batch.list
  for(j in seq_along(batch.list)){
    tmp <- batch.list[[j]]
    total <- matrix(0, nrow(tmp[[1]]), ncol(tmp[[1]]))
    for(l in seq_along(tmp)){
      total <- total+tmp[[l]]
    }
    batch2[[j]] <- total
  }
  batch <- batch2
  ##batch <- lapply(batch.list, function(x) do.call("+", x))
  component <- dens[["component"]]
  component.list <- split(component, km)
  length(component.list) == nclusters
  component <- lapply(component.list, function(x) rowSums(do.call(cbind, x)))
  overall <- rowSums(do.call(cbind, component))
  modes <- findModes(dens$quantiles, overall) ## should be the same
  list(batch=batch, component=component, overall=overall, modes=modes,
       quantiles=dens$quantiles,
       clusters=km, data=dens$data)
})

setMethod("densitiesCluster", "MultiBatchModel", function(object){
  dens <- densities(object)
  modes <- dens[["modes"]]
  if(length(modes(object)) > 0){
    object <- useModes(object)
  }
  ix <- order(mu(object))
  thetas <- mu(object)[ix]
  names(thetas) <- ix
  if(length(modes) == 1) {
    km <- setNames(rep(1, length(thetas)), seq_along(thetas))
  } else {
    if(length(thetas) > length(modes)){
      km <- kmeans(thetas, centers=modes)$cluster
      names(km) <- seq_along(thetas)
    } else {
      ## length modes = length thetas
      km <- setNames(seq_along(thetas), seq_along(thetas))
    }
  }
  nclusters <- length(unique(km))
  batch <- dens[["batch"]]
  batch.list <- split(batch, km)
  batch2 <- batch.list
  for(j in seq_along(batch.list)){
    tmp <- batch.list[[j]]
    total <- matrix(0, nrow(tmp[[1]]), ncol(tmp[[1]]))
    for(l in seq_along(tmp)){
      total <- total+tmp[[l]]
    }
    batch2[[j]] <- total
  }
  batch <- batch2
  ##batch <- lapply(batch.list, function(x) do.call("+", x))
  component <- dens[["component"]]
  component.list <- split(component, km)
  length(component.list) == nclusters
  component <- lapply(component.list, function(x) rowSums(do.call(cbind, x)))
  overall <- rowSums(do.call(cbind, component))
  modes <- findModes(dens$quantiles, overall) ## should be the same
  list(batch=batch, component=component, overall=overall, modes=modes,
       quantiles=dens$quantiles,
       clusters=km, data=dens$data)
})


doKmeans <- function(thetas, modes) length(thetas) > length(modes) && length(modes) > 1

setMethod("densitiesCluster", "MarginalModel", function(object){
  .Deprecated("See SingleBatchModel")
  dens <- densities(object)
  modes <- dens[["modes"]]
  if(length(modes(object)) > 0){
    object <- useModes(object)
  }
  ix <- order(theta(object))
  thetas <- theta(object)[ix]
  names(thetas) <- ix
  if(length(modes) == 1) {
    km <- setNames(rep(1, length(thetas)), seq_along(thetas))
  } else {
    if(length(thetas) > length(modes)){
      km <- kmeans(thetas, centers=modes)$cluster
      names(km) <- seq_along(thetas)
    } else {
      ## length modes = length thetas
      km <- setNames(seq_along(thetas), seq_along(thetas))
    }
  }
  nclusters <- length(unique(km))
  component <- dens[["component"]]
  component.list <- split(component, km)
  length(component.list) == nclusters
  component <- lapply(component.list, function(x) rowSums(do.call(cbind, x)))
  overall <- rowSums(do.call(cbind, component))
  modes <- findModes(dens$quantiles, overall) ## should be the same
  list(component=component, overall=overall, modes=modes,
       clusters=km, quantiles=dens$quantiles, data=dens$data)
})

setMethod("densitiesCluster", "SingleBatchModel", function(object){
  dens <- densities(object)
  modes <- dens[["modes"]]
  if(length(modes(object)) > 0){
    object <- useModes(object)
  }
  ix <- order(theta(object))
  thetas <- theta(object)[ix]
  names(thetas) <- ix
  if(length(modes) == 1) {
    km <- setNames(rep(1, length(thetas)), seq_along(thetas))
  } else {
    if(length(thetas) > length(modes)){
      km <- kmeans(thetas, centers=modes)$cluster
      names(km) <- seq_along(thetas)
    } else {
      ## length modes = length thetas
      km <- setNames(seq_along(thetas), seq_along(thetas))
    }
  }
  nclusters <- length(unique(km))
  component <- dens[["component"]]
  component.list <- split(component, km)
  length(component.list) == nclusters
  component <- lapply(component.list, function(x) rowSums(do.call(cbind, x)))
  overall <- rowSums(do.call(cbind, component))
  modes <- findModes(dens$quantiles, overall) ## should be the same
  list(component=component, overall=overall, modes=modes,
       clusters=km, quantiles=dens$quantiles, data=dens$data)
})

#' Constructor for DensityModel class
#'
#' DensityModel constructor has been deprecated.  
#' @seealso See \code{\link{ggSingleBatch}} and \code{\link{ggMultiBatch}} for visualization
#' @param object see \code{showMethods(DensityModel)}
#' @param merge Logical.  Whether to use kmeans clustering to cluster
#' the component means using the estimated modes from the overall
#' density as the centers for the \code{kmeans} function.
#' @return An object of class 'DensityModel'
#' @export
DensityModel <- function(object, merge=FALSE){
 .Deprecated()
  if(!missing(object)){
    if(!merge){
      dens <- densities(object)
    } else{
      dens <- densitiesCluster(object)
    }
    quantiles <- dens$quantiles
    clusters <- dens$clusters
    component <- dens$component
    overall <- dens$overall
    modes <- dens$modes
    data <- dens$data
  } else{
    ## object not provided
    component <- list()
    overall <- list()
    modes <- numeric()
    clusters <- numeric()
    quantiles <- numeric()
    data <- numeric()
  }
  if(isMarginalModel(object)){
    obj <- new("DensityModel", component=component, overall=overall, modes=modes,
               clusters=clusters, data=data, quantiles=quantiles)
    return(obj)
  }
  if(!missing(object)){
    batch <- dens$batch
  } else batch <- list()
  new("DensityBatchModel", batch=batch, component=component, overall=overall,
      modes=modes, clusters=clusters, data=data, quantiles=quantiles)
}

setMethod("component", "DensityModel", function(object) object@component)

#' @rdname batch-method
#' @aliases batch,DensityModel-method
setMethod("batch", "DensityModel", function(object) object@batch)

setMethod("overall", "DensityModel", function(object) object@overall)

#' @rdname modes-method
#' @aliases modes,DensityModel-method
setMethod("modes", "DensityModel", function(object) object@modes)

#' @rdname k-method
#' @aliases k,DensityModel-method
setMethod("k", "DensityModel", function(object) length(component(object)))

setMethod("show", "DensityModel", function(object){
  cat("An object of class 'DensityModel'\n")
  cat("   component densities:  list of", length(component(object)), "vectors \n")
  cat("   overall density:  vector of length",  length(overall(object)), "\n")
  cat("   modes: ",  paste(sort(round(modes(object), 2)), collapse=", "), "\n")
})

setMethod("show", "DensityBatchModel", function(object){
  cat("An object of class 'DensityBatchModel'\n")
  cat("   batch: list of ", length(component(object)), "matrices \n")
  cat("   component densities:  list of", length(component(object)), "vectors \n")
  cat("   overall density:  vector of length",  length(overall(object)), "\n")
  cat("   modes: ",  paste(sort(round(modes(object), 2)), collapse=", "), "\n")
})


.plotMarginal <- function(x, y, ...){
  .Deprecated("See ggSingleBatch and ggMultiBatch")
  object <- x
  args <- list(...)
  if(!"breaks" %in% names(args)){
    L <- length(y)
    hist(y, breaks=L/10,
         col="gray",
         border="gray",  xlab="y",
         freq=FALSE, ...)
  } else {
    hist(y, col="gray", border="gray", xlab="y",
         freq=FALSE, ...)
  }
  comp <- component(x)
  K <- length(comp)
  cols <- brewer.pal(max(K, 3), "Set1")
  cols <- cols[seq_len(K)]
  drawdens <- function(y, x, col, lwd=1) lines(x, y, col=col, lwd=lwd)
  mapply(drawdens, y=comp, col=cols, MoreArgs=list(x=quantiles(object), lwd=2))
  marg <- overall(object)
  lines(quantiles(object), marg, col="black", lwd=1)
}

.plotBatch <- function(x, y, show.batch=TRUE, ...){
  object <- x
  args <- list(...)
  if(!"breaks" %in% names(args)){
    L <- length(y)
    hist(y, breaks=L/10,
         col="gray",
         border="gray",  xlab="y",
         freq=FALSE, ...)
  } else {
    hist(y, col="gray", border="gray", xlab="y",
         freq=FALSE, ...)
  }
  comp <- component(x)
  K <- length(comp)
  cols <- brewer.pal(max(K, 3), "Set1")
  cols <- cols[seq_len(K)]
  drawdens <- function(y, x, col, lwd=1) lines(x, y, col=col, lwd=lwd)
  if(show.batch){
    batches <- batch(object)
    for(j in seq_along(batches)){
      B <- batches[[j]]
      apply(B, 2, drawdens, col=cols[j], x=quantiles(object), lwd=1)
    }
  }
  mapply(drawdens, y=comp, col=cols, MoreArgs=list(x=quantiles(object), lwd=2))
  marg <- overall(object)
  lines(quantiles(object), marg, col="black", lwd=1)
}

dens <- function(x, mean, sd, p1, p2){
  p.x <- p1*p2*dnorm(x, mean, sd)
  p.x
}

##drawEachComponent <- function(x, batches, thetas, sds, P, batchPr,
##cols, plot=TRUE){



findModes <- function(quantiles, x){ ##quantiles, density
  signs <- sign(diff(x))
  signs <- c(signs[1], signs)
  inflect <- c(0, cumsum(diff(signs) < 0))
  indices <- tapply(seq_along(inflect), inflect, max)
  indices <- indices[-length(indices)]
  modes <- quantiles[indices]
  modes
}

#' @rdname y-method
#' @aliases y,DensityModel-method
setMethod("y", "DensityModel", function(object) object@data)

#' @rdname plot
#' @aliases plot,DensityModel,numeric-method
setMethod("plot", "DensityModel", function(x, y, ...){
  .plotMarginal(x, x@data, ...)
})

#' @rdname plot
#' @aliases plot,MarginalModel,ANY-method
setMethod("plot", "MarginalModel", function(x, y, ...){
  object <- DensityModel(x)
  .plotMarginal(x=object, y=object@data, ...)
  return(object)
})

#' @rdname plot
#' @aliases plot,BatchModel,ANY-method
setMethod("plot", "BatchModel", function(x, y, show.batch=TRUE, ...){
  object <- DensityModel(x)
  .plotBatch(object, oned(x), show.batch, ...)
  return(object)
})

#' @rdname plot
#' @aliases plot,BatchModel,ANY-method
setMethod("plot", "DensityBatchModel", function(x, show.batch=TRUE, ...){
  .plotBatch(x, x@data, show.batch, ...)
})


setMethod("hist", "MixtureModel", function(x, ...){
  op <- par(las=1)
  yy <- y(x)
  ##if(rsample > 0) yy <- sample(yy, rsample)
  args <- list(...)
  if(!"breaks" %in% names(args)){
    L <- length(yy)
    hist(yy, breaks=L/10,
         col="gray",
         border="gray",  xlab="y",
         freq=FALSE, ...)
  } else {
    hist(yy, col="gray", border="gray", xlab="y",
         freq=FALSE, ...)
  }
  par(op)
})

setMethod("densities", "BatchModel", function(object){
  quantiles <- seq(min(observed(object)), max(observed(object)),  length.out=250)
  ## use the modes if available
  if(length(modes(object)) > 0){
    object <- useModes(object)
  }
  data <- y(object)
  thetas <- theta(object)
  sds <- sigma(object)
  P <- p(object)
  P <- matrix(P, nBatch(object), k(object), byrow=TRUE)
  rownames(P) <- uniqueBatch(object)
  ix <- order(mu(object))
  thetas <- thetas[, ix, drop=FALSE]
  sds <- sds[, ix, drop=FALSE]
  P <- P[, ix, drop=FALSE]
  batchPr <- table(batch(object))/length(y(object))
  dens.list <- batchDensities(quantiles, uniqueBatch(object),
                              thetas, sds, P, batchPr)
  component <- lapply(dens.list, rowSums)
  overall <- rowSums(do.call(cbind, component))
  modes <- findModes(quantiles, overall)
  clusters <- seq_len(k(object))
  names(clusters) <- ix
  list(batch=dens.list, component=component, overall=overall, modes=modes,
       clusters=clusters, quantiles=quantiles, data=data)
})

setMethod("densities", "MultiBatchModel", function(object){
  quantiles <- seq(min(observed(object)), max(observed(object)),  length.out=250)
  ## use the modes if available
  if(length(modes(object)) > 0){
    object <- useModes(object)
  }
  data <- y(object)
  thetas <- theta(object)
  sds <- sigma(object)
  P <- p(object)
  P <- matrix(P, nBatch(object), k(object), byrow=TRUE)
  rownames(P) <- uniqueBatch(object)
  ix <- order(mu(object))
  thetas <- thetas[, ix, drop=FALSE]
  sds <- sds[, ix, drop=FALSE]
  P <- P[, ix, drop=FALSE]
  batchPr <- table(batch(object))/length(y(object))
  dens.list <- batchDensities(quantiles, uniqueBatch(object),
                              thetas, sds, P, batchPr)
  component <- lapply(dens.list, rowSums)
  overall <- rowSums(do.call(cbind, component))
  modes <- findModes(quantiles, overall)
  clusters <- seq_len(k(object))
  names(clusters) <- ix
  list(batch=dens.list, component=component, overall=overall, modes=modes,
       clusters=clusters, quantiles=quantiles, data=data)
})


setMethod("densities", "MarginalModel", function(object){
  quantiles <- seq(min(observed(object)), max(observed(object)),  length.out=250)
  ## use the modes if available
  if(length(modes(object)) > 0){
    object <- useModes(object)
  }
  data <- y(object)
  thetas <- theta(object)
  sds <- sigma(object)
  P <- p(object)
  ix <- order(thetas)
  thetas <- thetas[ix]
  sds <- sds[ix]
  P <- P[ix]
  dens.list <- vector("list", length(thetas))
  for(i in seq_along(dens.list)){
    dens.list[[i]] <- P[i]*dnorm(quantiles, mean=thetas[i], sd=sds[i])
  }
  ##component <- lapply(dens.list, rowSums)
  overall <- rowSums(do.call(cbind, dens.list))
  modes <- findModes(quantiles, overall)
  clusters <- seq_len(k(object))
  names(clusters) <- ix
  list(component=dens.list, overall=overall, modes=modes,
       clusters=clusters, quantiles=quantiles, data=data)
})

setMethod("densities", "SingleBatchModel", function(object){
  quantiles <- seq(min(observed(object)), max(observed(object)),  length.out=250)
  ## use the modes if available
  if(length(modes(object)) > 0){
    object <- useModes(object)
  }
  data <- y(object)
  thetas <- theta(object)
  sds <- sigma(object)
  P <- p(object)
  ix <- order(thetas)
  thetas <- thetas[ix]
  sds <- sds[ix]
  P <- P[ix]
  dens.list <- vector("list", length(thetas))
  for(i in seq_along(dens.list)){
    dens.list[[i]] <- P[i]*dnorm(quantiles, mean=thetas[i], sd=sds[i])
  }
  ##component <- lapply(dens.list, rowSums)
  overall <- rowSums(do.call(cbind, dens.list))
  modes <- findModes(quantiles, overall)
  clusters <- seq_len(k(object))
  names(clusters) <- ix
  list(component=dens.list, overall=overall, modes=modes,
       clusters=clusters, quantiles=quantiles, data=data)
})

setMethod("densities", "SingleBatchPooled", function(object){
  quantiles <- seq(min(observed(object)), max(observed(object)),  length.out=250)
  ## use the modes if available
  if(length(modes(object)) > 0){
    object <- useModes(object)
  }
  data <- y(object)
  thetas <- theta(object)
  sd.pooled <- sigma(object)
  P <- p(object)
  ix <- order(thetas)
  thetas <- thetas[ix]
  P <- P[ix]
  dens.list <- vector("list", length(thetas))
  for(i in seq_along(dens.list)){
    dens.list[[i]] <- P[i]*dnorm(quantiles, mean=thetas[i], sd=sd.pooled)
  }
  ##component <- lapply(dens.list, rowSums)
  overall <- rowSums(do.call(cbind, dens.list))
  modes <- findModes(quantiles, overall)
  clusters <- seq_len(k(object))
  names(clusters) <- ix
  list(component=dens.list, overall=overall, modes=modes,
       clusters=clusters, quantiles=quantiles, data=data)
})

#' @rdname clusters-method
#' @aliases clusters,DensityModel-method
#' @return k-means clustering of the component means using the modes as centers.
setMethod("clusters", "DensityModel", function(object) {
  .Deprecated()
  object@clusters
})

##setMethod("quantiles", "DensityModel", function(object){
##  object@quantiles
##})


#' Create an object of class 'HyperparametersBatch' for the
#' batch mixture model
#'
#' @param k  length-one integer vector specifying number of components
#' (typically 1 <= k <= 4)
#' @param mu.0 length-one numeric vector of the of the normal prior
#' for the component means.
#' @param tau2.0 length-one numeric vector of the variance for the normal
#' prior of the component means
#' @param eta.0 length-one numeric vector of the shape parameter for
#' the Inverse Gamma prior of the component variances, tau2_h.  The
#' shape parameter is parameterized as 1/2 * eta.0.  In the batch
#' model, tau2_h describes the inter-batch heterogeneity of means for
#' component h.
#' @param m2.0 length-one numeric vector of the rate parameter for the
#' Inverse Gamma prior of the component variances, tau2_h.  The rate
#' parameter is parameterized as 1/2 * eta.0 * m2.0.  In the batch
#' model, tau2_h describes the inter-batch heterogeneity of means for
#' component h.
#' @param alpha length-k numeric vector of the shape parameters for
#' the dirichlet prior on the mixture probabilities
#' @param beta length-one numeric vector for the parameter of the
#' geometric prior for nu.0 (nu.0 is the shape parameter of the
#' Inverse Gamma sampling distribution for the component-specific
#' variances. Together, nu.0 and sigma2.0 model inter-component
#' heterogeneity in variances.).  beta is a probability and must be
#' in the interval [0,1].
#' @param a length-one numeric vector of the shape parameter for the
#' Gamma prior used for sigma2.0 (sigma2.0 is the shape parameter of
#' the Inverse Gamma sampling distribution for the component-specific
#' variances).
#' @param b a length-one numeric vector of the rate parameter for the
#' Gamma prior used for sigma2.0 (sigma2.0 is the rate parameter of
#' the Inverse Gamma sampling distribution for the component-specific
#' variances)
#' @return An object of class HyperparametersBatch
#' @examples
#' HyperparametersBatch(k=3)
#'
#' @export
HyperparametersBatch <- function(k=3L,
                                 mu.0=0,
                                 tau2.0=100,
                                 eta.0=1800,
                                 m2.0=1/60,
                                 alpha,
                                 beta=0.1, ## mean is 1/10
                                 a=1.8,
                                 b=6){
  .Deprecated("See HyperparametersMultiBatch")
  if(missing(alpha)) alpha <- rep(1, k)
  new("HyperparametersBatch",
      k=as.integer(k),
      mu.0=mu.0,
      tau2.0=tau2.0,
      eta.0=eta.0,
      m2.0=m2.0,
      alpha=alpha,
      beta=beta,
      a=a,
      b=b)
}

#' Create an object of class 'HyperparametersMarginal' for the
#' marginal mixture model
#'
#' @param k  length-one integer vector specifying number of components
#' (typically 1 <= k <= 4)
#' @param mu.0  length-one numeric vector of the mean for the normal
#' prior of the component means
#' @param tau2.0 length-one numeric vector of the variance for the normal
#' prior of the component means
#' @param eta.0 length-one numeric vector of the shape parameter for
#' the Inverse Gamma prior of the component variances.  The shape
#' parameter is parameterized as 1/2 * eta.0.
#' @param m2.0 length-one numeric vector of the rate parameter for
#' the Inverse Gamma prior of the component variances.  The rate
#' parameter is parameterized as 1/2 * eta.0 * m2.0.
#' @param alpha length-k numeric vector of the shape parameters for
#' the dirichlet prior on the mixture probabilities
#' @param beta length-one numeric vector for the parameter of the
#' geometric prior for nu.0 (nu.0 is the shape parameter of the
#' Inverse Gamma sampling distribution for the component-specific
#' variances).  beta is a probability and must be in the interval
#' [0,1].
#' @param a length-one numeric vector of the shape parameter for the
#' Gamma prior used for sigma2.0 (sigma2.0 is the shape parameter of
#' the Inverse Gamma sampling distribution for the component-specific
#' variances)
#' @param b a length-one numeric vector of the rate parameter for the
#' Gamma prior used for sigma2.0 (sigma2.0 is the rate parameter of
#' the Inverse Gamma sampling distribution for the component-specific
#' variances)
#'
#' @return An object of class HyperparametersMarginal
#' @examples
#' HyperparametersMarginal(k=3)
#'
#' @export
HyperparametersMarginal <- function(k=0L,
                                    mu.0=0,
                                    tau2.0=100,
                                    eta.0=1,
                                    m2.0=0.1,
                                    alpha,
                                    beta=0.1, ## mean is 1/10
                                    a=1.8,
                                    b=6){
  .Deprecated()
  if(missing(alpha)) alpha <- rep(1, k)
  ##if(missing(tau2)) tau2 <- rep(1, k)
  new("HyperparametersMarginal",
      k=as.integer(k),
      mu.0=mu.0,
      tau2.0=tau2.0,
      eta.0=eta.0,
      m2.0=m2.0,
      alpha=alpha,
      beta=beta,
      a=a,
      b=b)
}
