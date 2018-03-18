#' @include Defunct-classes.R
NULL

#' @rdname mu-method
#' @aliases mu,MarginalModel-method
#' @export
#' @rdname Defunct-functions
setMethod("mu", "MarginalModel", function(object){
  .Defunct("see SingleBatchModel")
  object@mu
})

#' @aliases tau2,MarginalModel-method
#' @rdname Defunct-functions
setMethod("tau2", "MarginalModel", function(object) object@tau2)

setMethod("show", "MarginalModel", function(object) callNextMethod())

setMethod("computeMeans", "MarginalModel", function(object){
  .Defunct("See SingleBatchModel")
  compute_means(object)
})

setMethod("computeVars", "MarginalModel", function(object){
  .Defunct("See SingleBatchModel")
  compute_vars(object)
})

setReplaceMethod("tau2", "MarginalModel", function(object, value){
  .Defunct("See SingleBatchModel")
  object@tau2 <- value
  object
})

setReplaceMethod("mu", "MarginalModel", function(object, value){
  .Defunct("See SingleBatchModel")
  object@mu <- value
  object
})


#' @aliases bic,MarginalModel-method
#' @rdname Defunct-functions
setMethod("bic", "MarginalModel", function(object){
  .Defunct("SingleBatchModel")
  object <- useModes(object)
  ## K: number of free parameters to be estimated
  ##   - component-specific parameters:  theta, sigma2   (3 x k(model))
  ##   - mixing probabilities:  k-1
  ##   - length-one parameters: mu, tau2, sigma2.0, nu.0             +4
  K <- 2*k(object) + (k(object)-1) + 4
  n <- length(y(object))
  -2*(log_lik(object) + logPrior(object)) + K*(log(n) - log(2*pi))
})

#' @rdname Defunct-functions
#' @aliases theta,MarginalModel-method
setMethod("theta", "MarginalModel", function(object) object@theta)

#' @rdname Defunct-functions
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
  .computeModesSingleBatch(object)
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

#' @rdname Defunct-functions
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

#' @rdname Defunct-functions
#' @aliases tau2,BatchModel-method
setMethod("tau2", "BatchModel", function(object) object@tau2)

setReplaceMethod("tau2", "BatchModel", function(object, value){
  object@tau2 <- value
  object
})

#' @rdname Defunct-functions
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

#' @param params list of parameters for computing marginal likelihood
#' @param model MarginalModel
#' @rdname Defunct-functions
#' @aliases marginalLikelihood,MarginalModel-method marginalLikelihood,MarginalModel,ANY-method
setMethod("marginalLikelihood", "MarginalModel",
          function(model, params=mlParams()) {
            .ml_singlebatch(model, params)
          })

#' @rdname Defunct-functions
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


#' @param bins length-one numeric vector specifying number of bins for plotting
#' @rdname Defunct-functions
setMethod("ggMultiBatch", "BatchModel", function(model, bins){
  .gg_multibatch(model, bins)
})

setMethod("densitiesCluster", "BatchModel", function(object){
  .Defunct("See MultiBatchModel")
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
  .Defunct("See SingleBatchModel")
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

#' DensityModel constructor has been deprecated.  
#' @seealso See \code{\link{ggSingleBatch}} and \code{\link{ggMultiBatch}} for visualization
#' @return An object of class 'DensityModel'
#' @export
#' @rdname Defunct-functions
DensityModel <- function(object, merge=FALSE){
 .Defunct()
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
 isMarginalModel <- function(object) NULL
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

#' @rdname Defunct-functions
#' @aliases batch,DensityModel-method
setMethod("batch", "DensityModel", function(object) object@batch)

setMethod("overall", "DensityModel", function(object) object@overall)

#' @rdname Defunct-functions
#' @aliases modes,DensityModel-method
setMethod("modes", "DensityModel", function(object) object@modes)

#' @rdname Defunct-functions
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
  .Defunct("See ggSingleBatch and ggMultiBatch")
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

#' @rdname Defunct-functions
#' @aliases y,DensityModel-method
setMethod("y", "DensityModel", function(object) object@data)

#' @rdname Defunct-functions
#' @aliases plot,DensityModel,numeric-method
setMethod("plot", "DensityModel", function(x, y, ...){
  .plotMarginal(x, x@data, ...)
})

#' @rdname Defunct-functions
#' @aliases plot,MarginalModel,ANY-method
setMethod("plot", "MarginalModel", function(x, y, ...){
  object <- DensityModel(x)
  .plotMarginal(x=object, y=object@data, ...)
  return(object)
})

#' @rdname Defunct-functions
#' @aliases plot,BatchModel,ANY-method
setMethod("plot", "BatchModel", function(x, y, show.batch=TRUE, ...){
  object <- DensityModel(x)
  .plotBatch(object, oned(x), show.batch, ...)
  return(object)
})

#' @rdname Defunct-functions
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


setGeneric("clusters", function(object) standardGeneric("clusters"))

setMethod("clusters", "DensityModel", function(object) {
  .Defunct()
  object@clusters
})

##setMethod("quantiles", "DensityModel", function(object){
##  object@quantiles
##})


#' @return An object of class HyperparametersBatch
#' @export
#' @rdname Defunct-functions
HyperparametersBatch <- function(k=3L,
                                 mu.0=0,
                                 tau2.0=100,
                                 eta.0=1800,
                                 m2.0=1/60,
                                 alpha,
                                 beta=0.1, ## mean is 1/10
                                 a=1.8,
                                 b=6){
  .Defunct("See HyperparametersMultiBatch")
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
#' @export
#' @rdname Defunct-functions
HyperparametersMarginal <- function(k=0L,
                                    mu.0=0,
                                    tau2.0=100,
                                    eta.0=1,
                                    m2.0=0.1,
                                    alpha,
                                    beta=0.1, ## mean is 1/10
                                    a=1.8,
                                    b=6){
  .Defunct()
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


#' @param data numeric vector of average log R ratios
#' @param mcmc.params a \code{McmcParams} object
#' @param ... additional arguments to \code{HyperparametersBatch}
#' @return a list. Each element of the list is a \code{BatchModel}
#' @seealso \code{\link{BatchModel}}.  For single-batch data, use 
#' @export
#' @rdname Defunct-functions
BatchModelList <- function(data=numeric(),
                           k=numeric(),
                           batch,
                           mcmc.params=McmcParams(),
                           ...){
  .Defunct("See MultiBatchModelList")
  model.list <- vector("list", length(k))
  for(i in seq_along(k)){
    hypp <- HyperparametersBatch(k=k[i], ...)
    model.list[[i]] <- BatchModel(data=data, k=k[i], batch=batch,
                                  mcmc.params=mcmc.params,
                                  hypp=hypp)
  }
  model.list
}

#' @return An object of class `BatchModel`
#' @export
#' @rdname Defunct-functions
BatchModel <- function(data=numeric(),
                       k=3,
                       batch,
                       hypp,
                       mcmc.params){
  .Defunct("See MultiBatchModel2")
  if(missing(batch)) batch <- as.integer(factor(rep("a", length(data))))
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  if(missing(hypp)) hypp <- HyperparametersBatch(k=k)
  if(missing(k) & !missing(hypp)){
    k <- k(hypp)
  }
  mcmc.chains <- McmcChains()
  bf <- factor(batch)
  batch <- as.integer(bf)
  ub <- unique(batch)
  ##ix <- order(batch)
  ix <- seq_along(batch)
  nbatch <- setNames(as.integer(table(batch)), levels(bf))
  B <- length(ub)
  if(B==1 && length(data) > 0){
    if(missing(hypp)) hypp <- HyperparametersMarginal(k=k)
    zz <- as.integer(factor(numeric(k)))
    zfreq <- as.integer(table(zz))
    obj <- SingleBatchModel(data, k, hypp, mcmc.params)
    return(obj)
  }
  if(k == 1) {
    if(missing(hypp)) hypp <- HyperparametersBatch(k=1)
    obj <- UnivariateBatchModel(data, k, batch, hypp, mcmc.params)
    return(obj)
  }
  if(missing(hypp)) hypp <- HyperparametersBatch(k=k)
  zz <- integer(length(data))
  zfreq <- as.integer(table(zz))
  if(length(data) != length(batch)) {
    stop("batch vector must be the same length as data")
  }
  obj <- new("BatchModel",
             k=as.integer(k),
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
             label_switch=FALSE,
             .internal.constraint=5e-4,
             .internal.counter=0L)
  obj <- startingValues(obj)
  obj
}

#' @rdname Defunct-functions
setMethod("posteriorSimulation", c("MixtureModel", "integer"),
          function(object, k) {
            ##.Defunct("Method is deprecated for signature 'MixtureModel, integer'.  Use MarginalModelList or BatchModelList prior to posteriorSimulation")
            stop("Specifying k not allowed.  See MutliBatchModelList or SingleBatchModelList for creating a list object.")
        if (length(k) > 1) {
          mlist <- vector("list", length(k))
          for (i in seq_along(k)) {
            k(object) <- k[i]
            mlist[[i]] <- .posteriorSimulation(object)
          }
          mlist
        } else {
          k(object) <- k
          .posteriorSimulation(object)
        }
    }
)

#' @rdname Defunct-functions
setMethod("posteriorSimulation", c("MixtureModel", "numeric"),
          function(object, k) {
            stop("Specifying k not allowed.  See MultiBatchModelList or SingleBatchModelList for creating a list object.")
            posteriorSimulation(object, as.integer(k))
    })

setMethod("densitiesCluster", "BatchModel", function(object){
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


#' Instantiates an instance of 'DensityModel' (or 'DensityBatchModel')
#' from a MarginalModel or BatchModel object. See the corresponding
#' class for additional details and examples.
#' @seealso \code{\link{DensityModel-class}} \code{\link{kmeans}}
#' @param object see \code{showMethods(DensityModel)}
#' @param merge Logical.  Whether to use kmeans clustering to cluster
#' the component means using the estimated modes from the overall
#' density as the centers for the \code{kmeans} function.
#' @return An object of class 'DensityModel'
#' @export
#' @rdname Defunct-functions
DensityModel <- function(object, merge=FALSE){
  .Defunct("see ggMixture")
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
  isMarginalModel <- function(x) NULL
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

#' @rdname Defunct-functions
#' @aliases batch,DensityModel-method
setMethod("batch", "DensityModel", function(object) object@batch)

setMethod("overall", "DensityModel", function(object) object@overall)

#' @rdname Defunct-functions
#' @aliases modes,DensityModel-method
setMethod("modes", "DensityModel", function(object) object@modes)

#' @rdname Defunct-functions
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

batchDensities <- function(x, batches, thetas, sds, P, batchPr){
  mlist <- list()
  for(j in seq_len(ncol(thetas))){
    mlist[[j]] <- .batchdens(x, batches, thetas[, j], sds[, j], P[, j], batchPr)
  }
  ##marginal <- do.call(cbind, mlist)
  mlist
}

.batchdens <- function(x, batches, thetas, sds, p1, p2){
  marginal <- matrix(NA, length(x), length(batches))
  for(b in seq_along(batches)){
    marginal[, b] <- dens(x, thetas[b], sds[b], p1[b], p2[b])
  }
  marginal
}

findModes <- function(quantiles, x){ ##quantiles, density
  signs <- sign(diff(x))
  signs <- c(signs[1], signs)
  inflect <- c(0, cumsum(diff(signs) < 0))
  indices <- tapply(seq_along(inflect), inflect, max)
  indices <- indices[-length(indices)]
  modes <- quantiles[indices]
  modes
}

#' @rdname Defunct-functions
#' @aliases y,DensityModel-method
setMethod("y", "DensityModel", function(object) object@data)

#' @rdname Defunct-functions
#' @aliases plot,DensityModel,numeric-method
setMethod("plot", "DensityModel", function(x, y, ...){
  .plotMarginal(x, x@data, ...)
})

#' @rdname Defunct-functions
#' @aliases plot,MarginalModel,ANY-method
setMethod("plot", "MarginalModel", function(x, y, ...){
  object <- DensityModel(x)
  .plotMarginal(x=object, y=object@data, ...)
  return(object)
})

#' @rdname Defunct-functions
#' @aliases plot,BatchModel,ANY-method
setMethod("plot", "BatchModel", function(x, y, show.batch=TRUE, ...){
  object <- DensityModel(x)
  .plotBatch(object, oned(x), show.batch, ...)
  return(object)
})

#' @rdname Defunct-functions
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

setMethod("clusters", "DensityModel", function(object) object@clusters)
setMethod("quantiles", "DensityModel", function(object) object@quantiles)


#' @rdname Defunct-functions
setMethod("ggMultiBatch", "MultiBatchModel", function(model, bins){
  .Defunct("See ggMixture")
  .gg_multibatch(model, bins)
})

#' @rdname Defunct-functions
setMethod("ggMultiBatch", "MultiBatchPooled", function(model, bins){
  .Defunct("See ggMixture")
  .gg_multibatch_pooled(model, bins)
})

setGeneric("ggSingleBatch", function(model, bins) standardGeneric("ggSingleBatch"))

#' @rdname Defunct-functions
#' @aliases ggSingleBatch,MarginalModel-method
setMethod("ggSingleBatch", "MarginalModel", function(model, bins){
  .Defunct("See ggMixture")
  .gg_singlebatch(model, bins)
})

#' @rdname Defunct-functions
#' @aliases ggSingleBatch,SingleBatchModel-method
setMethod("ggSingleBatch", "SingleBatchModel", function(model, bins){
  .Defunct("See ggMixture")
  .gg_singlebatch(model, bins)
})

#' @rdname Defunct-functions
setMethod("ggMultiBatch", "MultiBatchCopyNumber", function(model, bins){
  .Defunct("See ggMixture")
  .gg_multibatch_copynumber(model, bins)
})

#' @rdname Defunct-functions
#' @aliases ggSingleBatch,SingleBatchCopyNumber-method
setMethod("ggSingleBatch", "SingleBatchCopyNumber", function(model, bins){
  .Defunct("See ggMixture")
  .gg_singlebatch_copynumber(model, bins)
})

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
  mns <- c("\n", paste0(t(cbind(thetas, "\n")), collapse="\t"))
  mns <- paste0("\t", mns[2])
  mns <- paste0("\n", mns[1])
  mns
})



#' @param batch.file the name of a file contaning RDS data to be read in.
#' @param plate a vector containing the labels  from which batch each observation came from.
#' @param y in memory data
#' @param ntiles number of tiles in a batch
#' @param THR threshold above which to merge batches in Kolmogorov-Smirnov test.
#' @return Tile labels for each observation
#' @export
#' @rdname Defunct-functions
downsample <- function(batch.file, plate, y, ntiles=250, THR=0.1){
  .Defunct(downsample)
  if(file.exists(batch.file)){
    batches <- readRDS(batch.file)
  } else {
    message("merging plates... ")
    if(missing(plate)) stop("batch.file does not exist.  Chemistry plate must be specified")
    batches <- collapseBatch(y, plate, THR=THR)
    saveRDS(batches, file=batch.file)
  }
  ds <- downSampleEachBatch(y, ntiles, batches)
  if(length(unique(ds$batch)) > 12){
    batches <- collapseBatch(ds$y, ds$batch, THR=1e-3)
    ds$batch <- batches
  }
  if(any(table(ds$batch) <= 25)){
    tab <- table(ds$batch)
    keep <- ds$batch %in% names(tab)[tab > 25]
    ds$y <- ds$y[keep]
    ds$batch <- ds$batch[keep]
  }
  ds
}



#' @param nt the number of observations per batch
#' @return Tile labels for each observation
#' @seealso \code{\link[dplyr]{ntile}}
#' @export
#' @rdname Defunct-functions
downSampleEachBatch <- function(y, nt, batch){
  .Defunct("see tileMedians")
  ## NULL out these two variables to avoid NOTE about
  ## no visible binding for global variable
  yy <- y
  x <- obs.index <- NULL
  logratio <- NULL
  ##
  ## split observations by batch
  ##
  yb <- split(y, batch)
  indices <- split(seq_along(y), batch)
  S <- vector("list", length(yb))

  for (i in 1:length(yb)) {
    x <- yb[[i]]
    batch.id <- names(yb)[i]
    obs.index <- indices[[i]]
    ##
    ## observations can be quite different within a tile (e.g., homozygous deletions)
    ##
    tiles <- ntile(x, nt) %>% as.tibble %>%
      mutate(x=x) %>%
      set_colnames(c("tile", "logratio"))

    tiles2 <- tiles %>%
      group_by(tile) %>%
      summarize(spread=abs(diff(range(logratio))),
                n=n()) %>%
      arrange(-spread)
    ## Split tiles with large spread into multiple tiles
    tiles.keep <- tiles2 %>%
      filter(spread < 0.05)
    tiles.drop <- tiles2 %>%
      filter(spread >= 0.05)
    newtiles <- tiles %>% filter(tile %in% tiles.drop$tile)
    newtiles$tile <- seq_len(nrow(newtiles)) + max(tiles$tile)
    tiles3 <- filter(tiles, tile %in% tiles.keep$tile) %>%
      bind_rows(newtiles)
    tiles3$tile <- paste(batch.id, tiles3$tile, sep="_")
    S[[i]] <- tiles3
  }
  tiles <- do.call(bind_rows, S)
  tiles
}


#' Create an object for running hierarchical MCMC simulations.
#' @param batch a vector of the different batch numbers (must be sorted)
#' @param hypp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @return An object of class `MultiBatchModel`
#' @export
#' @rdname Defunct-functions
MultiBatchModel <- function(data=numeric(),
                            k=3,
                            batch,
                            hypp,
                            mcmc.params){
  .Defunct("See MultiBatchModel2")
  if(missing(batch)) batch <- as.integer(factor(rep("a", length(data))))
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  if(missing(hypp)) hypp <- HyperparametersMultiBatch(k=k)
  if(missing(k) & !missing(hypp)){
    k <- k(hypp)
  }
  mcmc.chains <- McmcChains()
  bf <- factor(batch)
  batch <- as.integer(bf)
  ub <- unique(batch)
  ##ix <- order(batch)
  ix <- seq_along(batch)
  nbatch <- setNames(as.integer(table(batch)), levels(bf))
  B <- length(ub)
  if(B==1 && length(data) > 0){
    if(missing(hypp)) hypp <- HyperparametersMarginal(k=k)
    zz <- as.integer(factor(numeric(k)))
    zfreq <- as.integer(table(zz))
    obj <- SingleBatchModel(data, k, hypp, mcmc.params)
    return(obj)
  }
  if(k == 1) {
    if(missing(hypp)) hypp <- HyperparametersMultiBatch(k=1)
    obj <- UnivariateBatchModel(data, k, batch, hypp, mcmc.params)
    return(obj)
  }
  if(missing(hypp)) hypp <- HyperparametersMultiBatch(k=k)
  zz <- integer(length(data))
  zfreq <- as.integer(table(zz))
  if(length(data) != length(batch)) {
    stop("batch vector must be the same length as data")
  }
  obj <- new("MultiBatchModel",
             k=as.integer(k),
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
             label_switch=FALSE,
             .internal.constraint=5e-4,
             .internal.counter=0L)
  obj <- startingValues(obj)
  obj
}

setMethod("computePrior", "MarginalModel", function(object){
  .Defunct("See SingleBatchModel")
  compute_logprior(object)
})

setMethod("computePrior", "BatchModel", function(object){
  compute_logprior_batch(object)
})

#' @rdname Defunct-functions
#' @aliases posteriorSimulation,list-method
setMethod("posteriorSimulation", "list",
          function(object) {
            .Defunct("See gibbs")
            params <- psParams()
            results <- vector("list", length(object))
            for(i in seq_along(results)){
              results[[i]] <- .posteriorSimulation(object[[i]], params)
            }
            ncomp <- sapply(results, k)
            if(is(results[[1]], "SingleBatchModel")){
              label <- "SB"
            } else label <- "MB"
            names(results) <- paste0(label, ncomp)
            ##isnull <- sapply(results, is.null)
            ##results <- results[!isnull]
            results
          })

#' @param x a \code{DensityModel}-derived object, or a
#' \code{MixtureModel}-derived object.
#' numeric vector of the one-dimensional summaries for a given copy
#' number polymorphism. If \code{x} is a \code{MixtureModel}, \code{y}
#' is ignored.
#' @param show.batch a logical. If true, batch specific densities
#' will be plotted.
#' @return A plot showing the density estimate
#' @export
#' @rdname Defunct-functions
setGeneric("plot")

#' DensityModel constructor and methods are Defunct
#'
#' @export
#' @docType methods
#' @rdname Defunct-functions
setGeneric("clusters", function(object) standardGeneric("clusters"))

#' @return A single proportion for a \code{MarginalModel} or a vector of proportions, one for each batch for a \code{BatchModel}
#' @export
#' @rdname Defunct-functions
setGeneric("labelSwitching", function(object, merge=TRUE) standardGeneric("labelSwitching"))

.posteriorSimulation <- function(post, params=psParams()){
  .Defunct("see .posteriorSimulation2")
  if(nStarts(post) > 1){
    post <- multipleStarts2(post)
  }
  if(burnin(post) > 0 ){
    post <- runBurnin(post)
  }
  if(!isOrdered(post)) label_switch(post) <- TRUE
  post <- sortComponentLabels(post)
  if( iter(post) < 1 ) return(post)
  post <- runMcmc(post)
  modes(post) <- computeModes(post)
  if(isOrdered(post)){
    label_switch(post) <- FALSE
    return(post)
  }
  ## not ordered: try additional MCMC simulations
  label_switch(post) <- TRUE
  post <- sortComponentLabels(post)
  ## reset counter for posterior probabilities
  post@probz[] <- 0
  post <- runMcmc(post)
  modes(post) <- computeModes(post)
  ##mcmcParams(post) <- mp.orig
  if(isOrdered(post)){
    label_switch(post) <- FALSE
    return(post)
  }
  label_switch(post) <- TRUE
  if(params[["warnings"]]) {
    ##
    ## at this point, we've tried to run the twice after burnin and we still
    ## have mixing. Most likely, we are fitting a model with k too big
    warning("label switching: model k=", k(post))
  }
  post <- sortComponentLabels(post)
  post
}



#' @rdname Defunct-functions
#' @aliases labelSwitching,MixtureModel-method
#' @export
setMethod("labelSwitching", "MixtureModel", 
    function(object, merge=TRUE) {
        # put together a map indicating which component a component
        # is merged into, if merging happens
        if (merge) {
            k.orig <- k(object)
            merged <- DensityModel(object, merge=TRUE)
            k.merged <- k(merged)
            comp_map <- clusters(merged)
            message("Merged from ", k.orig,
                    " components to ", k.merged,
                    " components")
        } else {
            # if merge==FALSE then this is an identity map
            comp_map <- clusters(object)
        }

        # a vector showing the batch of each observation
        # note that for a MarginalModel, all observations have batch=1
        batches <- unique(batch(object))

        # get the number of components
        components <- k(object)

        # empty vector for storing proportion of relabeling for 
        # each batch
        prop_relabeled <- numeric(length(batches))

        # grab the thetas for components/batch at each MCMC iteration
        thetas_all <- theta(chains(object))

        for (batch in batches) {
            # get indices for the given batch
            ind <- (batch - 1) * components + 1:components

            # get thetas at each iteration corresponding to the batch
            thetas_batch <- thetas_all[, ind]

            # calculate the proportion of relabeling for a given batch
            prop_relabeled[batch] <- 1 - relabeling(thetas_batch, 
                                                    comp_map)
        }

        return(prop_relabeled)
    }
    )

MultiBatchModelList <- function(data=numeric(),
                           k=numeric(),
                           batch,
                           mcmc.params=McmcParams(),
                           ...){
  model.list <- vector("list", length(k))
  for(i in seq_along(k)){
    hypp <- HyperparametersMultiBatch(k=k[i], ...)
    model.list[[i]] <- MultiBatchModel(data=data, k=k[i], batch=batch,
                                  mcmc.params=mcmc.params,
                                  hypp=hypp)
  }
  model.list
}

#' @param drop Not used.
#' @return An object of class 'BatchModel'
#' @aliases [,BatchModel-method [,BatchModel,ANY-method [,BatchModel,ANY,ANY-method [,BatchModel,ANY,ANY,ANY-method
#' @param i integer
#' @param j integer
#' @docType methods
#' @rdname Defunct-functions
setMethod("[", "BatchModel", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    y(x) <- y(x)[i]
    z(x) <- z(x)[i]
    batch(x) <- batch(x)[i]
  }
  x
})

#' @export
#' @rdname Defunct-functions
SingleBatchModel <- function(data=numeric(), k=3, hypp, mcmc.params){
  .Defunct("SingleBatchModel2")
  if(missing(hypp)) hypp <- Hyperparameters(k=3)
  batch <- rep(1L, length(data))
  if(missing(k)){
    k <- k(hypp)
  } else{
    k(hypp) <- k
  }
  if(missing(mcmc.params)) mcmc.params <- McmcParams(iter=1000, burnin=100)
  if(missing(hypp)) hypp <- HyperparametersSingleBatch(k=k)
  nbatch <- setNames(as.integer(table(batch)), levels(batch))
  zz <- sample(seq_len(k), length(data), replace=TRUE)
  zfreq <- as.integer(table(zz))
  object <- new("SingleBatchModel",
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


newMarginalModel <- function(object){
  .Defunct("newSingleBatchModel")
  mp <- mcmcParams(object)
  object2 <- SingleBatchModel(y(object), k=k(object), mcmc.params=mp,
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
  log_lik(object2) <- log_lik(object)
  logPrior(object2) <- logPrior(object)
  modes(object2) <- modes(object)
  object2
}

#' @export
#' @rdname Defunct-functions
multiBatchDensities <- function(model){
  .Defunct("ggMixture")
  dnorm_poly_multibatch(model)
}

#' @rdname collapseBatch-method
#' @aliases collapseBatch,BatchModel-method
setMethod("collapseBatch", "BatchModel", function(object){
  collapseBatch(y(object), as.character(batch(object)))
})
