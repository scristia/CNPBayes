#' Create simulated batch data for testing.
#'
#' @param N number of observations
#' @param p a vector indicating probability of membership to each component
#' @param theta a matrix of means.  Columns are components and rows are batches.
#' @param sds a matrix of standard deviations.  Columns are components and rows are batches.
#' @param batch a vector of labels indication from which batch each simulation should come from
#' @param zz a vector indicating latent variable membership. Can be omitted.
#' @param df length-1 numeric vector for the t-distribution degrees of freedom
#' @return An object of class 'MultiBatchModel'
simulateBatchData <- function(N=2500, p, theta, sds, batch, zz, df=10){
  ## order ys by batch
  if(!is.matrix(p)) p <- matrix(p, nrow(theta), ncol(theta), byrow=TRUE)
  if(ncol(p) != ncol(theta)) stop("length of p must be same as ncol(theta)")
  if(!all(rowSums(p)==1)) stop("elements of p must sum to 1")
  if(missing(batch)) {
    batch <- rep(1L, N)
  } else {
    batch <- as.integer(factor(batch))
  }
  batch <- sort(batch)
  if(missing(zz)) {
    zz <- simulateZ(N, p)
  }
  yy <- rep(NA, N)
  ub <- unique(batch)
  rownames(theta) <- rownames(sds) <- ub
  for(b in ub){
    index <- which(batch==b)
    nn <- length(index)
    cn <- zz[index]
    mu <- theta[b, ]
    s <- sds[b, ]
    yy[index] <- rnorm(nn, mu[cn], s[cn])/sqrt(rchisq(nn, df)/df)
  }
  ##ix <- order(batch)
  ##object <- MultiBatchModel(yy, batch=batch, k=ncol(theta))
  object <- MultiBatchModel2(dat=yy, batches=batch,
                             hpList(k=ncol(theta))[["MB"]])
  z(object) <- as.integer(factor(zz))
  ##
  ## Must initialize the slots independently (must point to different
  ## locations in memory)
  ##
  theta(object) <- computeMeans(object)
  dataMean(object) <- computeMeans(object)
  mu(object) <- colMeans(dataMean(object))
  sigma2(object) <- computeVars(object)
  dataPrec(object) <- 1/computeVars(object)
  counts <- tableBatchZ(object)
  P <- counts/rowSums(counts)
  p(object) <- P
  log_lik(object) <- computeLoglik(object)
  object
}

#' Create simulated data for testing.
#'
#' @param N number of observations
#' @param p a vector indicating probability of membership to each component
#' @param theta a vector of means, one per component
#' @param sds a vector of standard deviations, one per component
#' @param df length-1 numeric vector for the t-distribution degrees of freedom
#' @return An object of class 'SingleBatchModel'
simulateData <- function(N, p, theta, sds, df=10){
  theta <- matrix(theta, nrow=1)
  sds <- matrix(sds, nrow=1)
  p <- matrix(p, nrow=1)
  object <- simulateBatchData(N=N,
                              p=p,
                              theta=theta,
                              sds=sds,
                              df=df,
                              batch=rep(1L, N))
  return(object)
}

simulateZ <- function(N, p){
  P <- rdirichlet(N, p[1, ])
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

