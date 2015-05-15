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
##
#### compute distance for a given permutation of columns
##.computeDistanceOnePerm <- function(mc, column.permutation, modes){
##  param.sds <- sapply(modes, sd)
##  .computeDistOneBatch(th=theta(mc)[, column.permutation],
##                       s2=sigma2(mc)[, column.permutation],
##                       mus=mu(mc)[, column.permutation],
##                       tau2s=tau2(mc)[, column.permutation],
##                       P=p(mc)[, column.permutation],
##                       modes=modes,
##                       param.sds=param.sds)
##}
##
##setMethod("computeDistance", "BatchModel", function(object){
##  modal.params <- modes(object)
##  ##param.sds <- sapply(modal.params, function(x) sd(x[1,]))
##  mc <- mcmcChains(object)
##  th <- theta(mc)
##  s2 <- sigma2(mc)
##  ## ix is the number of possible orderings for each batch
##  ix <- permutations(k(object), k(object))##
##  nr <- nrow(th)
##  nc <- nrow(ix)
##  Dlist <- vector("list", nBatch(bmodel))
##  ## iterate over batches
##  for(b in seq_len(nBatch(object))){
##    mc2 <- mc
##    batch.index <- seq(b, nBatch(object)*k(object), by=nBatch(object))
##    theta(mc2) <- th[, batch.index]
##    sigma2(mc2) <- s2[, batch.index]
##    D <- matrix(NA, nr, nc)
##    m.params <- list(theta=modal.params[["theta"]][b,],
##                     sigma2=modal.params[["sigma2"]][b,],
##                     mu=modal.params[["mu"]],
##                     tau2=modal.params[["tau2"]],
##                     mixprob=modal.params[["mixprob"]])
##    ## iterate over all possible permutations
##    for(j in 1:nrow(ix)){
##      J <- ix[j, ]
##      D[, j] <- .computeDistanceOnePerm(mc=mc2, column.permutation=J,
##                                        modes=m.params)
##    }
##    Dlist[[b]] <- D
##  }
##  Dlist
##})
##
##setMethod("switchLabels", "BatchModel", function(object){
##  Dlist <- computeDistance(object)
##  mc <- mcmcChains(object)
##  warn <- FALSE
##  for(b in seq_along(Dlist)){
##    D <- Dlist[[b]]
##    ordering_index <- apply(D, 1, which.min)
##    if(all(ordering_index == 1)) next()
##    warn <- TRUE
##    batch.index <- seq(b, nBatch(object)*k(object), by=nBatch(object))
##    mc2 <- mc
##    theta(mc2) <- theta(mc)[, batch.index]
##    sigma2(mc2) <- sigma2(mc)[, batch.index]
##    perms <- permutations(k(object), k(object))
##    tab <- as.integer(names(table(ordering_index)))
##    tab <- tab[tab!=1]
##    for(i in seq_along(tab)){
##      mcmc.index <- which(ordering_index == tab[i])
##      j <- perms[tab[i], ]
##      ## rewrite batch.index in mc from the permuted index in mc2
##      theta(mc)[mcmc.index, batch.index] <- theta(mc2)[mcmc.index, j]
##      sigma2(mc)[mcmc.index, batch.index] <- sigma2(mc2)[mcmc.index, j]
##      p(mc)[mcmc.index, ] <- p(mc)[mcmc.index, j]
##      mu(mc)[mcmc.index,] <- mu(mc)[mcmc.index, j]
##      tau2(mc)[mcmc.index,] <- tau2(mc)[mcmc.index, j]
##    }
##    mcmcChains(object) <- mc
##  }
##  if(warn) warning("Label switching occurred. Posterior probabilities for z may be incorrect")
##  object
##})
##



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

orderCols <- function(th, x){
  ix <- t(apply(th, 1, order))
  for(j in 1:nrow(th)){
    k <- ix[j, ]
    x[j, ] <- x[j, k]
  }
  x
}

# Berkhof et al., Statistica Sinica 2003
setMethod("reorderComponents", "MixtureModel", function(object, new_levels){
  ##
  ## Choose another ordering, then force the other to switch by fixing z
  ##
##  mcmcp <- McmcParams(iter=20, burnin=0)
##  mcmcChains(object) <- McmcChains(object, mcmcp)
  zz <- z(object)
  zz <- factor(as.integer(factor(zz, levels=new_levels)), levels=seq_len(k(object)))
  z(object) <- zz
  dataMean(object) <- computeMeans(object)
  dataPrec(object) <- 1/computeVars(object)
##  for(s in 1:iter(mcmcp)){
##    object <- reducedGibbsZ(object, TRUE)
##    object <- moveChain(object, s)
##  }
  object
})
##
##setMethod("switchLabels", "MarginalModel", function(object){
##  D <- computeDistance(object)
##  ordering_index <- apply(D, 1, which.min)
##  if(all(ordering_index == 1)) return(object)
##  warning("Label switching occurred. Posterior probabilities for z may be incorrect")
##  mc <- mcmcChains(object)
##  perms <- permutations(k(object), k(object))
##  tab <- as.integer(names(table(ordering_index)))
##  tab <- tab[tab!=1]
##  for(i in seq_along(tab)){
##    mcmc.index <- which(ordering_index == tab[i])
##    j <- perms[tab[i], ]
##    theta(mc)[mcmc.index, ] <- theta(mc)[mcmc.index, j]
##    sigma2(mc)[mcmc.index, ] <- sigma2(mc)[mcmc.index, j]
##    p(mc)[mcmc.index, ] <- p(mc)[mcmc.index, j]
##  }
##  mcmcChains(object) <- mc
##  object
##})

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

##.updateLabels <- function(object){
##  modal.params <- modes(object)
##  param.sds <- sapply(modal.params, sd)
##  th <- theta(object)
##  s2 <- sigma2(object)
##  P <- p(object)
##  ix <- permutations(k(object), k(object))##
##  nc <- nrow(ix)
##  D <- rep(NA, nc)
##  ##
##  ## Compute distance to the modes for the current ordering (1,2,3)
##  ## at each iteration of the chain.
##  ##
##  ## subtrace a vector from each row of chain matrix
##  D[1] <- .computeDistanceMarginal(th, s2, P, modal.params, param.sds)
##  if(FALSE) plot.ts(D[,1], col="gray")
##  for(j in 2:nc){
##    J <- ix[j, ]
##    browser()
##    D[j] <- .computeDistanceMarginal(th[J], s2[J], P[J], modal.params, param.sds)
##  }
##  reordering <- ix[which.min(D)]
##  theta(object) <- th[reordering]
##  sigma2(object) <- s2[reordering]
##  p(object) <- P[reordering]
##  z(object) <- factor(z(object), levels=reordering)
##  dataMean(object) <- dataMean(object)[reordering]
##  dataPrec(object) <- dataPrec(object)[reordering]
##  object
##}


##
##setMethod("updateLabels", "MarginalModel", function(object){
##  .updateLabels(object)
##})
##
##
##setMethod("computeDistance", "MarginalModel", function(object){
##  modal.params <- modes(object)
##  param.sds <- sapply(modal.params, sd)
##  mc <- mcmcChains(object)
##  th <- theta(mc)
##  s2 <- sigma2(mc)
##  P <- p(mc)
##  ix <- permutations(k(object), k(object))##
##  nr <- nrow(th)
##  nc <- nrow(ix)
##  D <- matrix(NA, nr, nc)
##  ##
##  ## Compute distance to the modes for the current ordering (1,2,3)
##  ## at each iteration of the chain.
##  ##
##  ## subtrace a vector from each row of chain matrix
##  D[, 1] <- .computeDistanceMarginal(th, s2, P, modal.params, param.sds)
##  if(FALSE) plot.ts(D[,1], col="gray")
##  for(j in 2:nrow(ix)){
##    J <- ix[j, ]
##    D[, j] <- .computeDistanceMarginal(th[, J], s2[, J], P[, J], modal.params, param.sds)
##  }
##  D
##})
##
