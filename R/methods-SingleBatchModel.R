#' @include methods-MixtureModel.R
NULL

#' Constructors for SB and SBP models
#'
#' Create objects of class SingleBatchModel or SingleBatchPooled
#'
#' @param dat numeric vector of average log R ratios
#' @param mp an object of class \code{McmcParams}
#' @param hp an object of class \code{Hyperparameters}
#' @seealso \code{\link{MultiBatchModel2}}
#' @return An instance of \code{MultiBatchModel}
#' @examples
#' SingleBatchModel2()
#' SingleBatchModel2(dat=rnorm(100), hpList(k=2)[["SB"]])
#' SingleBatchPooled()
#' @export
SingleBatchModel2 <- function(dat=numeric(),
                              hp=Hyperparameters(),
                              mp=McmcParams(iter=1000, burnin=1000,
                                            thin=10, nStarts=4)){
  sb <- MB(dat=dat, hp=hp, mp=mp, batches=rep(1L, length(dat)))
##  if(length(dat) == 0){
##    return(.SingleBatchModel2(dat, hp, mp))
##  }
##  iter <- 0
##  validZ <- FALSE
##  mp.tmp <- McmcParams(iter=0, burnin=5, thin=1, nStarts=1)
##  while(!validZ){
##    sb <- .SingleBatchModel2(dat, hp, mp.tmp)
##    sb <- runBurnin(sb)
##    tabz <- table(z(sb))
##    if(length(tabz) == k(hp)) validZ <- TRUE
##    iter <- iter + 1
##    if(iter > 50) stop("Trouble initializing valid model")
##  }
##  mcmcParams(sb) <- mp
  sb
}

#' @export
SB <- SingleBatchModel2


getK <- function(object){
  hypp <- hyperParams(object)
  getK(hypp)
}

#' @rdname mu-method
#' @aliases mu,SingleBatchModel-method
#' @export
setMethod("mu", "SingleBatchModel", function(object) object@mu)



#' @rdname tau2-method
#' @aliases tau2,SingleBatchModel-method
setMethod("tau2", "SingleBatchModel", function(object) object@tau2)

setMethod("show", "SingleBatchModel", function(object) callNextMethod())

setMethod("computeMeans", "SingleBatchModel", function(object){
  compute_means(object)
})

setMethod("computeVars", "SingleBatchModel", function(object){
  compute_vars(object)
})



setReplaceMethod("tau2", "SingleBatchModel", function(object, value){
  object@tau2 <- value
  object
})


setReplaceMethod("mu", "SingleBatchModel", function(object, value){
  object@mu <- value
  object
})

#' @rdname bic-method
#' @aliases bic,SingleBatchModel-method
setMethod("bic", "SingleBatchModel", function(object){
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
#' @aliases theta,SingleBatchModel-method
setMethod("theta", "SingleBatchModel", function(object) object@theta)

#' @rdname sigma2-method
#' @aliases sigma2,SingleBatchModel-method
setMethod("sigma2", "SingleBatchModel", function(object) object@sigma2)



newSingleBatchModel <- function(object){
  mp <- mcmcParams(object)
  object2 <- SingleBatchModel2(dat=y(object),
                               mp=mp,
                               hp=hyperParams(object))
  theta(object2) <- theta(object)
  sigma2(object2) <- sigma2(object)
  p(object2) <- p(object)
  z(object2) <- z(object)
  object2@u <- u(object)
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

setMethod("relabel", "SingleBatchModel", function(object, zindex){
  object <- newSingleBatchModel(object)
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

.computeModesSingleBatch <- function(object){
  i <- argMax(object)
  if(length(i) == 0) i <- iter(object)
  mc <- chains(object)
  thetamax <- theta(mc)[i, ]
  sigma2max <- sigma2(mc)[i, ]
  pmax <- p(mc)[i, ]
  modes <- list(theta=thetamax,
                sigma2=sigma2max,
                mixprob=pmax,
                mu=mu(mc)[i],
                tau2=tau2(mc)[i],
                nu0=nu.0(mc)[i],
                sigma2.0=sigma2.0(mc)[i],
                zfreq=zFreq(mc)[i, ],
                loglik=log_lik(mc)[i],
                logprior=logPrior(mc)[i])
  modes
}

setMethod("computeModes", "SingleBatchModel", function(object){
  .computeModesSingleBatch(object)
})

setMethod("showMeans", "SingleBatchModel", function(object){
  paste(round(theta(object), 2), collapse=", ")
})


setMethod("showSigmas", "SingleBatchModel", function(object){
  paste(round(sqrt(sigma2(object)), 2), collapse=", ")
})

setMethod("tablez", "SingleBatchModel", function(object) table(z(object)))

modelOtherModes <- function(model, maxperm=5){
  kperm <- permnK(k(model), maxperm)
  model.list <- vector("list", nrow(kperm))
  for(i in seq_along(model.list)){
    model.list[[i]] <- relabel(model, kperm[i, ])
  }
  model.list
}

#' Compute the log bayes factor between models.
#'
#' Models of varying component sizes are compared. The log bayes factor is 
#' calculated comparing each set of two models by marginal likelihood, 
#' as computed by \code{marginalLikelihood}.
#' @param x the result of a call to \code{computeMarginalLik}.
#' @return Log Bayes factor comparing the two models with highest likelihood.
#' @export
logBayesFactor <- function(x) {
    k <- length(x)
    mat <- matrix(0, nrow=k, ncol=k, dimnames=list(names(x), names(x)))
    for (i in seq_len(length(x))) {
        x_i <- x[c(i:k, 0:(i-1))]
        diff <- x_i[1] - x_i
        mat[i, names(diff)] <- diff
    }

    return(mat)
}

setMethod("updateMultinomialProb", "SingleBatchModel", function(object){
  update_multinomialPr(object)
})

setMethod("computeLoglik", "SingleBatchModel", function(object){
  loglik(object)
})

reorderSingleBatchChains <- function(model){
  K <- k(model)
  if(K < 2) return(model)
  ##pstar=.blockUpdates(model, mcmcParams(mlist2[[2]]))
  ch <- chains(model)
  theta.ch <- thetac(model)
  ix <- apply(theta.ch, 1, function(x) paste0(order(x), collapse=","))
  tab.ix <- table(ix)
  ord <- names(tab.ix)[which.max(tab.ix)]
  ord <- as.integer(unlist(strsplit(ord, ",")))
  if(identical(ord, seq_len(K))){
    return(model)
  }
  theta.ch <- theta.ch[, ord]
  sigma.ch <- sigmac(model)
  sigma.ch <- sigma.ch[, ord]
  p.ch <- pic(model)
  p.ch <- p.ch[, ord]
  zfreq.ch <- zFreq(ch)

  ch@theta <- theta.ch
  ch@sigma2 <- sigma.ch^2
  ch@pi <- p.ch
  ch@zfreq <- zfreq.ch[, ord]
  chains(model) <- ch
  model
}

setMethod("updateZ", "SingleBatchModel", function(object){
  update_z(object)
})
