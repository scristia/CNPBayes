#' @include methods-SingleBatchModel.R
NULL

setMethod("computeLoglik", "SBPt", function(object){
  loglik_pooled(object)
})

## ad-hoc constructor
##.SingleBatchPooled <- function(dat=numeric(),
##                              hp=Hyperparameters(),
##                              mp=McmcParams(iter=1000, burnin=1000,
##                                            thin=10, nStarts=4)){
##  model <- .SingleBatchPooled2(dat=dat, hp=hp, mp=mp)
##  model
##}


.SBPt <- function(dat=numeric(),
                               hp=Hyperparameters(),
                               mp=McmcParams(iter=1000, burnin=1000,
                                             thin=10, nStarts=4)){
  if(length(dat) == 0){
    sb <- .empty_singlebatch_model(hp, mp)
    sbp <- as(sb, "SBPt")
    return(sbp)
  }
  K <- k(hp)
  mu <- rnorm(1, median(dat), sd(dat))
  tau2 <- 1/rgamma(1, 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
  p <- rdirichlet(1, alpha(hp))[1, ]
  theta <- sort(rnorm(k(hp), mu, sqrt(tau2)))
  ##nu.0 <- 3.5
  ##sigma2.0 <- 0.25
  ##sigma2 <- 1/rgamma(k(hp), 0.5 * nu.0, 0.5 * nu.0 * sigma2.0)
  ##sigma2 <- abs(rnorm(1, 0.1, 0.1))^2
  nu.0 <- rgeom(1, betas(hp))
  sigma2.0 <- rgamma(1, 0.5*a(hp), 0.5*a(hp)*b(hp))
  sigma2 <- 1/rgamma(1, 0.5 * nu.0, 0.5 * nu.0 * sigma2.0)
  object <- new("SBPt",
                k=as.integer(K),
                hyperparams=hp,
                theta=theta,
                sigma2=sigma2,
                mu=mu,
                tau2=tau2,
                nu.0=nu.0,
                sigma2.0=sigma2.0,
                pi=p,
                data=dat,
                data.mean=numeric(K),
                data.prec=numeric(K),
                z=integer(length(dat)),
                zfreq=integer(K),
                probz=matrix(0, length(dat), K),
                logprior=numeric(1),
                loglik=numeric(1),
                mcmc.chains=McmcChains(),
                batch=rep(1L, length(dat)),
                batchElements=1L,
                modes=list(),
                mcmc.params=mp,
                label_switch=FALSE,
                marginal_lik=as.numeric(NA),
                .internal.constraint=5e-4,
                .internal.counter=0L)
  chains(object) <- McmcChains(object)
  object
}

.computeModesSBPt <- function(object){
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

setMethod("computeModes", "SBPt", function(object){
  .computeModesSBPt(object)
})

## We initialize the z's to be all zeros, then z is updated as the first step of
## the mcmc. However, if the first update of z results in some components not
## being observed, then z will not be updated and will stay at zero. This
## creates NaNs for the thetas and several other parameters. To try to
## circumvent this issue, we have a while loop that simulates a model and runs 5
## iterations burnin. If the parameter values are valid, we stop. If not, we
## simulate a new model.
#' @export
#' @rdname SingleBatchModel2
SBPt <- function(dat=numeric(),
                              hp=Hyperparameters(),
                              mp=McmcParams(iter=1000, burnin=1000,
                                            thin=10, nStarts=4)){
  if(length(dat) == 0){
    return(.SBPt(dat, hp, mp))
  }
  iter <- 0
  validZ <- FALSE
  ##
  ## Burnin with the more flexible SingleBatch model to obtain starting values,
  ## then convert to SingleBatchPooled class
  ##
  mp.tmp <- McmcParams(iter=0, burnin=burnin(mp), thin=1, nStarts=1)
  while(!validZ){
    ##
    ## Burnin with SB model
    ##
    sb <- .SingleBatchModel2(dat, hp, mp.tmp)
    sb <- runBurnin(sb)
    tabz <- table(z(sb))
    if(length(tabz) == k(hp)) validZ <- TRUE
    iter <- iter + 1
    if(iter > 50) stop("Trouble initializing valid model. Try increasing the burnin")
  }
  sbp <- as(sb, "SBPt")
  dfr(sbp) <- hp@dfr
  mcmcParams(sbp) <- mp
  sbp <- sortComponentLabels(sbp)
  log_lik(sbp) <- loglik_pooled(sbp)
  sbp
}

setValidity("SingleBatchPooled", function(object){
  s2 <- sigma2(object)
  if(length(s2) != 1){
    return("sigma2 slot should be length-one numeric vector")
  }
  TRUE
})

setMethod("u", "SBPt", function(object) object@u )

#' @rdname dfr-method
#' @aliases dfr,SBPt-method
setMethod("dfr", "SBPt", function(object) object@hyperparams@dfr )

setReplaceMethod("dfr", "SBPt", function(object, value) {
                     object@hyperparams@dfr <- value
                     object@u <- rchisq(length(y(object)), value)
                     object
})
