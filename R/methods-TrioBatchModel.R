.empty_batch_model <- function(hp, mp){
  K <- k(hp)
  B <- 0
  N <- 0
  obj <- new("TrioBatchModel",
             k=as.integer(K),
             hyperparams=hp,
             theta=matrix(NA, 0, K),
             sigma2=matrix(NA, 0, K),
             mu=numeric(K),
             tau2=numeric(K),
             nu.0=numeric(1),
             sigma2.0=numeric(1),
             pi=numeric(K),
             data=numeric(0),
             triodata=as_tibble(0),
             data.mean=matrix(NA, B, K),
             data.prec=matrix(NA, B, K),
             z=integer(0),
             zfreq=integer(K),
             probz=matrix(0, N, K),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=McmcChains(),
             mcmc.params=mp,
             batch=integer(0),
             batchElements=integer(0),
             label_switch=FALSE,
             marginal_lik=as.numeric(NA),
             .internal.constraint=5e-4,
             .internal.counter=0L)
  chains(obj) <- McmcChains(obj)
  obj
}

.TBM <- function(dat=numeric(),
                hp=HyperparametersTrios(),
                mp=McmcParams(iter=1000, thin=10,
                              burnin=1000, nStarts=4),
                batches=integer()){
  ## If the data is not ordered by batch,
  ## its a little harder to sort component labels
  dat2 <- tibble(y=dat$data, batch=batches)
  ub <- unique(batches)
  nbatch <- setNames(as.integer(table(batches)), ub)
  B <- length(ub)
  N <- nrow(dat$data)
  ## move to setValidity
  if(nrow(dat$data) != length(batches)) {
    stop("batch vector must be the same length as data")
  }
  K <- k(hp)
  ## mu_k is the average across batches of the thetas for component k
  ## tau_k is the sd of the batch means for component k
  mu <- sort(rnorm(k(hp), mu.0(hp), sqrt(tau2.0(hp))))
  tau2 <- 1/rgamma(k(hp), 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
  p <- rdirichlet(1, alpha(hp))[1, ]
  sim_theta <- function(mu, tau, B) sort(rnorm(B, mu, tau))
  . <- NULL
  thetas <- map2(mu, sqrt(tau2), sim_theta, B) %>%
    do.call(cbind, .) %>%
    apply(., 1, sort) %>%
    t
  if(K == 1) thetas <- t(thetas)
  nu.0 <- 3.5
  sigma2.0 <- 0.25
  sigma2s <- 1/rgamma(k(hp) * B, 0.5 * nu.0, 0.5 * nu.0 * sigma2.0) %>%
    matrix(B, k(hp))
  u <- rchisq(length(dat), hp@dfr)
  obj <- new("TrioBatchModel",
             k=as.integer(K),
             hyperparams=hp,
             theta=thetas,
             sigma2=sigma2s,
             mu=mu,
             tau2=tau2,
             nu.0=nu.0,
             sigma2.0=sigma2.0,
             pi=p,
             data=dat,
             triodata=dat2,
             u=u,
             data.mean=matrix(NA, B, K),
             data.prec=matrix(NA, B, K),
             z=sample(seq_len(K), N, replace=TRUE),
             zfreq=integer(K),
             probz=matrix(0, N, K),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=McmcChains(),
             mcmc.params=mp,
             batch=batches,
             batchElements=nbatch,
             label_switch=FALSE,
             marginal_lik=as.numeric(NA),
             .internal.constraint=5e-4,
             .internal.counter=0L)
  obj
}

#' Constructor for TrioBatchModel
#'
#' Initializes a TrioBatchModel, a container for storing data, parameters, and MCMC output for mixture models with batch- and component-specific means and variances.
#'
#' @param dat the data for the simulation.
#' @param batches an integer-vector of the different batches
#' @param hp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mp An object of class 'McmcParams'
#' @return An object of class `TrioBatchModel`
#' @export
TrioBatchModel <- function(dat=numeric(),
                             hp=HyperparametersTrios(),
                             mp=McmcParams(iter=1000, thin=10,
                                           burnin=1000, nStarts=4),
                             batches=integer()){
  if(length(dat) == 0){
    return(.empty_batch_model(hp, mp))
  }
  iter <- 0
  validZ <- FALSE
  mp.tmp <- McmcParams(iter=0, burnin=burnin(mp), thin=1, nStarts=1)
  while(!validZ){
    ##
    ## Burnin with MB model
    ##
    tbm <- .TBM(dat, hp, mp.tmp, batches)
    tbm <- runBurnin(tbm)
    tabz1 <- table(batch(tbm), z(tbm))
    tabz2 <- table(z(tbm))
    validZ <- length(tabz2) == k(hp) && all(tabz1 > 1)
    iter <- iter + 1
    if(iter == 50) {
      message("Trouble initializing a valid model. The number of components is likely too large")
      return(NULL)
    }
  }
  tbm2 <- sortComponentLabels(tbm)
  mcmcParams(tbm2) <- mp
  chains(tbm2) <- McmcChains(tbm2)
  tbm2
}

#' @export
TBM <- TrioBatchModel
