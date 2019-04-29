.empty_trio_model <- function(hp, mp){
  K <- k(hp)
  B <- 0L
  S <- iter(mp)
  ch <- initialize_mcmcT(K, S, B, 0)
  obj <- new("TrioBatchModel",
             k=as.integer(K),
             hyperparams=hp,
             theta=matrix(NA, 0, K),
             sigma2=matrix(NA, 0, K),
             mu=numeric(K),
             tau2=numeric(K),
             nu.0=numeric(1),
             sigma2.0=numeric(1),
             pi=matrix(NA, 0, K),
             pi_parents=matrix(NA, 0, K),
             #data=numeric(K),
             triodata=as_tibble(0),
             mprob=matrix(NA, 0, 0),
             #maplabel=numeric(K),
             predictive=numeric(K*B),
             zstar=integer(K*B),
             data.mean=matrix(NA, B, K),
             data.prec=matrix(NA, B, K),
             z=integer(0),
             zfreq=integer(K),
             zfreq_parents=integer(K),
             probz=matrix(0, S, K),
             # note this assumes parents are 2/3 of N
             probz_par=matrix(0, 2/3*S, K),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=ch,
             mcmc.params=mp,
             batch=integer(0),
             batchElements=integer(0),
             label_switch=FALSE,
             marginal_lik=as.numeric(NA),
             .internal.constraint=5e-4,
             .internal.counter=0L)
  chains(obj) <- McmcChainsTrios(obj)
  obj
}

.TBM <- function(triodata=as_tibble(),
                 hp=HyperparametersTrios(),
                 mp=McmcParams(iter=1000, thin=10,
                               burnin=1000, nStarts=1),
                 mprob=matrix(), 
                 maplabel=maplabel){
  ## If the data is not ordered by batch,
  ## its a little harder to sort component labels
  log_ratio <- triodata$log_ratio
  batches <- triodata$batches
  triodata <- select(triodata, -c(log_ratio, batches))
  ub <- unique(batches)
  nbatch <- setNames(as.integer(table(batches)), ub)
  B <- length(ub)
  N <- nrow(triodata)
  ## move to setValidity
  if(nrow(triodata) != length(batches)) {
    stop("batch vector must be the same length as data")
  }
  K <- k(hp)
  ##mprob <- mprob
  ##maplabel <- maplabel
  ## mu_k is the average across batches of the thetas for component k
  ## tau_k is the sd of the batch means for component k
  mu <- sort(rnorm(k(hp), mu.0(hp), sqrt(tau2.0(hp))))
  tau2 <- 1/rgamma(k(hp), 1/2*eta.0(hp), 1/2*eta.0(hp) * m2.0(hp))
  pp <- p <- matrix(rdirichlet(1, alpha(hp))[1, ], B, K, byrow=TRUE)
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
  u <- rchisq(nrow(triodata), hp@dfr)
  S <- iter(mp)
  T <- length(unique(triodata$id))
  ch <- initialize_mcmcT(K, S, B, T)
  index <- match(c("father", "mother"), colnames(mprob))
  mprob2 <- mprob[, -index]
  father <- mprob[, "father"]
  mother <- mprob[, "mother"]
  is_mendel <- rep(1L, T)
  ##zamp <- sample(seq_len(K), N, replace=TRUE)
  ##is_offspring <- model@triodata$family_member=="o"
  # write R update z module
  #zo <- update
  ##zamp[is_offspring] <- zo
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
             pi_parents=pp,
             data=log_ratio,
             batch=batches,
             triodata=triodata,
             mprob=mprob2,
             maplabel=maplabel,
             u=u,
             predictive=numeric(K*B),
             zstar=integer(K*B),
             data.mean=matrix(0, B, K),
             data.prec=matrix(0, B, K),
             z=sample(seq_len(K), N, replace=TRUE),
             zfreq=integer(K),
             zfreq_parents=integer(K),
             probz=matrix(0, N, K),
             # note this assumes parents are 2/3 of N
             probz_par=matrix(0, 2/3*N, K),
             logprior=numeric(1),
             loglik=numeric(1),
             mcmc.chains=ch,
             mcmc.params=mp,
             batchElements=nbatch,
             label_switch=FALSE,
             marginal_lik=as.numeric(NA),
             father=as.integer(father),
             mother=as.integer(mother),
             is_mendelian=is_mendel,
             .internal.constraint=5e-4,
             .internal.counter=0L)
  obj
}

triodata <- function(object, wide=FALSE, map=FALSE){
  dat <- object@triodata
  if(map){
    dat$z <- map_z(object)
  } else {
    dat$z <- z(object)
  }
  dat$is_mendelian <- isMendelian(object)
  dat$log_ratio <- y(object)
  dat$batch <- batch(object)
  if(wide){
    dat2 <- dat %>%
      select(c("id", "family_member", "z")) %>%
      spread("family_member", "z")
    dat3 <- dat %>%
      select(c("id", "family_member", "log_ratio")) %>%
      spread("family_member", "log_ratio") %>%
      set_colnames(c("id", paste0("logr_", c("m", "f", "o"))))
    ## assumes that family members in a trio were processed in same batch
    dat4 <- group_by(dat, id) %>%
      summarize(batch=unique(batch))
    dat5 <- left_join(dat2, dat3, by="id") %>%
      left_join(dat4, by="id")
    dat <- dat5
  }
  dat
}

setReplaceMethod("probzpar", "TrioBatchModel", function(object, value){
  object@probz_par <- value
  object
})

## TODO Dangerous to have accessor do something more than return the value of it
## slot.  Further
## probzpar(object) <- probzpar(object)
## will not behave as expected
#' @rdname probzpar-method
#' @aliases probzpar,TrioBatchModel-method
setMethod("probzpar", "TrioBatchModel", function(object) {
  ## because first iteration not saved
  object@probz_par/(iter(object)-1)
})

combine_batchTrios <- function(model.list, batches){
  . <- NULL
  ch.list <- map(model.list, chains)
  th <- map(ch.list, theta) %>% do.call(rbind, .)
  s2 <- map(ch.list, sigma2) %>% do.call(rbind, .)
  ll <- map(ch.list, log_lik) %>% unlist
  ppi <- map(ch.list, p) %>% do.call(rbind, .)
  pp.par <- map(ch.list, pp) %>% do.call(rbind, .)
  n0 <- map(ch.list, nu.0) %>% unlist
  s2.0 <- map(ch.list, sigma2.0) %>% unlist
  logp <- map(ch.list, logPrior) %>% unlist
  .mu <- map(ch.list, mu) %>% do.call(rbind, .)
  .tau2 <- map(ch.list, tau2) %>% do.call(rbind, .)
  zfreq <- map(ch.list, zFreq) %>% do.call(rbind, .)
  zfreqpar <- map(ch.list, zFreqPar) %>% do.call(rbind, .)
  pred <- map(ch.list, predictive) %>% do.call(rbind, .)
  zsta <- map(ch.list, zstar) %>% do.call(rbind, .)
  ismendel <- map(ch.list, isMendelian) %>% Reduce("+", .)
  mc <- new("McmcChainsTrios",
            theta=data.matrix(th),
            sigma2=data.matrix(s2),
            pi=data.matrix(ppi),
            pi_parents=data.matrix(pp.par),
            mu=data.matrix(.mu),
            tau2=data.matrix(.tau2),
            nu.0=n0,
            sigma2.0=s2.0,
            zfreq=zfreq,
            zfreq_parents=zfreqpar,
            logprior=logp,
            loglik=ll,
            iter=nrow(th),
            predictive=pred,
            zstar=zsta,
            k=k(model.list[[1]]),
            B=length(unique(batches)),
            is_mendelian=ismendel)
  hp <- hyperParams(model.list[[1]])
  mp <- mcmcParams(model.list[[1]])
  triodata <- model.list[[1]]@triodata
  mprob <- model.list[[1]]@mprob
  father <- model.list[[1]]@father
  mother <- model.list[[1]]@mother
  maplabel <- model.list[[1]]@maplabel
  iter(mp) <- nrow(th)
  B <- length(unique(batches))
  K <- k(model.list[[1]])
  pm.th <- matrix(colMeans(th), B, K)
  pm.s2 <- matrix(colMeans(s2), B, K)
  pm.p <- matrix(colMeans(ppi), B, K)
  pm.par <- matrix(colMeans(pp.par), B, K)
  pm.n0 <- median(n0)
  pm.mu <- colMeans(.mu)
  pm.tau2 <- colMeans(.tau2)
  pm.s20 <- mean(s2.0)
  pz <- map(model.list, probz) %>% Reduce("+", .)
  pz <- pz/length(model.list)
  ## the accessor divides by number of iterations, so rescale
  ## RS: this is super confusing
  pz <- pz * (iter(mp) - 1)
  zz <- max.col(pz)
  yy <- y(model.list[[1]])
  y_mns <- as.numeric(tapply(yy, zz, mean))
  y_prec <- as.numeric(1/tapply(yy, zz, var))
  zfreq <- as.integer(table(zz))
  pzpar <- map(model.list, probzpar) %>% Reduce("+", .)
  pzpar <- pzpar/length(model.list)
  pzpar <- pzpar * (iter(mp) - 1)
  zzpar <- max.col(pzpar)
  #fam.class <- model.list[[1]]@triodata$family_member
  #zzpar <- zzpar[fam.class!="o"]
  zfreq_parents <- as.integer(table(zzpar))
  any_label_swap <- any(map_lgl(model.list, label_switch))
  ## use mean marginal likelihood in combined model,
  ## or NA if marginal likelihood has not been estimated
  ml <- map_dbl(model.list, marginal_lik)
  if(all(is.na(ml))) {
    ml <- as.numeric(NA)
  } else ml <- mean(ml, na.rm=TRUE)
  nbatch <- as.integer(table(batch(model.list[[1]])))
  N <- sapply(model.list, iter) %>% sum
  prob.mendel <- ismendel/N
  is_mendel <- rbinom(1000, 1, prob.mendel)
  model <- new(class(model.list[[1]]),
               triodata=triodata,
               mprob=mprob,
               father=father,
               mother=mother,
               maplabel=maplabel,
               k=k(hp),
               hyperparams=hp,
               theta=pm.th,
               sigma2=pm.s2,
               mu=pm.mu,
               tau2=pm.tau2,
               nu.0=pm.n0,
               sigma2.0=pm.s20,
               pi=pm.p,
               pi_parents=pm.par,
               data=y(model.list[[1]]),
               u=u(model.list[[1]]),
               data.mean=y_mns,
               data.prec=y_prec,
               z=zz,
               zfreq=zfreq,
               zfreq_parents=zfreq_parents,
               probz=pz,
               probz_par=pzpar,
               predictive=predictive(mc)[nrow(th), ],
               zstar=zstar(mc)[nrow(th), ],
               logprior=numeric(1),
               loglik=numeric(1),
               mcmc.chains=mc,
               batch=batch(model.list[[1]]),
               batchElements=nbatch,
               modes=list(),
               mcmc.params=mp,
               label_switch=any_label_swap,
               marginal_lik=ml,
               is_mendelian=is_mendel,
               .internal.constraint=5e-4,
               .internal.counter=0L)
  modes(model) <- computeModes(model)
  log_lik(model) <- computeLoglik(model)
  logPrior(model) <- computePrior(model)
  model
}

#' Constructor for TrioBatchModel
#'
#' Initializes a TrioBatchModel, a container for storing data, parameters, and MCMC output for mixture models with batch- and component-specific means and variances.
#'
#' @param triodata the data for the simulation.
#' @param batches an integer-vector of the different batches
#' @param hp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mp An object of class 'McmcParams'
#' @return An object of class `TrioBatchModel`
..TrioBatchModel <- function(triodata=tibble(),
                           hp=HyperparametersTrios(),
                           mp=McmcParams(iter=1000, thin=10,
                                         burnin=1000, nStarts=4),
                           mprob=mprob,
                           maplabel=maplabel
                           ){
  if(length(triodata) == 0){
    return(.empty_trio_model(hp, mp))
  }
  iter <- 0
  validZ <- FALSE
  mp.tmp <- McmcParams(iter=0, burnin=burnin(mp), thin=1, nStarts=1)
  while(!validZ){
    ##
    ## Burnin with TBM model
    ##
    tbm <- .TBM(triodata, hp, mp.tmp, mprob, maplabel)
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
  chains(tbm2) <- McmcChainsTrios(tbm2)
  tbm2
}


#' Initializes a TrioBatchModel, a container for storing data, parameters, and MCMC output for mixture models with batch- and component-specific means and variances.
#'
#' @param triodata the data for the simulation.
#' @param batches an integer-vector of the different batches
#' @param hp An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @param mp An object of class 'McmcParams'
#' @return An object of class `TrioBatchModel`
#' @export
TrioBatchModel <- function(triodata=tibble(),
                           hp=HyperparametersTrios(),
                           mp=McmcParams(iter=1000, thin=10,
                                         burnin=1000, nStarts=4),
                           mprob=mprob,
                           maplabel=maplabel
                           ){
  if(length(triodata) == 0){
    return(.empty_trio_model(hp, mp))
  }
  tbm <- .TBM(triodata, hp, mp, mprob, maplabel)
  z(tbm) <- update_z(tbm)
  chains(tbm) <- McmcChainsTrios(tbm)
  tbm
}

#' @export
TBM <- TrioBatchModel




theta.fix <- function(model, values) {
  # values is in the form of a [1,K] matrix
  model@theta <- values
  model
}

theta.fix2 <- function(model) {
  modes(model) <- computeModes(model)
  model@theta <- as.matrix(model@modes$theta)
  model
}

mp.fix <- function(model, mp){
  mcmcParams(model) <- mp
  model
}

sigma2.fix <- function(model){
  modes(model) <- computeModes(model)
  model@sigma2 <- as.matrix(model@modes$sigma2)
  model
}

pi.fix <- function(model){
  modes(model) <- computeModes(model)
  model@pi <- model@modes$mixprob
  model
}



gibbs_trios_K2 <- function(hp,
                          mp,
                          k_range=c(1, 4),
                          dat,
                          batches,
                          maplabel,
                          mprob,
                          max_burnin=32000,
                          reduce_size=TRUE,
                          min_effsize=500, 
                          theta.feed){
  K <- seq(k_range[1], k_range[2])
  hp.list <- map(K, updateK, hp)
  model.list <- map(hp.list,
                    .gibbs_trios_mcmc2,
                    mp=mp,
                    dat=dat,
                    batches=batches,
                    maplabel=maplabel,
                    mprob=mprob,
                    max_burnin=max_burnin,
                    min_effsize=min_effsize,
                    theta.feed=theta.feed)
  names(model.list) <- paste0("TBM", map_dbl(model.list, k))
  ## sort by marginal likelihood
  ##
  ## if(reduce_size) TODO:  remove z chain, keep y in one object
  ##
  ix <- order(map_dbl(model.list, marginal_lik), decreasing=TRUE)
  models <- model.list[ix]
}

gibbs_trios2 <- function(model="TBM",
                        dat,
                        mp,
                        hp.list,
                        batches,
                        k_range=c(1, 4),
                        max_burnin=32e3,
                        top=2,
                        df=100,
                        min_effsize=500, 
                        theta.feed){
  #model <- unique(model)
  max_burnin <- max(max_burnin, burnin(mp)) + 1
  if(missing(hp.list)){
    #note this line will have to be incorporated in gibbs.R when generalised
    hp.list <- hpList(df=df)[["TBM"]]
  }
  message("Fitting TBM models")
  tbm <- gibbs_trios_K2(hp.list,
                       k_range=k_range,
                       mp=mp,
                       dat=dat,
                       batches=rep(1L, length(dat)),
                       maplabel=maplabel,
                       mprob=mprob,
                       max_burnin=max_burnin,
                       min_effsize=min_effsize,
                       theta.feed=theta.feed)
}

par.fix <- function(model, values) {
  # values is in the form of a model@z input
  model@z <- values
  model
}


gibbs_trios_K3 <- function(hp,
                           mp,
                           k_range=c(1, 4),
                           dat,
                           batches,
                           maplabel,
                           mprob,
                           max_burnin=32000,
                           reduce_size=TRUE,
                           min_effsize=500, 
                           theta.feed, par.feed){
  K <- seq(k_range[1], k_range[2])
  hp.list <- map(K, updateK, hp)
  model.list <- map(hp.list,
                    .gibbs_trios_mcmc3,
                    mp=mp,
                    dat=dat,
                    batches=batches,
                    maplabel=maplabel,
                    mprob=mprob,
                    max_burnin=max_burnin,
                    min_effsize=min_effsize,
                    theta.feed=theta.feed,
                    par.feed=par.feed)
  names(model.list) <- paste0("TBM", map_dbl(model.list, k))
  ## sort by marginal likelihood
  ##
  ## if(reduce_size) TODO:  remove z chain, keep y in one object
  ##
  ix <- order(map_dbl(model.list, marginal_lik), decreasing=TRUE)
  models <- model.list[ix]
}

gibbs_trios3 <- function(model="TBM",
                         dat,
                         mp,
                         hp.list,
                         batches,
                         k_range=c(1, 4),
                         max_burnin=32e3,
                         top=2,
                         df=100,
                         min_effsize=500, 
                         theta.feed,
                         par.feed){
  #model <- unique(model)
  max_burnin <- max(max_burnin, burnin(mp)) + 1
  if(missing(hp.list)){
    #note this line will have to be incorporated in gibbs.R when generalised
    hp.list <- hpList(df=df)[["TBM"]]
  }
  message("Fitting TBM models")
  tbm <- gibbs_trios_K3(hp.list,
                        k_range=k_range,
                        mp=mp,
                        dat=dat,
                        batches=rep(1L, length(dat)),
                        maplabel=maplabel,
                        mprob=mprob,
                        max_burnin=max_burnin,
                        min_effsize=min_effsize,
                        theta.feed=theta.feed,
                        par.feed=par.feed)
}


gMendelian.multi <- function(tau=c(0.5, 0.5, 0.5)){
  tau1 <- tau[1]
  tau2 <- tau[2]
  tau3 <- tau[3]
  mendelian.probs <- array(dim=c(5, 5, 5))
  genotypes <- c("CN0", "CN1", "CN2", "CN3", "CN4")
  dimnames(mendelian.probs) <- list(paste0("O_", genotypes),
                                    paste0("M_", genotypes),
                                    paste0("F_", genotypes))
  mendelian.probs[, 1, 1] <- c(1,0,0,0,0)
  mendelian.probs[, 1, 2] <- c(tau1, 1 - tau1, 0,0,0)
  mendelian.probs[, 1, 3] <- c(tau2 / 2, 0.5, (1 - tau2) / 2, 0,0)
  mendelian.probs[, 1, 4] <- c(0, tau3, 1 - tau3, 0,0)
  mendelian.probs[, 1, 5] <- c(0,0,1,0,0)
  mendelian.probs[, 2, 1] <- c(tau1, 1 - tau1, 0,0,0)
  mendelian.probs[, 2, 2] <- c((tau1^2), 2 * (tau1 * (1 - tau1)), ((1 - tau1)^2), 0,0)
  mendelian.probs[, 2, 3] <- c((tau1 * tau2) / 2, (tau2 * (1 - tau1) + tau1) / 2, (tau1 * (1-tau2) + (1 - tau1)) / 2, ((1 - tau1) * (1 - tau2)) / 2, 0)
  mendelian.probs[, 2, 4] <- c(0, tau1 * tau3, tau1 * (1 - tau3) + (1 - tau1) * tau3, (1- tau1) * (1 - tau3), 0)
  mendelian.probs[, 2, 5] <- c(0, 0, tau1, (1 - tau1), 0)
  mendelian.probs[, 3, 1] <- c(tau2 / 2, 0.5, (1 - tau2) / 2, 0, 0)
  mendelian.probs[, 3, 2] <- c((tau1 * tau2) / 2, (tau1 + tau2 * (1 - tau1)) / 2, ((1 - tau1) + (tau1 * (1-tau2)) ) / 2, (1 - tau1) * (1 - tau2) / 2, 0)
  mendelian.probs[, 3, 3] <- c(tau2^2 / 4, tau2 / 2, (0.5 + tau2 * (1 - tau2)) / 2, (1 - tau2) / 2, (1 - tau2)^2 / 4)
  mendelian.probs[, 3, 4] <- c(0, tau2 * tau3 / 2, (tau3 + tau2 * (1 - tau3)) / 2, (((1 - tau3) + (1 - tau2) * tau3) / 2), (1 - tau2) * (1 - tau1) /2)
  mendelian.probs[, 3, 5] <- c(0, 0, tau2 / 2, 0.5, (1 - tau2) / 2)
  mendelian.probs[, 4, 1] <- c(0, tau3, (1-tau3), 0, 0)
  mendelian.probs[, 4, 2] <- c(0, tau1 * tau3, tau1 * (1 - tau3) + (1 - tau1) * tau3, (1 - tau1) * (1 - tau3), 0)
  mendelian.probs[, 4, 3] <- c(0, tau2 * tau3 / 2, (tau3 + tau2 * (1 - tau3)) / 2, ((1 - tau3) + (1 - tau2) * tau3) / 2, (1 - tau2) * (1 - tau3) / 2)
  mendelian.probs[, 4, 4] <- c(0,0, tau3^2, 2 * tau3 * (1 - tau3), (1 - tau3)^2)
  mendelian.probs[, 4, 5] <- c(0,0,0, tau3, 1-tau3)
  mendelian.probs[, 5, 1] <- c(0,0,1,0,0)
  mendelian.probs[, 5, 2] <- c(0,0, tau1, 1 - tau1, 0)
  mendelian.probs[, 5, 3] <- c(0,0, tau2 / 2, 0.5, (1 - tau2) / 2)
  mendelian.probs[, 5, 4] <- c(0,0,0, tau3, 1 - tau3)
  mendelian.probs[, 5, 5] <- c(0,0,0,0,1)
  mendelian.probs
}

mprob.matrix <-  function(tau=c(0.5, 0.5, 0.5), maplabel, error){
  
  # to avoid confusion, maplabel1 will always be unique applied throughout all functions
  # and maplabel will not
  # so line below will repeat as needed
  maplabel1 <- unique(maplabel)

  # these conditionals means if maplabel has only one component = 2,
  # then it defaults to biallelic matrix
  if (all(any(maplabel1 < 2) & any(maplabel1 > 2))) {
    mprob.mat <- mprob.matrix.mallelic(tau, maplabel)
  } else {
    mprob.mat <- mprob.matrix.biallelic(tau, maplabel, error)
  }
 
  setDT(mprob.mat)[, c("father","mother") := tstrsplit(parents, "")]
  
  mprob.mat <- mprob.mat[, -1] %>%
    as_tibble() %>%
    mutate(father=as.numeric(father),
           mother=as.numeric(mother)) %>%
    as.matrix

# makes sure transmission matrix same size as maplabel, including repeats CN states
mprob.mat2 <- mprob.mat.resize(mprob.mat, maplabel)  

# makes sure each row adds to 1
mprob.mat3 <- mprob.balance(mprob.mat2, maplabel)

mprob.mat3
  
}    

mprob.matrix.mallelic <- function (tau=c(0.5, 0.5, 0.5), maplabel){
  tau1 <- tau[1]
  tau2 <- tau[2]
  tau3 <- tau[3]
  mendelian.probs <- array(dim=c(25,6))
  colnames(mendelian.probs) <- c("parents", "p(0|f,m)", "p(1|f,m)", "p(2|f,m)", "p(3|f,m)", "p(4|f,m)")
  
  mendelian.probs[1, 2:6] <- c(1,0,0,0,0)
  mendelian.probs[2, 2:6] <- c(tau1, 1 - tau1, 0,0,0)
  mendelian.probs[3, 2:6] <- c(tau2 / 2, 0.5, (1 - tau2) / 2, 0,0)
  mendelian.probs[4, 2:6] <- c(0, tau3, 1 - tau3, 0,0)
  mendelian.probs[5, 2:6] <- c(0,0,1,0,0)
  mendelian.probs[6, 2:6] <- c(tau1, 1 - tau1, 0,0,0)
  mendelian.probs[7, 2:6] <- c((tau1^2), 2 * (tau1 * (1 - tau1)), ((1 - tau1)^2), 0,0)
  mendelian.probs[8, 2:6] <- c((tau1 * tau2) / 2, (tau2 * (1 - tau1) + tau1) / 2, (tau1 * (1-tau2) + (1 - tau1)) / 2, ((1 - tau1) * (1 - tau2)) / 2, 0)
  mendelian.probs[9, 2:6] <- c(0, tau1 * tau3, tau1 * (1 - tau3) + (1 - tau1) * tau3, (1- tau1) * (1 - tau3), 0)
  mendelian.probs[10, 2:6] <- c(0, 0, tau1, (1 - tau1), 0)
  mendelian.probs[11, 2:6] <- c(tau2 / 2, 0.5, (1 - tau2) / 2, 0, 0)
  mendelian.probs[12, 2:6] <- c((tau1 * tau2) / 2, (tau1 + tau2 * (1 - tau1)) / 2, ((1 - tau1) + (tau1 * (1-tau2)) ) / 2, (1 - tau1) * (1 - tau2) / 2, 0)
  mendelian.probs[13, 2:6] <- c(tau2^2 / 4, tau2 / 2, (0.5 + tau2 * (1 - tau2)) / 2, (1 - tau2) / 2, (1 - tau2)^2 / 4)
  mendelian.probs[14, 2:6] <- c(0, tau2 * tau3 / 2, (tau3 + tau2 * (1 - tau3)) / 2, (((1 - tau3) + (1 - tau2) * tau3) / 2), (1 - tau2) * (1 - tau1) /2)
  mendelian.probs[15, 2:6] <- c(0, 0, tau2 / 2, 0.5, (1 - tau2) / 2)
  mendelian.probs[16, 2:6] <- c(0, tau3, (1-tau3), 0, 0)
  mendelian.probs[17, 2:6] <- c(0, tau1 * tau3, tau1 * (1 - tau3) + (1 - tau1) * tau3, (1 - tau1) * (1 - tau3), 0)
  mendelian.probs[18, 2:6] <- c(0, tau2 * tau3 / 2, (tau3 + tau2 * (1 - tau3)) / 2, ((1 - tau3) + (1 - tau2) * tau3) / 2, (1 - tau2) * (1 - tau3) / 2)
  mendelian.probs[19, 2:6] <- c(0,0, tau3^2, 2 * tau3 * (1 - tau3), (1 - tau3)^2)
  mendelian.probs[20, 2:6] <- c(0,0,0, tau3, 1-tau3)
  mendelian.probs[21, 2:6] <- c(0,0,1,0,0)
  mendelian.probs[22, 2:6] <- c(0,0, tau1, 1 - tau1, 0)
  mendelian.probs[23, 2:6] <- c(0,0, tau2 / 2, 0.5, (1 - tau2) / 2)
  mendelian.probs[24, 2:6] <- c(0,0,0, tau3, 1 - tau3)
  mendelian.probs[25, 2:6] <- c(0,0,0,0,1)
  
  if(all((rowSums(mendelian.probs, na.rm=T))==1)==F) stop("mendelian matrix is incorrect")
 
  # must do as_tibble step as cbind converts all cells into one type i.e. characters otherwise
  mprob.mat <- as_tibble(mendelian.probs)
  ref.geno <- reference.genotype(maplabel)
  mprob.mat[, 1] <- ref.geno
  
  mprob.mat <- mprob.subset2(mprob.mat, maplabel)
  mprob.mat
}

# error term is hard coded into bi-allelic matrix for now
mprob.matrix.biallelic <- function (tau=c(0.5, 0.5, 0.5), maplabel, error = 0.001){
  # note in biallelic, only 1 tau value required
  tau1 <- tau[1]
  allele1 <- c(1 - error, 1-tau1-error, 0 + error)
  allele2 <- 1 - allele1
  
  # kronecker's product
  a1a1 <- allele1 %x% allele1
  a1a2 <- allele1 %x% allele2 + allele2 %x% allele1
  a2a2 <- allele2 %x% allele2
  ref.geno <- reference.genotype(maplabel)
  
  # put the three vectors together for the biallelic transmission matrix
  #for single autosomal locus
  mendelian.probs <- array(dim=c(9,4))
  offspring.geno <- cbind(a1a1, a1a2, a2a2)
  mendelian.probs[,2:4] <- offspring.geno
  mprob.mat <- as_tibble(mendelian.probs)
  mprob.mat[,1] <- ref.geno

  # correctly label the columns 
  maplabel1 <- unique(maplabel)
  
  if (all(maplabel1 < 3)){
    maplabel2 <- c(0,1,2)
  } else {
    maplabel2 <- c(2,3,4)
  }
  
  colnames.label <- vector(length = 3)
  for (i in 1:3) {
    cn <- maplabel2[i]
    label <- paste0("p(",cn,"|f,m)")
    colnames.label[i] <- label
  }
  
  colnames(mprob.mat)[1] <- "parents"
  colnames(mprob.mat)[2:4] <- colnames.label
  
  mprob.mat <- mprob.subset2(mprob.mat, maplabel)
  mprob.mat
}

reference.genotype <- function(maplabel){
 
   maplabel1 <- unique(maplabel)
  
  if (all(any(maplabel1 < 2) & any(maplabel1 > 2))) {
    ref.geno <- c("00", "01", "02", "03", "04", 
                  "10", "11", "12", "13", "14",
                  "20", "21", "22", "23", "24",
                  "30", "31", "32", "33", "34",
                  "40", "41", "42", "43", "44")
  } else {
    if (all(maplabel1 < 3)){ # this is deletion matrix reference genotypes
      ref.geno <- c("00", "01", "02", 
                    "10", "11", "12",
                    "20", "21", "22")
    } else { # this would be duplication matrix reference genotypes
      ref.geno <- c("22", "23", "24", 
                    "32", "33", "34",
                    "42", "43", "44")
    }
  }
}

mprob.subset2 <- function(mprob.mat, maplabel) {
  maplabel1 <- unique(maplabel)
  K <- length(maplabel1)
  M <- ifelse(all(maplabel1 > 1), 0, 2)
  col.a <- maplabel1[1] + M
  col.b <- maplabel1[K] + M
  
  # reference genotype is an index to enable the relevant rows of the 
  #transmission matrix to be subset 
  ref.geno <- reference.genotype(maplabel)
  index <- mprob.label2(maplabel)
  rows <- match(index, ref.geno)
  
  mprob.subset <- mprob.mat[rows, c(col.a:col.b)]
  mprob.rows <- mprob.mat[rows, 1]
  mprob.subset <- cbind(mprob.rows, mprob.subset)
  mprob.subset
}

mprob.label2 <- function(maplabel){
  maplabel1 <- unique(maplabel)
  n <- length(maplabel1)
  combo <- permutations(n=n, r=2, v=maplabel1, repeats.allowed=T)
  geno.combo <- paste0(combo[,1], combo[,2])
  geno.combo
}

mprob.mat.resize <- function(mprob.matrix, maplabel){
  maplabel1 <- unique(maplabel)
  K <- ncol(mprob.matrix)
  mprob.mat2 <- mprob.matrix[,1:(K-2)]
  mprob.mat2t <- t(mprob.mat2)
  mprob.mat3 <- cbind(mprob.mat2t, maplabel1)
  mprob.mat3 <- data.frame(mprob.mat3)
  
  maplabel.sort <- data.frame(sort(maplabel))
  colnames(maplabel.sort) <- "maplabels"
  
  mprob.sized <- left_join(maplabel.sort, mprob.mat3, by=c("maplabels"="maplabel1"))
  sized.matrix <- as.matrix(t(mprob.sized))
  sized.matrix <- sized.matrix[-1,]
  rownames(sized.matrix) <- NULL

  L <- length(maplabel)
  colnames.label <- vector(length = L)
  for (i in 1:L) {
    cn <- maplabel[i]
    label <- paste0("p(",cn,"|f,m)")
    colnames.label[i] <- label
  }
  
  colnames(sized.matrix) <- colnames.label
  mprob.parents <- mprob.matrix[,(K-1):K]
  mprob.mat.sized <- cbind(sized.matrix, mprob.parents)
  
  mprob.mat.sized
}

mprob.balance <- function (mprob.matrix, maplabel){
    cols <- length(maplabel)
    mprob2 <- mprob.matrix[,1:cols]
    mprob2 <- (mprob2)/(rowSums(mprob2))
    mprob.matrix[, 1:cols] <- mprob2
    mprob.matrix
}


component_stats <- function(tbl){
  tbl2 <- tbl %>% group_by(batches, copy_number) %>%
    summarize(mean=mean(log_ratio),
              sd=sd(log_ratio),
              n=n())
  k <- length(unique(tbl2$batches))
  denom <- sum(tbl2$n)/k
  tbl3 <- tbl2 %>% group_by(batches, copy_number) %>% mutate(p = n / denom)
  tbl3
  }

# sensitivity and specificity function

# deprecate the following simulation functions:
.simulate_data_multi2 <- function(params, n, mprob, maplabel){
  K <- nrow(params)
  p <- params$p
  theta <- params$theta
  sigma2 <- params$sigma2
  sigma <- sqrt(sigma2)
  c_mk <- sample(1:K, size = n, replace = TRUE, prob = p)
  c_fk <- sample(1:K, size = n, replace = TRUE, prob = p)
  c_m <- maplabel[c_mk]
  c_f <- maplabel[c_fk]
  c_ok <- rep(NA, length = n)
  C <- ncol(mprob)
  for(i in 1:n){
    cn_ms <- c_m[i]
    cn_fs <- c_f[i]
    mprob2 <- mprob[mprob[,"mother"]==cn_ms,]
    mprob2 <- mprob2[mprob2[,"father"]==cn_fs,]
    mprob3 <- mprob2[1:(C-2)]
    c_ok[i] <- sample(1:K, size = 1, prob = mprob3)
  }
  c_o <- maplabel[c_ok]
  id.index <- formatC(seq_len(n), flag="0", width=3)
  logr.tbl <- tibble(m=rnorm(n, mean = theta[c_mk], sd = sigma[c_mk]),
                     f=rnorm(n, mean = theta[c_fk], sd = sigma[c_fk]),
                     o=rnorm(n, mean = theta[c_ok], sd = sigma[c_ok]),
                     id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="log_ratio", -id) 
  cn.mat <- cbind(c_m, c_f, c_o)
  colnames(cn.mat) <- c("m", "f", "o")
  cn.tbl <- as_tibble(cn.mat) %>%
    mutate(id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="copy_number", -id)
  tbl <- left_join(logr.tbl, cn.tbl, by=c("id", "family_member")) %>%
    mutate(family_member=factor(family_member, levels=c("m", "f", "o"))) %>%
    arrange(id, family_member)
  tbl
}

# this function ensures the simulation generates the requisite number of components regardless of MAF parameter (p)
simulate_data_multi2 <- function(params, N, batches,
                                 error=0, mendelian.probs,
                                 maplabel){
  tbl <- .simulate_data_multi2(params, N, mendelian.probs, maplabel)
  # append batch info (assumes ordering by trios and id)
  if(missing(batches)) {
    batches <- rep(1L, 3*N)
  } else {
    batches <- as.integer(factor(batches))
  }
  batches <- sort(batches)
  tbl2 <- cbind(tbl, batches)
  tbl3 <- as_tibble(tbl2)
  ##
  ## update parameters to be same as empirical values
  ##
  stats <- component_stats(tbl3)
  while(nrow(stats)!=length(maplabel)){
    tbl <- .simulate_data_multi2(params, N, mendelian.probs, maplabel)
    batches <- sort(batches)
    tbl2 <- cbind(tbl, batches)
    tbl3 <- as_tibble(tbl2)
    stats <- component_stats(tbl3)
  }
  p <- stats %>% group_by(copy_number) %>% summarize(p=mean(p))
  sd <- stats %>% group_by(copy_number) %>% summarize(sd=mean(sd))
  mean <- stats %>% group_by(copy_number) %>% summarize(mean=mean(mean))
  params$p <- p$p
  params$sigma2 <- (sd$sd)^2
  params$theta <- mean$mean
  komponent <- frank(tbl$copy_number, ties.method=c("dense"))
  loglik <- sum(dnorm(tbl$log_ratio, params$theta[komponent],
                      sqrt(params$sigma2[komponent]), log=TRUE))
  truth <- list(data=tbl3, params=params, loglik=loglik)
  truth
}

simulate_data_multi3 <- function(params, N, batches, error=0, mendelian.probs, maplabel){
  tbl <- .simulate_data_multi2(params, N, mendelian.probs, maplabel)
  
  # append batch info (assumes ordering by trios and id)
  if(missing(batches)) {
    batches <- rep(1L, 3*N)
  } else {
    batches <- as.integer(factor(batches))
  }
  batches <- sort(batches)
  tbl2 <- cbind(tbl, batches)
  tbl3 <- as_tibble(tbl2)
  ##
  ## update parameters to be same as empirical values
  ##
  stats <- component_stats(tbl3)
  
  p <- stats %>% group_by(copy_number) %>% summarize(p=mean(p))
  sd <- stats %>% group_by(copy_number) %>% summarize(sd=mean(sd))
  mean <- stats %>% group_by(copy_number) %>% summarize(mean=mean(mean))
  params$p <- p$p
  params$sigma2 <- (sd$sd)^2
  params$theta <- mean$mean
  komponent <- frank(tbl$copy_number, ties.method=c("dense"))
  loglik <- sum(dnorm(tbl$log_ratio, params$theta[komponent],
                      sqrt(params$sigma2[komponent]), log=TRUE))
  truth <- list(data=tbl3, params=params, loglik=loglik)
  truth
}

#
startAtTrueValues2 <- function(model, truth_stats, data){
  rows <- length(unique(truth_stats$batches))
  cols <- length(unique(truth_stats$copy_number))
  theta(model) <- matrix(truth_stats$mean, nrow=rows, ncol=cols, byrow=T)
  sigma2(model) <- matrix((truth_stats$sd)^2, nrow=rows, ncol=cols, byrow=T)
  p(model) <- data$params$p
  z(model) <- as.integer(data$data$copy_number)
  zFreq1 <- truth_stats %>% group_by(copy_number) %>% summarize(zfreq=sum(n))
  zFreq(model) <- zFreq1$zfreq
  model
}

setMethod("triodata_lrr", "TrioBatchModel", function(object){
  object@triodata$log_ratio
}
)

setMethod("updateZ", "TrioBatchModel", function(object){
  update_z(object)
})


setMethod("computeMeans", "TrioBatchModel", function(object){
  compute_means(object)
})


setMethod("computePrec", "TrioBatchModel", function(object){
  compute_prec(object)
})

setMethod("computeModes", "TrioBatchModel", function(object){
  .computeModesBatch(object)
})

componentVariances <- function(y, z)  v <- sapply(split(y, z), var)


setMethod("computeVars", "TrioBatchModel", function(object){
  compute_vars(object)
})


setMethod("show", "TrioBatchModel", function(object){
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

#' @rdname k-method
#' @aliases k,TrioBatchModel-method
setMethod("k", "TrioBatchModel", function(object) object@k)

#' @rdname k-method
#' @aliases k<-,TrioBatchModel-method
setReplaceMethod("k", "TrioBatchModel",
                 function(object, value) {
                   k <- as.integer(value)
                   hypp <- hyperParams(object)
                   hypp@k <- k
                   hypp@alpha <- rep(1, k)
                   hyperParams(object) <- hypp
                   object@k <- k
                   object@pi <- matrix(rep(1/k, k), numBatch(object), k, byrow=TRUE)
                   object@probz <- matrix(0, length(y(object)), k)
                   object <- startingValues(object)
                   object
                 }
)

#' Retrieve mixture proportions of parents.
#'
#' @examples
#'   \dontrun{
#'      pp(TrioBatchModelExample)
#'    } 
#' @param object an object of class TrioBatchModel
#' @return A vector of length the number of components
#' @export
pp <- function(object) object@pi_parents

setReplaceMethod("pp", "TrioBatchModel", function(object, value){
  object@pi_parents <- value
  object
})

# adapted from old SNR script
# MC: will modify variables to be more self-explanatory later
snr.calc.median<-function(count.num, sd.num, deltas, lrr.list){
  loop.length<-ifelse(length(lrr.list)!=1, length(lrr.list)-1, 1)
  variance<-rep(NA,loop.length)
  for (i in 1:loop.length){
    var.calc<-((count.num[[i]]-1)*(sd.num[[i]]^2)+(count.num[[i+1]]-1)*(sd.num[[i+1]]^2))/(count.num[i]+count.num[i+1]-2)
    variance[i]<-var.calc
  }
  snr<-deltas/(variance^0.5)
  median.snr<-median(snr)
  return(median.snr)
}

# wrapper function 
snr.calc <- function(model){
  lrr<-as.numeric(y(model))
  calls<-as.numeric(z(model))
  avglrr.list <- split(lrr, calls)
  comp.means <- sapply(avglrr.list, mean)
  comp.means <- sort(comp.means)
  deltas <- diff(comp.means)
  sds <- sapply(avglrr.list, sd)
  count<-sapply(avglrr.list,length)
  SNR<-snr.calc.median(count, sds, deltas, avglrr.list)
  SNR
}

snr.calc.mean<-function(count.num, sd.num, deltas, lrr.list){
  loop.length<-ifelse(length(lrr.list)!=1, length(lrr.list)-1, 1)
  variance<-rep(NA,loop.length)
  for (i in 1:loop.length){
    var.calc<-((count.num[[i]]-1)*(sd.num[[i]]^2)+(count.num[[i+1]]-1)*(sd.num[[i+1]]^2))/(count.num[i]+count.num[i+1]-2)
    variance[i]<-var.calc
  }
  snr<-deltas/(variance^0.5)
  mean.snr<-mean(snr)
  return(mean.snr)
}

snr.calc2 <- function(model){
  lrr<-as.numeric(y(model))
  calls<-as.numeric(z(model))
  avglrr.list <- split(lrr, calls)
  comp.means <- sapply(avglrr.list, mean)
  comp.means <- sort(comp.means)
  deltas <- diff(comp.means)
  sds <- sapply(avglrr.list, sd)
  count<-sapply(avglrr.list,length)
  SNR<-snr.calc.mean(count, sds, deltas, avglrr.list)
  SNR
}


# cross table in the same format as caret package
# specifically for 3 state cross table
modified.sens.spec.calc <- function(cross.table, cn.type=c("DEL", "DUP")){
  # sensitivity = TP / TP+ FN
  # specificity = TN / TN + FP
  if(cn.type == "DEL"){
    tp <- cross.table[2,2] + cross.table[1,1]
    fn <- cross.table[3,1] + cross.table[3,2]
    tn <- cross.table[3,3]
    fp <- cross.table[1,2] + cross.table[1,3] + cross.table[2,1] + cross.table[2,3]
    sens <- tp / (tp + fn)
    spec <- tn / (tn + fp)
  } else{
    tp <- cross.table[2,2] + cross.table[3,3]
    fn <- cross.table[1,2] + cross.table[1,3]
    tn <- cross.table[1,1]
    fp <- cross.table[2,1] + cross.table[2,3] + cross.table[3,1] + cross.table[3,2]
    sens <- tp / (tp + fn)
    spec <- tn / (tn + fp)
  }
  sens.spec <- list("sensitivity" = sens, "specificity" = spec)
  return(sens.spec)
}

ci95.wrap <- function(vector) {
  test <- t.test(vector)
  ci95 <- test$conf.int
  ci95
}

setMethod("isMendelian", "TrioBatchModel", function(object) {
  object@is_mendelian
})

nTrios <- function(object) {
  dat <- triodata(object)
  length(unique(dat$id))
}

simulateTrioData <- function(theta=c(-1.2, 0.3, 1.7), maplabel=c(0,1,2),
                             error=0){
  set.seed(123)
  p <- c(0.24, 0.43, 0.33)
  sigma2 <- c(0.05, 0.05, 0.05)
  params <- data.frame(cbind(p, theta, sigma2))
  mp <- McmcParams(iter=50, burnin=5)

  nbatch <- 1
  N <- 300
  maplabel <- c(0,1,2)
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=error)
  dat2 <- simulate_data_multi2(params,
                               N=N,
                               batches = rep(c(1:nbatch),
                                             length.out = 3*N),
                               error=0, mprob, maplabel)
  hp <- HyperparametersTrios(k = 3)
  model <- TBM(triodata=dat2$data,
                hp=hp,
                mp=mp,
                mprob=mprob,
               maplabel=maplabel)
  dat <- dat2$data
  list(model=model, true_z=as.integer(dat$copy_number) + 1L)
}


#' @rdname mcmcParams-method
#' @aliases mcmcParams,MixtureModel-method
setReplaceMethod("mcmcParams", "TrioBatchModel", function(object, value){
  S <- iter(value)
  B <- numBatch(object)
  K <- k(object)
  T <- nTrios(object)
  if(iter(object) != S){
    if(S > iter(object)){
      object@mcmc.params <- value
      ## create a new chain
      mcmc_chains <- initialize_mcmcT(K, S, B, T)
    } else {
      object@mcmc.params <- value
      index <- seq_len(S)
      mcmc_chains <- chains(object)[index, ]
      mcmc_chains@iter <- S
      mcmc_chains@B <- B
      mcmc_chains@k <- K
    }
    chains(object) <- mcmc_chains
    return(object)
  }
  ## if we've got to this point, it must be safe to update mcmc.params
  ## (i.e., size of chains is not effected)
  object@mcmc.params <- value
  object
})

probMendelian <- function(model){
  isMendelian(chains(model))/(iter(model))
}

setMethod("sortComponentLabels", "TrioBatchModel", function(model){
  is_ordered <- .ordered_thetas_multibatch(model)
  if(is_ordered) return(model)
  ## thetas are not all ordered
  thetas <- theta(model)
  s2s <- sigma2(model)
  K <- k(model)
  ix <- order(thetas[1, ])
  B <- nBatch(model)
  . <- NULL
  tab <- tibble(z_orig=z(model),
                z=z(model),
                batch=batch(model)) %>%
    mutate(index=seq_len(nrow(.)))
  z_relabel <- NULL
  for(i in seq_len(B)){
    ix.next <- order(thetas[i, ])
    thetas[i, ] <- thetas[i, ix.next]
    s2s[i, ] <- s2s[i, ix]
    index <- which(tab$batch == i)
    tab2 <- tab[index, ] %>%
      mutate(z_relabel=factor(z, levels=ix.next)) %>%
      mutate(z_relabel=as.integer(z_relabel))
    tab$z[index] <- tab2$z_relabel
  }
  ps <- p(model)[ix]
  mu(model) <- mu(model)[ix]
  tau2(model) <- tau2(model)[ix]
  model@pi_parents <- model@pi_parents[, ix, drop=FALSE]
  model@zfreq_parents <- model@zfreq_parents[ix]
  model@probz_par <- model@probz_par[, ix, drop=FALSE]
  sigma2(model) <- s2s
  theta(model) <- thetas
  p(model) <- ps
  z(model) <- tab$z
  log_lik(model) <- computeLoglik(model)
  model
})
