.empty_trio_model <- function(hp, mp){
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
             #pi_child=numeric(K),
             #data=numeric(K),
             triodata=as_tibble(0),
             mprob=matrix(NA, 0, 0),
             #maplabel=numeric(K),
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

.TBM <- function(triodata=as_tibble(),
                 hp=HyperparametersTrios(),
                 mp=McmcParams(iter=1000, thin=10,
                               burnin=1000, nStarts=4),
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
  p <- rdirichlet(1, alpha(hp))[1, ]
  p.child <- rdirichlet(1, alpha(hp))[1, ]
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
  index <- match(c("father", "mother"), colnames(mprob))
  mprob2 <- mprob[, -index]
  father <- mprob[, "father"]
  mother <- mprob[, "mother"]
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
             #pi_child=p.child,
             data=log_ratio,
             batch=batches,
             triodata=triodata,
             mprob=mprob2,
             maplabel=maplabel,
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
             batchElements=nbatch,
             label_switch=FALSE,
             marginal_lik=as.numeric(NA),
             father=as.integer(father),
             mother=as.integer(mother),
             .internal.constraint=5e-4,
             .internal.counter=0L)
  obj
}

triodata <- function(object){
  dat <- object@triodata
  dat$log_ratio <- y(object)
  dat$batches <- batch(object)
  dat
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
#' @export
TrioBatchModel <- function(triodata=tibble(),
                           hp=HyperparametersTrios(),
                           mp=McmcParams(iter=1000, thin=10,
                                         burnin=1000, nStarts=4),
                           mprob=mprob,
                           maplabel=maplabel
                           ){
  browser()
  if(nrow(triodata) == 0){
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
  chains(tbm2) <- McmcChains(tbm2)
  tbm2
}

#' @export
TBM <- TrioBatchModel


.gibbs_trios_mcmc <- function(hp, mp, dat, max_burnin=32000, batches, min_effsize=500){
  nchains <- nStarts(mp)
  nStarts(mp) <- 1L ## because posteriorsimulation uses nStarts in a different way
  if(iter(mp) < 500){
    warning("Require at least 500 Monte Carlo simulations")
    MIN_EFF <- ceiling(iter(mp) * 0.5)
  } else MIN_EFF <- min_effsize
  MIN_CHAINS <- 3
  MIN_GR <- 1.2
  neff <- 0; r <- 2
  while(burnin(mp) < max_burnin && thin(mp) < 100){
    message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
    mod.list <- replicate(nchains, TBM(dat=dat,
                                       hp=hp,
                                       mp=mp,
                                       batches=batches))

    ## MC: Define a method for runBurnin / runMcmc  :  define corresponding CPP code
    mod.list <- suppressWarnings(map(mod.list, posteriorSimulation))
    no_label_swap <- !map_lgl(mod.list, label_switch)
    if(sum(no_label_swap) < MIN_CHAINS){
      burnin(mp) <- as.integer(burnin(mp) * 2)
      mp@thin <- as.integer(thin(mp) + 2)
      nStarts(mp) <- nStarts(mp) + 1
      next()
    }
    mod.list <- mod.list[ no_label_swap ]
    ## MC:  check what selectModels does?
    mod.list <- mod.list[ selectModels(mod.list) ]
    ## MC:  need mcmcList method defined for this class, or new function
    mlist <- mcmcList(mod.list)
    neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
    if(is.null(neff)) neff <- 0
    r <- gelman_rubin(mlist, hp)
    message("     Gelman-Rubin: ", round(r$mpsrf, 2))
    message("     eff size (median): ", round(min(neff), 1))
    message("     eff size (mean): ", round(mean(neff), 1))
    if((mean(neff) > min_effsize) && r$mpsrf < MIN_GR) break()
    burnin(mp) <- as.integer(burnin(mp) * 2)
    mp@thin <- as.integer(thin(mp) + 2)
    nStarts(mp) <- nStarts(mp) + 1
    ##mp@thin <- as.integer(thin(mp) * 2)
  }
  ## MC: combine_batch / new function needed for TrioBatchModel
  model <- combine_batch(mod.list, batches)
  meets_conditions <- (mean(neff) > min_effsize) &&
    r$mpsrf < MIN_GR &&
    !label_switch(model)
  if(meets_conditions){
    ## Commented by Rob for now
    ## model <- compute_marginal_lik(model)
  }
  model
}

gibbs_trios_K <- function(hp,
                          mp,
                          k_range=c(1, 4),
                          dat,
                          batches,
                          max_burnin=32000,
                          reduce_size=TRUE,
                          min_effsize=500){
  K <- seq(k_range[1], k_range[2])
  hp.list <- map(K, updateK, hp)
  model.list <- map(hp.list,
                    .gibbs_trio_mcmc,
                    mp=mp,
                    dat=dat,
                    batches=batches,
                    max_burnin=max_burnin,
                    min_effsize=min_effsize)
  names(model.list) <- paste0("MB", map_dbl(model.list, k))
  ## sort by marginal likelihood
  ##
  ## if(reduce_size) TODO:  remove z chain, keep y in one object
  ##
  ix <- order(map_dbl(model.list, marginal_lik), decreasing=TRUE)
  models <- model.list[ix]
}

gibbs_trios <- function(model="trio",
                        dat,
                        mp,
                        hp.list,
                        batches,
                        k_range=c(1, 4),
                        max_burnin=32e3,
                        top=2,
                        df=100,
                        min_effsize=500){
  model <- unique(model)
  max_burnin <- max(max_burnin, burnin(mp)) + 1
  if(missing(hp.list)){
    ## MC:  implement
    hp.list <- hpTrioList(df=df)
  }
  message("Fitting SB models")
  sb <- gibbs_trios_K(hp.list,
                      k_range=k_range,
                      mp=mp,
                      dat=dat,
                      batches=rep(1L, length(dat)),
                      max_burnin=max_burnin,
                      min_effsize=min_effsize)
}


.simulate_data_multi <- function(params, n, mendelian.probs, GP){
  gp <- GP
  states <- gp$states + 1
  K <- nrow(params)
  p <- params$p
  theta <- params$theta
  sigma <- params$sigma
  c_mk <- sample(1:K, size = n, replace = TRUE, prob = p)
  c_fk <- sample(1:K, size = n, replace = TRUE, prob = p)
  c_o <- rep(NA, length = n)
  M <- cn_adjust2(gp)
  c_m <- c_mk + M
  c_f <- c_fk + M
  for(i in 1:n){
    cn_m <- c_m[i] + 1
    cn_f <- c_f[i] + 1
    p.offspring <- mendelian.probs[, cn_m, cn_f]
    p.offspring <- p.offspring[states]
    c_o[i] <- sample(1:K, size = 1, prob = p.offspring)
  }
  c_o <- c_o + M
  c_ok <- c_o - M
  id.index <- formatC(seq_len(n), flag="0", width=3)
  logr.tbl <- tibble(m=rnorm(n, mean = theta[c_mk], sd = sigma[c_mk]),
                     f=rnorm(n, mean = theta[c_fk], sd = sigma[c_fk]),
                     o=rnorm(n, mean = theta[c_ok], sd = sigma[c_ok]),
                     id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="log_ratio", -id) 
  cn.mat <- cbind(c_m, c_f, c_o)
  colnames(cn.mat) <- c("m", "f", "o")
  cn.tbl <- as.tibble(cn.mat) %>%
    mutate(id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="copy_number", -id)
  tbl <- left_join(logr.tbl, cn.tbl, by=c("id", "family_member")) %>%
    mutate(family_member=factor(family_member, levels=c("m", "f", "o"))) %>%
    arrange(id, family_member)
  tbl
}

geneticParams <- function(K=5,
                          states=0:4,
                          tau=c(1, 0.5, 0),
                          xi=c(1.5, 1, 1, 1, 1), 
                          mu=c(-3, -0.5, 0, 1, 2), ## theoretical
                          nu=1,
                          sigma2.0=0.001,
                          a=rep(1, K),
                          eta=c(0.5, 0.5),
                          error=5e-4,
                          ncp=30,
                          model="Genetic"){
  ##m.probs <- gMendelian(tau)
  list(K=K,
       states=states,
       tau=tau,
       mu=mu,
       xi=xi,
       nu=nu,
       sigma2.0=sigma2.0,
       a=a,
       eta=eta,
       error=error,
       ncp=ncp,
       ##m.probs=m.probs,
       model=model)
}

cn_adjust2 <- function(gp){
  states <- gp$states
  start.state <-states[1]
  M <- (-1)
  M <- ifelse(start.state==1, M <- 0, M)
  M <- ifelse(start.state==2, M <- 1, M)
  M <- ifelse(start.state==3, M <- 2, M)
  M <- as.numeric(M)
  M
}

simulate_data_multi <- function(params, N, batches, error=0, GP){
  gp <- GP
  mendelian.probs <- gMendelian.multi()
  tbl <- .simulate_data_multi(params, N, mendelian.probs, GP)
  ## standardize
  #tbl <- tbl %>%
   # mutate(log_ratio=(log_ratio-median(log_ratio))/sd(log_ratio))
  M <- cn_adjust2(gp)
  z <- tbl$copy_number - M
 
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
  stats <- component_stats(tbl)
  params$p <- stats$n/sum(stats$n)
  params$sigma <- stats$sd
  params$theta <- stats$mean
  loglik <- sum(dnorm(tbl$log_ratio, params$theta[z],
                      params$sigma[z], log=TRUE))
  truth <- list(data=tbl3, params=params, loglik=loglik)
  truth
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

mprob.matrix <-  function(tau=c(0.5, 0.5, 0.5), maplabel){
  
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
  
  mendelian.probs[, 1] <- c(00, 01, 02, 03, 04, 
                            10, 11, 12, 13, 14,
                            20, 21, 22, 23, 24,
                            30, 31, 32, 33, 34,
                            40, 41, 42, 43, 44)
  
  mprob.mat <- as.tibble(mendelian.probs)
  
  
  mprob.mat[, 1] <- c("00", "01", "02", "03", "04", 
                      "10", "11", "12", "13", "14",
                      "20", "21", "22", "23", "24",
                      "30", "31", "32", "33", "34",
                      "40", "41", "42", "43", "44")
  
  mprob.mat <- mprob.subset2(mprob.mat, maplabel)
  setDT(mprob.mat)[, c("father","mother") := tstrsplit(parents, "")]
  
  mprob.mat <- mprob.mat[, -1] %>%
    as.tibble() %>%
    mutate(father=as.numeric(father),
           mother=as.numeric(mother)) %>%
    as.matrix
  
  #K <- gp$K
  #ST <- gp$states[1]
  mprob.mat
  
  #extdata <- system.file("extdata", package="marimba2")
  #filename <- paste0("mendelian_probs2_",K,"_",ST,".rds")
  #saveRDS(mprob.mat, file.path(extdata, filename))
}

# assumes contiguous CN states but not 1:1 K:states correspondence
mprob.subset2 <- function(mprob.mat, maplabel) {
  states <- unique(maplabel)
  K <- length(states)
  col.a <- states[1] + 2
  col.b <- states[K] + 2
  
  ref.geno <- c("00", "01", "02", "03", "04", 
                "10", "11", "12", "13", "14",
                "20", "21", "22", "23", "24",
                "30", "31", "32", "33", "34",
                "40", "41", "42", "43", "44")
  index <- mprob.label2(states)
  rows <- match(index, ref.geno)
  
  mprob.subset <- mprob.mat[rows, c(col.a:col.b)]
  mprob.rows <- mprob.mat[rows, 1]
  mprob.subset <- cbind(mprob.rows, mprob.subset)
  mprob.subset
}

mprob.label2 <- function(states){
  n <- length(states)
  v <- states
  combo <- permutations(n=n, r=2, v=v, repeats.allowed=T)
  geno.combo <- paste0(combo[,1], combo[,2])
  geno.combo
}

# deprecate this subsetting matrix function based on gp 
# and 1:1 K:CN correspondence
mprob.subset <- function(mprob.mat, gp) {
  K <- gp$K
  states <- gp$states
  col.a <- states[1] + 2
  col.b <- states[K] + 2
  
  ref.geno <- c("00", "01", "02", "03", "04", 
                "10", "11", "12", "13", "14",
                "20", "21", "22", "23", "24",
                "30", "31", "32", "33", "34",
                "40", "41", "42", "43", "44")
  index <- mprob.label(gp)
  rows <- match(index, ref.geno)
  
  mprob.subset <- mprob.mat[rows, c(col.a:col.b)]
  mprob.rows <- mprob.mat[rows, 1]
  mprob.subset <- cbind(mprob.rows, mprob.subset)
  mprob.subset
}

# deprecate this as well as associated with mprob.subset
mprob.label <- function(gp){
  n <- gp$K
  v <- gp$states
  combo <- permutations(n=n, r=2, v=v, repeats.allowed=T)
  geno.combo <- paste0(combo[,1], combo[,2])
  geno.combo
}

component_stats <- function(tbl){
  tbl %>% group_by(copy_number) %>%
    summarize(mean=mean(log_ratio),
              sd=sd(log_ratio),
              n=n())
}

setMethod("triodata_lrr", "TrioBatchModel", function(object){
  object@triodata$log_ratio
}
)

setMethod("updateZ", "TrioBatchModel", function(object){
  update_ztrio(object)
})
