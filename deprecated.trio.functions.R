##################################
#deprecated functions
###############################

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

# deprecate the following simulation functions:
.simulate_data_multi <- function(params, n, mendelian.probs, GP){
  gp <- GP
  states <- gp$states + 1
  K <- nrow(params)
  p <- params$p
  theta <- params$theta
  sigma2 <- params$sigma2
  sigma <- sqrt(sigma2)
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
  stats <- component_stats(tbl3)
  p <- stats %>% group_by(copy_number) %>% summarise(p=mean(p))
  sd <- stats %>% group_by(copy_number) %>% summarise(sd=mean(sd))
  mean <- stats %>% group_by(copy_number) %>% summarise(mean=mean(mean))
  params$p <- p$p
  params$sigma2 <- (sd$sd)^2
  params$theta <- mean$mean
  loglik <- sum(dnorm(tbl$log_ratio, params$theta[z],
                      sqrt(params$sigma2[z]), log=TRUE))
  truth <- list(data=tbl3, params=params, loglik=loglik)
  truth
}
