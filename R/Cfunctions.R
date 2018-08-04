## R functions for checking C code

## R version of marginal_theta
## r_marginal_theta <- function(model){
##   model2 <- useModes(model)
##   thetastar <- theta(model2)
##   S <- nrow(theta(chains(model)))
##   Z <- z(chains(model))
##   ##U <- u(chains(model))
##   K <- k(model)
##   N <- length(y(model))
##   df <- dfr(model)
##   tau2c <- tau2(chains(model))
##   sigma2 <- sigma2(chains(model))
##   muc <- mu(chains(model))
##   p_theta <- rep(NA, S)
##   for(s in seq_len(S)){
##     zz <- Z[s, ]
##     ##zz = Z(s, _) ;
##     ##uu = U(s, _) ;
##     uu <-  rchisq(N, df) ;
##     ##uu <- U[s, ]
##     z(model) <- zz
##     model@u <- uu
##     ##model.slot("z") = zz ;
##     ##model.slot("u") = uu ;
##     counts <- tableZ(K, zz) ;
##     ##NumericVector data_mean =  compute_heavy_means(xmod) ;
##     data_mean <- compute_heavy_means(model)/df
##     ##data_mean =  data_mean/df ;
##     ##NumericVector sumu = compute_u_sums(xmod) ;
##     sumu <- compute_u_sums(model)/df
##     sumu = sumu/df ;
##     ##nn = as<NumericVector>(counts) * sumu ;
##     nn <- counts
##     ##nn <- counts*sumu
##     nn <- sumu
##     tau2_tilde = 1/tau2c[s]
##     ##CHECK: length-k vector?
##     s2 <- sigma2[s, ]
##     ##sigma2_tilde = 1.0/sigma2[s, ] ;
##     ##CHECK: length-k vector?
##     sigma2_tilde <- 1.0/s2
##     ##tmp = dnorm(thetastar, muc[s], tauc[s]) ;
##     ##double prod = 1.0;
##     prod <- 1.0
##     ##for(int k = 0; k < K; ++k) {
##     for(k in seq_len(K)){
##       post_prec <- tau2_tilde + sigma2_tilde[k] * nn[k];
##       tau_n = sqrt(1/post_prec);
##       w1 = tau2_tilde/post_prec;
##       w2 = nn[k]*sigma2_tilde[k]/post_prec;
##       mu_n = w1*muc[s] + w2*data_mean[k];
##       tmp = dnorm(thetastar[k], mu_n, tau_n) ;
##       prod = prod * tmp[1] ;
##     }
##     p_theta[s] = prod ;
##   }
## }

##p_theta_batch2 <- function(model, thetastar){
##  prod <- 1.0
##  B <- length(unique(batch(model)))
##  K <- k(model)
##  n_hb <- tableBatchZ(model)
##  data_mean = compute_heavy_sums_batch(model);
##  sumu = compute_u_sums_batch(model) ;
##  df <- dfr(model)
##  ##invs2 <- 1.0 / sigma2(s, Rcpp::_);    // this is a vector of length B*K
##  invs2 <- 1.0/sigma2(model)
##  ##sigma2_tilde <- toMatrix(invs2, B, K);
##  sigma2_tilde <- invs2
##  for (b in seq_len(B)){
##    for(k in seq_len(K)){
##      heavyn <- n_hb[b, k] * sumu[b, k] / df;
##      heavy_mean <- data_mean[b, k] / df;
##      post_prec <- 1.0/tau2[k] + heavyn*1.0 * sigma2_tilde[b, k] ;
##      tau_n <- sqrt(1.0/post_prec) ;
##      w1 <- (1.0/tau2[k])/post_prec ;
##      w2 <- (n_hb[b, k] * sigma2_tilde[b, k])/post_prec ;
##      mu_n <- w1*mu[k] + w2*heavy_mean ;
##      theta <- thetastar[b, k];
##      tmp <- dnorm(theta, mu_n, tau_n);
##      prod <- prod*tmp[0];
##    }
##  }
##  prod
##}
##
##r_marginal_theta_batch <- function(model){
##  thetastar <- modes(model)[["theta"]];
##  xmod <- model
##  S <- iter(model)
##  p_theta <- rep(NA, S)
##  for(i in seq_len(S)){
##    ##model.slot("z") = update_z_batch(xmod) ;
##    z(model) <- update_z_batch(xmod) ;
##    ##model.slot("zfreq") = tableZ(K, model.slot("z")) ;
##    zFreq(model) <- tableZ(K, z(model)) 
##    ##model.slot("data.mean") = compute_means_batch(xmod) ;
##    dataMean(model) <- compute_means_batch(xmod) ;
##    ##model.slot("data.prec") = compute_prec_batch(xmod) ;
##    ##dataPrec(model) <- compute_prec_batch(xmod) ;
##    ##model.slot("theta") = update_theta_batch(xmod) ;
##    theta(model) <- update_theta_batch(xmod) ;
##    model.slot("sigma2") = update_sigma2_batch(xmod) ;
##    model.slot("mu") = update_mu_batch(xmod) ;
##    model.slot("tau2") = update_tau2_batch(xmod) ;
##    model.slot("sigma2.0") = update_sigma20_batch(xmod) ;
##    model.slot("nu.0") = update_nu0_batch(xmod) ;
##    model.slot("pi") = update_p_batch(xmod) ;
##    model.slot("u") = Rcpp::rchisq(N, df) ;
##    p_theta[i] <- p_theta_batch(model, thetastar)
##  }
##  p_theta
##
##  for(s in seq_len(S)){
##    tauc <- sqrt(tau2c[s, ])
##    ##tauc = sqrt(tau2c(s, Rcpp::_));
##    ##// extract z from the chain
##    ##zz = Z(s, Rcpp::_);
##    zz <- Z[s, ]
##    ##// add z to the current slot in the model
##    model@z <- zz
##    ## make a B x k matrix of the counts
##    n_hb = tableBatchZ(model);
##    data_mean <- compute_means_batch(model)
##    ##data_mean = compute_heavy_sums_batch(model);
##    ##sumu = compute_u_sums_batch(model) ;
##    mu = muc(s, Rcpp::_);
##    tau2 = tau2c(s, Rcpp::_);
##    tau2_tilde = 1.0 / tau2;
##    invs2 = 1.0 / sigma2(s, Rcpp::_);    // this is a vector of length B*K
##    sigma2_tilde = Rcpp::as<Rcpp::NumericVector>(toMatrix(invs2, B, K));
##    double prod = 1.0;
##    double heavyn = 0.0;
##    double heavy_mean = 0.0;
##    for (int b = 0; b < B; ++b) {
##      for(int k = 0; k < K; ++k){
##        heavyn = n_hb(b, k) * sumu(b, k) / df;
##        heavy_mean = data_mean(b, k) / df;
##        post_prec = 1.0/tau2[k] + heavyn*1.0 * sigma2_tilde(b, k) ;
##        if (post_prec == R_PosInf) {
##          throw std::runtime_error("Bad simulation. Run again with different start.");
##        }
##        tau_n = sqrt(1.0/post_prec) ;
##        w1 = (1.0/tau2[k])/post_prec ;
##        w2 = (n_hb(b, k) * sigma2_tilde(b, k))/post_prec ;
##        // mu_n = w1*mu[k] + w2*data_mean(b, k) ;
##        mu_n = w1*mu[k] + w2*heavy_mean ;
##        //theta_new(b, k) = as<double>(rnorm(1, mu_n, tau_n)) ;
##        theta[0] = thetastar(b, k);
##        tmp = dnorm(theta, mu_n, tau_n);
##        prod *= tmp[0];
##      }
##    }
##    p_theta[s] = prod;
##  }
##  return p_theta;
##}

##probSigma2 <- function(model, sigma2star){
##  prec <- 1/sigma2star
##  B <- batch(model)
##  K <- k(model)
##  ub <- unique(b)
##  zz <- z(model)
##  uu <- u(model)
##  ss <- matrix(NA, nrow(prec), ncol(prec))
##  thetastar <- modes(model)[["theta"]]
##  yy <- y(model)
##  for(b in seq_along(ub)){
##    for(k in seq_len(K)){
##      index <- B == ub[i] & zz == k
##      ss[b, k] <- sum(uu[index] * (yy[index] - thetastar[b, k])^2)
##    }
##  }
##  nn <- matrix(NA, nrow(thetastar), ncol(thetastar))
##  nn <- tableBatchZ(model)
##  nu0 <- nu.0(model)
##  s20 <- sigma2.0(model)
##  df <- dfr(model)
##  total <- 0
##  for(b in seq_along(ub)){
##    for(k in seq_len(K)){
##      ## calculate nu_n and sigma2_n
##      nu_n = nu0 + nn[b, k];
##      sigma2_n = 1.0 / nu_n * (nu0 * s20 + ss[b, k]/df);
##      ## calculate shape and rate
##      shape = 0.5 * nu_n;
##      rate = shape * sigma2_n;
##      ## calculate probability
##      prec_typed = prec[b, k];
##      ##tmp = Rcpp::dgamma(prec_typed, shape, 1/rate, TRUE);
##      tmp <- dgamma(prec_typed, shape=shape, rate=rate, log=TRUE)
##      ##tmp2 <- rgamma(1000, shape=shape, rate=rate)
##      ##hist(tmp2, breaks=250)
##      ##abline(v=prec[1])
##      total <- total+tmp;
##    }
##  }
##  total
##}

probTheta <- function(model){
  thetastar <- modes(model)[["theta"]]
  zz <- z(model)
  K <- k(model)
  nn = tableZ(K, zz);
  tau2 <- tau2(model)
  tauc = sqrt(tau2);
  data_sum = compute_heavy_sums_batch(model);
  sumu = compute_u_sums_batch(model) ;
  sigma2 <- sigma2(model)
  ##data_mean = compute_means_batch(model);
  tau2_tilde = 1.0 / tau2;
  invs2 = 1.0 / sigma2;    ##this is a vector of length B
  sigma2_tilde = matrix(invs2, nrow=1)
  total = 0.0;
  for (k in seq_len(K)) {
    for (b in seq_len(nrow(thetastar))) {
      heavyn = sumu[b, k] / df;
      post_prec = tau2_tilde[k] + heavyn*1.0 * sigma2_tilde[b] ;
      tau_n = sqrt(1/post_prec);
      w1 = tau2_tilde[k]/post_prec;
      w2 = heavyn * sigma2_tilde[b, k]/post_prec;
      heavy_mean = data_sum[b, k] / heavyn / df;
      mu_n = w1*mu[k] + w2*heavy_mean;
      theta = thetastar[b, k];
      tmp = dnorm(theta, mu_n, tau_n, log=TRUE);
      total = total + tmp;
    }
  }
  total
}

updateThetaPooled <- function(model){
  n_hb <- tableBatchZ(model)
  sigma2 <- sigma2(model)
  mu <- mu(model)
  tau2 <- tau2(model)
  thetas <- theta(model)
  B <- nrow(thetas)
  K <- ncol(thetas)
  df <- dfr(model)
  data_sum =  compute_heavy_sums_batch(model) ;
  sumu = compute_u_sums_batch(model) ;
  theta_new <- thetas
  heavy_mean = 0.0;
  for(b in seq_len(B)){
    for(k in seq_len(K)){
      heavyn = sumu[b, k] / df;
      post_prec = 1.0/tau2[k] + heavyn*1.0/sigma2[b] ;
      tau_n = sqrt(1.0/post_prec) ;
      w1 = (1.0/tau2[k])/post_prec ;
      w2 = (heavyn * 1.0/sigma2[b])/post_prec ;
      heavy_mean = data_sum[b, k] / heavyn / df;
      mu_n = w1*mu[k] + w2*heavy_mean ;
      theta_new[b, k] = rnorm(1, mu_n, tau_n) ;
    }
  }
  theta_new
}

updateThetaBatch <- function(model){
  n_hb <- tableBatchZ(model)
  sigma2 <- sigma2(model)
  mu <- mu(model)
  tau2 <- tau2(model)
  thetas <- theta(model)
  B <- nrow(thetas)
  K <- ncol(thetas)
  df <- dfr(model)
  data_sum =  compute_heavy_sums_batch(model) ;
  sumu = compute_u_sums_batch(model) ;
  theta_new <- thetas
  heavy_mean = 0.0;
  for(b in seq_len(B)){
    for(k in seq_len(K)){
      heavyn = sumu[b, k] / df;
      post_prec = 1.0/tau2[k] + heavyn*1.0/sigma2[b,k] ;
      tau_n = sqrt(1.0/post_prec) ;
      w1 = (1.0/tau2[k])/post_prec ;
      w2 = (heavyn * 1.0/sigma2[b,k])/post_prec ;
      heavy_mean = data_sum[b, k] / heavyn / df;
      mu_n = w1*mu[k] + w2*heavy_mean ;
      theta_new[b, k] = rnorm(1, mu_n, tau_n) ;
    }
  }
  theta_new
}

updateZbatch <- function(model) {
  K = k(model)
  x = y(model)
  nn <- length(x)
  theta <- theta(model)
  batch=batch(model)
  B = nrow(theta)
  p <- matrix(NA, nn, K)
  p = update_multinomialPr_batch(model) ;
  u = runif(nn) ;
  zz=rep(NA, length(nn))
  freq <- matrix(0, B, K)
  for (i in seq_len(nn)){
    ##initialize accumulator ;
    acc = 0 ;
    for(k in seq_len(K)){
      acc = acc + p[i, k] ;
      if( u[i] < acc ) {
        zz[i] = k + 1 ;
        b = batch[i] ;
        freq[b, k] = freq[b, k] + 1 ;
        break ;
      }
    }
  }
  if(all(freq > 1)){
    return (zz) ;
  }
  return(z(model))
}
