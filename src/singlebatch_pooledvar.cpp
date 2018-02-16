#include "miscfunctions.h" // for rdirichlet, tableZ, ...
#include "singlebatch.h" // getK
#include <Rmath.h>

using namespace Rcpp ;


// update theta for pooled variance model
// [[Rcpp::export]]
Rcpp::NumericVector theta_pooled(Rcpp::S4 xmod) {
    RNGScope scope ;
    Rcpp::S4 model(xmod) ;
    NumericVector theta = model.slot("theta") ;
    double tau2 = model.slot("tau2") ;
    double tau2_tilde = 1/tau2 ;
    NumericVector sigma2 = model.slot("sigma2") ;
    NumericVector data_mean = model.slot("data.mean") ;
    NumericVector sigma2_tilde = 1.0/sigma2 ;
    IntegerVector z = model.slot("z");
    int K = getK(model.slot("hyperparams"));
    double mu = model.slot("mu");
    //
    // Initialize nn, vector of component-wise sample size
    //
    IntegerVector nn = tableZ(K, z) ;
    double post_prec = 1.0;
    double tau_n = 0.0;
    double mu_n = 0.0;
    double w1 = 0.0;
    double w2 = 0.0;
    NumericVector x = model.slot("data") ;
    //int n = x.size() ;
    NumericVector thetas(K);

    for(int k = 0; k < K; ++k) {
      post_prec = tau2_tilde + sigma2_tilde[0] * nn[k];
      if (post_prec == R_PosInf) {
        throw std::runtime_error("Bad simulation. Run again with different start.");
      }
      tau_n = sqrt(1/post_prec);
      w1 = tau2_tilde/post_prec;
      w2 = nn[k]*sigma2_tilde[0]/post_prec;
      mu_n = w1*mu + w2*data_mean[k];
      thetas[k] = as<double>(rnorm(1, mu_n, tau_n));
    }
    return thetas;
}


// [[Rcpp::export]]
Rcpp::NumericVector theta_heavy(Rcpp::S4 xmod) {
    RNGScope scope ;
    Rcpp::S4 model(xmod) ;
    NumericVector theta = model.slot("theta") ;
    double tau2 = model.slot("tau2") ;
    double tau2_tilde = 1/tau2 ;
    NumericVector sigma2 = model.slot("sigma2") ;
    //NumericVector data_mean = model.slot("data.mean") ;
    NumericVector sigma2_tilde = 1.0/sigma2 ;
    IntegerVector z = model.slot("z");
    int K = getK(model.slot("hyperparams"));
    double mu = model.slot("mu");
    //
    // Initialize nn, vector of component-wise sample size
    //
    IntegerVector counts = tableZ(K, z) ;
    double post_prec = 1.0;
    double tau_n = 0.0;
    double mu_n = 0.0;
    double w1 = 0.0;
    double w2 = 0.0;
    NumericVector x = model.slot("data") ;
    // int n = x.size() ;
    NumericVector thetas(K);
    NumericVector u = model.slot("u") ;
    double df = getDf(model.slot("hyperparams")) ;
    NumericVector data_mean =  compute_heavy_means(xmod) ;
    data_mean =  data_mean/df ;
    NumericVector sumu = compute_u_sums(xmod) ;
    sumu = sumu/df ;
    NumericVector nn(K) ;
    nn = as<NumericVector>(counts)  * sumu ;
    // NumericVector sumu = sum(u);
    // return u;
    //sigma2_tilde[0] = sigma_tilde[0] / sumu[0] * 10.0 ;
    for(int k = 0; k < K; ++k) {
      post_prec = tau2_tilde + sigma2_tilde[0] * nn[k];
      if (post_prec == R_PosInf) {
        throw std::runtime_error("Bad simulation. Run again with different start.");
      }
      tau_n = sqrt(1/post_prec);
      w1 = tau2_tilde/post_prec;
      w2 = nn[k]*sigma2_tilde[0]/post_prec;
      mu_n = w1*mu + w2*data_mean[k];
      thetas[k] = as<double>(rnorm(1, mu_n, tau_n));
    }
    return thetas;
}

// [[Rcpp::export]]
Rcpp::NumericVector loglik_pooled(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  NumericVector p = model.slot("pi") ;
  int K = getK(model.slot("hyperparams")) ;
  NumericVector theta = model.slot("theta") ;
  NumericVector sigma2 = model.slot("sigma2") ;
  NumericVector sigma = sqrt(sigma2) ;
  int n = x.size() ;
  double df = getDf(model.slot("hyperparams")) ;
  //double lik;
  NumericVector loglik(1) ;
  loglik[0] = 0.0 ;
  NumericVector y(1);
  NumericVector lik(n);
  // Below is equivalent to rowSums(lik) in .loglikMarginal
  for(int k = 0; k < K; k++) {
    //lik += p[k]*dnorm(x, theta[k], sigma[0]);
    lik += p[k]*dlocScale_t(x, df, theta[k], sigma[0]);
  }
  for(int i = 0; i < n; i++){
    loglik[0] += log(lik[i]);
  }
  return loglik;
}

// [[Rcpp::export]]
Rcpp::NumericVector stageTwoLogLik_pooled(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  int K = getK(model.slot("hyperparams")) ;
  NumericVector theta = model.slot("theta") ;
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector mu = model.slot("mu") ;
  NumericVector nu0 = model.slot("nu.0") ;
  NumericVector s20 = model.slot("sigma2.0") ;
  NumericVector sigma2=model.slot("sigma2") ;
  NumericVector sigma2_tilde = 1.0/sigma2 ;
  NumericVector loglik(1) ;
  double tau = sqrt(tau2[0]) ;
  NumericVector liknorm(K) ;
  NumericVector likprec(K) ;
  liknorm = dnorm(theta, mu[0], tau) ;
  likprec = dgamma(sigma2_tilde, 0.5*nu0[0], 1.0/(0.5 * nu0[0] * s20[0])) ;
  NumericVector LL(K) ;
  LL = log(liknorm * likprec[0]);
  //  for(int k=0; k<K; ++k){
  //  LL[K] = log(liknorm[k] * likprec[0]) ;
  //}
  double tmp = 0.0 ;
  for(int k = 0; k < K; k++) {
    tmp += LL[k] ;
  }
  loglik[0] = tmp ;
  return loglik ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix multinomialPr_pooled(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector theta = model.slot("theta") ;
  NumericVector sigma2 = model.slot("sigma2") ;
  NumericVector sigma = sqrt(sigma2) ;
  NumericVector p = model.slot("pi") ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  NumericMatrix lik(n, K) ;
  NumericMatrix probs(n, K) ;
  NumericVector tmp(n) ;
  NumericVector total(n) ;
  for(int k = 0; k < K; k++) {
    tmp = p[k]*dnorm(x, theta[k], sigma[0]) ;
    for(int i = 0; i < n; i++){
      lik(i, k) = tmp[i] ;
    }
    total += tmp ;
  }
  for(int k = 0; k < K; k++){
    for(int i = 0; i < n; i++){
      probs(i, k) = lik(i,k) / total[i] ;
    }
  }
  return probs ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix multinomialPr_heavy(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector theta = model.slot("theta") ;
  NumericVector sigma2 = model.slot("sigma2") ;
  NumericVector sigma = sqrt(sigma2) ;
  NumericVector p = model.slot("pi") ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  NumericMatrix lik(n, K) ;
  NumericMatrix probs(n, K) ;
  NumericVector tmp(n) ;
  NumericVector total(n) ;
  double df = getDf(hypp) ;

  for(int k = 0; k < K; k++) {
    tmp = p[k]*dlocScale_t(x, df, theta[k], sigma[0]) ;
    for(int i = 0; i < n; i++){
      lik(i, k) = tmp[i] ;
    }
    total += tmp ;
  }
  for(int k = 0; k < K; k++){
    for(int i = 0; i < n; i++){
      probs(i, k) = lik(i,k) / total[i] ;
    }
  }
  return probs ;
}

// Be smarter about assigning samples to a component with zero or one
// observations
// [[Rcpp::export]]
Rcpp::IntegerVector z_pooled(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  NumericMatrix p(n, K);
  // p = multinomialPr_pooled(xmod) ;
  p = multinomialPr_heavy(xmod) ;
  NumericMatrix cumP(n, K);
  for(int i=0; i < n; i++){
    for(int k = 0; k < K; k++){
      if(k > 0){
        cumP(i, k) = cumP(i, k-1) + p(i, k) ;
      } else {
        cumP(i, k) = p(i, k) ;
      }
    }
    cumP(i, K-1) = 1.000001 ;
  }
  NumericVector u = runif(n) ;
  IntegerVector zz(n) ;
  IntegerVector freq(K) ;
  for(int k = 0; k < K; k++) {
    freq[k]=0;
  }
  for(int i = 0; i < n; i++){
    int k = 0 ;
    while(k < K) {
      if( u[i] < cumP(i, k)){
        zz[i] = k + 1 ;
        freq[k] += 1 ;
        break ;
      }
      k += 1 ;
    }
    cumP(i, K-1) = 1.00001 ;  // just to be certain
  }
  if(is_true(all(freq > 1))){
    return zz ;
  }
  //
  // Don't update z if there are states with zero frequency
  //
  int counter = model.slot(".internal.counter");
  counter++;
  model.slot(".internal.counter") = counter;
  return model.slot("z") ;  
  //  return zz ;
}


// [[Rcpp::export]]
Rcpp::NumericVector nu0_pooled(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  double sigma2_0 = model.slot("sigma2.0") ;
  NumericVector sigma2 = model.slot("sigma2") ;
  double betas = hypp.slot("beta") ;
  //
  // sample nu0 from an unnormalized probability distribution
  //
  int MAX = 1000 ;
  NumericVector x(MAX) ;  // 1000 is the maximum allowed value for nu_0
  NumericVector lpnu0(MAX);
  double prec = 0.0 ;
  double lprec = 0.0 ;
  prec = 1.0/sigma2[0] ;
  //for(int k = 0; k < K; k++) prec += 1.0/sigma2[k] ;
  lprec = log(prec) ;
  //for(int k = 0; k < K; k++) lprec += log(1.0/sigma2[k]) ;
  x = seq_len(MAX) ;
  NumericVector y1(MAX) ;
  NumericVector y2(MAX) ;
  NumericVector y3(MAX) ;
  y1 = K*(0.5*x*log(sigma2_0*0.5*x) - lgamma(x*0.5)) ;
  y2 = (0.5*x - 1.0) * lprec ;
  y3 = x*(betas + 0.5*sigma2_0*prec) ;
  lpnu0 =  y1 + y2 - y3 ;
  NumericVector prob(MAX) ;
  prob = exp(lpnu0) ; // - maxprob) ;
  prob = prob/sum(prob) ;
  //double maxprob = max(lpnu0) ;
  NumericVector nu0(1) ;
  //int u ;
  NumericVector u(1) ;
  double cumprob = 0.0 ;
  // sample x with probability prob
  for(int i = 0; i < MAX; i++){
    cumprob += prob[i] ;
    u = runif(1) ;
    if (u[0] < cumprob){
      nu0[0] = x[i] ;
      break ;
    }
  }
  if(nu0[0] < 1) nu0[0] = 1 ;
  return nu0 ;
}

// [[Rcpp::export]]
Rcpp::NumericVector sigma2_0_pooled(Rcpp::S4 xmod) {
    RNGScope scope ;
    Rcpp::S4 model(xmod) ;
    Rcpp::S4 hypp(model.slot("hyperparams")) ;
    double a = hypp.slot("a") ;
    double b = hypp.slot("b") ;
    double nu_0 = model.slot("nu.0") ;
    NumericVector sigma2 = model.slot("sigma2") ;
    Rcpp::NumericVector sigma2_0_old = model.slot("sigma2.0");
    int K = getK(hypp) ;
    double a_k = a + 0.5*K*nu_0 ;
    double b_k = b + 0.5*nu_0/sigma2[0];

    NumericVector sigma2_0(1);
    sigma2_0[0] = as<double>(rgamma(1, a_k, 1.0/b_k)) ;
    double constraint = model.slot(".internal.constraint");
    if (constraint > 0) {
        if (sigma2_0[0] < constraint) {
            return sigma2_0_old;
        }
    }
    return sigma2_0;
}

// [[Rcpp::export]]
Rcpp::NumericVector sigma2_pooled(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;
    // get model
    Rcpp::S4 model(xmod);

    // get parameter estimates
    Rcpp::NumericVector theta = model.slot("theta");
    Rcpp::IntegerVector z = model.slot("z");
    double nu_0 = model.slot("nu.0");
    double sigma2_0 = model.slot("sigma2.0");

    // get data and size attributes
    Rcpp::NumericVector x = model.slot("data");
    Rcpp::S4 hypp(model.slot("hyperparams"));
    int K = theta.size();
    int n = x.size();

    Rcpp::NumericVector nu_n(1) ;
    nu_n[0] = 0.5*(nu_0 + n);
    Rcpp::NumericVector ss(1);
    ss[0] = 0.0 ;
    for(int i = 0; i < n; i++){
        int k = 0;
        while (k < K) {
            if ( z[i] == k + 1 ) {
              ss[0] += pow(x[i] - theta[k], 2.0);
              break;
            }
            k++;
        }
    }
    Rcpp::NumericVector sigma2_new(1);
    double sigma2_n = 1.0 / (0.5*(nu_0 * sigma2_0 + ss[0])) ;
    sigma2_new[0] = 1.0 / Rcpp::as<double>(rgamma(1, nu_n[0], sigma2_n)) ;
    return sigma2_new ;
}

// [[Rcpp::export]]
Rcpp::NumericVector sigma2_heavy(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;
    // get model
    Rcpp::S4 model(xmod);

    // get parameter estimates
    Rcpp::NumericVector theta = model.slot("theta");
    Rcpp::IntegerVector z = model.slot("z");
    double nu_0 = model.slot("nu.0");
    double sigma2_0 = model.slot("sigma2.0");

    // get data and size attributes
    Rcpp::NumericVector x = model.slot("data");
    Rcpp::S4 hypp(model.slot("hyperparams"));
    int K = theta.size();
    int n = x.size();
    NumericVector u = model.slot("u") ;
    double df = getDf(model.slot("hyperparams")) ;

    Rcpp::NumericVector nu_n(1) ;
    nu_n[0] = 0.5*(nu_0 + n);
    Rcpp::NumericVector ss(1);
    ss[0] = 0.0 ;
    for(int i = 0; i < n; i++){
        int k = 0;
        while (k < K) {
            if ( z[i] == k + 1 ) {
              ss[0] += u[i] * pow(x[i] - theta[k], 2.0);
              break;
            }
            k++;
        }
    }
    Rcpp::NumericVector sigma2_new(1);
    double sigma2_n = 1.0 / (0.5*(nu_0 * sigma2_0 + ss[0]/df)) ;
    sigma2_new[0] = 1.0 / Rcpp::as<double>(rgamma(1, nu_n[0], sigma2_n)) ;
    return sigma2_new ;
}


//
// BURNIN: Don't move chains, no thinning, no need to compute loglik at each iteration
//
// [[Rcpp::export]]
Rcpp::S4 burnin_singlebatch_pooled(Rcpp::S4 object, Rcpp::S4 mcmcp) {
  RNGScope scope ;
  Rcpp::S4 xmod = clone(object) ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  Rcpp::S4 params(mcmcp) ;
  IntegerVector up = params.slot("param_updates") ;
  int S = params.slot("burnin") ;
  double df = getDf(model.slot("hyperparams")) ;
  NumericVector x = model.slot("data") ;
  int N = x.size() ;
  if( S < 1 ){
    return xmod ;
  }
  for(int s = 0; s < S; ++s){
    if(up[7] > 0){
      model.slot("z") = z_pooled(xmod) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    }
    model.slot("data.mean") = compute_means(xmod) ;
    model.slot("data.prec") = compute_prec(xmod) ;
    if(up[0] > 0)
      model.slot("theta") = theta_heavy(xmod) ;
      //model.slot("theta") = theta_pooled(xmod) ;
    if(up[1] > 0)
      model.slot("sigma2") = sigma2_heavy(xmod) ;
    // model.slot("sigma2") = sigma2_pooled(xmod) ;
    if(up[2] > 0)
      model.slot("pi") = update_p(xmod) ;
    if(up[3] > 0)
      model.slot("mu") = update_mu(xmod) ;
    if(up[4] > 0)
      model.slot("tau2") = update_tau2(xmod) ;
    if(up[5] > 0)
      model.slot("nu.0") = nu0_pooled(xmod) ;
    if(up[6] > 0)
      model.slot("sigma2.0") = sigma2_0_pooled(xmod) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
  }
  // compute log prior probability from last iteration of burnin
  // compute log likelihood from last iteration of burnin
  NumericVector ll = loglik_pooled(xmod) ;
  //NumericVector lls2 = stageTwoLogLik(xmod) ;
  NumericVector lls2 = stageTwoLogLik_pooled(xmod) ;
  model.slot("loglik") = ll + lls2 ;
  model.slot("logprior") = compute_logprior(xmod) ;
  return xmod ;
}


// Function has gotten pretty long. Might be useful to have a separate
// function whose sole job is to move the chains.
//
// [[Rcpp::export]]
Rcpp::S4 mcmc_singlebatch_pooled(Rcpp::S4 object, Rcpp::S4 mcmcp) {
  RNGScope scope ;
  Rcpp::S4 xmod = clone(object) ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 chain(model.slot("mcmc.chains")) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  Rcpp::S4 params(mcmcp) ;
  IntegerVector up = params.slot("param_updates") ;
  int K = getK(hypp) ;
  int T = params.slot("thin") ;
  int S = params.slot("iter") ;
  if( S < 1 ) return xmod ;
  NumericVector x = model.slot("data") ;
  double df = getDf(model.slot("hyperparams")) ;
  //NumericVector u = model.slot("u") ;
  int N = x.size() ;
  NumericMatrix theta = chain.slot("theta") ;
  NumericMatrix sigma2 = chain.slot("sigma2") ;
  NumericMatrix pmix = chain.slot("pi") ;
  NumericMatrix zfreq = chain.slot("zfreq") ;
  IntegerMatrix Z = chain.slot("z") ;
  NumericMatrix U = chain.slot("u") ;
  NumericVector mu = chain.slot("mu") ;
  NumericVector tau2 = chain.slot("tau2") ;
  NumericVector nu0 = chain.slot("nu.0") ;
  NumericVector sigma2_0 = chain.slot("sigma2.0") ;
  NumericVector loglik_ = chain.slot("loglik") ;
  NumericVector logprior_ = chain.slot("logprior") ;
  IntegerMatrix pz(N, K) ;
  // initialize probabilities to zero
  model.slot("probz") = pz ;
  NumericVector th(K) ;
  NumericVector s2(1) ;
  NumericVector p(K) ;
  NumericVector m(1) ; //mu
  NumericVector t2(1) ;//tau2
  NumericVector n0(1) ;//nu0
  IntegerVector z(N) ;
  NumericVector u(N) ;
  NumericVector s20(1) ; //sigma2_0
  NumericVector mns(1) ;
  NumericVector precs(1) ;
  NumericVector ll(1) ;
  NumericVector lls2(1) ; // stage 2 log likelihood
  NumericVector lp(1) ;
  NumericVector stage2(1) ;
  IntegerVector tmp(K) ;
  IntegerVector zf(K) ;
  // Initial values
  th = model.slot("theta") ;
  s2 = model.slot("sigma2") ;
  p = model.slot("pi") ;
  m = model.slot("mu") ;
  t2 = model.slot("tau2") ;
  n0 = model.slot("nu.0") ;
  s20 = model.slot("sigma2.0") ;
  zf = model.slot("zfreq") ;
  z = model.slot("z") ;
  u = model.slot("u") ;
  ll = model.slot("loglik") ;
  lp = model.slot("logprior") ;
  // Record initial values in chains
  mu[0] = m[0] ;
  tau2[0] = t2[0] ;
  nu0[0] = n0[0] ;
  sigma2_0[0] = s20[0] ;
  loglik_[0] = ll[0] ;
  logprior_[0] = lp[0] ;
  //  stagetwo = stage2[0] ;
  theta(0, _) = th ;
  sigma2(0, 0) = s2[0] ;
  pmix(0, _) = p ;
  zfreq(0, _) = zf ;
  Z(0, _) = z ;
  U(0, _) = u ;
  // start at 1 instead of zero. Initial values are as above
  for(int s = 1; s < S; ++s){
    if(up[7] > 0){
      z = z_pooled(xmod) ;
      model.slot("z") = z ;
      tmp = tableZ(K, z) ;
      model.slot("zfreq") = tmp ;
      model.slot("probz") = compute_probz(xmod) ;
    } else {
      z = model.slot("z") ;
      tmp = model.slot("zfreq") ;
    }
    Z(s, _) = z ;
    zfreq(s, _) = tmp ;
    model.slot("data.mean") = compute_means(xmod) ;
    model.slot("data.prec") = compute_prec(xmod) ;
    if(up[0] > 0) {
      try {
          th = theta_heavy(xmod) ;
      } catch(std::runtime_error &ex) {
          forward_exception_to_r(ex);
      } catch(...) {
          ::Rf_error("c++ exception (unknown reason)");
      }
      model.slot("theta") = th ;
    } else {
      th = model.slot("theta") ;
    }
    theta(s, _) = th ;
    if(up[1] > 0){
      s2 = sigma2_heavy(xmod) ;
      model.slot("sigma2") = s2 ;
    } else { s2= model.slot("sigma2") ; }
    sigma2(s, 0) = s2[0] ;
    if(up[2] > 0){
      p = update_p(xmod) ;
      model.slot("pi") = p ;
    } else {
      p = model.slot("pi") ;
    }
    pmix(s, _) = p ;
    if(up[3] > 0){
      m = update_mu(xmod) ;
      model.slot("mu") = m ;
    } else {
      m = model.slot("mu") ;
    }
    mu[s] = m[0] ;
    if(up[4] > 0){
      t2 = update_tau2(xmod) ;
      model.slot("tau2") = t2 ;
    } else {
      t2 = model.slot("tau2") ;
    }
    tau2[s] = t2[0] ;
    if(up[5] > 0){
      n0 = nu0_pooled(xmod) ;
      model.slot("nu.0") = n0 ;
    } else {
      n0 = model.slot("nu.0") ;
    }
    nu0[s] = n0[0] ;
    if(up[6] > 0){
      s20 = sigma2_0_pooled(xmod) ;
      model.slot("sigma2.0") = s20 ;
    } else {
      s20 = model.slot("sigma2.0") ;
    }
    sigma2_0[s] = s20[0] ;
    ll = loglik_pooled(xmod) ;
    lls2 = stageTwoLogLik_pooled(xmod) ;
    ll = ll + lls2 ;
    loglik_[s] = ll[0] ;
    model.slot("loglik") = ll ;
    lp = compute_logprior(xmod) ;
    logprior_[s] = lp[0] ;
    model.slot("logprior") = lp ;
    u = Rcpp::rchisq(N, df) ;
    model.slot("u") = u;
    U(s, _) = u ;
    // Thinning
    for(int t = 0; t < T; ++t){
      if(up[7] > 0){
        model.slot("z") = z_pooled(xmod) ;
        model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      }
      model.slot("data.mean") = compute_means(xmod) ;
      model.slot("data.prec") = compute_prec(xmod) ;
      if(up[0] > 0)
        model.slot("theta") = theta_heavy(xmod) ;
      if(up[1] > 0)
        model.slot("sigma2") = sigma2_heavy(xmod) ;
      if(up[2] > 0)
        model.slot("pi") = update_p(xmod) ;
      if(up[3] > 0)
        model.slot("mu") = update_mu(xmod) ;
      if(up[4] > 0)
        model.slot("tau2") = update_tau2(xmod) ;
      if(up[5] > 0)
        model.slot("nu.0") = nu0_pooled(xmod) ;
      if(up[6] > 0)
        model.slot("sigma2.0") = sigma2_0_pooled(xmod) ;
      model.slot("u") = Rcpp::rchisq(N, df) ;
    }
  }
  //
  // assign chains back to object
  //
  chain.slot("theta") = theta ;
  chain.slot("sigma2") = sigma2 ;
  chain.slot("pi") = pmix ;
  chain.slot("mu") = mu ;
  chain.slot("tau2") = tau2 ;
  chain.slot("nu.0") = nu0 ;
  chain.slot("sigma2.0") = sigma2_0 ;
  chain.slot("zfreq") = zfreq ;
  chain.slot("loglik") = loglik_ ;
  chain.slot("logprior") = logprior_ ;
  chain.slot("z") = Z ;
  chain.slot("u") = U ;
  model.slot("mcmc.chains") = chain ;
  return xmod ;
}
