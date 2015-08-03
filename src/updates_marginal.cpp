#include "update.h" // getK
#include "miscfunctions.h" // for rdirichlet, tableZ, ...
#include <Rmath.h>
#include <Rcpp.h>

using namespace Rcpp ;

// [[Rcpp::export]]
Rcpp::NumericVector loglik(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  NumericVector p = model.slot("pi") ;
  int K = getK(model.slot("hyperparams")) ;
  NumericVector theta = model.slot("theta") ;
  NumericVector sigma2 = model.slot("sigma2") ;
  NumericVector sigma = sqrt(sigma2) ;
  int n = x.size() ;
  //double lik;
  NumericVector loglik(1) ;
  NumericVector y(1);    
  NumericVector lik(n);
  // Below is equivalent to rowSums(lik) in .loglikMarginal
  for(int k = 0; k < K; k++) {
    lik += p[k]*dnorm(x, theta[k], sigma[k]);
  }
  for(int i = 0; i < n; i++){
    loglik[0] += log(lik[i]);
  }
  return loglik;
}

// [[Rcpp::export]]
Rcpp::NumericVector log_ddirichlet_(Rcpp::NumericVector x_, 
                                    Rcpp::NumericVector alpha_) {
  // NumericVector x = as<NumericVector>(x_) ;
  NumericVector x = clone(x_) ;
  int K = x.size() ;
  NumericVector alpha = clone(alpha_) ;
  // NumericVector alpha = as<NumericVector>(alpha_) ;
  NumericVector total_lg(1) ;
  NumericVector tmp(1);
  NumericVector total_lalpha(1) ;
  NumericVector logD(1) ;
  NumericVector result(1) ;
  double s = 0.0 ;
  for(int k=0; k < K; ++k){
    total_lg[0] += lgamma(alpha[k]) ;
    tmp[0] += alpha[k] ;
    s += (alpha[k] - 1.0) * log(x[k]) ;
  }
  total_lalpha = lgamma(tmp[0]) ;
  logD[0] = total_lg[0] - total_lalpha[0] ;
  result[0] = s - logD[0] ;
  return result ;
}

// [[Rcpp::export]]
Rcpp::NumericVector stageTwoLogLik(Rcpp::S4 xmod) {
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
  LL = log(liknorm * likprec) ;
  double tmp = 0.0 ;
  for(int k = 0; k < K; k++) {
    tmp += LL[k] ;
  }
  loglik[0] = tmp ;
  return loglik ;
}

//
// This function does not reproduce the R update .updateMu when the
// same seed is used...
//


// [[Rcpp::export]]
Rcpp::NumericVector update_mu(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;  
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  double tau2_0 = hypp.slot("tau2.0") ;
  double tau20_tilde = 1.0/tau2_0 ;
  //double tau2 = model.slot("tau2") ;
  NumericVector tau2 = model.slot("tau2") ;
  // double tau2_tilde = 1.0/tau2 ;
  NumericVector tau2_tilde = 1.0/tau2 ;
  //if(is_true(any(is_nan(tau2_tilde)))) tau2_tilde[0] = 1.0 ;
  double mu_0 = hypp.slot("mu.0") ;
  double K = getK(hypp) ;
  NumericVector theta = model.slot("theta") ;
  //NumericVector nn = model.slot("zfreq") ;
  IntegerVector z = model.slot("z") ;
  IntegerVector nn = tableZ(K, z) ;
  double thetabar ;
  double total = 0.0 ;
  for(int k = 0; k < K; k++) total += nn[k] ;
  for(int k = 0; k < K; k++) thetabar += nn[k] * theta[k] / total ;
  double mu_K ;
  double post_prec = 1.0/tau2_0 + K*tau2_tilde[0] ;
  double w1 ;
  double w2 ;
  w1 = tau20_tilde/post_prec ;
  w2 = K*tau2_tilde[0]/post_prec ;
  mu_K =  w1*mu_0 +  w2*thetabar ;
  NumericVector mu_new(1);
  double tau_k = sqrt(1.0/post_prec) ;
  mu_new[0] = as<double>(rnorm(1, mu_K, tau_k)) ;
  if(is_true(any(is_nan(mu_new)))){
    mu_new[0] = as<double>(rnorm(1, mu_0, sqrt(tau2_0))) ;
  }
  return mu_new ;
  //return mu_K ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_tau2(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;  
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  double m2_0 = hypp.slot("m2.0") ;
  int K = getK(hypp) ;
  double eta_0 = hypp.slot("eta.0") ;
  // should eta_k be double or integer?
  double eta_k = eta_0 + K ;
  double mu = model.slot("mu") ;
  NumericVector theta = model.slot("theta") ;
  NumericVector s2_k(1) ;
  for(int k = 0; k < K; k++) s2_k[0] += (theta[k] - mu) * (theta[k] - mu) ;
  NumericVector m2_k(1) ;
  m2_k[0] = 1/eta_k*(eta_0*m2_0 + s2_k[0]) ;

  NumericVector tau2(1) ;
  //  rgamma is parameterized by scale.  In R, I've parameterized by rate
  tau2[0] = 1.0/as<double>(rgamma(1, 0.5*eta_k, 1.0/(0.5*eta_k*m2_k[0]))) ;
  if(is_true(any(is_nan(tau2)))){
    tau2[0] = 1.0/as<double>(rgamma(1, 0.5*eta_0, 1.0/(0.5*eta_0*m2_0))) ;
  }
  //   tau2[0] = 1/as<double>(rgamma(1, 1/2*eta_k, (1/2*eta_k*m2_k[0]))) ;
  // In R, I check that this is valid and simulate from the prior if not
  return tau2 ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_sigma2_0(Rcpp::S4 xmod) {
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
    double b_k ;

    for (int k=0; k < K; k++) {
      b_k += 0.5*nu_0/sigma2[k];
    }

    b_k += b ;
    NumericVector sigma2_0(1);
    sigma2_0[0] = as<double>(rgamma(1, a_k, 1.0/b_k)) ;
    double constraint = model.slot(".internal.constraint");
    if (constraint > 0) {
        if (sigma2_0[0] < constraint) {
            return sigma2_0_old;
        }
        else {
            return sigma2_0;
        }
    }
    else {
        return sigma2_0;
    }
}

// [[Rcpp::export]]
Rcpp::NumericVector update_nu0(Rcpp::S4 xmod) {
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
  NumericVector x(100) ;  // 100 is the maximum allowed value for nu_0
  NumericVector lpnu0(100);
  double prec = 0.0 ;
  double lprec = 0.0 ;
  for(int k = 0; k < K; k++) prec += 1.0/sigma2[k] ;
  for(int k = 0; k < K; k++) lprec += log(1.0/sigma2[k]) ;
  x = seq_len(100) ;
  NumericVector y1(100) ;
  NumericVector y2(100) ;
  NumericVector y3(100) ;
  y1 = K*(0.5*x*log(sigma2_0*0.5*x) - lgamma(x*0.5)) ;
  y2 = (0.5*x - 1.0) * lprec ;
  y3 = x*(betas + 0.5*sigma2_0*prec) ;
  lpnu0 =  y1 + y2 - y3 ;
  NumericVector prob(100) ;
  prob = exp(lpnu0) ; // - maxprob) ;
  prob = prob/sum(prob) ;
  //double maxprob = max(lpnu0) ;
  NumericVector nu0(1) ;
  //int u ;
  NumericVector u(1) ;
  double cumprob = 0.0 ;
  // sample x with probability prob
  for(int i = 0; i < 100; i++){
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
Rcpp::NumericVector update_p(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;  
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  // IntegerVector nn = model.slot("zfreq");
  IntegerVector z = model.slot("z") ;
  IntegerVector nn = tableZ(K, z) ;  
  IntegerVector alpha = hypp.slot("alpha") ;
  NumericVector alpha_n(K) ;  // really an integer vector, but rdirichlet expects numeric
  for(int k=0; k < K; k++) alpha_n[k] = alpha[k] + nn[k] ;
  NumericVector p(K) ;
  // pass by reference
  rdirichlet(alpha_n, p) ;
  return p ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_multinomialPr(Rcpp::S4 xmod) {
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
    tmp = p[k]*dnorm(x, theta[k], sigma[k]) ;
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
Rcpp::IntegerVector update_z(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;  
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  NumericMatrix p(n, K);
  p = update_multinomialPr(xmod) ;
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
  //return cumP ;
  NumericVector u = runif(n) ;
  IntegerVector zz(n) ;
  IntegerVector freq(K) ;
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
  return model.slot("z") ;  
  // To prevent 0 frequencies, arbitrarily switch the label
  //while(!is_true(all(freq > 0))){
  for(int k = 0; k < K; ++k){
    if( freq[k] >= 2 ) continue ;
    NumericVector r(2) ;
    IntegerVector i(2) ;
    r = runif(2, 0, 1) * n ;
    for(int j = 0; j < 2; ++j){
      // cast as integer
      i = (int) r[j] ;
      zz[i] = k + 1 ;
    }
    freq[k] = sum(zz == (k+1)) ;
  }
  return zz ;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_means(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  IntegerVector z = model.slot("z") ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;  
  // IntegerVector nn = model.slot("zfreq") ;
  IntegerVector nn = tableZ(K, z) ;
  NumericVector means( K ) ;
  for(int i = 0; i < n; i++){
    for(int k = 0; k < K; k++){
      if(z[i] == k+1){
        means[k] += x[i] ;
      }
    }
  }
  for(int k = 0; k < K; k++){
    means[k] /= nn[k] ;
  }
  return means ;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_vars(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  IntegerVector z = model.slot("z") ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;  
  // IntegerVector nn = model.slot("zfreq") ;
  IntegerVector nn ;
  nn = tableZ(K, z) ;
  //  return nn ;
  NumericVector mn = model.slot("theta") ;
  if(is_true(any(is_nan(mn)))){
    mn = compute_means(xmod) ;
  }
  //NumericVector mn = model.slot("data.mean") ;
  NumericVector vars(K) ;
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector is_z(n) ;
  //  for(int i = 0; i < n; i++){
  for(int k = 0; k < K; k++){
    is_z = z == (k + 1) ;
    if(nn[k] <= 1){
      vars[k] = tau2[0] ;
    } else {
      vars[k] = sum(pow(x - mn[k], 2.0) * is_z) / (nn[k]-1) ;
    }      
  }
  return vars ;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_prec(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;  
  NumericVector vars(K) ;
  NumericVector prec(K) ;
  vars = compute_vars(xmod) ;
  for(int k = 0; k < K; ++k){
    prec[k] = 1.0/vars[k] ;
  }
  return prec ;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_logprior(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector mu = model.slot("mu") ;
  NumericVector sigma2_0 = model.slot("sigma2.0") ;
  double a = hypp.slot("a") ;
  double b = hypp.slot("b") ;
  NumericVector nu_0 = model.slot("nu.0") ;
  double mu_0 = hypp.slot("mu.0") ;
  double tau2_0 = hypp.slot("tau2.0") ;
  double tau_0 = sqrt(tau2_0) ;
  double betas = hypp.slot("beta") ;
  NumericVector pmix = model.slot("pi") ;

  double eta = hypp.slot("eta.0") ;
  double m2 = hypp.slot("m2.0") ;
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector alpha = hypp.slot("alpha") ;
  
  NumericVector p_sigma2_0(1) ;
  NumericVector p_mu(1) ;
  NumericVector p_nu_0(1) ;
  NumericVector p_tau2(1) ;
  NumericVector logp_pmix(1) ;

  logp_pmix = log_ddirichlet_(pmix, alpha) ;
  p_tau2 = dgamma(1.0/tau2, 0.5*eta, 2.0/(eta * m2)) ;
  p_sigma2_0 = dgamma(sigma2_0, a, 1/b) ;
  p_nu_0 = dgeom(nu_0, betas) ;
  p_mu = dnorm(mu, mu_0, tau_0) ;
  NumericVector prior_prob(1) ;
  prior_prob = log(p_sigma2_0) + log(p_nu_0) + log(p_mu) + log(p_tau2) + logp_pmix ;
  return prior_prob ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_sigma2(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector x = model.slot("data") ;
  NumericVector theta = model.slot("theta") ;
  double nu_0 = model.slot("nu.0") ;
  double sigma2_0 = model.slot("sigma2.0") ;
  int n = x.size() ;
  IntegerVector z = model.slot("z") ;
  NumericVector nu_n(K) ;
  IntegerVector nn = model.slot("zfreq") ;
  
  for(int k = 0; k < K; ++k){
    nu_n[k] = nu_0 + nn[k] ;
  }
  
  //  return nn ;
  NumericVector ss(K) ;
  for(int i = 0; i < n; i++){
    int k = 0 ;
    while(k <= K) {
      if( z[i] == k + 1 ){
        ss[k] += pow(x[i] - theta[k], 2.0) ;
        break ;
      }
      k++ ;
    }
  }
  double sigma2_n ;
  NumericVector sigma2_new(K) ;
  for (int k = 0; k < K; k++){
    sigma2_n = 1.0/nu_n[k]*(nu_0*sigma2_0 + ss[k]) ;
    sigma2_new[k] = 1.0/as<double>(rgamma(1, 0.5*nu_n[k], 1.0/(0.5*nu_n[k]*sigma2_n))) ;
    //if(sigma2_new[k] < 0.001) sigma2_new[k] = 0.001 ;
  }
  return sigma2_new ;
}

// From stackoverflow http://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder

// [[Rcpp::export]]
Rcpp::IntegerVector ordertheta_(Rcpp::NumericVector x) {
  NumericVector sorted = clone(x).sort();
  //return match(sorted, x);
  return match(x, sorted) ;
}

//
// For the posterior probability of a copy number state, order the z
// labels by theta
//
// [[Rcpp::export]]
Rcpp::IntegerMatrix compute_probz(Rcpp::S4 xmod){
    Rcpp::S4 model(xmod);
    Rcpp::S4 hypp(model.slot("hyperparams"));
    int K = hypp.slot("k");
    IntegerVector z = model.slot("z");
    int N = z.size();
    IntegerMatrix pZ = model.slot("probz");
    NumericVector theta = model.slot("theta");
    NumericVector cn(K);
    NumericVector ordering(K);
    cn = ordertheta_(theta);

    for (int k = 0; k < K; ++k) {
        cn[k] = cn[k] - 1;
    }
    
    NumericVector is_z(N);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < K; ++k) {
            if (z[i] == (k + 1)) {
                pZ(i, cn[k]) += 1;
            }
        }
    }

    return pZ;
}
  


//
// BURNIN: Don't move chains, no thinning, no need to compute loglik at each iteration
//
// [[Rcpp::export]]
Rcpp::S4 mcmc_marginal_burnin(Rcpp::S4 xmod, Rcpp::S4 mcmcp) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;  
  Rcpp::S4 params(mcmcp) ;
  IntegerVector up = params.slot("param_updates") ;
  int S = params.slot("burnin") ;
  if( S < 1 ){
    return xmod ;
  }
  for(int s = 0; s < S; ++s){
    if(up[0] > 0)
      model.slot("theta") = update_theta(xmod) ;
    if(up[1] > 0)
      model.slot("sigma2") = update_sigma2(xmod) ;
    if(up[2] > 0)
      model.slot("pi") = update_p(xmod) ;
    if(up[3] > 0)
      model.slot("mu") = update_mu(xmod) ;
    if(up[4] > 0)    
      model.slot("tau2") = update_tau2(xmod) ;
    if(up[5] > 0)    
      model.slot("nu.0") = update_nu0(xmod) ;
    if(up[6] > 0)        
      model.slot("sigma2.0") = update_sigma2_0(xmod) ;
    if(up[7] > 0){        
      model.slot("z") = update_z(xmod) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    }
    model.slot("data.mean") = compute_means(xmod) ;
    model.slot("data.prec") = compute_prec(xmod) ;
  }
  // compute log prior probability from last iteration of burnin
  // compute log likelihood from last iteration of burnin
  NumericVector ll = loglik(xmod) ;
  NumericVector lls2 = stageTwoLogLik(xmod) ;
  model.slot("loglik") = ll + lls2 ;
  model.slot("logprior") = compute_logprior(xmod) ;    
  return xmod ;
}


// Function has gotten pretty long. Might be useful to have a separate
// function whose sole job is to move the chains.  
//
// [[Rcpp::export]]
Rcpp::S4 mcmc_marginal(Rcpp::S4 object, Rcpp::S4 mcmcp) {
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
  int N = x.size() ;
  NumericMatrix theta = chain.slot("theta") ;
  NumericMatrix sigma2 = chain.slot("sigma2") ;
  NumericMatrix pmix = chain.slot("pi") ;
  NumericMatrix zfreq = chain.slot("zfreq") ;
  IntegerMatrix Z = chain.slot("z") ;
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
  NumericVector s2(K) ;
  NumericVector p(K) ;
  NumericVector m(1) ; //mu
  NumericVector t2(1) ;//tau2
  NumericVector n0(1) ;//nu0
  IntegerVector z(N) ;
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
  sigma2(0, _) = s2 ;
  pmix(0, _) = p ;
  zfreq(0, _) = zf ;
  Z(0, _) = z ;
  // start at 1 instead of zero. Initial values are as above
  for(int s = 1; s < S; ++s){
    if(up[0] > 0) {
      try {
          th = update_theta(xmod) ;
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
      s2 = update_sigma2(xmod) ;
      model.slot("sigma2") = s2 ;
    } else { s2= model.slot("sigma2") ; }
    sigma2(s, _) = s2 ;
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
      n0 = update_nu0(xmod) ;
      model.slot("nu.0") = n0 ;
    } else {
      n0 = model.slot("nu.0") ;
    }
    nu0[s] = n0[0] ;
    if(up[6] > 0){        
      s20 = update_sigma2_0(xmod) ;
      model.slot("sigma2.0") = s20 ;
    } else {
      s20 = model.slot("sigma2.0") ;
    }
    sigma2_0[s] = s20[0] ;
    if(up[7] > 0){
      z = update_z(xmod) ;
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
    ll = loglik(xmod) ;
    lls2 = stageTwoLogLik(xmod) ;
    // ll = ll + lls2 ;
    loglik_[s] = ll[0] ;
    model.slot("loglik") = ll ;
    lp = compute_logprior(xmod) ;
    logprior_[s] = lp[0] ;
    model.slot("logprior") = lp ;
    // Thinning
    for(int t = 0; t < T; ++t){
      if(up[0] > 0)
        model.slot("theta") = update_theta(xmod) ;
      if(up[1] > 0)      
        model.slot("sigma2") = update_sigma2(xmod) ;
      if(up[2] > 0)
        model.slot("pi") = update_p(xmod) ;
      if(up[3] > 0)      
        model.slot("mu") = update_mu(xmod) ;
      if(up[4] > 0)      
        model.slot("tau2") = update_tau2(xmod) ;
      if(up[5] > 0)
        model.slot("nu.0") = update_nu0(xmod) ;
      if(up[6] > 0)
        model.slot("sigma2.0") = update_sigma2_0(xmod) ;
      if(up[7] > 0){
        model.slot("z") = update_z(xmod) ;
        model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      } 
      model.slot("data.mean") = compute_means(xmod) ;
      model.slot("data.prec") = compute_prec(xmod) ;
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
  model.slot("mcmc.chains") = chain ;
  return xmod ;
}


//
// Estimating marginals
//

// [[Rcpp::export]]
Rcpp::NumericVector marginal_theta(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains(model.slot("mcmc.chains")) ;  
  int S = params.slot("iter") ;
  // Assume current values are the modes (in R, useModes(object) ensures this)
  // List modes = model.slot("modes") ;
  // NumericVector thetastar = as<NumericVector>(modes["theta"]) ;
  List modes = model.slot("modes") ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector thetastar=clone(theta_) ;
  int K = thetastar.size() ;
  NumericVector p_theta(S) ;
  NumericVector muc = chains.slot("mu") ;
  NumericVector tau2c = chains.slot("tau2") ;
  NumericMatrix sigma2 = chains.slot("sigma2") ;
  NumericVector tauc = sqrt(tau2c) ;
  NumericVector tmp(K) ;

  IntegerMatrix Z = chains.slot("z") ;
  IntegerVector zz ;

  double tau2_tilde ;
  NumericVector sigma2_tilde(K) ;

  // this should be updated for each iteration
  NumericVector data_mean(K) ;
  IntegerVector nn(K) ;
  double post_prec;
  double tau_n;
  double mu_n;
  double w1;
  double w2;
  double prod ;

  for(int s=0; s < S; ++s){
    zz = Z(s, _) ;
    model.slot("z") = zz ;
    nn = tableZ(K, zz) ;
    data_mean = compute_means(model) ;
    tau2_tilde = 1/tau2c[s] ;
    sigma2_tilde = 1.0/sigma2(s, _) ;
    //tmp = dnorm(thetastar, muc[s], tauc[s]) ;
    double prod = 1.0;
    for(int k = 0; k < K; ++k) {
      post_prec = tau2_tilde + sigma2_tilde[k] * nn[k];
      tau_n = sqrt(1/post_prec);
      w1 = tau2_tilde/post_prec;
      w2 = nn[k]*sigma2_tilde[k]/post_prec;
      mu_n = w1*muc[s] + w2*data_mean[k];
      tmp = dnorm(thetastar, mu_n, tau_n) ;
      prod = prod * tmp[k] ;
    }
    p_theta[s] = prod ;
  }
  return p_theta ;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_theta_zpermuted(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 mcmcp = model.slot("mcmc.params") ;
  //Rcpp::S4 params(mcmcp) ;
  int S = mcmcp.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  int K = thetastar.size() ;
  NumericVector logp_theta(S) ;
  Rcpp::S4 chains(model.slot("mcmc.chains")) ;
  double mu ;
  NumericVector tau(1) ;
  NumericVector tmp(K) ;
  IntegerMatrix Z = chains.slot("z") ;
  int N = Z.ncol() ;
  IntegerVector h (N) ;
  NumericVector tau2(1) ;
  for(int s=0; s < S; ++s){
    h = Z(s, _ ) ;
    model.slot("z") = h ;
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    model.slot("theta") = update_theta(model) ;
    model.slot("sigma2") = update_sigma2(model) ;
    model.slot("pi") = update_p(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
    mu = model.slot("mu") ;
    tau2 = model.slot("tau2") ;
    tau = sqrt(tau2) ;
    tmp = dnorm(thetastar, mu, tau[0]) ;
    double prod = 0.0;
    for(int k = 0; k < K; ++k) {
      prod += log(tmp[k]) ;
    }
    logp_theta[s] = prod ;
  }
  return logp_theta ;
}

Rcpp::NumericVector marginal_sigma2(Rcpp::S4 xmod, Rcpp::S4 mcmcp) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;  
  Rcpp::S4 params(mcmcp) ;
  int S = params.slot("iter") ;
  // Assume current values are the modes (in R, useModes(object) ensures this)
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  //NumericVector sigma2_ = model.slot("sigma2") ;
  NumericVector sigma2star = clone(sigma2_) ;
  NumericVector thetastar = clone(theta_) ;
  NumericVector prec = pow(sigma2star, -1.0) ;
  int K = prec.size() ;
  NumericVector logp_prec(S) ;
  //
  // Run reduced Gibbs  -- theta is fixed at modal ordinate
  //
  Rcpp::S4 chains(model.slot("mcmc.chains")) ;
  NumericVector tmp(K) ;
  NumericVector nu0 = chains.slot("nu.0") ;
  NumericVector s20 = chains.slot("sigma2.0") ;
  for(int s=0; s < S; ++s){
    model.slot("z") = update_z(model) ;
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    // model.slot("theta") = update_theta(model) ;  Do not update theta!
    model.slot("sigma2") = update_sigma2(model) ;
    model.slot("pi") = update_p(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
    nu0 = model.slot("nu.0") ;
    s20 = model.slot("sigma2.0") ;
    tmp = dgamma(prec, 0.5*nu0[0], 2.0 / (nu0[0]*s20[0])) ;
    double total = 0.0 ;
    for(int k = 0; k < K; ++k){
      total += log(tmp[k]) ;
    }
    logp_prec[s] = total ;
  }
  return logp_prec ;
}

// [[Rcpp::export]]
Rcpp::S4 simulate_z_reduced1(Rcpp::S4 object) {
  RNGScope scope ;
  Rcpp::S4 model_(object) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector logp_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector y = model.slot("data") ;
  int N = y.size() ;
  NumericVector tau2(1) ;
  NumericVector nu0 (1) ;
  NumericVector s20 (1) ;
  NumericVector s2 (1) ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  NumericVector nu0chain = chains.slot("nu.0") ;
  NumericVector s20chain = chains.slot("sigma2.0") ;
  IntegerVector h(N) ;
  model.slot("theta") = thetastar ;
  //
  // Run reduced Gibbs  -- theta is fixed at modal ordinate
  //  
  for(int s=0; s < S; ++s){
    model.slot("z") = update_z(model) ;
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    //model.slot("theta") = update_theta(model) ; Do not update theta !
    model.slot("sigma2") = update_sigma2(model) ;
    model.slot("pi") = update_p(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
    nu0chain[s] = model.slot("nu.0") ;
    s20chain[s] = model.slot("sigma2.0") ;
    h = model.slot("z") ;
    Z(s, _) = h ;
  }
  //return logp_prec ;
  chains.slot("z") = Z ;
  chains.slot("nu.0") = nu0chain ;
  chains.slot("sigma2.0") = s20chain ;
  model.slot("mcmc.chains") = chains ;
  //return logp_prec ;
  return model ;
}

// [[Rcpp::export]]
Rcpp::S4 simulate_z_reduced2(Rcpp::S4 object) {
  RNGScope scope ;
  Rcpp::S4 model_(object) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector logp_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector y = model.slot("data") ;
  int N = y.size() ;
  NumericVector tau2(1) ;
  NumericVector nu0 (1) ;
  NumericVector s20 (1) ;
  NumericVector s2 (1) ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  NumericVector nu0chain = chains.slot("nu.0") ;
  NumericVector s20chain = chains.slot("sigma2.0") ;
  IntegerVector h(N) ;
  model.slot("theta") = thetastar ;
  model.slot("sigma2") = sigma2star ;  
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  for(int s=0; s < S; ++s){
    model.slot("z") = update_z(model) ;
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    // model.slot("theta") = update_theta(model) ; Do not update theta !
    // model.slot("sigma2") = update_sigma2(model) ;
    model.slot("pi") = update_p(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
    h = model.slot("z") ;
    Z(s, _) = h ;
  }
  //return logp_prec ;
  chains.slot("z") = Z ;
  // chains.slot("nu.0") = nu0chain ;
  // chains.slot("sigma2.0") = s20chain ;
  model.slot("mcmc.chains") = chains ;
  //return logp_prec ;
  return model ;
}

// [[Rcpp::export]]
Rcpp::S4 permutedz_reduced1(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector logp_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector y = model.slot("data") ;
  int N = y.size() ;
  NumericVector tau2(1) ;
  NumericVector nu0 (1) ;
  NumericVector s20 (1) ;
  NumericVector s2 (1) ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  NumericVector nu0chain = chains.slot("nu.0") ;
  NumericVector s20chain = chains.slot("sigma2.0") ;
  IntegerVector h(N) ;
  model.slot("theta") = thetastar ;
  //
  // Run reduced Gibbs  -- theta is fixed at modal ordinate
  //  
  for(int s=0; s < S; ++s){
    h = Z(s, _) ;
    //model.slot("z") = update_z(xmod) ;
    model.slot("z") = h ;
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    //model.slot("theta") = update_theta(model) ; Do not update theta !
    model.slot("sigma2") = update_sigma2(model) ;
    model.slot("pi") = update_p(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
    nu0chain[s] = model.slot("nu.0") ;
    s20chain[s] = model.slot("sigma2.0") ;
  }
  //return logp_prec ;
  // chains.slot("z") = Z ;
  chains.slot("nu.0") = nu0chain ;
  chains.slot("sigma2.0") = s20chain ;
  model.slot("mcmc.chains") = chains ;
  //return logp_prec ;
  return model ;
}


// [[Rcpp::export]]
Rcpp::S4 permutedz_reduced2(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector logp_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector y = model.slot("data") ;
  int N = y.size() ;
  NumericVector tau2(1) ;
  NumericVector nu0 (1) ;
  NumericVector s20 (1) ;
  NumericVector s2 (1) ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  NumericVector nu0chain = chains.slot("nu.0") ;
  NumericVector s20chain = chains.slot("sigma2.0") ;
  IntegerVector h(N) ;
  model.slot("theta") = thetastar ;
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  //  
  for(int s=0; s < S; ++s){
    h = Z(s, _) ;
    //model.slot("z") = update_z(xmod) ;
    model.slot("z") = h ;
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    // model.slot("theta") = update_theta(model) ; Do not update theta !
    // model.slot("sigma2") = update_sigma2(model) ;
    model.slot("pi") = update_p(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
    nu0chain[s] = model.slot("nu.0") ;
    s20chain[s] = model.slot("sigma2.0") ;
  }
  //return logp_prec ;
  // chains.slot("z") = Z ;
  chains.slot("nu.0") = nu0chain ;
  chains.slot("sigma2.0") = s20chain ;
  model.slot("mcmc.chains") = chains ;
  //return logp_prec ;
  return model ;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_pmix_reduced(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 mcmcp = model.slot("mcmc.params") ;
  Rcpp::S4 chains = model.slot("mcmc.chains") ;
  Rcpp::S4 hypp = model.slot("hyperparams") ;
  List modes = model.slot("modes") ;
  //
  //
  NumericVector x = model.slot("data") ;    
  int K = hypp.slot("k") ;
  int S = mcmcp.slot("iter") ;  
  int N = x.size() ;
  //
  NumericVector p_=as<NumericVector>(modes["mixprob"]) ;
  NumericVector pstar = clone(p_) ;
  NumericMatrix Z = chains.slot("z") ;
  NumericVector alpha = hypp.slot("alpha") ;
  NumericVector pmix(S) ;
  //
  // Run reduced Gibbs  -- theta,sigma2 fixed at modal ordinates
  //
  NumericVector h(N) ;
  NumericVector alpha_n(K) ;
  NumericVector tmp(1) ;
  for(int s=0; s < S; ++s){
    h = Z(s, _ ) ;    
    for(int k = 0 ; k < K; ++k){
      alpha_n[k] = sum(h == k+1) + alpha[k] ;
    }
    tmp = log_ddirichlet_(pstar, alpha_n) ;
    pmix[s] = exp(tmp[0]) ;
  }
  // return tmp ;
  return pmix ;
}


// [[Rcpp::export]]
Rcpp::S4 reduced_sigma(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector logp_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector y = model.slot("data") ;
  int N = y.size() ;
  NumericVector tau2(1) ;
  NumericVector nu0 (1) ;
  NumericVector s20 (1) ;
  NumericVector s2 (1) ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  NumericVector nu0chain = chains.slot("nu.0") ;
  NumericVector s20chain = chains.slot("sigma2.0") ;
  NumericVector muchain = chains.slot("mu") ;
  NumericVector tauchain = chains.slot("tau2") ;
  NumericMatrix pichain = chains.slot("pi") ;
  NumericMatrix sigmachain = chains.slot("sigma2") ;
  
  NumericVector sigma2 = model.slot("sigma2") ;
  NumericVector pi = model.slot("pi") ;
  NumericVector tau = model.slot("tau2") ;
  NumericVector mu = model.slot("mu") ;

  
  IntegerVector h(N) ;
  model.slot("theta") = thetastar ;
  IntegerVector zz ;
  //
  // Run reduced Gibbs  -- theta is fixed at modal ordinate
  //  
  for(int s=0; s < S; ++s){
    zz = update_z(model) ;
    model.slot("z") = zz ;
    Z(s, _) = zz ;
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    //model.slot("theta") = update_theta(model) ; Do not update theta !
    model.slot("sigma2") = update_sigma2(model) ;
    model.slot("pi") = update_p(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
    nu0chain[s] = model.slot("nu.0") ;
    s20chain[s] = model.slot("sigma2.0") ;
    // update the following chains for debugging small sigma2.0 values
    sigma2 = model.slot("sigma2") ;
    sigmachain(s, _) = sigma2 ;
    pi = model.slot("pi") ;
    pichain(s, _) = pi ;
    tau = model.slot("tau2") ;
    tauchain[s] = tau[0] ;
    mu = model.slot("mu") ;
    muchain[s] = mu[0] ;
  }
  //return logp_prec ;
  chains.slot("z") = Z ;
  chains.slot("nu.0") = nu0chain ;
  chains.slot("sigma2.0") = s20chain ;

  // update the following chains for debugging
  chains.slot("pi") = pichain ;
  chains.slot("sigma2") = sigmachain ;
  chains.slot("tau2") = tauchain ;
  chains.slot("mu") = muchain ;

  model.slot("mcmc.chains") = chains ;
  //return logp_prec ;
  return model ;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_sigma_reduced(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector p_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector nu0 (1) ;
  NumericVector s20 (1) ;
  IntegerMatrix Z = chains.slot("z") ;
  IntegerVector zz ;
  IntegerVector nn(K) ;
  NumericVector nu_n(K) ;
  //
  // We need to keep the Z|y,theta* chain
  //
  NumericVector nu0chain = chains.slot("nu.0") ;
  NumericVector s20chain = chains.slot("sigma2.0") ;
  //
  NumericVector sigma2_n(1) ;
  NumericVector sigma2_new(K) ;  
  //
  // Run reduced Gibbs  -- theta is fixed at modal ordinate
  //
  prec = 1.0/sigma2star ;
  for(int s=0; s < S; ++s){
    zz = Z(s, _) ;
    nn = tableZ(K, zz) ;
    s20 = s20chain[s] ;
    nu0 = nu0chain[s] ;

    NumericVector ss(K) ;
    for(int i = 0; i < n; i++){
      int k = 0 ;
      while(k <= K) {
        if( zz[i] == k + 1 ){
          ss[k] += pow(x[i] - thetastar[k], 2.0) ;
          break ;
        }
        k++ ;
      }
    }
    
    double total = 1.0 ;    
    for(int k = 0; k < K; ++k){
      nu_n = nu0 + nn[k] ;
      sigma2_n = 1.0/nu_n[0]*(nu0*s20 + ss[k]) ;
      tmp = dgamma(prec, 0.5*nu_n[0], 2.0 / (nu_n[0]*sigma2_n[0])) ;
      total = total*tmp[k] ;
    }
    p_prec[s] = total ;
  }
  return p_prec ;
}


// [[Rcpp::export]]
Rcpp::S4 reduced_pi(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector pi_ = as<NumericVector>(modes["mixprob"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector logp_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector y = model.slot("data") ;
  int N = y.size() ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  NumericVector zz(N) ;
  model.slot("theta") = thetastar ;
  model.slot("sigma2") = sigma2star ;
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  //  
  for(int s=0; s < S; ++s){
    zz = update_z(model) ;
    model.slot("z") = zz ;
    Z(s, _) = zz ;
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    // model.slot("theta") = update_theta(model) ; Do not update theta !
    // model.slot("sigma2") = update_sigma2(model) ;
    model.slot("pi") = update_p(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
  }
  chains.slot("z") = Z ;
  model.slot("mcmc.chains") = chains ;
  return model ;
}

// [[Rcpp::export]]
Rcpp::S4 reduced_mu(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector pi_ = as<NumericVector>(modes["mixprob"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  NumericVector pistar=clone(pi_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector logp_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector y = model.slot("data") ;
  NumericVector tau2chain(S) ;
  int N = y.size() ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  IntegerVector zz(N) ;
  model.slot("theta") = thetastar ;
  model.slot("sigma2") = sigma2star ;
  model.slot("pi") = pistar ;

  // keep chains for debugging
  NumericVector nu0chain = chains.slot("nu.0") ;
  NumericVector s20chain = chains.slot("sigma2.0") ;
  NumericVector muchain = chains.slot("mu") ;
  NumericVector tauchain = chains.slot("tau2") ;
  NumericVector mu = model.slot("mu") ;
  NumericVector nu0 = model.slot("nu.0") ;
  NumericVector s20 = model.slot("sigma2.0") ;
  NumericVector tau2 = model.slot("tau2") ;  
  
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  //  
  for(int s=0; s < S; ++s){
    zz = update_z(model) ;
    model.slot("z") = zz ;
    Z(s, _) = zz ;    
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    // model.slot("theta") = update_theta(model) ; Do not update theta !
    // model.slot("sigma2") = update_sigma2(model) ;
    // model.slot("pi") = update_p(model) ;
    mu = update_mu(model) ;
    tau2 = update_tau2(model) ;
    nu0 = update_nu0(model) ;
    s20 = update_sigma2_0(model) ;
    
    model.slot("mu") = mu ;
    model.slot("tau2") = tau2 ;
    model.slot("nu.0") = nu0 ;
    model.slot("sigma2.0") = s20 ;
    nu0chain[s] = nu0[0] ;
    s20chain[s] = s20[0] ;
    muchain[s] = mu[0] ;
    tau2chain[s] = tau2[0] ;
  }
  chains.slot("tau2") = tau2chain ;
  chains.slot("mu") = muchain ;
  chains.slot("nu.0") = nu0chain ;
  chains.slot("sigma2.0") = s20chain ;  
  chains.slot("z") = Z ;
  model.slot("mcmc.chains") = chains ;
  return model ;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_mu_reduced(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 mcmcp = model.slot("mcmc.params") ;
  Rcpp::S4 chains = model.slot("mcmc.chains") ;
  Rcpp::S4 hypp = model.slot("hyperparams") ;
  List modes = model.slot("modes") ;
  //
  //
  NumericVector x = model.slot("data") ;    
  int K = hypp.slot("k") ;
  int S = mcmcp.slot("iter") ;  
  int N = x.size() ;
  //
  NumericVector p_=as<NumericVector>(modes["mixprob"]) ;
  NumericVector theta_=as<NumericVector>(modes["theta"]) ;
  NumericVector mu_=as<NumericVector>(modes["mu"]) ;
  NumericVector pstar = clone(p_) ;
  NumericVector mustar = clone(mu_) ;
  NumericVector thetastar = clone(theta_) ;
  IntegerMatrix Z = chains.slot("z") ;
  NumericVector tau2chain = chains.slot("tau2") ;
  IntegerVector zz(N) ;

  IntegerVector nn(K) ;
  NumericVector mu0=hypp.slot("mu.0") ;
  double mu_0 = mu0[0] ;
  NumericVector tau2_0 = hypp.slot("tau2.0") ;
  double tau20_tilde = 1.0/tau2_0[0] ;

  double mu_k ;
  NumericVector tau2_tilde = 1.0/tau2chain ;
  double w1 ;
  double w2 ;
  double tau_k ;
  NumericVector p_mu(S) ;
  NumericVector tmp(1) ;
  for(int s = 0; s < S; ++s){
    zz = Z(s, _) ;
    nn = tableZ(K, zz) ;
    double total = 0.0 ;
    double thetabar = 0.0 ;    
    for(int k = 0; k < K; k++) total += nn[k] ;    
    for(int k = 0; k < K; k++) thetabar += nn[k] * thetastar[k] / total ;
    double post_prec = tau20_tilde + K*tau2_tilde[s] ;
    w1 = tau20_tilde/post_prec ;
    w2 = K*tau2_tilde[s]/post_prec ;
    mu_k =  w1*mu_0 +  w2*thetabar ;
    tau_k = sqrt(1.0/post_prec) ;
    tmp  = dnorm(mustar, mu_k, tau_k) ;
    p_mu[s] = tmp[0] ;
  }
  return p_mu ;
}

// [[Rcpp::export]]
Rcpp::S4 reduced_tau(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector pi_ = as<NumericVector>(modes["mixprob"]) ;
  NumericVector mu_ = as<NumericVector>(modes["mu"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  NumericVector pistar=clone(pi_) ;
  NumericVector mustar=clone(mu_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector logp_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector y = model.slot("data") ;
  //NumericVector s20chain(S) ;
  int N = y.size() ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  IntegerVector zz(N) ;
  model.slot("theta") = thetastar ;
  model.slot("sigma2") = sigma2star ;
  model.slot("pi") = pistar ;
  model.slot("mu") = mustar ;
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  //  
  for(int s=0; s < S; ++s){
    zz = update_z(model) ;
    model.slot("z") = zz ;
    Z(s, _) = zz ;    
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    // model.slot("theta") = update_theta(model) ; Do not update theta !
    // model.slot("sigma2") = update_sigma2(model) ;
    // model.slot("pi") = update_p(model) ;
    // model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
    //s20chain[s] = model.slot("sigma2.0") ;
  }
  chains.slot("z") = Z ;
  //chains.slot("sigma2.0") = s20chain ;
  model.slot("mcmc.chains") = chains ;
  return model ;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_tau_reduced(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 mcmcp = model.slot("mcmc.params") ;
  Rcpp::S4 chains = model.slot("mcmc.chains") ;
  Rcpp::S4 hypp = model.slot("hyperparams") ;
  List modes = model.slot("modes") ;
  //
  //
  NumericVector x = model.slot("data") ;    
  int K = hypp.slot("k") ;
  int S = mcmcp.slot("iter") ;  
  int N = x.size() ;
  //
  NumericVector p_=as<NumericVector>(modes["mixprob"]) ;
  NumericVector theta_=as<NumericVector>(modes["theta"]) ;
  NumericVector mu_=as<NumericVector>(modes["mu"]) ;
  NumericVector tau2_=as<NumericVector>(modes["tau2"]) ;
  NumericVector pstar = clone(p_) ;
  NumericVector mustar = clone(mu_) ;
  NumericVector tau2star = clone(tau2_) ;
  NumericVector thetastar = clone(theta_) ;
  IntegerMatrix Z = chains.slot("z") ;
  IntegerVector zz(N) ;
  NumericVector p_tau(S) ;

  double m2_0 = hypp.slot("m2.0") ;
  double eta_0 = hypp.slot("eta.0") ;
  double eta_k = eta_0 + K ;
  NumericVector s2_k(1) ;
  for(int k = 0; k < K; k++) s2_k[0] += pow(thetastar[k] - mustar[0], 2.0) ;
  NumericVector m2_k(1) ;
  m2_k[0] = 1/eta_k * (eta_0 * m2_0 + s2_k[0]) ;
  //  NumericVector tmp(1) ;
  //  for(int s = 0; s < S; ++s) {
  //    tmp = dgamma(1.0/tau2star, 0.5 * eta_k, 1.0/(0.5 * eta_k * m2_k[0])) ;
  //    p_tau[s] = tmp ;
  //  }
  p_tau = dgamma(1.0/tau2star, 0.5 * eta_k, 1.0/(0.5 * eta_k * m2_k[0])) ;
  return p_tau ;
}

// [[Rcpp::export]]
Rcpp::S4 reduced_nu0(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector pi_ = as<NumericVector>(modes["mixprob"]) ;
  NumericVector mu_ = as<NumericVector>(modes["mu"]) ;
  NumericVector tau2_ = as<NumericVector>(modes["tau2"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  NumericVector pistar=clone(pi_) ;
  NumericVector mustar=clone(mu_) ;
  NumericVector tau2star=clone(tau2_) ;
  int K = thetastar.size() ;
  NumericVector prec(K) ;
  NumericVector logp_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector y = model.slot("data") ;
  NumericVector s20chain(S) ;
  int N = y.size() ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  IntegerVector zz(N) ;
  model.slot("theta") = thetastar ;
  model.slot("sigma2") = sigma2star ;
  model.slot("pi") = pistar ;
  model.slot("mu") = mustar ;
  model.slot("tau2") = tau2star ;

  for(int s=0; s < S; ++s){
    zz = update_z(model) ;
    model.slot("z") = zz ;
    Z(s, _) = zz ;    
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    // model.slot("theta") = update_theta(model) ; Do not update theta !
    // model.slot("sigma2") = update_sigma2(model) ;
    // model.slot("pi") = update_p(model) ;
    // model.slot("mu") = update_mu(model) ;
    // model.slot("tau2") = update_tau2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
    s20chain[s] = model.slot("sigma2.0") ;
  }
  chains.slot("z") = Z ;
  chains.slot("sigma2.0") = s20chain ;
  model.slot("mcmc.chains") = chains ;
  return model ;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_nu0_reduced(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 mcmcp = model.slot("mcmc.params") ;
  Rcpp::S4 chains = model.slot("mcmc.chains") ;
  Rcpp::S4 hypp = model.slot("hyperparams") ;
  List modes = model.slot("modes") ;
  //
  //
  NumericVector x = model.slot("data") ;    
  int K = hypp.slot("k") ;
  int S = mcmcp.slot("iter") ;  
  int N = x.size() ;
  //
  NumericVector p_=as<NumericVector>(modes["mixprob"]) ;
  NumericVector theta_=as<NumericVector>(modes["theta"]) ;
  NumericVector mu_=as<NumericVector>(modes["mu"]) ;
  NumericVector tau2_=as<NumericVector>(modes["tau2"]) ;
  IntegerVector nu0_=as<IntegerVector>(modes["nu0"]) ;
  NumericVector sigma2_=as<NumericVector>(modes["sigma2"]) ;
  NumericVector pstar = clone(p_) ;
  NumericVector mustar = clone(mu_) ;
  NumericVector tau2star = clone(tau2_) ;
  NumericVector thetastar = clone(theta_) ;
  NumericVector sigma2star = clone(sigma2_) ;
  IntegerVector nu0=clone(nu0_) ;
  NumericVector p_nu0(S) ;
  int nu0star = nu0[0] ;

  NumericVector s20chain = chains.slot("sigma2.0") ;
  double betas = hypp.slot("beta") ;

  //
  // compute p(nu0*, ) from *normalized* probabilities
  //
  NumericVector d(100) ;  // 100 is the maximum allowed value for nu_0
  NumericVector lpnu0(100);
  double prec = 0.0 ;
  double lprec = 0.0 ;
  for(int k = 0; k < K; k++) prec += 1.0/sigma2star[k] ;
  for(int k = 0; k < K; k++) lprec += log(1.0/sigma2star[k]) ;
  d = seq_len(100) ;
  NumericVector y1(100) ;
  NumericVector y2(100) ;
  NumericVector y3(100) ;
  for(int s = 0; s < S; ++s) {
    y1 = K*(0.5*d*log(s20chain[s]*0.5*d) - lgamma(d*0.5)) ;
    y2 = (0.5*d - 1.0) * lprec ;
    y3 = d*(betas + 0.5*s20chain[s]*prec) ;
    lpnu0 =  y1 + y2 - y3 ;
    NumericVector prob(100) ;
    prob = exp(lpnu0) ; // - maxprob) ;
    prob = prob/sum(prob) ;  // this is now normalized
    p_nu0[s] = prob[nu0star] ;
  }
  return p_nu0 ;
}

// [[Rcpp::export]]
Rcpp::S4 reduced_s20(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericVector sigma2_ = as<NumericVector>(modes["sigma2"]) ;
  NumericVector theta_ = as<NumericVector>(modes["theta"]) ;
  NumericVector pi_ = as<NumericVector>(modes["mixprob"]) ;
  NumericVector mu_ = as<NumericVector>(modes["mu"]) ;
  NumericVector tau2_ = as<NumericVector>(modes["tau2"]) ;
  IntegerVector nu0_ = as<IntegerVector>(modes["nu0"]) ;
  NumericVector sigma2star=clone(sigma2_) ;
  NumericVector thetastar=clone(theta_) ;
  NumericVector pistar=clone(pi_) ;
  NumericVector mustar=clone(mu_) ;
  NumericVector tau2star=clone(tau2_) ;
  IntegerVector nu0star=clone(nu0_) ;
  int K = thetastar.size() ;
  NumericVector y = model.slot("data") ;
  int N = y.size() ;
  //
  // We need to keep the Z|y,theta* chain
  //
  IntegerMatrix Z = chains.slot("z") ;
  IntegerVector zz(N) ;
  model.slot("theta") = thetastar ;
  model.slot("sigma2") = sigma2star ;
  model.slot("pi") = pistar ;
  model.slot("mu") = mustar ;
  model.slot("tau2") = tau2star ;
  model.slot("nu.0") = nu0star ;

  for(int s=0; s < S; ++s){
    zz = update_z(model) ;
    model.slot("z") = zz ;
    Z(s, _) = zz ;    
    model.slot("data.mean") = compute_means(model) ;
    model.slot("data.prec") = compute_prec(model) ;
    // model.slot("theta") = update_theta(model) ; Do not update theta !
    // model.slot("sigma2") = update_sigma2(model) ;
    // model.slot("pi") = update_p(model) ;
    // model.slot("mu") = update_mu(model) ;
    // model.slot("tau2") = update_tau2(model) ;
    // model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma2_0(model) ;
  }
  chains.slot("z") = Z ;
  model.slot("mcmc.chains") = chains ;
  return model ;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_s20_reduced(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 mcmcp = model.slot("mcmc.params") ;
  Rcpp::S4 chains = model.slot("mcmc.chains") ;
  Rcpp::S4 hypp = model.slot("hyperparams") ;
  List modes = model.slot("modes") ;
  //
  //
  NumericVector x = model.slot("data") ;    
  int K = hypp.slot("k") ;
  int S = mcmcp.slot("iter") ;  
  int N = x.size() ;
  //
  NumericVector p_=as<NumericVector>(modes["mixprob"]) ;
  NumericVector theta_=as<NumericVector>(modes["theta"]) ;
  NumericVector mu_=as<NumericVector>(modes["mu"]) ;
  NumericVector tau2_=as<NumericVector>(modes["tau2"]) ;
  IntegerVector nu0_=as<IntegerVector>(modes["nu0"]) ;
  NumericVector sigma2_=as<NumericVector>(modes["sigma2"]) ;
  NumericVector s20_=as<NumericVector>(modes["sigma2.0"]) ;
  NumericVector pstar = clone(p_) ;
  NumericVector mustar = clone(mu_) ;
  NumericVector tau2star = clone(tau2_) ;
  NumericVector thetastar = clone(theta_) ;
  NumericVector sigma2star = clone(sigma2_) ;
  NumericVector s20star = clone(s20_) ;
  IntegerVector nu0=clone(nu0_) ;
  int nu0star = nu0[0] ;
  
  NumericVector p_s20(S) ;

  double a = hypp.slot("a") ;
  double b = hypp.slot("b") ;  
  double a_k = a + 0.5*K*nu0star ;
  double b_k ;

  for (int k=0; k < K; k++) {
    b_k += 0.5*nu0star/sigma2star[k];
  }
  b_k += b ;
  p_s20 = dgamma(s20star, a_k, 1.0/b_k) ;
  return p_s20 ;
}

