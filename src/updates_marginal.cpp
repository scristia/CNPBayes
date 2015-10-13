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
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector tau2_tilde = 1.0/tau2 ;
  double mu_0 = hypp.slot("mu.0") ;
  double K = getK(hypp) ;
  NumericVector theta = model.slot("theta") ;
  IntegerVector z = model.slot("z") ;
  IntegerVector nn = tableZ(K, z) ;
  double thetabar = 0.0;
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
    double b_k = 0.0;

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
  // return model.slot("z") ;  
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
    // set RNG
    Rcpp::RNGScope scope;
    
    // Get model/accessories
    Rcpp::S4 model(xmod);
    Rcpp::S4 hypp(model.slot("hyperparams"));

    // get hyperparameters
    double a = hypp.slot("a");
    double b = hypp.slot("b");
    double mu_0 = hypp.slot("mu.0");
    double betas = hypp.slot("beta");
    double eta = hypp.slot("eta.0");
    double m2 = hypp.slot("m2.0");
    Rcpp::NumericVector alpha = hypp.slot("alpha");
    double tau2_0 = hypp.slot("tau2.0");
    double tau_0 = sqrt(tau2_0);

    // get parameter estimates
    Rcpp::NumericVector mu = model.slot("mu");
    Rcpp::NumericVector sigma2_0 = model.slot("sigma2.0");
    Rcpp::NumericVector nu_0 = model.slot("nu.0");
    Rcpp::NumericVector pmix = model.slot("pi");
    Rcpp::NumericVector tau2 = model.slot("tau2");
    
    // calculate probabilities
    Rcpp::NumericVector logp_pmix = log_ddirichlet_(pmix, alpha);
    Rcpp::NumericVector p_tau2 = dgamma(1.0 / tau2, 0.5 * eta, 2.0 / (eta * m2));
    Rcpp::NumericVector p_sigma2_0 = dgamma(sigma2_0, a, 1 / b);
    Rcpp::NumericVector p_nu_0 = dgeom(nu_0, betas);
    Rcpp::NumericVector p_mu = dnorm(mu, mu_0, tau_0);

    Rcpp::NumericVector prior_prob = log(p_sigma2_0) + log(p_nu_0) + log(p_mu) + log(p_tau2) + logp_pmix;
    return prior_prob;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_sigma2(Rcpp::S4 xmod) {
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

    Rcpp::NumericVector nu_n(K);
    Rcpp::IntegerVector nn = model.slot("zfreq");
    
    for (int k = 0; k < K; ++k) {
        nu_n[k] = nu_0 + nn[k];
    }
    
    Rcpp::NumericVector ss(K);

    for(int i = 0; i < n; i++){
        int k = 0;

        while (k <= K) {
            if ( z[i] == k + 1 ) {
                ss[k] += pow(x[i] - theta[k], 2.0);
                break;
            }

            k++;
        }
    }

    double sigma2_n;
    NumericVector sigma2_new(K);

    for (int k = 0; k < K; k++) {
        sigma2_n = 1.0 / nu_n[k] * (nu_0 * sigma2_0 + ss[k]);
        sigma2_new[k] = 1.0 / Rcpp::as<double>(Rcpp::rgamma(1, 0.5 * nu_n[k], 1.0 / (0.5 * nu_n[k] * sigma2_n)));
    }

    return sigma2_new;
}

// From stackoverflow http://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder

Rcpp::IntegerVector ordertheta_(Rcpp::NumericVector x) {
  NumericVector sorted = clone(x).sort();
  //return match(sorted, x);
  return match(x, sorted) ;
}

//
// For the posterior probability of a copy number state, order the z
// labels by theta
//
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

