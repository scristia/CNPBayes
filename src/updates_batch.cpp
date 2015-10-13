#include "update.h" // getK
#include "miscfunctions.h" // for rdirichlet
#include <Rmath.h>
#include <Rcpp.h>

using namespace Rcpp ;

// [[Rcpp::export]]
Rcpp::IntegerVector uniqueBatch(Rcpp::IntegerVector x) {
  IntegerVector tmp = unique(x) ;
  IntegerVector b = clone(tmp) ;
  std::sort(b.begin(), b.end()) ;
  return b ;
}

Rcpp::NumericMatrix tableBatchZ(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  int K = getK(model.slot("hyperparams")) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  int B = ub.size() ;
  IntegerVector z = model.slot("z") ;
  NumericMatrix nn(B, K) ;
  for(int j = 0; j < B; ++j){
    for(int k = 0; k < K; k++){
      nn(j, k) = sum((z == (k+1)) & (batch == ub[j]));
    }
  }  
  return nn ;  
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_loglik_batch(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  int K = getK(model.slot("hyperparams")) ;
  NumericVector x = model.slot("data") ;
  int N = x.size() ;  
  NumericVector p = model.slot("pi") ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch);
  NumericMatrix theta = model.slot("theta") ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  IntegerVector batch_freq = model.slot("batchElements") ;
  NumericVector loglik_(1) ;
  NumericMatrix tabz = tableBatchZ(xmod) ;
  int B = ub.size() ;
  NumericMatrix P(B, K) ;
  NumericMatrix sigma(B, K) ;
  // component probabilities for each batch
  for(int b = 0; b < B; ++b){
    int rowsum = 0 ;
    for(int k = 0; k < K; ++k){
      rowsum += tabz(b, k) ;
      sigma(b, k) = sqrt(sigma2(b, k)) ;
    }
    for(int k = 0; k < K; ++k){
      P(b, k) = tabz(b, k)/rowsum ;
    }
  }
  NumericMatrix lik(N, K) ;
  // double y ;
  NumericVector tmp(1) ;
  NumericVector this_batch(N) ;
  for(int k = 0; k < K; ++k){
    NumericVector dens(N) ;
    for(int b = 0; b < B; ++b){
      this_batch = batch == ub[b] ;
      tmp = P(b, k) * dnorm(x, theta(b, k), sigma(b, k)) * this_batch ;
      dens = dens + tmp ;
    }
    lik(_, k) = dens ;
  }
  NumericVector marginal_prob(N) ;
  for(int i = 0; i < N; ++i){
    for(int k = 0; k < K; ++k) {
      marginal_prob[i] += lik(i, k) ;
    }
    loglik_[0] += log(marginal_prob[i]) ;
  }
  return loglik_ ;
}

Rcpp::NumericVector update_mu_batch(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  double tau2_0 = hypp.slot("tau2.0") ;
  double tau2_0_tilde = 1/tau2_0 ;
  double mu_0 = hypp.slot("mu.0") ;
  
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector tau2_tilde = 1/tau2 ;
  IntegerVector z = model.slot("z") ;
  NumericMatrix theta = model.slot("theta") ;
  IntegerVector nn = model.slot("zfreq") ;

  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  int B = ub.size() ;
  
  NumericVector tau2_B_tilde(K) ;;
  for(int k = 0; k < K; ++k) tau2_B_tilde[k] = tau2_0_tilde + B*tau2_tilde[k] ;
  
  NumericVector w1(K) ;
  NumericVector w2(K) ;
  for(int k = 0; k < K; ++k){
    w1[k] = tau2_0_tilde/(tau2_0_tilde + B*tau2_tilde[k]) ;
    w2[k] = B*tau2_tilde[k]/(tau2_0_tilde + B*tau2_tilde[k]) ;
  }
  NumericMatrix n_b = tableBatchZ(xmod) ;
  NumericVector theta_bar(K) ;
  NumericVector th(K) ;

  for(int k = 0; k < K; ++k){
    double n_k = 0.0 ; // number of observations for component k
    double colsumtheta = 0.0;
    for(int i = 0; i < B; ++i){
      colsumtheta += n_b(i, k)*theta(i, k) ;
      n_k += n_b(i, k) ;
    }
    theta_bar[k] = colsumtheta/n_k ;
  }
  NumericVector mu_n(K) ;
  NumericVector mu_new(K) ;
  double post_prec ;
  for(int k=0; k<K; ++k){
    post_prec = sqrt(1.0/tau2_B_tilde[k]) ;
    mu_n[k] = w1[k]*mu_0 + w2[k]*theta_bar[k] ;
    mu_new[k] = as<double>(rnorm(1, mu_n[k], post_prec)) ;
  }
  // simulate from prior if NAs
  LogicalVector isnan = is_nan(mu_new) ;
  if(!is_true(any(isnan)))
    return mu_new ;
  
  for(int k = 0; k < K; ++k){
    if(isnan[k])
      mu_new[k] = as<double>(rnorm(1, mu_0, sqrt(tau2_0))) ;
  }
  return mu_new ;
}

Rcpp::NumericVector update_tau2_batch(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  double m2_0 = hypp.slot("m2.0") ;
  int K = getK(hypp) ;
  double eta_0 = hypp.slot("eta.0") ;

  NumericVector mu = model.slot("mu") ;
  NumericMatrix theta = model.slot("theta") ;
  
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  int B = ub.size() ;
  double eta_B = eta_0 + B ;
  
  NumericVector s2_k(K) ;
  for(int k = 0; k < K; ++k){
    for(int i = 0; i < B; ++i){
      s2_k[k] += pow(theta(i, k) - mu[k], 2) ;
    }
  }
  NumericVector tau2(K) ;
  for(int k = 0; k < K; ++k) {
    double m2_k = 0.0 ;
    m2_k = 1.0/eta_B*(eta_0*m2_0 + s2_k[k]) ;
    tau2[k] = 1.0/as<double>(rgamma(1, 0.5*eta_B, 2.0/(eta_B*m2_k))) ;
  }
  return tau2 ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_sigma20_batch(Rcpp::S4 xmod){
    RNGScope scope ;
    Rcpp::S4 model(xmod) ;
    Rcpp::S4 hypp(model.slot("hyperparams")) ;
    int K = getK(hypp) ;
    IntegerVector batch = model.slot("batch") ;
    IntegerVector ub = uniqueBatch(batch) ;
    int B = ub.size() ;
    NumericVector a = hypp.slot("a") ;
    NumericVector b = hypp.slot("b") ;
    NumericVector nu_0 = model.slot("nu.0") ;  
    NumericMatrix sigma2 = model.slot("sigma2") ;
    Rcpp::NumericVector sigma2_0_old = model.slot("sigma2.0");
    NumericVector prec(1) ;

    for(int i = 0; i < B; ++i){
        for(int k = 0; k < K; ++k){
            prec[0] += 1.0/sigma2(i, k) ;
        }
    }

    NumericVector a_k(1) ;
    NumericVector b_k(1) ;
    a_k[0] = a[0] + 0.5*(K * B)*nu_0[0] ;
    b_k[0] = b[0] + 0.5*nu_0[0]*prec[0] ;
    double rate ;
    rate = 1.0/b_k[0] ;
    //return b_k ;
    NumericVector sigma2_0(1) ;
    sigma2_0[0] = as<double>(rgamma(1, a_k[0], rate)) ;
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

Rcpp::NumericVector update_nu0_batch(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  int B = sigma2.nrow() ;  
  double sigma2_0 = model.slot("sigma2.0") ;  
  double prec = 0.0;
  double lprec = 0.0 ;
  double betas = hypp.slot("beta") ;
  for(int i = 0; i < B; ++i){
    for(int k = 0; k < K; ++k){
      prec += 1.0/sigma2(i, k) ;
      lprec += log(1.0/sigma2(i, k)) ;
    }
  }
  NumericVector x(100) ;  
  for(int i = 0; i < 100; i++)  x[i] = i+1 ;
  NumericVector lpnu0(100);
  NumericVector y1(100) ;
  NumericVector y2(100) ;
  NumericVector y3(100) ;
  NumericVector prob(100) ;
  y1 = (B*K)*(0.5*x*log(sigma2_0*0.5*x) - lgamma(x*0.5)) ;
  y2 = (0.5*x - 1.0) * lprec ;
  y3 = x*(betas + 0.5*sigma2_0*prec) ;
  lpnu0 =  y1 + y2 - y3 ;
  prob = exp(lpnu0) ; // - maxprob) ;
  prob = prob/sum(prob) ;
  NumericVector nu0(1) ;
  //int u ;
  NumericVector u(1) ;
  double cumprob = 0.0;
  // sample x with probability prob
  for(int i = 0; i < 100; i++){
    cumprob += prob[i] ;
    u = runif(1) ;
    if (u[0] < cumprob){
      nu0[0] = x[i] ;
      break ;
    }
  }
  return nu0 ;  
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_multinomialPr_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  NumericVector p = model.slot("pi") ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  NumericMatrix theta = model.slot("theta") ;
  int B = sigma2.nrow() ;
  NumericVector x = model.slot("data") ;
  IntegerVector nb = model.slot("batchElements") ;
  int N = x.size() ;
  NumericMatrix lik(N, K) ;
  NumericVector this_batch(N) ;
  NumericVector tmp(N) ;
  NumericVector rowtotal(N) ;
  for(int k = 0; k < K; ++k){
    NumericVector dens(N) ;
    for(int b = 0; b < B; ++b){
      this_batch = batch == ub[b] ;
      tmp = p[k] * dnorm(x, theta(b, k), sqrt(sigma2(b, k))) * this_batch ;
      // if(is_true(any(tmp < 1e-10))) tmp[tmp < 1e-10] = 1e-10 ;
      dens += tmp ;
    }
    lik(_, k) = dens ;
    rowtotal += dens ;
  }
  //return lik ;
  NumericMatrix P(N, K) ;
  for(int k=0; k<K; ++k){
    P(_, k) = lik(_, k)/rowtotal ;
  }
  return P ;
}

//
// Note, this is currently the same as the marginal model
//
Rcpp::NumericVector update_p_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;  
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector z = model.slot("z") ;  
  IntegerVector nn = model.slot("zfreq");
  IntegerVector alpha = hypp.slot("alpha") ;
  NumericVector alpha_n(K) ;  // really an integer vector, but rdirichlet expects numeric
  for(int k=0; k < K; k++) alpha_n[k] = alpha[k] + nn[k] ;
  NumericVector p(K) ;
  // pass by reference
  rdirichlet(alpha_n, p) ;
  return p ;
  //  return alpha_n ;
}

// [[Rcpp::export]]
Rcpp::IntegerVector update_z_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;  
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector x = model.slot("data") ;
  NumericMatrix theta = model.slot("theta") ;
  IntegerVector batch = model.slot("batch") ;
  int B = theta.nrow() ;
  int n = x.size() ;
  NumericMatrix p(n, K);
  p = update_multinomialPr_batch(xmod) ;
  //NumericMatrix cumP(n, K) ;
  //  Make more efficient?
  //return cumP ;
  NumericVector u = runif(n) ;
  IntegerVector zz_(n) ;
  IntegerVector zz = clone(zz_) ;
  IntegerMatrix freq(B, K) ;
  int b ;  
  for(int i=0; i < n; i++){
    //initialize accumulator ;
    double acc = 0 ;
    for(int k = 0; k < K; k++){
      acc += p(i, k) ;
      if( u[i] < acc ) {
        zz[i] = k + 1 ;
        b = batch[i] - 1 ;
        freq(b, k) += 1 ;
        break ;
      }
    }
  }
  if(is_true(all(freq > 1))){
    return zz ;
  }
  //
  // Don't update z if there are states with zero frequency.  
  //
  return model.slot("z") ;
  //
  // Code below is no longer evaluated
  //
  IntegerVector index = seq_len(n) ;
  NumericVector r(2) ;
  LogicalVector is_batch (n) ;
  IntegerVector S = model.slot("batchElements") ;
  int w ;
  int start = 0 ;
  for(int k = 0; k < K; ++k){
    //
    // for zero frequency states, switch label of a random sample
    //
    for(int b = 0; b < B; ++b){
      if(freq(b, k) > 1) continue ;
      if(b > 0) start = S[b-1] ;
      IntegerVector J(S[b]) ;
      //int b = 1 ;
      is_batch = batch == (b + 1) ;
      J = zz [ is_batch ] ;
      r = runif(2, 0, 1) * (S[b] - 1) ;
      for(int j = 0; j < 2; ++j){
        w = (int) r[j] ;
        J[w] = (k + 1) ;
        zz[is_batch] = J ;
      }
    }
  }
  return zz ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_means_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  NumericVector mu = model.slot("mu") ;
  int n = x.size() ;
  IntegerVector z = model.slot("z") ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;  
  IntegerVector nn = model.slot("zfreq") ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  int B = ub.size() ;
  NumericMatrix means(B, K) ;
  NumericVector is_z(n) ;
  NumericVector this_batch(n) ;
  IntegerVector total(1) ;
  for(int b = 0; b < B; ++b){
    this_batch = batch == ub[b] ;
    for(int k = 0; k < K; ++k){
      is_z = z == (k + 1) ;
      total[0] = sum(is_z * this_batch) ;
      if(total[0] == 0){
        means(b, k) = mu[k] ;
      } else {
        means(b, k) = sum(x * this_batch * is_z)/total[0] ;
      }
    }
  }
  return means ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_vars_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  IntegerVector z = model.slot("z") ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;  
  IntegerVector nn = model.slot("zfreq") ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  int B = ub.size() ;
  NumericMatrix vars(B, K) ;
  NumericMatrix tabz = tableBatchZ(xmod) ;
  NumericMatrix mn = model.slot("data.mean") ;
  
  NumericVector this_batch(n) ;
  NumericVector is_z(n) ;
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector total(1) ;
  for(int b = 0; b < B; ++b){
    this_batch = batch == ub[b] ;
    for(int k = 0; k < K; ++k){
      is_z = z == (k + 1) ;
      total[0] = sum(is_z * this_batch) ;
      if(total[0] <= 1){
        vars(b, k) = tau2[k] ;
      } else {
        vars(b, k) = sum(pow(x - mn(b,k), 2.0) * this_batch * is_z) / (tabz(b,k)-1) ;
      }
    }
  }
  return vars ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_prec_batch(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericMatrix vars = compute_vars_batch(xmod) ;
  int B = vars.nrow() ;
  int K = vars.ncol() ;
  NumericMatrix prec(B, K) ;
  for(int b = 0; b < B; ++b){
    for(int k = 0; k < K; ++k){
      prec(b, k) = 1.0/vars(b, k) ;
    }
  }
  return prec ;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_logprior_batch(Rcpp::S4 xmod) {
    // set RNG
    RNGScope scope;
    
    // Get model/accessories
    Rcpp::S4 model(xmod);
    Rcpp::S4 hypp(model.slot("hyperparams"));

    // get hyperparameters
    double mu_0 = hypp.slot("mu.0");
    double tau2_0 = hypp.slot("tau2.0");
    double a = hypp.slot("a");
    double b = hypp.slot("b");
    double beta = hypp.slot("beta");

    // get parameter estimates
    Rcpp::NumericVector tau2 = model.slot("tau2");
    Rcpp::NumericVector mu = model.slot("mu");
    Rcpp::NumericVector sigma2_0 = model.slot("sigma2.0");
    Rcpp::IntegerVector nu0 = model.slot("nu.0");

    Rcpp::NumericVector p_mu = dnorm(mu, mu_0, sqrt(tau2_0));
    Rcpp::NumericVector p_sigma2_0 = dgamma(sigma2_0, a, 1.0/b);
    Rcpp::NumericVector p_nu0 = dgeom(nu0, beta);
    Rcpp::NumericVector logprior = sum(log(p_mu)) + log(p_sigma2_0) + log(p_nu0);

    return logprior;
}

// [[Rcpp::export]]
Rcpp::NumericVector stageTwoLogLikBatch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericMatrix theta = model.slot("theta") ;
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector mu = model.slot("mu") ;
  NumericVector nu0 = model.slot("nu.0") ;
  NumericVector s20 = model.slot("sigma2.0") ;
  NumericMatrix sigma2=model.slot("sigma2") ;
  NumericVector sigma2_tilde(1) ; // = 1.0/sigma2 ;
  NumericVector loglik(1) ;
  NumericVector tau = sqrt(tau2) ;
  int K = theta.ncol() ;  
  NumericVector liknorm(1) ;
  NumericVector likprec(1) ;
  NumericVector th(1) ;
  NumericVector LL(1) ;
  LL[0] = 0.0 ;
  for(int k = 0; k < K; ++k){
    liknorm[0] = sum(log(dnorm(theta(_, k), mu[k], tau[k]))) ;
    likprec[0] = sum(log(dgamma(1.0/sigma2(_, k), 0.5*nu0[0], 1.0/(0.5*nu0[0]*s20[0])))) ;
    LL[0] = LL[0] + liknorm[0] + likprec[0] ;
  }
  return LL ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix update_theta_batch(Rcpp::S4 xmod){
    RNGScope scope ;
    Rcpp::S4 model(xmod) ;
    Rcpp::S4 hypp(model.slot("hyperparams")) ;
    int K = getK(hypp) ;
    NumericVector x = model.slot("data") ;
    // NumericMatrix theta = model.slot("theta") ;
    NumericVector tau2 = model.slot("tau2") ;
    NumericMatrix sigma2 = model.slot("sigma2") ;
    NumericMatrix n_hb = tableBatchZ(xmod) ;
    NumericVector mu = model.slot("mu") ;
    int B = n_hb.nrow() ;
    NumericMatrix ybar = model.slot("data.mean") ;
    NumericMatrix theta_new(B, K) ;
    double w1 = 0.0 ;
    double w2 = 0.0 ;
    double mu_n = 0.0 ;
    double tau_n = 0.0 ;
    double post_prec = 0.0 ;

    for (int b = 0; b < B; ++b) {
        for(int k = 0; k < K; ++k){
            post_prec = 1.0/tau2[k] + n_hb(b, k)*1.0/sigma2(b, k) ;
            if (post_prec == R_PosInf) {
                throw std::runtime_error("Bad simulation. Run again with different start.");
            }
            tau_n = sqrt(1.0/post_prec) ;
            w1 = (1.0/tau2[k])/post_prec ;
            w2 = (n_hb(b, k) * 1.0/sigma2(b, k))/post_prec ;
            mu_n = w1*mu[k] + w2*ybar(b, k) ;
            theta_new(b, k) = as<double>(rnorm(1, mu_n, tau_n)) ;
        }
    }
    return theta_new ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_sigma2_batch(Rcpp::S4 xmod){
    Rcpp::RNGScope scope;

    // get model
    Rcpp::S4 model(xmod);

    // get parameter estimates
    Rcpp::NumericMatrix theta = model.slot("theta");
    Rcpp::IntegerVector z = model.slot("z");
    double nu_0 = model.slot("nu.0");
    double sigma2_0 = model.slot("sigma2.0");

    // get data and size attributes
    Rcpp::NumericVector x = model.slot("data");
    int n = x.size();
    int K = theta.ncol();
    int B = theta.nrow();

    //IntegerVector nn = model.slot("zfreq");
    // get batch info
    Rcpp::NumericMatrix tabz = tableBatchZ(xmod);
    Rcpp::IntegerVector batch = model.slot("batch");
    Rcpp::IntegerVector ub = uniqueBatch(batch);
    Rcpp::NumericMatrix ss(B, K);

    for (int i = 0; i < n; ++i) {
        for (int b = 0; b < B; ++b) {
            if (batch[i] != ub[b]) {
                continue;
            }

            for (int k = 0; k < K; ++k){
                if (z[i] == k+1){
                    ss(b, k) += pow(x[i] - theta(b, k), 2);
                }
            }
        }
    }

    //NumericMatrix sigma2_nh(B, K);
    double shape;
    double rate;
    double sigma2_nh;
    double nu_n;
    Rcpp::NumericMatrix sigma2_tilde(B, K);
    Rcpp::NumericMatrix sigma2_(B, K);

    for (int b = 0; b < B; ++b) {
        for (int k = 0; k < K; ++k) {
            nu_n = nu_0 + tabz(b, k);
            sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss(b, k));
            shape = 0.5 * nu_n;
            rate = shape * sigma2_nh;
            sigma2_tilde(b, k) = Rcpp::as<double>(rgamma(1, shape, 1.0/rate));
            sigma2_(b, k) = 1.0 / sigma2_tilde(b, k);
        }
    }

    return sigma2_;
}

// From stackoverflow http://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder

Rcpp::IntegerVector order_(Rcpp::NumericVector x) {
  NumericVector sorted = clone(x).sort();
  return match(x, sorted);
}

Rcpp::IntegerMatrix update_probz_batch(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = hypp.slot("k") ;
  IntegerVector z = model.slot("z") ;
  NumericMatrix theta = model.slot("theta") ;
  int N = z.size() ;
  IntegerMatrix pZ = model.slot("probz") ;
  //
  // update probz such that the z value corresponding to the lowest
  // mean is 1, the second lowest mean is 2, etc.
  //
  // Assume that all batches have the same ordering and so here we
  // just use the first batch
  //
  NumericVector means(K) ;
  means = theta(0, _) ;
  NumericVector cn(K) ;
  cn = order_(means) ;
  for(int k = 0; k < K; ++k) cn[k] = cn[k] - 1 ;
  
  for(int i = 0; i < N; ++i){
    for(int k = 0; k < K; ++k){
      if(z[i] == (k + 1)){
        pZ(i, cn[k]) += 1;
      }
    }
  }
  return pZ ;
}



// [[Rcpp::export]]
Rcpp::S4 mcmc_batch_burnin(Rcpp::S4 xmod, Rcpp::S4 mcmcp) {
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
      try {
        model.slot("theta") = update_theta_batch(xmod) ;
      } catch(std::runtime_error &ex) {
          forward_exception_to_r(ex);
      } catch(...) {
          ::Rf_error("c++ exception (unknown reason)");
      }
    if(up[1] > 0)
      model.slot("sigma2") = update_sigma2_batch(xmod) ;
    if(up[3] > 0)
      model.slot("mu") = update_mu_batch(xmod) ;    
    if(up[4] > 0)    
      model.slot("tau2") = update_tau2_batch(xmod) ;
    if(up[6] > 0)        
      model.slot("sigma2.0") = update_sigma20_batch(xmod) ;    
    if(up[5] > 0)    
      model.slot("nu.0") = update_nu0_batch(xmod) ;
    if(up[2] > 0)
      model.slot("pi") = update_p_batch(xmod) ;
    if(up[7] > 0){        
      model.slot("z") = update_z_batch(xmod) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    }
    model.slot("data.mean") = compute_means_batch(xmod) ;
    model.slot("data.prec") = compute_prec_batch(xmod) ;
  }
  // compute log prior probability from last iteration of burnin
  // compute log likelihood from last iteration of burnin
  model.slot("loglik") = compute_loglik_batch(xmod) ;
  model.slot("logprior") = compute_logprior_batch(xmod) ;    
  return xmod ;
  // return vars ;
}

// [[Rcpp::export]]
Rcpp::S4 mcmc_batch(Rcpp::S4 object, Rcpp::S4 mcmcp) {
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
  NumericMatrix mu = chain.slot("mu") ;  
  NumericMatrix tau2 = chain.slot("tau2") ;
  NumericVector nu0 = chain.slot("nu.0") ;
  NumericVector sigma2_0 = chain.slot("sigma2.0") ;
  NumericVector loglik_ = chain.slot("loglik") ;
  NumericVector logprior_ = chain.slot("logprior") ;
  NumericVector th(K) ;
  NumericVector s2(K) ;
  NumericVector p(K) ;
  NumericVector m(K) ; //mu
  NumericVector t2(K) ;//tau2
  NumericVector n0(1) ;//nu0
  IntegerVector z(N) ;
  NumericVector s20(1) ; //sigma2_0
  //  NumericVector mns(1) ;   
  // NumericVector precs(1) ;
  NumericVector ll(1) ;
  NumericVector lp(1) ;
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
  mu(0, _) = m ;
  nu0[0] = n0[0] ;
  sigma2_0[0] = s20[0] ;
  loglik_[0] = ll[0] ;
  logprior_[0] = lp[0] ;
  theta(0, _) = th ;
  tau2(0, _) = t2 ;
  sigma2(0, _) = s2 ;
  pmix(0, _) = p ;
  zfreq(0, _) = zf ;
  Z(0, _) = z ;

  // Is accessing a slot in an object expensive?
  
  // Currently, there is no alternative as the current values are
  // stored in the object.  Hence, the entire object has to be passed
  // to the updating functions.
  
  // start at 1 instead of zero. Initial values are as above
  for(int s = 1; s < S; ++s){
    if(up[0] > 0) {
      th = as<Rcpp::NumericVector>(update_theta_batch(xmod));
      model.slot("theta") = th ;
    } else {
      th = model.slot("theta") ;
    }
    theta(s, _) = th ;
    if(up[1] > 0){
      s2 = as<Rcpp::NumericVector>(update_sigma2_batch(xmod));
      model.slot("sigma2") = s2 ;
    } else {
      s2 = model.slot("sigma2") ;
    }
    sigma2(s, _) = s2 ;
    if(up[2] > 0){
      p = update_p_batch(xmod) ;
      model.slot("pi") = p ;      
    } else {
      p = model.slot("pi") ;
    }
    pmix(s, _) = p ;
    if(up[3] > 0){
      m = update_mu_batch(xmod) ;
      model.slot("mu") = m ;
    } else {
      m = model.slot("mu") ;
    }
    mu(s, _) = m ;
    if(up[4] > 0){    
      t2 = update_tau2_batch(xmod) ;
      model.slot("tau2") = t2 ;
    } else {
      t2 = model.slot("tau2") ;
    }
    tau2(s, _) = t2 ;
    if(up[5] > 0){        
      n0 = update_nu0_batch(xmod) ;
      model.slot("nu.0") = n0 ;
    } else {
      n0 = model.slot("nu.0") ;
    }
    nu0[s] = n0[0] ;
    if(up[6] > 0){        
      s20 = update_sigma20_batch(xmod) ;
      model.slot("sigma2.0") = s20 ;
    } else {
      s20 = model.slot("sigma2.0") ;
    }
    sigma2_0[s] = s20[0] ;
    if(up[7] > 0){
      z = update_z_batch(xmod) ;
      model.slot("z") = z ;
      tmp = tableZ(K, z) ;
      model.slot("probz") = update_probz_batch(xmod) ;    
      model.slot("zfreq") = tmp ;
      // mean and prec only change if z is updated
      model.slot("data.mean") = compute_means_batch(xmod) ;
      model.slot("data.prec") = compute_prec_batch(xmod) ;      
    } else {
      tmp = model.slot("zfreq") ;
    }
    Z(s, _) = z ;
    zfreq(s, _) = tmp ;    
    ll = compute_loglik_batch(xmod) ;
    loglik_[s] = ll[0] ;
    model.slot("loglik") = ll ;
    lp = compute_logprior_batch(xmod) ;
    logprior_[s] = lp[0] ;
    model.slot("logprior") = lp ;
    // Thinning
    for(int t = 0; t < T; ++t){
      if(up[0] > 0)
        model.slot("theta") = update_theta_batch(xmod) ;
      if(up[1] > 0)      
        model.slot("sigma2") = update_sigma2_batch(xmod) ;
      if(up[2] > 0)
        model.slot("pi") = update_p_batch(xmod) ;
      if(up[3] > 0)      
        model.slot("mu") = update_mu_batch(xmod) ;
      if(up[4] > 0)      
        model.slot("tau2") = update_tau2_batch(xmod) ;
      if(up[5] > 0)
        model.slot("nu.0") = update_nu0_batch(xmod) ;
     if(up[6] > 0)
        model.slot("sigma2.0") = update_sigma20_batch(xmod) ;
      if(up[7] > 0){
        model.slot("z") = update_z_batch(xmod) ;
        model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      }
      model.slot("data.mean") = compute_means_batch(xmod) ;
      model.slot("data.prec") = compute_prec_batch(xmod) ;
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
  chain.slot("z") = Z ;
  chain.slot("loglik") = loglik_ ;
  chain.slot("logprior") = logprior_ ;
  model.slot("mcmc.chains") = chain ;
  return xmod ;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_theta_batch(Rcpp::S4 xmod) {
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
  NumericMatrix theta_ = as<NumericMatrix>(modes["theta"]) ;
  NumericMatrix thetastar=clone(theta_) ;
  int K = thetastar.ncol() ;
  int B = thetastar.nrow() ;
  NumericVector p_theta(S) ;
  NumericMatrix muc = chains.slot("mu") ;
  NumericMatrix tau2c = chains.slot("tau2") ;
  NumericVector theta(1) ;
  NumericVector tmp(1) ;
  NumericVector d(1) ;
  double tau ;
  double mu ;  
  for(int s = 0; s < S; ++s){
    d = 1.0 ;
    for(int k = 0; k < K; ++k){
      tau = sqrt(tau2c(s, k)) ;
      mu = muc(s, k) ;
      for(int b = 0; b < B; ++b){
        theta[0] = thetastar(b, k) ;
        tmp = dnorm(theta, mu, tau) ;      
        d = d * tmp[0] ;
      }
    }
    p_theta[s] = d[0] ;
  }
  return p_theta ;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_theta_zfixed_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 mcmcp = model.slot("mcmc.params") ;
  //Rcpp::S4 params(mcmcp) ;
  int S = mcmcp.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericMatrix sigma2_ = as<NumericMatrix>(modes["sigma2"]) ;
  NumericMatrix theta_ = as<NumericMatrix>(modes["theta"]) ;
  NumericMatrix sigma2star=clone(sigma2_) ;
  NumericMatrix thetastar=clone(theta_) ;
  int K = thetastar.ncol() ;
  int B = thetastar.nrow() ;
  NumericVector p_theta(S) ;
  Rcpp::S4 chains(model.slot("mcmc.chains")) ;
  //double mu ;
  NumericVector tau(K) ;
  NumericVector tmp(K) ;
  IntegerMatrix Z = chains.slot("z") ;
  int N = Z.ncol() ;
  IntegerVector h (N) ;
  NumericVector tau2(K) ;
  NumericVector mu(K) ;
  NumericVector theta(1) ;
  NumericVector d(1) ;
  for(int s=0; s < S; ++s){
    h = Z(s, _ ) ;
    model.slot("z") = h ;
    model.slot("data.mean") = compute_means_batch(model) ;
    model.slot("data.prec") = compute_prec_batch(model) ;
    model.slot("theta") = update_theta_batch(model) ;
    model.slot("sigma2") = update_sigma2_batch(model) ;
    model.slot("pi") = update_p_batch(model) ;
    model.slot("mu") = update_mu_batch(model) ;
    model.slot("tau2") = update_tau2_batch(model) ;
    model.slot("nu.0") = update_nu0_batch(model) ;
    model.slot("sigma2.0") = update_sigma20_batch(model) ;
    mu = model.slot("mu") ;
    tau2 = model.slot("tau2") ;
    tau = sqrt(tau2) ;
    // tmp = dnorm(thetastar, mu, tau[0]) ;
    // double prod = 0.0;
    // for(int k = 0; k < K; ++k) {
    // prod += log(tmp[k]) ;
    // }
    d[0] = 1.0 ;
    for(int k = 0; k < K; ++k){
      for(int b = 0; b < B; ++b){
        theta[0] = thetastar(b, k) ;
        tmp = dnorm(theta, mu[k], tau[k]) ;      
        d = d * tmp[0] ;
      }
    }
    p_theta[s] = d[0] ;
    //logp_theta[s] = prod ;
  }
  return p_theta ;
}

// [[Rcpp::export]]
Rcpp::S4 reduced_z_theta_fixed(Rcpp::S4 object) {
  RNGScope scope ;
  Rcpp::S4 model_(object) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericMatrix theta_ = as<NumericMatrix>(modes["theta"]) ;
  NumericMatrix thetastar=clone(theta_) ;
  NumericVector y = model.slot("data") ;
  int N = y.size() ;
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
    model.slot("z") = h ;
    model.slot("data.mean") = compute_means_batch(model) ;
    model.slot("data.prec") = compute_prec_batch(model) ;
    //model.slot("theta") = update_theta(model) ; Do not update theta !
    model.slot("sigma2") = update_sigma2_batch(model) ;
    model.slot("pi") = update_p_batch(model) ;
    model.slot("mu") = update_mu_batch(model) ;
    model.slot("tau2") = update_tau2_batch(model) ;
    model.slot("nu.0") = update_nu0_batch(model) ;
    model.slot("sigma2.0") = update_sigma20_batch(model) ;
    nu0chain[s] = model.slot("nu.0") ;
    s20chain[s] = model.slot("sigma2.0") ;
  }
  //return logp_prec ;
  chains.slot("z") = Z ;
  chains.slot("nu.0") = nu0chain ;
  chains.slot("sigma2.0") = s20chain ;
  model.slot("mcmc.chains") = chains ;
  return model ;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_sigma2_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericMatrix sigma2_ = as<NumericMatrix>(modes["sigma2"]) ;
  NumericMatrix theta_ = as<NumericMatrix>(modes["theta"]) ;
  NumericMatrix sigma2star=clone(sigma2_) ;
  NumericMatrix thetastar=clone(theta_) ;
  int K = thetastar.ncol() ;
  int B = thetastar.nrow() ;
  NumericMatrix prec(B, K) ;
  NumericVector p_prec(S) ;
  NumericVector tmp(K) ;
  NumericVector nu0 (1) ;
  NumericVector s20 (1) ;
  //
  // These chains have already been updated
  //
  NumericVector nu0chain = chains.slot("nu.0") ;
  NumericVector s20chain = chains.slot("sigma2.0") ;
  NumericVector d(1) ;
  //
  // Run reduced Gibbs  -- theta is fixed at modal ordinate
  //
  for(int k=0; k < K; ++k){
    prec(_, k) = 1.0/sigma2star(_, k) ;
  }  
  for(int s=0; s < S; ++s){
    s20 = s20chain[s] ;
    nu0 = nu0chain[s] ;
    double total = 1.0 ;
    for(int b=0; b < B; ++b) {
      tmp = dgamma(prec(b, _), 0.5*nu0[0], 2.0 / (nu0[0]*s20[0])) ;
      for(int k = 0; k < K; ++k){
        //total = log(tmp[k]) ;
        total = total * tmp[k] ;
      }
    }
    p_prec[s] = total ;
  }
  return p_prec ;
}

