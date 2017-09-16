#include "update.h" // for rdirichlet
#include "miscfunctions.h" // for rdirichlet
#include "multibatch.h" // getK
#include <Rmath.h>
#include <Rcpp.h>

using namespace Rcpp ;

// [[Rcpp::export]]
Rcpp::NumericVector loglik_multibatch_pvar(Rcpp::S4 xmod){
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

Rcpp::NumericVector mu_multibatch_pvar(Rcpp::S4 xmod){
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

Rcpp::NumericVector tau2_multibatch_pvar(Rcpp::S4 xmod){
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
Rcpp::NumericVector sigma20_multibatch_pvar(Rcpp::S4 xmod){
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

Rcpp::NumericVector nu0_multibatch_pvar(Rcpp::S4 xmod){
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
Rcpp::NumericMatrix multinomialPr_multibatch_pvar(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  NumericVector p = model.slot("pi") ;
  NumericVector sigma2 = model.slot("sigma2") ;
  NumericMatrix theta = model.slot("theta") ;
  int B = theta.nrow() ;
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
      //tmp = p[k] * dnorm(x, theta(b, k), sqrt(sigma2(b, k))) * this_batch ;
      tmp = p[k] * dnorm(x, theta(b, k), sqrt(sigma2[b])) * this_batch ;
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
Rcpp::NumericVector p_multibatch_pvar(Rcpp::S4 xmod) {
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
Rcpp::IntegerVector z_multibatch_pvar(Rcpp::S4 xmod) {
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
  p = multinomialPr_multibatch_pvar(xmod) ;
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
  int counter = model.slot(".internal.counter");
  counter++;
  model.slot(".internal.counter") = counter;
  return model.slot("z") ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix means_multibatch_pvar(Rcpp::S4 xmod) {
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
Rcpp::NumericMatrix vars_multibatch_pvar(Rcpp::S4 xmod) {
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
Rcpp::NumericMatrix prec_multibatch_pvar(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericMatrix vars = vars_multibatch_pvar(xmod) ;
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
Rcpp::NumericVector logprior_multibatch_pvar(Rcpp::S4 xmod) {
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
Rcpp::NumericVector stagetwo_multibatch_pvar(Rcpp::S4 xmod) {
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
Rcpp::NumericMatrix theta_multibatch_pvar(Rcpp::S4 xmod){
    RNGScope scope ;
    Rcpp::S4 model(xmod) ;
    Rcpp::S4 hypp(model.slot("hyperparams")) ;
    int K = getK(hypp) ;
    NumericVector x = model.slot("data") ;
    // NumericMatrix theta = model.slot("theta") ;
    NumericVector tau2 = model.slot("tau2") ;
    NumericVector sigma2 = model.slot("sigma2") ;
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
            post_prec = 1.0/tau2[k] + n_hb(b, k)*1.0/sigma2[b] ;
            if (post_prec == R_PosInf) {
                throw std::runtime_error("Bad simulation. Run again with different start.");
            }
            tau_n = sqrt(1.0/post_prec) ;
            w1 = (1.0/tau2[k])/post_prec ;
            w2 = (n_hb(b, k) * 1.0/sigma2[b])/post_prec ;
            mu_n = w1*mu[k] + w2*ybar(b, k) ;
            theta_new(b, k) = as<double>(rnorm(1, mu_n, tau_n)) ;
        }
    }
    return theta_new ;
}

// [[Rcpp::export]]
Rcpp::NumericVector sigma2_multibatch_pvar(Rcpp::S4 xmod){
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

    // get batch info
    Rcpp::NumericMatrix tabz = tableBatchZ(xmod);
    Rcpp::IntegerVector batch = model.slot("batch");
    Rcpp::IntegerVector ub = uniqueBatch(batch);
    //Rcpp::NumericMatrix ss(B, K);
    Rcpp::NumericVector ss(B);

    for (int i = 0; i < n; ++i) {
      for (int b = 0; b < B; ++b) {
        if (batch[i] != ub[b]) {
          continue;
        }
        for (int k = 0; k < K; ++k){
          if (z[i] == k+1){
            //ss(b, k) += pow(x[i] - theta(b, k), 2);
            ss[b] += pow(x[i] - theta(b, k), 2);
          }
        }
      }
    }
    //NumericMatrix sigma2_nh(B, K);
    double shape;
    double rate;
    double sigma2_nh;
    double nu_n;
    //Rcpp::NumericMatrix sigma2_tilde(B, K);
    //Rcpp::NumericMatrix sigma2_(B, K);
    Rcpp::NumericVector sigma2_tilde(B);
    Rcpp::NumericVector sigma2_(B);

    for (int b = 0; b < B; ++b) {
      //for (int k = 0; k < K; ++k) {
        //nu_n = nu_0 + tabz(b, k);
      nu_n = nu_0 + sum(tabz(b, _));
        //sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss(b, k));
      sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss[b]);
      shape = 0.5 * nu_n;
      rate = shape * sigma2_nh;
        //sigma2_tilde(b, k) = Rcpp::as<double>(rgamma(1, shape, 1.0/rate));
      sigma2_tilde[b] = Rcpp::as<double>(rgamma(1, shape, 1.0/rate));
        //sigma2_(b, k) = 1.0 / sigma2_tilde(b, k);
      sigma2_[b] = 1.0 / sigma2_tilde[b];
        //}
    }
    return sigma2_;
}

Rcpp::IntegerMatrix probz_multibatch_pvar(Rcpp::S4 xmod){
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
Rcpp::S4 burnin_multibatch_pvar(Rcpp::S4 xmod, Rcpp::S4 mcmcp) {
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
    if(up[7] > 0){
      model.slot("z") = z_multibatch_pvar(xmod) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    }
    model.slot("data.mean") = means_multibatch_pvar(xmod) ;
    model.slot("data.prec") = prec_multibatch_pvar(xmod) ;
    if(up[0] > 0)
      try {
        model.slot("theta") = theta_multibatch_pvar(xmod) ;
      } catch(std::runtime_error &ex) {
          forward_exception_to_r(ex);
      } catch(...) {
          ::Rf_error("c++ exception (unknown reason)");
      }
    if(up[1] > 0)
      model.slot("sigma2") = sigma2_multibatch_pvar(xmod) ;
    if(up[3] > 0)
      model.slot("mu") = mu_multibatch_pvar(xmod) ;
    if(up[4] > 0)
      model.slot("tau2") = tau2_multibatch_pvar(xmod) ;
    if(up[6] > 0)
      model.slot("sigma2.0") = sigma20_multibatch_pvar(xmod) ;
    if(up[5] > 0)
      model.slot("nu.0") = nu0_multibatch_pvar(xmod) ;
    if(up[2] > 0)
      model.slot("pi") = p_multibatch_pvar(xmod) ;
  }
  // compute log prior probability from last iteration of burnin
  // compute log likelihood from last iteration of burnin
  model.slot("loglik") = loglik_multibatch_pvar(xmod) ;
  model.slot("logprior") = logprior_multibatch_pvar(xmod) ;
  return xmod ;
  // return vars ;
}

// [[Rcpp::export]]
Rcpp::S4 mcmc_multibatch_pvar(Rcpp::S4 object, Rcpp::S4 mcmcp) {
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
      th = as<Rcpp::NumericVector>(theta_multibatch_pvar(xmod));
      model.slot("theta") = th ;
    } else {
      th = model.slot("theta") ;
    }
    theta(s, _) = th ;
    if(up[1] > 0){
      s2 = as<Rcpp::NumericVector>(sigma2_multibatch_pvar(xmod));
      model.slot("sigma2") = s2 ;
    } else {
      s2 = model.slot("sigma2") ;
    }
    sigma2(s, _) = s2 ;
    if(up[2] > 0){
      p = p_multibatch_pvar(xmod) ;
      model.slot("pi") = p ;
    } else {
      p = model.slot("pi") ;
    }
    pmix(s, _) = p ;
    if(up[3] > 0){
      m = mu_multibatch_pvar(xmod) ;
      model.slot("mu") = m ;
    } else {
      m = model.slot("mu") ;
    }
    mu(s, _) = m ;
    if(up[4] > 0){
      t2 = tau2_multibatch_pvar(xmod) ;
      model.slot("tau2") = t2 ;
    } else {
      t2 = model.slot("tau2") ;
    }
    tau2(s, _) = t2 ;
    if(up[5] > 0){
      n0 = nu0_multibatch_pvar(xmod) ;
      model.slot("nu.0") = n0 ;
    } else {
      n0 = model.slot("nu.0") ;
    }
    nu0[s] = n0[0] ;
    if(up[6] > 0){
      s20 = sigma20_multibatch_pvar(xmod) ;
      model.slot("sigma2.0") = s20 ;
    } else {
      s20 = model.slot("sigma2.0") ;
    }
    sigma2_0[s] = s20[0] ;
    if(up[7] > 0){
      z = z_multibatch_pvar(xmod) ;
      model.slot("z") = z ;
      tmp = tableZ(K, z) ;
      model.slot("probz") = probz_multibatch_pvar(xmod) ;
      model.slot("zfreq") = tmp ;
      // mean and prec only change if z is updated
      model.slot("data.mean") = means_multibatch_pvar(xmod) ;
      model.slot("data.prec") = prec_multibatch_pvar(xmod) ;
    } else {
      tmp = model.slot("zfreq") ;
    }
    Z(s, _) = z ;
    zfreq(s, _) = tmp ;
    ll = loglik_multibatch_pvar(xmod) ;
    loglik_[s] = ll[0] ;
    model.slot("loglik") = ll ;
    lp = logprior_multibatch_pvar(xmod) ;
    logprior_[s] = lp[0] ;
    model.slot("logprior") = lp ;
    // Thinning
    for(int t = 0; t < T; ++t){
      if(up[0] > 0)
        model.slot("theta") = theta_multibatch_pvar(xmod) ;
      if(up[1] > 0)
        model.slot("sigma2") = sigma2_multibatch_pvar(xmod) ;
      if(up[2] > 0)
        model.slot("pi") = p_multibatch_pvar(xmod) ;
      if(up[3] > 0)
        model.slot("mu") = mu_multibatch_pvar(xmod) ;
      if(up[4] > 0)
        model.slot("tau2") = tau2_multibatch_pvar(xmod) ;
      if(up[5] > 0)
        model.slot("nu.0") = nu0_multibatch_pvar(xmod) ;
     if(up[6] > 0)
        model.slot("sigma2.0") = sigma20_multibatch_pvar(xmod) ;
      if(up[7] > 0){
        model.slot("z") = z_multibatch_pvar(xmod) ;
        model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      }
      model.slot("data.mean") = means_multibatch_pvar(xmod) ;
      model.slot("data.prec") = prec_multibatch_pvar(xmod) ;
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
