#include "miscfunctions.h" // for rdirichlet
#include <Rmath.h>
#include <Rcpp.h>

using namespace Rcpp ;


// [[Rcpp::export]]
Rcpp::NumericVector compute_loglik(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
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
  NumericMatrix tabz = tableBatchZ(model) ;
  int B = ub.size() ;
  NumericMatrix P(B, K) ;
  NumericMatrix sigma(B, K) ;
  double df = getDf(model.slot("hyperparams")) ;
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
      tmp = P(b, k) * dlocScale_t(x, df, theta(b, k), sigma(b, k)) * this_batch ;
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

// [[Rcpp::export]]
Rcpp::NumericVector update_mu(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
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
  NumericMatrix n_b = tableBatchZ(model) ;
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

// [[Rcpp::export]]
Rcpp::NumericVector update_tau2(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod));
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
Rcpp::NumericVector update_sigma20(Rcpp::S4 xmod){
    RNGScope scope ;
    Rcpp::S4 model(clone(xmod));
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

// [[Rcpp::export]]
Rcpp::NumericVector update_nu0(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
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
Rcpp::NumericMatrix update_multinomialPr(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
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
  double df = getDf(hypp) ;
  for(int k = 0; k < K; ++k){
    NumericVector dens(N) ;
    for(int b = 0; b < B; ++b){
      this_batch = batch == ub[b] ;
      double sigma = sqrt(sigma2(b, k));
      tmp = p[k] * dlocScale_t(x, df, theta(b, k), sigma) * this_batch ;
      dens += tmp ;
    }
    lik(_, k) = dens ;
    rowtotal += dens ;
  }
  NumericMatrix P(N, K) ;
  for(int k=0; k<K; ++k){
    P(_, k) = lik(_, k)/rowtotal ;
  }
  return P ;
}

//
// Note, this is currently the same as the marginal model
//
// [[Rcpp::export]]
Rcpp::NumericVector update_p(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
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
Rcpp::IntegerVector update_z(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector x = model.slot("data") ;
  NumericMatrix theta = model.slot("theta") ;
  IntegerVector batch = model.slot("batch") ;
  int B = theta.nrow() ;
  int n = x.size() ;
  NumericMatrix p(n, K);
  p = update_multinomialPr(xmod) ;
  //NumericMatrix cumP(n, K) ;
  //  Make more efficient
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
Rcpp::NumericMatrix compute_means(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
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
Rcpp::NumericMatrix compute_vars(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
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
  NumericMatrix tabz = tableBatchZ(model) ;
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
        vars(b, k) = sum(pow(x - mn(b,k), 2.0) * this_batch * is_z) / (tabz(b, k) - 1) ;
      }
    }
  }
  return vars ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_prec(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  NumericMatrix vars = compute_vars(model) ;
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
Rcpp::NumericVector compute_logprior(Rcpp::S4 xmod) {
    // set RNG
    RNGScope scope;
    
    // Get model/accessories
    Rcpp::S4 model(clone(xmod)) ;
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
  Rcpp::S4 model(clone(xmod)) ;
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
Rcpp::NumericMatrix update_theta(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector x = model.slot("data") ;
  NumericVector tau2 = model.slot("tau2") ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  NumericMatrix n_hb = tableBatchZ(model) ;
  NumericVector mu = model.slot("mu") ;
  int B = n_hb.nrow() ;
  NumericMatrix theta_new(B, K) ;
  double w1 = 0.0 ;
  double w2 = 0.0 ;
  double mu_n = 0.0 ;
  double tau_n = 0.0 ;
  double post_prec = 0.0 ;
  NumericVector u = model.slot("u") ;
  double df = getDf(hypp) ;
  // find heavy means by batch
  NumericMatrix data_sum =  compute_heavy_sums_batch(model) ;
  NumericMatrix sumu = compute_u_sums_batch(model) ;
  double heavyn = 0.0;
  double heavy_mean = 0.0;
  for (int b = 0; b < B; ++b) {
    for(int k = 0; k < K; ++k){
      heavyn = sumu(b, k) / df;
      post_prec = 1.0/tau2[k] + heavyn*1.0/sigma2(b, k) ;
      if (post_prec == R_PosInf) {
        throw std::runtime_error("Bad simulation. Run again with different start.");
      }
      w1 = (1.0/tau2[k])/post_prec ;
      w2 = (heavyn * 1.0/sigma2(b, k))/post_prec ;
      heavy_mean = data_sum(b, k) / heavyn / df;
      mu_n = w1*mu[k] + w2*heavy_mean ;
      tau_n = sqrt(1.0/post_prec) ;
      theta_new(b, k) = as<double>(rnorm(1, mu_n, tau_n)) ;
    }
  }
  return theta_new ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_sigma2(Rcpp::S4 xmod){
    Rcpp::RNGScope scope;
    Rcpp::S4 model(clone(xmod)) ;
    // get model
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
    NumericVector u = model.slot("u") ;
    double df = getDf(model.slot("hyperparams")) ;

    //IntegerVector nn = model.slot("zfreq");
    // get batch info
    Rcpp::NumericMatrix tabz = tableBatchZ(model);
    Rcpp::IntegerVector batch = model.slot("batch");
    Rcpp::IntegerVector ub = uniqueBatch(batch);
    Rcpp::NumericMatrix ss(B, K);

    for (int i = 0; i < n; ++i) {
        for (int b = 0; b < B; ++b) {
            if (batch[i] != ub[b]) {
                continue;
            }

            for (int k = 0; k < K; ++k){
              if((z[i] == k+1) && (batch[i] == b+1)){
                    ss(b, k) += u[i] * pow(x[i] - theta(b, k), 2);
                    //ss(b, k) += pow(x[i] - theta(b, k), 2);
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
            sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss(b, k)/df);
           // sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss(b, k));
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

// [[Rcpp::export]]
Rcpp::IntegerMatrix update_probz(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
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
Rcpp::S4 cpp_burnin(Rcpp::S4 object, Rcpp::S4 mcmcp) {
  RNGScope scope ;
  Rcpp::S4 model(clone(object)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  Rcpp::S4 params(mcmcp) ;
  IntegerVector up = params.slot("param_updates") ;
  int S = params.slot("burnin") ;
  NumericVector x = model.slot("data") ;
  int N = x.size() ;
  double df = getDf(model.slot("hyperparams")) ;
  if( S < 1 ){
    return model ;
  }
  for(int s = 0; s < S; ++s){
    if(up[7] > 0){
      model.slot("z") = update_z(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    }
    if(up[0] > 0)
      try {
        model.slot("theta") = update_theta(model) ;
      } catch(std::runtime_error &ex) {
          forward_exception_to_r(ex);
      } catch(...) {
          ::Rf_error("c++ exception (unknown reason)");
      }
    if(up[1] > 0)
      model.slot("sigma2") = update_sigma2(model) ;
    if(up[3] > 0)
      model.slot("mu") = update_mu(model) ;
    if(up[4] > 0)
      model.slot("tau2") = update_tau2(model) ;
    if(up[6] > 0)
      model.slot("sigma2.0") = update_sigma20(model) ;
    if(up[5] > 0)
      model.slot("nu.0") = update_nu0(model) ;
    if(up[2] > 0)
      model.slot("pi") = update_p(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
  }
  // compute log prior probability from last iteration of burnin
  // compute log likelihood from last iteration of burnin
  NumericVector lls2(1);
  NumericVector ll(1);
  lls2 = stageTwoLogLikBatch(model);
  ll = compute_loglik(model);
  ll = ll + lls2;
  model.slot("loglik") = ll;
  model.slot("logprior") = compute_logprior(model) ;
  return model ;
  // return vars ;
}

// [[Rcpp::export]]
Rcpp::S4 cpp_mcmc(Rcpp::S4 object, Rcpp::S4 mcmcp) {
  RNGScope scope ;
  Rcpp::S4 model(clone(object)) ;
  Rcpp::S4 chain(model.slot("mcmc.chains")) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  Rcpp::S4 params(mcmcp) ;
  IntegerVector up = params.slot("param_updates") ;
  int K = getK(hypp) ;
  int T = params.slot("thin") ;
  int S = params.slot("iter") ;
  if( S < 1 ) return model ;
  NumericVector x = model.slot("data") ;
  int N = x.size() ;
  double df = getDf(model.slot("hyperparams")) ;
  NumericMatrix thetac = chain.slot("theta") ;
  NumericMatrix sigma2c = chain.slot("sigma2") ;
  NumericMatrix th = model.slot("theta");
  int B = th.nrow();
  NumericMatrix s2(B, K);
  NumericMatrix pmix = chain.slot("pi") ;
  NumericMatrix zfreq = chain.slot("zfreq") ;
  NumericMatrix mu = chain.slot("mu") ;
  NumericMatrix tau2 = chain.slot("tau2") ;
  NumericVector nu0 = chain.slot("nu.0") ;
  NumericVector sigma2_0 = chain.slot("sigma2.0") ;
  NumericVector loglik_ = chain.slot("loglik") ;
  NumericVector logprior_ = chain.slot("logprior") ;
  //NumericVector th(K) ;
  //NumericVector s2(K) ;
  NumericVector p(K) ;
  NumericVector m(K) ; //mu
  NumericVector t2(K) ;//tau2
  NumericVector n0(1) ;//nu0
  IntegerVector z(N) ;
  NumericVector u(N) ;
  NumericVector s20(1) ; //sigma2_0
  NumericVector lls2(1) ;  // stage 2 log lik
  NumericVector ll(1) ;
  NumericVector lp(1) ;
  IntegerVector tmp(K) ;
  IntegerVector zf(K) ;
  // Initial values
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
  mu(0, _) = m ;
  nu0[0] = n0[0] ;
  sigma2_0[0] = s20[0] ;
  loglik_[0] = ll[0] ;
  logprior_[0] = lp[0] ;
  thetac(0, _) = as<Rcpp::NumericVector>(model.slot("theta")) ;
  tau2(0, _) = t2 ;
  sigma2c(0, _) = as<Rcpp::NumericVector>(model.slot("sigma2")) ;
  pmix(0, _) = p ;
  zfreq(0, _) = zf ;

  // Is accessing a slot in an object expensive?
  // Currently, there is no alternative as the current values are
  // stored in the object.  Hence, the entire object has to be passed
  // to the updating functions.
  // start at 1 instead of zero. Initial values are as above
  //up[7] = 0;
  //up[6] = 0;
  //up[5] = 0;
  //up[4] = 0;
  //up[3] = 0;
  //up[2] = 0;
  //up[1] = 0;
  //up[0] = 0;
  for(int s = 1; s < S; ++s){
    if(up[7] > 0){
      z = update_z(model) ;
      model.slot("z") = z ;
      tmp = tableZ(K, z) ;
      model.slot("probz") = update_probz(model) ;
      model.slot("zfreq") = tmp ;
    } else {
      tmp = model.slot("zfreq") ;
    }
    zfreq(s, _) = tmp ;
    if(up[0] > 0) {
      model.slot("theta") = update_theta(model) ;
    }
    thetac(s, _) = as<Rcpp::NumericVector>(model.slot("theta")) ;
    if(up[1] > 0){
      model.slot("sigma2") = update_sigma2(model) ;
    }
    sigma2c(s, _) = as<Rcpp::NumericVector>(model.slot("sigma2"));
    if(up[2] > 0){
      p = update_p(model) ;
      model.slot("pi") = p ;
    } else {
      p = model.slot("pi") ;
    }
    pmix(s, _) = p ;
    if(up[3] > 0){
      m = update_mu(model) ;
      model.slot("mu") = m ;
    } else {
      m = model.slot("mu") ;
    }
    mu(s, _) = m ;
    if(up[4] > 0){
      t2 = update_tau2(model) ;
      model.slot("tau2") = t2 ;
    } else {
      t2 = model.slot("tau2") ;
    }
    tau2(s, _) = t2 ;
    if(up[5] > 0){
      n0 = update_nu0(model) ;
      model.slot("nu.0") = n0 ;
    } else {
      n0 = model.slot("nu.0") ;
    }
    nu0[s] = n0[0] ;
    if(up[6] > 0){
      s20 = update_sigma20(model) ;
      model.slot("sigma2.0") = s20 ;
    } else {
      s20 = model.slot("sigma2.0") ;
    }
    sigma2_0[s] = s20[0] ;
    ll = compute_loglik(model) ;
    lls2 = stageTwoLogLikBatch(model) ;
    ll = ll + lls2 ;
    loglik_[s] = ll[0] ;
    model.slot("loglik") = ll ;
    lp = compute_logprior(model) ;
    logprior_[s] = lp[0] ;
    model.slot("logprior") = lp ;
    u = Rcpp::rchisq(N, df) ;
    model.slot("u") = u;
    // Thinning
    for(int t = 0; t < T; ++t){
      if(up[7] > 0){
        model.slot("z") = update_z(model) ;
        model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      }
      if(up[0] > 0)
        model.slot("theta") = update_theta(model) ;
      if(up[1] > 0)
        model.slot("sigma2") = update_sigma2(model) ;
      if(up[2] > 0)
        model.slot("pi") = update_p(model) ;
      if(up[3] > 0)
        model.slot("mu") = update_mu(model) ;
      if(up[4] > 0)
        model.slot("tau2") = update_tau2(model) ;
      if(up[5] > 0)
        model.slot("nu.0") = update_nu0(model) ;
     if(up[6] > 0)
       model.slot("sigma2.0") = update_sigma20(model) ;
     model.slot("u") = Rcpp::rchisq(N, df) ;
    }
  }
  //
  // assign chains back to object
  //
  chain.slot("theta") = thetac ;
  chain.slot("sigma2") = sigma2c ;
  chain.slot("pi") = pmix ;
  chain.slot("mu") = mu ;
  chain.slot("tau2") = tau2 ;
  chain.slot("nu.0") = nu0 ;
  chain.slot("sigma2.0") = sigma2_0 ;
  chain.slot("zfreq") = zfreq ;
  chain.slot("loglik") = loglik_ ;
  chain.slot("logprior") = logprior_ ;
  model.slot("mcmc.chains") = chain ;
  return model ;
}
