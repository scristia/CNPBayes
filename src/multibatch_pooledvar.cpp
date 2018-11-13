#include "miscfunctions.h" // for rdirichlet
#include "multibatch.h" // getK
#include <Rmath.h>
#include <Rcpp.h>

using namespace Rcpp ;

//[[Rcpp::export]]
Rcpp::IntegerVector sample_componentsP(Rcpp::IntegerVector x, int size, Rcpp::NumericVector prob){
  int n = x.size() ;
  Rcpp::IntegerVector z = clone(x);
  for(int i=0; i < n; i++){
    //initialize accumulator ;
    double accept = 0 ;
    Rcpp::NumericVector u=runif(1);
    for(int j = 0; j < n; j++){
      accept += prob[j] ;
      if( u[0] < accept ) {
        z[i] = x[j] ;
        break ;
      }
    }
  }
  return(z);
}

//[[Rcpp::export]]
Rcpp::S4 update_predictiveP(Rcpp::S4 xmod){
  Rcpp::RNGScope scope;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::NumericMatrix theta = model.slot("theta");
  // sigma2 is a vector of length B (B = number batches)
  Rcpp::NumericVector sigma2 = model.slot("sigma2");
  Rcpp::NumericVector prob = model.slot("pi");
  Rcpp::IntegerVector batch = model.slot("batch") ;
  Rcpp::IntegerVector ub = uniqueBatch(batch);
  int K = theta.ncol();
  int nb = theta.nrow();
  int B = ub.size() ;
  Rcpp::IntegerVector components=seq_len(K);
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::IntegerVector z(nb);
  Rcpp::NumericMatrix ystar(B, K) ;
  Rcpp::IntegerMatrix zstar(B, K) ;
  Rcpp::NumericVector u=Rcpp::rchisq(K*B, df) ;;
  int n=1;
  int j=0;
  // sample components according to mixture probabilities
  // mixture probabilities are assumed to be the same for each batch
  z=sample_componentsP(components, K, prob);
  z=z-1;  // Need to subtract 1 for this to index the right column
  Rcpp::NumericVector yy(K*B);
  Rcpp::IntegerVector zz(K*B);
  for(int k=0; k < K; k++){
    for(int b = 0; b < B; b++){
      int index = z[k];
      //double sigma=pow(sigma2[index], 0.5) ;
      double sigma=pow(sigma2[b], 0.5) ;
      zstar(b, k) = index ;
      ystar(b, k) = (rlocScale_t(n, theta(b, index), sigma, df, u[j]))[0];
      yy[j] = ystar(b, k);
      zz[j] = zstar(b, k);
      j++;
    }
  }
  model.slot("predictive") = yy;
  model.slot("zstar") = zz;
  return model ;
}


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
  NumericVector sigma2 = model.slot("sigma2") ;
  IntegerVector batch_freq = model.slot("batchElements") ;
  NumericVector loglik_(1) ;
  NumericMatrix tabz = tableBatchZ(xmod) ;
  int B = ub.size() ;
  NumericMatrix P(B, K) ;
  //NumericMatrix sigma(B, K) ;
  NumericVector sigma(B);
  double df = getDf(model.slot("hyperparams")) ;
  // component probabilities for each batch
  for(int b = 0; b < B; ++b){
    int rowsum = 0 ;
    rowsum = sum(tabz(b, _));
    //for(int k = 0; k < K; ++k){
    //  rowsum += tabz(b, k) ;
    //sigma(b, k) = sqrt(sigma2(b, k)) ;
    sigma[b] = sqrt(sigma2[b]) ;
    //}
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
      tmp = P(b, k) * dlocScale_t(x, df, theta(b, k), sigma[b]) * this_batch ;
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
Rcpp::NumericVector sigma20_multibatch_pvar(Rcpp::S4 xmod){
    RNGScope scope ;
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 hypp(model.slot("hyperparams")) ;
    int K = getK(hypp) ;
    IntegerVector batch = model.slot("batch") ;
    IntegerVector ub = uniqueBatch(batch) ;
    int B = ub.size() ;
    NumericVector a = hypp.slot("a") ;
    NumericVector b = hypp.slot("b") ;
    NumericVector nu_0 = model.slot("nu.0") ;
    NumericVector sigma2 = model.slot("sigma2") ;
    Rcpp::NumericVector sigma2_0_old = model.slot("sigma2.0");
    NumericVector prec(1) ;

    for(int i = 0; i < B; ++i){
      //for(int k = 0; k < K; ++k){
      //prec[0] += 1.0/sigma2(i, k) ;
      prec[0] += 1.0/sigma2[i] ;
      //  }
    }

    NumericVector a_k(1) ;
    NumericVector b_k(1) ;
    // is this right?
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
Rcpp::NumericVector nu0_multibatch_pvar(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector sigma2 = model.slot("sigma2") ;
  int B = sigma2.size() ;
  double sigma2_0 = model.slot("sigma2.0") ;
  double prec = 0.0;
  double lprec = 0.0 ;
  double betas = hypp.slot("beta") ;
  for(int i = 0; i < B; ++i){
    //for(int k = 0; k < K; ++k){
    //prec += 1.0/sigma2(i, k) ;
    prec += 1.0/sigma2[i] ;
    //lprec += log(1.0/sigma2(i, k)) ;
    lprec += log(1.0/sigma2[i]);
    //}
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
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
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
  double df = getDf(hypp) ;
  for(int k = 0; k < K; ++k){
    NumericVector dens(N) ;
    for(int b = 0; b < B; ++b){
      this_batch = batch == ub[b] ;
      //tmp = p[k] * dnorm(x, theta(b, k), sqrt(sigma2(b, k))) * this_batch ;
      tmp = p[k] * dlocScale_t(x, df, theta(b, k), sqrt(sigma2[b])) * this_batch ;
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

// [[Rcpp::export]]
Rcpp::IntegerVector z_multibatch_pvar(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
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
Rcpp::NumericVector stagetwo_multibatch_pvar(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericMatrix theta = model.slot("theta") ;
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector mu = model.slot("mu") ;
  NumericVector nu0 = model.slot("nu.0") ;
  NumericVector s20 = model.slot("sigma2.0") ;
  NumericVector sigma2=model.slot("sigma2") ;
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
    //likprec[0] = sum(log(dgamma(1.0/sigma2(_, k), 0.5*nu0[0], 1.0/(0.5*nu0[0]*s20[0])))) ;
    likprec[0] = sum(log(dgamma(1.0/sigma2, 0.5*nu0[0], 1.0/(0.5*nu0[0]*s20[0])))) ;
    LL[0] += liknorm[0] + likprec[0] ;
  }
  return LL ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix theta_multibatch_pvar(Rcpp::S4 xmod){
    RNGScope scope ;
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 hypp(model.slot("hyperparams")) ;
    int K = getK(hypp) ;
    NumericVector x = model.slot("data") ;
    NumericVector tau2 = model.slot("tau2") ;
    NumericVector sigma2 = model.slot("sigma2") ;
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
        post_prec = 1.0/tau2[k] + heavyn*1.0/sigma2[b] ;
        if (post_prec == R_PosInf) {
          throw std::runtime_error("Bad simulation. Run again with different start.");
        }
        w1 = (1.0/tau2[k])/post_prec ;
        w2 = (heavyn * 1.0/sigma2[b])/post_prec ;
        heavy_mean = data_sum(b, k) / heavyn / df;
        mu_n = w1*mu[k] + w2*heavy_mean ;
        tau_n = sqrt(1.0/post_prec) ;
        theta_new(b, k) = as<double>(rnorm(1, mu_n, tau_n)) ;
      }
    }
    return theta_new ;
}

// [[Rcpp::export]]
Rcpp::NumericVector sigma2_multibatch_pvar(Rcpp::S4 xmod){
    Rcpp::RNGScope scope;
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
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

    // get batch info
    Rcpp::NumericMatrix tabz = tableBatchZ(model);
    Rcpp::IntegerVector batch = model.slot("batch");
    Rcpp::IntegerVector ub = uniqueBatch(batch);
    Rcpp::NumericVector ss(B);

    for (int i = 0; i < n; ++i) {
      for (int b = 0; b < B; ++b) {
        if (batch[i] != ub[b]) {
          continue;
        }
        for (int k = 0; k < K; ++k){
          if (z[i] == k+1){
            //ss(b, k) += pow(x[i] - theta(b, k), 2);
            ss[b] += u[i] * pow(x[i] - theta(b, k), 2);
          }
        }
      }
    }
    double shape;
    double rate;
    double sigma2_nh;
    double nu_n;
    Rcpp::NumericVector sigma2_tilde(B);
    Rcpp::NumericVector sigma2_(B);

    for (int b = 0; b < B; ++b) {
      nu_n = nu_0 + sum(tabz(b, _));
      sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss[b]/df);
      shape = 0.5 * nu_n;
      rate = shape * sigma2_nh;
      sigma2_tilde[b] = Rcpp::as<double>(rgamma(1, shape, 1.0/rate));
      sigma2_[b] = 1.0 / sigma2_tilde[b];
    }
    return sigma2_;
}

// [[Rcpp::export]]
Rcpp::S4 burnin_multibatch_pvar(Rcpp::S4 object, Rcpp::S4 mcmcp) {
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
    model.slot("z") = z_multibatch_pvar(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    model.slot("theta") = theta_multibatch_pvar(model) ;
    model.slot("sigma2") = sigma2_multibatch_pvar(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("sigma2.0") = sigma20_multibatch_pvar(model) ;
    model.slot("nu.0") = nu0_multibatch_pvar(model) ;
    model.slot("pi") = update_p(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
  }
  // compute log prior probability from last iteration of burnin
  // compute log likelihood from last iteration of burnin
  NumericVector ll(1);
  NumericVector lls2(1);
  ll = loglik_multibatch_pvar(model);
  lls2 = stagetwo_multibatch_pvar(model);
  ll = ll + lls2;
  model.slot("loglik") = ll ;
  model.slot("logprior") = compute_logprior(model) ;
  return model ;
}

// [[Rcpp::export]]
Rcpp::S4 mcmc_multibatch_pvar(Rcpp::S4 object, Rcpp::S4 mcmcp) {
  RNGScope scope ;
  Rcpp::S4 model(clone(object)) ;
  Rcpp::S4 chain(model.slot("mcmc.chains")) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  Rcpp::S4 params(mcmcp) ;
  IntegerVector up = params.slot("param_updates") ;
  int K = getK(hypp) ;
  int T = params.slot("thin") ;
  int S = params.slot("iter") ;
  S--; // decrement so that S is a zero-based index
  T--;
  NumericVector x = model.slot("data") ;
  int N = x.size() ;
  double df = getDf(model.slot("hyperparams")) ;
  NumericMatrix theta = chain.slot("theta") ;
  int B = theta.nrow();
  NumericMatrix sigma2 = chain.slot("sigma2") ;
  NumericMatrix pmix = chain.slot("pi") ;
  NumericMatrix zfreq = chain.slot("zfreq") ;
  NumericMatrix mu = chain.slot("mu") ;
  NumericMatrix tau2 = chain.slot("tau2") ;
  NumericVector nu0 = chain.slot("nu.0") ;
  NumericVector sigma2_0 = chain.slot("sigma2.0") ;
  NumericVector loglik_ = chain.slot("loglik") ;
  NumericVector logprior_ = chain.slot("logprior") ;
  NumericMatrix predictive_ = chain.slot("predictive") ;
  IntegerMatrix zstar_ = chain.slot("zstar") ;
  NumericVector th(K) ;
  NumericVector s2(K) ;
  NumericVector p(K) ;
  NumericVector m(K) ; //mu
  NumericVector t2(K) ;//tau2
  NumericVector n0(1) ;//nu0
  IntegerVector z(N) ;
  NumericVector u(N) ;
  NumericVector s20(1) ; //sigma2_0
  NumericVector ll(1) ;
  NumericVector lls2(1) ;
  NumericVector lp(1) ;
  IntegerVector tmp(K) ;
  IntegerVector zf(K) ;
  Rcpp::NumericVector ystar = NumericVector(B*K);
  Rcpp::IntegerVector zstar = IntegerVector(B*K);
  //
  // This for-loop uses a zero-based index
  //
  // For each parameter update
  //       a.  update values in 'current value' slots
  //       b.  update chain
  //
  // The current value at simulation s is the s-1 row of the chain.
  //
  for(int s = 0; s < (S + 1); ++s){
    z = z_multibatch_pvar(model) ;
    model.slot("z") = z ;
    tmp = tableZ(K, z) ;
    model.slot("probz") = update_probz(model) ;
    model.slot("zfreq") = tmp ;
    // mean and prec only change if z is updated
    zfreq(s, _) = tmp ;
    th = as<Rcpp::NumericVector>(theta_multibatch_pvar(model));
    model.slot("theta") = th ;
    theta(s, _) = th ;
    s2 = as<Rcpp::NumericVector>(sigma2_multibatch_pvar(model));
    model.slot("sigma2") = s2 ;
    sigma2(s, _) = s2 ;
    p = update_p(model) ;
    model.slot("pi") = p ;
    pmix(s, _) = p ;
    m = update_mu(model) ;
    model.slot("mu") = m ;
    mu(s, _) = m ;
    t2 = update_tau2(model) ;
    model.slot("tau2") = t2 ;
    tau2(s, _) = t2 ;
    n0 = nu0_multibatch_pvar(model) ;
    model.slot("nu.0") = n0 ;
    nu0[s] = n0[0] ;
    s20 = sigma20_multibatch_pvar(model) ;
    model.slot("sigma2.0") = s20 ;
    sigma2_0[s] = s20[0] ;
    ll = loglik_multibatch_pvar(model) ;
    lls2 = stagetwo_multibatch_pvar(model) ;
    ll = ll + lls2 ;
    loglik_[s] = ll[0] ;
    model.slot("loglik") = ll ;
    lp = compute_logprior(model) ;
    logprior_[s] = lp[0] ;
    model.slot("logprior") = lp ;
    u = Rcpp::rchisq(N, df) ;
    model.slot("u") = u ;
    model = update_predictiveP(model);
    ystar = model.slot("predictive");
    zstar = model.slot("zstar");
    predictive_(s, _) = ystar ;
    zstar_(s, _) = zstar ;
    // Thinning
    for(int t = 0; t < T; ++t){
      model.slot("z") = z_multibatch_pvar(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      model.slot("theta") = theta_multibatch_pvar(model) ;
      model.slot("sigma2") = sigma2_multibatch_pvar(model) ;
      model.slot("pi") = update_p(model) ;
      model.slot("mu") = update_mu(model) ;
      model.slot("tau2") = update_tau2(model) ;
      model.slot("nu.0") = nu0_multibatch_pvar(model) ;
      model.slot("sigma2.0") = sigma20_multibatch_pvar(model) ;
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
  model.slot("mcmc.chains") = chain ;
  return model ;
}
