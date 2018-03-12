#include "miscfunctions.h" // for rdirichlet, tableZ, ...
#include "multibatch.h"
#include "multibatch_reduced.h"
#include "multibatch_pooledvar.h"
#include <Rcpp.h>
#include <Rmath.h>

Rcpp::NumericMatrix toMatrix_pvar(Rcpp::NumericVector x, int NR, int NC) {
  // x is length NR
  Rcpp::NumericMatrix Y(NR, NC);
  int iter = 0;
  for(int i = 0; i < NR; ++i) {
    for(int j = 0; j < NC; ++j) {
      // every component in batch i has the same variance
      Y(i, j) = x[iter];
    }
    iter++;
  }
  return Y;
}

double log_prob_thetap(Rcpp::S4 xmod, Rcpp::NumericMatrix thetastar){
  Rcpp::RNGScope scope;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 params=model.slot("mcmc.params");
  Rcpp::List modes = model.slot("modes");
  int K = thetastar.ncol();
  Rcpp::NumericVector mu = model.slot("mu");
  Rcpp::NumericVector tau2 = model.slot("tau2");
  Rcpp::NumericVector sigma2 = model.slot("sigma2");
  int B = thetastar.nrow();
  Rcpp::NumericVector tmp(1);
  Rcpp::IntegerVector zz = model.slot("z");
  Rcpp::NumericVector tau2_tilde(K);
  Rcpp::NumericVector sigma2_tilde(K);
  Rcpp::NumericVector u = model.slot("u") ;
  // this should be updated for each iteration
  Rcpp::NumericMatrix data_sum(B, K);
  Rcpp::NumericMatrix sumu(B, K);
  Rcpp::IntegerVector nn(K);
  double post_prec;
  double tau_n;
  double mu_n;
  double w1;
  double w2;
  Rcpp::NumericVector tauc(K);
  Rcpp::NumericMatrix iSigma2(B, K);
  Rcpp::NumericVector invs2;
  Rcpp::NumericVector theta(1);
  double df = getDf(model.slot("hyperparams")) ;
  double heavyn = 0.0;
  double heavy_mean = 0.0;
  nn = tableZ(K, zz);
  tauc = sqrt(tau2);
  data_sum = compute_heavy_sums_batch(model);
  sumu = compute_u_sums_batch(model) ;
  //data_mean = compute_means_batch(model);
  tau2_tilde = 1.0 / tau2;
  invs2 = 1.0 / sigma2;    // this is a vector of length B
  sigma2_tilde = Rcpp::as<Rcpp::NumericVector>(toMatrix_pvar(invs2, B, K));
  double total = 0.0;
  for (int k = 0; k < K; ++k) {
    for (int b = 0; b < B; ++b) {
      heavyn = sumu(b, k) / df;
      post_prec = tau2_tilde[k] + heavyn*1.0 * sigma2_tilde[b] ;
      tau_n = sqrt(1/post_prec);
      w1 = tau2_tilde[k]/post_prec;
      w2 = heavyn * sigma2_tilde(b, k)/post_prec;
      heavy_mean = data_sum(b, k) / heavyn / df;
      mu_n = w1*mu[k] + w2*heavy_mean;
      theta[0] = thetastar(b, k);
      tmp = dnorm(theta, mu_n, tau_n, true);
      total += tmp[0];
    }
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::NumericVector marginal_theta_pooled(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params=model.slot("mcmc.params");
    int S = params.slot("iter");
    Rcpp::NumericVector y = model.slot("data");
    int N=y.size();
    Rcpp::List modes = model.slot("modes");
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix thetastar=clone(theta_);
    Rcpp::NumericVector logp(S);
    int K = thetastar.ncol();
    double df = getDf(model.slot("hyperparams")) ;
    for (int s=0; s < S; ++s) {
      model.slot("z") = z_multibatch_pvar(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      model.slot("theta") = theta_multibatch_pvar(model) ;
      model.slot("sigma2") = sigma2_multibatch_pvar(model) ;
      model.slot("mu") = update_mu_batch(model) ;
      model.slot("tau2") = update_tau2_batch(model) ;
      model.slot("sigma2.0") = sigma20_multibatch_pvar(model) ;
      model.slot("nu.0") = nu0_multibatch_pvar(model) ;
      model.slot("pi") = update_p_batch(model) ;
      model.slot("u") = Rcpp::rchisq(N, df) ;
      logp[s]=log_prob_thetap(model, thetastar) ;
    }
    return logp;
}

// [[Rcpp::export]]
double log_prob_sigmap(Rcpp::S4 xmod, Rcpp::NumericVector sigma2star) {
  Rcpp::RNGScope scope;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 params = model.slot("mcmc.params");
  Rcpp::List modes = model.slot("modes");
  Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
  Rcpp::NumericMatrix thetastar = clone(theta_);
  Rcpp::NumericVector x = model.slot("data");
  Rcpp::IntegerVector batch = model.slot("batch") ;
  Rcpp::IntegerVector ub = uniqueBatch(batch) ;
  int n = x.size();
  int K = thetastar.ncol();
  int B = thetastar.nrow();
  Rcpp::NumericVector prec(B);
  Rcpp::NumericVector tmp(K);
  Rcpp::IntegerVector zz = model.slot("z");
  Rcpp::NumericVector nu0 = model.slot("nu.0");
  Rcpp::NumericVector s20 = model.slot("sigma2.0");
  Rcpp::NumericVector nu_n(1);
  Rcpp::NumericVector sigma2_n(1);
  Rcpp::NumericVector shape(1);
  Rcpp::NumericVector rate(1);
  Rcpp::NumericVector sigma2_new(1) ;
  Rcpp::NumericVector u = model.slot("u") ;
  prec = 1.0 / sigma2star;
  Rcpp::NumericMatrix nn(B, K);
  Rcpp::NumericVector ss(B);
  double df = getDf(model.slot("hyperparams")) ;
  for (int i = 0; i < n; i++) {
    for (int b = 0; b < B; ++b) {
      if (batch[i] != ub[b]) {
        continue;
      }
      for (int k = 0; k < K; ++k) {
        if (zz[i] == k + 1) {
          ss[b] += u[i] * pow(x[i] - thetastar(b, k), 2);
        }
      }
    }
  }
  double total = 0.0;
  Rcpp::NumericVector prec_typed(1);
  Rcpp::NumericMatrix tabz = tableBatchZ(xmod);
  for (int b = 0; b < B; ++b) {
    nu_n = nu0 + sum(tabz(b, Rcpp::_));
    sigma2_n = 1.0 / nu_n * (nu0 * s20 + ss[b]/df);
    // calculate shape and rate
    shape = 0.5 * nu_n;
    rate = shape * sigma2_n;
    // calculate probability
    prec_typed[0] = prec[b];
    tmp = Rcpp::dgamma(prec_typed, shape[0], 1.0 / rate[0], true);
    total += tmp[0];
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduced_sigma_pooled(Rcpp::S4 xmod) {
  Rcpp::RNGScope scope;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 params=model.slot("mcmc.params");
  int S = params.slot("iter");
  Rcpp::NumericVector y = model.slot("data");
  int N=y.size();
  Rcpp::List modes = model.slot("modes");
  Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
  Rcpp::NumericMatrix thetastar=clone(theta_);
  Rcpp::NumericVector sigma2_ = Rcpp::as<Rcpp::NumericVector>(modes["sigma2"]);
  Rcpp::NumericVector sigma2star = clone(sigma2_);
  int K = thetastar.ncol();
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::NumericVector logp(S) ;
  //
  // Run reduced Gibbs    -- theta is fixed at modal ordinate
  //
  for (int s = 0; s < S; ++s) {
    model.slot("z") = z_multibatch_pvar(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    //model.slot("theta") = theta_multibatch_pvar(model) ;
    model.slot("sigma2") = sigma2_multibatch_pvar(model) ;
    model.slot("mu") = update_mu_batch(model) ;
    model.slot("tau2") = update_tau2_batch(model) ;
    model.slot("sigma2.0") = sigma20_multibatch_pvar(model) ;
    model.slot("nu.0") = nu0_multibatch_pvar(model) ;
    model.slot("pi") = update_p_batch(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
    logp[s]=log_prob_sigmap(model, sigma2star) ;
  }
  return logp ;
}


// [[Rcpp::export]]
Rcpp::NumericVector reduced_pi_pooled(Rcpp::S4 xmod) {
  Rcpp::RNGScope scope;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 params=model.slot("mcmc.params");
  int S = params.slot("iter");
  Rcpp::NumericVector y = model.slot("data");
  int N=y.size();
  Rcpp::List modes = model.slot("modes");
  Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
  Rcpp::NumericMatrix thetastar=clone(theta_);
  Rcpp::NumericVector p_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
  Rcpp::NumericVector pstar = clone(p_);
  int K = thetastar.ncol();
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::NumericVector logp(S) ;
  //
  // Run reduced Gibbs    -- theta is fixed at modal ordinate
  //
  for (int s = 0; s < S; ++s) {
    model.slot("z") = z_multibatch_pvar(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    //model.slot("theta") = theta_multibatch_pvar(model) ;
    //model.slot("sigma2") = sigma2_multibatch_pvar(model) ;
    model.slot("mu") = update_mu_batch(model) ;
    model.slot("tau2") = update_tau2_batch(model) ;
    model.slot("sigma2.0") = sigma20_multibatch_pvar(model) ;
    model.slot("nu.0") = nu0_multibatch_pvar(model) ;
    model.slot("pi") = update_p_batch(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
    logp[s]=log_prob_pmix(model, pstar) ;
  }
  return logp ;
}


// [[Rcpp::export]]
Rcpp::NumericVector reduced_mu_pooled(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;
    // model and accessories
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::NumericVector y = model.slot("data");
    int N=y.size();
    Rcpp::S4 params=model.slot("mcmc.params");
    int S = params.slot("iter");
    Rcpp::List modes = model.slot("modes");
    // get modal ordinates
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix thetastar=clone(theta_);
    Rcpp::NumericVector mu_ = Rcpp::as<Rcpp::NumericVector>(modes["mu"]);
    Rcpp::NumericVector mustar=clone(mu_);
    Rcpp::NumericVector logp(S) ;
    int K = thetastar.ncol();
    double df = getDf(model.slot("hyperparams")) ;
    for (int s = 0; s < S; ++s) {
      model.slot("z") = z_multibatch_pvar(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      //model.slot("theta") = theta_multibatch_pvar(model) ;
      //model.slot("sigma2") = sigma2_multibatch_pvar(model) ;
      //model.slot("pi") = update_p_batch(model) ;
      model.slot("mu") = update_mu_batch(model) ;
      model.slot("tau2") = update_tau2_batch(model) ;
      model.slot("sigma2.0") = sigma20_multibatch_pvar(model) ;
      model.slot("nu.0") = nu0_multibatch_pvar(model) ;
      model.slot("u") = Rcpp::rchisq(N, df) ;
      logp[s]=log_prob_mu(model, mustar) ;
    }
    return logp ;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduced_tau2_pooled(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;
    // model and accessories
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::NumericVector y = model.slot("data");
    int N=y.size();
    Rcpp::S4 params=model.slot("mcmc.params");
    int S = params.slot("iter");
    Rcpp::List modes = model.slot("modes");
    // get modal ordinates
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix thetastar=clone(theta_);
    Rcpp::NumericVector logp(S) ;
    int K = thetastar.ncol();
    double df = getDf(model.slot("hyperparams")) ;
    for (int s = 0; s < S; ++s) {
      model.slot("z") = z_multibatch_pvar(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      //model.slot("theta") = theta_multibatch_pvar(model) ;
      //model.slot("sigma2") = sigma2_multibatch_pvar(model) ;
      //model.slot("pi") = update_p_batch(model) ;
      //model.slot("mu") = update_mu_batch(model) ;
      model.slot("tau2") = update_tau2_batch(model) ;
      model.slot("sigma2.0") = sigma20_multibatch_pvar(model) ;
      model.slot("nu.0") = nu0_multibatch_pvar(model) ;
      model.slot("u") = Rcpp::rchisq(N, df) ;
      logp[s]=log_prob_tau2(model) ;
    }
    return logp ;
}

// [[Rcpp::export]]
double log_prob_nu0p(Rcpp::S4 xmod, int nu0star) {
    Rcpp::RNGScope scope;
    // get model and accessories
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 mcmcp = model.slot("mcmc.params");
    Rcpp::S4 hypp = model.slot("hyperparams");
    Rcpp::List modes = model.slot("modes");
    // get modal ordinates
    Rcpp::NumericVector sigma2_ = Rcpp::as<Rcpp::NumericVector>(modes["sigma2"]);
    Rcpp::NumericVector sigma2star = clone(sigma2_);
    // hyperparameters
    int K = hypp.slot("k");
    double betas = hypp.slot("beta");
    int B = sigma2star.size() ;
    Rcpp::NumericVector s20 = model.slot("sigma2.0");
    Rcpp::NumericVector d(100) ;  // 100 is the maximum allowed value for nu_0
    Rcpp::NumericVector lpnu0(100);
    double prec = 0.0;
    double lprec = 0.0;
    d = Rcpp::seq_len(100);
    for (int b = 0; b < B; ++b) {
      prec += 1.0/sigma2star[b];
      lprec += log(1.0 / sigma2star[b]);
    }
    Rcpp::NumericVector y1(100);
    Rcpp::NumericVector y2(100);
    Rcpp::NumericVector y3(100);
    y1 = B * K * (0.5 * d * log(s20[0] * 0.5 * d) - lgamma(d * 0.5));
    y2 = (0.5 * d - 1.0) * lprec;
    y3 = d * (betas + 0.5 * s20[0] * prec);
    lpnu0 = y1 + y2 - y3;
    Rcpp::NumericVector prob(100);
    prob = exp(lpnu0); //
    prob = prob / sum(prob);  // this is now normalized
    double log_p = log(prob[nu0star]) ;
    return log_p ;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduced_nu0_pooled(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;
    // model and accessories
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params=model.slot("mcmc.params");
    Rcpp::S4 hypp=model.slot("hyperparams");
    int S = params.slot("iter");
    int K = hypp.slot("k");
    Rcpp::NumericVector y = model.slot("data");
    int N=y.size();
    Rcpp::List modes = model.slot("modes");
    // get modal ordinates
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix thetastar=clone(theta_);
    Rcpp::IntegerVector nu0_ = Rcpp::as<Rcpp::IntegerVector>(modes["nu0"]);
    Rcpp::IntegerVector nu0star=clone(nu0_);
    Rcpp::NumericVector logp(S) ;
    double df = getDf(model.slot("hyperparams")) ;
    for (int s = 0; s < S; ++s) {
      model.slot("z") = z_multibatch_pvar(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      //model.slot("theta") = theta_multibatch_pvar(model) ;
      //model.slot("sigma2") = sigma2_multibatch_pvar(model) ;
      //model.slot("pi") = update_p_batch(model) ;
      //model.slot("mu") = update_mu_batch(model) ;
      //model.slot("tau2") = update_tau2_batch(model) ;
      model.slot("sigma2.0") = sigma20_multibatch_pvar(model) ;
      model.slot("nu.0") = nu0_multibatch_pvar(model) ;
      model.slot("u") = Rcpp::rchisq(N, df) ;
      logp[s]=log_prob_nu0p(model, nu0star[0]) ;
    }
    return logp ;
}


// [[Rcpp::export]]
double log_prob_s20p(Rcpp::S4 xmod) {
  Rcpp::RNGScope scope;
  // get model and accessories
  Rcpp::S4 model(xmod);
  Rcpp::S4 mcmcp = model.slot("mcmc.params");
  Rcpp::S4 hypp = model.slot("hyperparams");
  Rcpp::List modes = model.slot("modes");
  Rcpp::IntegerVector batch = model.slot("batch");
  Rcpp::IntegerVector ub = uniqueBatch(batch);
  int B = ub.size();
  // get ordinal modes.
  Rcpp::IntegerVector nu0_ = Rcpp::as<Rcpp::IntegerVector>(modes["nu0"]);
  Rcpp::NumericVector sigma2_ = Rcpp::as<Rcpp::NumericVector>(modes["sigma2"]);
  Rcpp::NumericVector s20_ = Rcpp::as<Rcpp::NumericVector>(modes["sigma2.0"]);
  Rcpp::NumericVector sigma2star = clone(sigma2_);
  Rcpp::NumericVector s20star = clone(s20_);
  double nu0star = clone(nu0_)[0];
  // get hyperparameters
  int K = hypp.slot("k");
  double a = hypp.slot("a");
  double b = hypp.slot("b");
  // calculate a_k, b_k
  double prec = 0.0;
  for (int b = 0; b < B; ++b) {
    prec += 1.0 / sigma2star[b];
  }
  double a_k = a + 0.5 * K * B * nu0star;
  double b_k = b + 0.5 * nu0star * prec;
  Rcpp::NumericVector p_s20 = dgamma(s20star, a_k, 1.0 / b_k, true);
  double logp = p_s20[0];
  return logp ;
}
