#include "miscfunctions.h" // for rdirichlet, tableZ, ...
#include "multibatch.h"
#include "singlebatch.h" // for log_ddirichlet_
#include <Rcpp.h>
#include <Rmath.h>

Rcpp::NumericMatrix toMatrix(Rcpp::NumericVector x, int NR, int NC) {
    Rcpp::NumericMatrix Y(NR, NC);
    int iter = 0;
    for(int j = 0; j < NC; ++j) {
        for(int i = 0; i < NR; ++i) {
            Y(i, j) = x[iter];
            iter++;
        }
    }
    return Y;
}

// [[Rcpp::export]]
double log_prob_theta(Rcpp::S4 xmod, Rcpp::NumericMatrix thetastar) {
  Rcpp::RNGScope scope;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  int K = thetastar.ncol();
  int B = thetastar.nrow();
  Rcpp::NumericVector tau2 = model.slot("tau2") ;
  Rcpp::NumericVector tau2_tilde(K);
  Rcpp::NumericMatrix sigma2_tilde(K);
  Rcpp::NumericMatrix sigma2 = model.slot("sigma2");
  Rcpp::NumericVector u = model.slot("u") ;
  Rcpp::NumericMatrix data_sum(B, K);
  Rcpp::NumericMatrix n_hb(B, K);
  Rcpp::NumericMatrix sumu(B, K);
  Rcpp::NumericVector tauc(K);
  Rcpp::NumericMatrix iSigma2(B, K);
  Rcpp::NumericVector invs2;
  Rcpp::NumericVector theta(1);
  Rcpp::NumericVector zz;
  Rcpp::NumericVector mu=model.slot("mu");
  double post_prec;
  double tau_n;
  double mu_n;
  double w1;
  double w2;
  Rcpp::NumericVector tau=sqrt(tau2) ;
  Rcpp::NumericVector tmp(1);
  zz = model.slot("z");
  // make a B x k matrix of the counts
  n_hb = tableBatchZ(model);
  data_sum = compute_heavy_sums_batch(model);
  sumu = compute_u_sums_batch(model) ;
  mu=model.slot("mu") ;
  tau2_tilde = 1.0 / tau2;
  sigma2_tilde = 1.0 / sigma2;
  double df = getDf(model.slot("hyperparams")) ;
  double total = 0.0;
  double heavyn = 0.0;
  double heavy_mean = 0.0;
  for (int b = 0; b < B; ++b) {
    for(int k = 0; k < K; ++k) {
      //heavyn = n_hb(b, k) * sumu(b, k) / df;
      heavyn = sumu(b, k) / df;
      post_prec = tau2_tilde[k] + heavyn*1.0 * sigma2_tilde(b, k) ;
      if (post_prec == R_PosInf) {
        throw std::runtime_error("Bad simulation. Run again with different start.");
      }
      tau_n = sqrt(1.0/post_prec) ;
      w1 = (tau2_tilde[k])/post_prec ;
      w2 = (n_hb(b, k) * sigma2_tilde(b, k))/post_prec ;
      w1 = w1/(w1 + w2);
      w2 = 1-w1; 
      heavy_mean = data_sum(b, k) / heavyn / df;
      mu_n = w1*mu[k] + w2*heavy_mean ;
      theta[0] = thetastar(b, k);
      tmp = dnorm(theta, mu_n, tau_n, true);
      total += tmp[0];
    }
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::NumericVector marginal_theta_batch(Rcpp::S4 xmod) {
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
      model.slot("z") = update_z_batch(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      model.slot("theta") = update_theta_batch(model) ;
      model.slot("sigma2") = update_sigma2_batch(model) ;
      model.slot("mu") = update_mu_batch(model) ;
      model.slot("tau2") = update_tau2_batch(model) ;
      model.slot("sigma2.0") = update_sigma20_batch(model) ;
      model.slot("nu.0") = update_nu0_batch(model) ;
      model.slot("pi") = update_p_batch(model) ;
      model.slot("u") = Rcpp::rchisq(N, df) ;
      logp[s]=log_prob_theta(model, thetastar) ;
    }
    return logp;
}

// [[Rcpp::export]]
double log_prob_sigma2(Rcpp::S4 model, Rcpp::NumericMatrix sigma2star){
  Rcpp::RNGScope scope;
  int K = sigma2star.ncol();
  int B = sigma2star.nrow();
  Rcpp::NumericVector tmp(K);
  Rcpp::NumericMatrix prec(B, K);
  for (int k = 0; k < K; ++k) {
    prec(Rcpp::_, k) = 1.0 / sigma2star(Rcpp::_, k);
  }
  // calculate probability of sigma2
  Rcpp::NumericVector nu0 = model.slot("nu.0");
  Rcpp::NumericVector s20 = model.slot("sigma2.0");
  Rcpp::IntegerVector batch = model.slot("batch") ;
  Rcpp::IntegerVector ub = uniqueBatch(batch) ;
  Rcpp::IntegerVector zz = model.slot("z") ;
  Rcpp::NumericVector u = model.slot("u") ;
  Rcpp::NumericVector x = model.slot("data") ;
  Rcpp::NumericMatrix ss(B, K);
  Rcpp::NumericVector sigma2_n(1);
  Rcpp::NumericVector nu_n(1);
  Rcpp::NumericMatrix thetastar = model.slot("theta");
  int n = x.size();
  for (int i = 0; i < n; i++) {
    for (int b = 0; b < B; ++b) {
      if (batch[i] != ub[b]) {
        continue;
      }
      for (int k = 0; k < K; ++k) {
        if (zz[i] == k + 1) {
          ss(b, k) += u[i] * pow(x[i] - thetastar(b, k), 2);
        }
      }
    }
  }
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::NumericMatrix nn = tableBatchZ(model) ;
  Rcpp::NumericVector prec_typed(1);
  Rcpp::NumericVector shape(1);
  Rcpp::NumericVector rate(1);
  double total = 0.0;
  for (int b = 0; b < B; ++b) {
    for (int k = 0; k < K; ++k) {
      // calculate nu_n and sigma2_n
      nu_n = nu0 + nn(b, k);
      sigma2_n = 1.0 / nu_n * (nu0 * s20 + ss(b, k)/df);
      // calculate shape and rate
      shape = 0.5 * nu_n;
      rate = shape * sigma2_n;
      // calculate probability
      prec_typed[0] = prec(b, k);
      tmp = Rcpp::dgamma(prec_typed, shape[0], 1.0 / rate[0], true);
      total += tmp[0];
    }
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduced_sigma_batch(Rcpp::S4 xmod) {
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
  Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
  Rcpp::NumericMatrix sigma2star = clone(sigma2_);
  int K = thetastar.ncol();
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::NumericVector logp(S) ;
  //
  // Run reduced Gibbs    -- theta is fixed at modal ordinate
  //
  for (int s = 0; s < S; ++s) {
    model.slot("z") = update_z_batch(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    // THETA FIXED AT MODAL ORDINATE
    //model.slot("theta") = update_theta_batch(model) ;
    model.slot("sigma2") = update_sigma2_batch(model) ;
    model.slot("mu") = update_mu_batch(model) ;
    model.slot("tau2") = update_tau2_batch(model) ;
    model.slot("sigma2.0") = update_sigma20_batch(model) ;
    model.slot("nu.0") = update_nu0_batch(model) ;
    model.slot("pi") = update_p_batch(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
    logp[s]=log_prob_sigma2(model, sigma2star) ;
  }
  return logp ;
}

// [[Rcpp::export]]
double log_prob_pmix(Rcpp::S4 xmod, Rcpp::NumericVector pstar) {
  Rcpp::RNGScope scope;
  Rcpp::S4 model(xmod);
  Rcpp::S4 hypp = model.slot("hyperparams");
  Rcpp::NumericVector x = model.slot("data");
  Rcpp::S4 mp=model.slot("mcmc.params");
  int K = hypp.slot("k");
  Rcpp::NumericVector alpha = hypp.slot("alpha");
  Rcpp::NumericVector alpha_n(K);
  Rcpp::NumericVector logp(1);
  Rcpp::NumericVector z=model.slot("z") ;
  for (int k = 0 ; k < K; ++k) {
    alpha_n[k] = sum(z == k+1) + alpha[k];
  }
  logp = log_ddirichlet_(pstar, alpha_n);
  return logp[0] ;
}


// [[Rcpp::export]]
Rcpp::NumericVector reduced_pi_batch(Rcpp::S4 xmod) {
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
  Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
  Rcpp::NumericMatrix sigma2star = clone(sigma2_);
  Rcpp::NumericVector p_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
  Rcpp::NumericVector pstar = clone(p_);
  int K = thetastar.ncol();
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::NumericVector logp(S) ;
  model.slot("theta") = thetastar;
  model.slot("sigma2") = sigma2star;
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  for (int s = 0; s < S; ++s) {
    model.slot("z") = update_z_batch(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    // theta, sigma2 FIXED AT MODAL ORDINATES
    //model.slot("theta") = update_theta_batch(model) ;
    //model.slot("sigma2") = update_sigma2_batch(model) ;
    model.slot("mu") = update_mu_batch(model) ;
    model.slot("tau2") = update_tau2_batch(model) ;
    model.slot("sigma2.0") = update_sigma20_batch(model) ;
    model.slot("nu.0") = update_nu0_batch(model) ;
    model.slot("pi") = update_p_batch(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
    logp[s]=log_prob_pmix(model, pstar) ;
  }
  return logp;
}


// [[Rcpp::export]]
double log_prob_mu(Rcpp::S4 xmod, Rcpp::NumericVector mustar) {
  Rcpp::RNGScope scope;
  Rcpp::S4 model(xmod);
  Rcpp::S4 mcmcp = model.slot("mcmc.params");
  Rcpp::S4 chains = model.slot("mcmc.chains");
  Rcpp::S4 hypp = model.slot("hyperparams");
  Rcpp::List modes = model.slot("modes");
  Rcpp::NumericVector x = model.slot("data");
  Rcpp::IntegerVector batch = model.slot("batch") ;
  Rcpp::IntegerVector ub = uniqueBatch(batch) ;
  int S = mcmcp.slot("iter");
  int B = ub.size() ;
  // get hyperparameters
  int K = hypp.slot("k");
  double mu_0 = hypp.slot("mu.0");
  double tau2_0 = hypp.slot("tau2.0");
  double tau2_0_tilde = 1.0 / tau2_0;
  // get ordinal modes
  Rcpp::NumericVector p_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
  Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
  Rcpp::NumericVector mu_ = Rcpp::as<Rcpp::NumericVector>(modes["mu"]);
  Rcpp::NumericVector pstar = clone(p_);
  Rcpp::NumericMatrix thetastar = clone(theta_);
  // tau2
  Rcpp::NumericVector tau2=model.slot("tau2");
  Rcpp::NumericVector tau2_tilde(K);
  Rcpp::NumericVector tau2_B_tilde(K);

  for (int k = 0; k < K; ++k) {
    tau2_tilde[k] = 1.0 / tau2[k];
    tau2_B_tilde[k] = tau2_0_tilde + B * tau2_tilde[k];
  }
  // initalize variables for calculation of p_mu
  Rcpp::NumericMatrix n_b = tableBatchZ(model);
  Rcpp::NumericVector p_mu(S);
  // initialize variables
  Rcpp::NumericVector w1(K);
  Rcpp::NumericVector w2(K);
  Rcpp::NumericVector thetabar(K);
  Rcpp::NumericVector post_prec(K);
  Rcpp::NumericVector mu_k(K);
  Rcpp::NumericVector tau_k(K);
  double total = 0.0;
  for (int k = 0; k < K; ++k) {
    // calculate weights
    w1[k] = tau2_0_tilde / (tau2_0_tilde + B * tau2_tilde[k]);
    w2[k] = B * tau2_tilde[k] / (tau2_0_tilde + B * tau2_tilde[k]);
    // calculate thetabar
    double n_k = 0.0 ; // number of observations for component k
    double colsumtheta = 0.0;
    for (int b = 0; b < B; ++b) {
      colsumtheta += n_b(b, k) * thetastar(b, k);
      n_k += n_b(b, k);
    }
    thetabar[k] = colsumtheta / n_k;
    // calculate post_prec
    post_prec[k] = sqrt(1.0 / tau2_B_tilde[k]);
    // calculate mu_k
    mu_k[k] = w1[k] * mu_0 + w2[k] * thetabar[k];
    // calculate p_mu[s]
    Rcpp::NumericVector mu(1);
    mu[0] = mustar[k];
    total += Rcpp::dnorm(mu, mu_k[k], post_prec[k], true)[0];
  }
  return total ;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduced_mu_batch(Rcpp::S4 xmod) {
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
  Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
  Rcpp::NumericMatrix sigma2star = clone(sigma2_);
  Rcpp::NumericVector p_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
  Rcpp::NumericVector pstar = clone(p_);
  Rcpp::NumericVector mu_ = Rcpp::as<Rcpp::NumericVector>(modes["mu"]);
  Rcpp::NumericVector mustar = clone(mu_);
  int K = thetastar.ncol();
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::NumericVector logp(S) ;
  model.slot("theta") = thetastar;
  model.slot("sigma2") = sigma2star;
  model.slot("pi") = pstar;
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  for (int s = 0; s < S; ++s) {
    model.slot("z") = update_z_batch(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    // FIXED AT MODAL ORDINATES
    //model.slot("theta") = update_theta_batch(model) ;
    //model.slot("sigma2") = update_sigma2_batch(model) ;
    //model.slot("pi") = update_p_batch(model) ;
    model.slot("mu") = update_mu_batch(model) ;
    model.slot("tau2") = update_tau2_batch(model) ;
    model.slot("sigma2.0") = update_sigma20_batch(model) ;
    model.slot("nu.0") = update_nu0_batch(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
    logp[s]=log_prob_mu(model, mustar) ;
  }
  return logp;
}

// [[Rcpp::export]]
double log_prob_tau2(Rcpp::S4 xmod, Rcpp::NumericVector tau2star) {
  Rcpp::RNGScope scope;
  // get model and accessories
  Rcpp::S4 model(xmod);
  Rcpp::S4 mcmcp = model.slot("mcmc.params");
  Rcpp::S4 hypp = model.slot("hyperparams");
  Rcpp::List modes = model.slot("modes");
  // get modal ordinates
  Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
  Rcpp::NumericVector mu_ = Rcpp::as<Rcpp::NumericVector>(modes["mu"]);
  Rcpp::NumericVector mustar = clone(mu_);
  Rcpp::NumericMatrix thetastar = clone(theta_);

  // hyperparameters
  int K = hypp.slot("k");
  double m2_0 = hypp.slot("m2.0");
  double eta_0 = hypp.slot("eta.0");
  // batch and component parameters
  Rcpp::IntegerVector batch = model.slot("batch");
  Rcpp::IntegerVector ub = uniqueBatch(batch);
  int B = ub.size();
  // updated parameters
  double eta_B = eta_0 + B;
  Rcpp::NumericVector s2_k(K);
  Rcpp::NumericVector m2_k(K);
  double total = 0.0;
  for (int k = 0; k < K; ++k) {
    for (int b = 0; b < B; ++b) {
      s2_k[k] += pow(thetastar(b, k) - mustar[k], 2) ;
    }
    m2_k[k] = 1.0 / eta_B * (eta_0 * m2_0 + s2_k[k]);
    Rcpp::NumericVector tau2star_typed(1);
    tau2star_typed[0] = 1.0 / tau2star[k];
    total += Rcpp::dgamma(tau2star_typed, 0.5 * eta_B,
                          1.0 / (0.5 * eta_B * m2_k[k]), true)[0];
  }
  return total ;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduced_tau_batch(Rcpp::S4 xmod) {
  Rcpp::RNGScope scope;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 params=model.slot("mcmc.params");
  Rcpp::S4 hypp=model.slot("hyperparams");
  int S = params.slot("iter");
  int K = hypp.slot("k");
  Rcpp::NumericVector y = model.slot("data");
  int N=y.size();
  Rcpp::List modes = model.slot("modes");
  Rcpp::NumericVector tau2_ = Rcpp::as<Rcpp::NumericVector>(modes["tau2"]);
  Rcpp::NumericVector tau2star = clone(tau2_);
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::NumericVector logp(S) ;
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  for (int s = 0; s < S; ++s) {
    model.slot("z") = update_z_batch(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    // FIXED AT MODAL ORDINATES
    //model.slot("theta") = update_theta_batch(model) ;
    //model.slot("sigma2") = update_sigma2_batch(model) ;
    //model.slot("pi") = update_p_batch(model) ;
    //model.slot("mu") = update_mu_batch(model) ;
    model.slot("tau2") = update_tau2_batch(model) ;
    model.slot("sigma2.0") = update_sigma20_batch(model) ;
    model.slot("nu.0") = update_nu0_batch(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
    logp[s]=log_prob_tau2(model, tau2star) ;
  }
  return logp;
}

// [[Rcpp::export]]
double log_prob_nu0(Rcpp::S4 xmod, int nu0star) {
  Rcpp::RNGScope scope;
  // get model and accessories
  Rcpp::S4 model(xmod);
  Rcpp::S4 mcmcp = model.slot("mcmc.params");
  Rcpp::S4 hypp = model.slot("hyperparams");
  Rcpp::List modes = model.slot("modes");
  // get modal ordinates
  Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
  Rcpp::NumericMatrix sigma2star = clone(sigma2_);
  // hyperparameters
  int K = hypp.slot("k");
  double betas = hypp.slot("beta");
  int B = sigma2star.nrow() ;
  Rcpp::NumericVector s20=model.slot("sigma2.0") ;
  //
  // compute p(nu0*, ) from *normalized* probabilities
  //
  Rcpp::NumericVector d(100) ;  // 100 is the maximum allowed value for nu_0
  Rcpp::NumericVector lpnu0(100);
  double prec = 0.0;
  double lprec = 0.0;
  d = Rcpp::seq_len(100);
  for (int b = 0; b < B; ++b) {
    for (int k = 0; k < K; ++k) {
      prec += 1.0 / sigma2star(b, k);
      lprec += log(1.0 / sigma2star(b, k));
    }
  }
  Rcpp::NumericVector y1(100);
  Rcpp::NumericVector y2(100);
  Rcpp::NumericVector y3(100);
  y1 = B * K * (0.5 * d * log(s20[0] * 0.5 * d) - lgamma(d * 0.5));
  y2 = (0.5 * d - 1.0) * lprec;
  y3 = d * (betas + 0.5 * s20[0] * prec);
  lpnu0 = y1 + y2 - y3;
  Rcpp::NumericVector prob(100);
  prob = exp(lpnu0); // - maxprob);
  prob = prob / sum(prob);  // this is now normalized
  //p_nu0[s] = prob[nu0star];
  double log_p = log(prob[nu0star]) ;
  return log_p ;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduced_nu0_batch(Rcpp::S4 xmod) {
  Rcpp::RNGScope scope;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 params=model.slot("mcmc.params");
  Rcpp::S4 hypp=model.slot("hyperparams");
  int S = params.slot("iter");
  int K = hypp.slot("k");
  Rcpp::NumericVector y = model.slot("data");
  int N=y.size();
  Rcpp::List modes = model.slot("modes");
  Rcpp::IntegerVector nu0_ = Rcpp::as<Rcpp::IntegerVector>(modes["nu0"]);
  Rcpp::IntegerVector nu0star = clone(nu0_);
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::NumericVector logp(S) ;
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  for (int s = 0; s < S; ++s) {
    model.slot("z") = update_z_batch(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    // FIXED AT MODAL ORDINATES
    //model.slot("theta") = update_theta_batch(model) ;
    //model.slot("sigma2") = update_sigma2_batch(model) ;
    //model.slot("pi") = update_p_batch(model) ;
    //model.slot("mu") = update_mu_batch(model) ;
    //model.slot("tau2") = update_tau2_batch(model) ;
    model.slot("sigma2.0") = update_sigma20_batch(model) ;
    model.slot("nu.0") = update_nu0_batch(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
    logp[s]=log_prob_nu0(model, nu0star[0]) ;
  }
  return logp;
}

// [[Rcpp::export]]
double log_prob_sigma2_0(Rcpp::S4 xmod, Rcpp::NumericVector s20star) {
  Rcpp::RNGScope scope;
  // get model and accessories
  Rcpp::S4 model(xmod);
  Rcpp::S4 mcmcp = model.slot("mcmc.params");
  Rcpp::S4 chains = model.slot("mcmc.chains");
  Rcpp::S4 hypp = model.slot("hyperparams");
  Rcpp::List modes = model.slot("modes");
  Rcpp::IntegerVector batch = model.slot("batch");
  Rcpp::IntegerVector ub = uniqueBatch(batch);
  int B = ub.size();
  // get ordinal modes.
  Rcpp::IntegerVector nu0_ = Rcpp::as<Rcpp::IntegerVector>(modes["nu0"]);
  Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
  Rcpp::NumericMatrix sigma2star = clone(sigma2_);
  double nu0star = clone(nu0_)[0];
  // get hyperparameters
  int K = hypp.slot("k");
  double a = hypp.slot("a");
  double b = hypp.slot("b");
  // calculate a_k, b_k
  double prec = 0.0;
  for (int b = 0; b < B; ++b) {
    for (int k = 0; k < K; ++k) {
      prec += 1.0 / sigma2star(b, k);
    }
  }
  double a_k = a + 0.5 * K * B * nu0star;
  double b_k = b + 0.5 * nu0star * prec;
  Rcpp::NumericVector logp_s20 = dgamma(s20star, a_k, 1.0 / b_k, true);
  return logp_s20[0] ;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduced_s20_batch(Rcpp::S4 xmod) {
  Rcpp::RNGScope scope;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 params=model.slot("mcmc.params");
  Rcpp::S4 hypp=model.slot("hyperparams");
  int S = params.slot("iter");
  int K = hypp.slot("k");
  Rcpp::NumericVector y = model.slot("data");
  int N=y.size();
  Rcpp::List modes = model.slot("modes");
  Rcpp::NumericVector s20_ = Rcpp::as<Rcpp::NumericVector>(modes["sigma2.0"]);
  Rcpp::NumericVector s20star = clone(s20_);
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::NumericVector logp(S) ;
  //
  // Run reduced Gibbs:
  //   -- theta is fixed at modal ordinate
  //   -- sigma2 is fixed at modal ordinate
  for (int s = 0; s < S; ++s) {
    model.slot("z") = update_z_batch(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    // FIXED AT MODAL ORDINATES
    //model.slot("theta") = update_theta_batch(model) ;
    //model.slot("sigma2") = update_sigma2_batch(model) ;
    //model.slot("pi") = update_p_batch(model) ;
    //model.slot("mu") = update_mu_batch(model) ;
    //model.slot("tau2") = update_tau2_batch(model) ;
    //model.slot("nu.0") = update_nu0_batch(model) ;
    model.slot("sigma2.0") = update_sigma20_batch(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
    logp[s]=log_prob_sigma2_0(model, s20star) ;
  }
  return logp;
}
