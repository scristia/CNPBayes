#include "miscfunctions.h" // for rdirichlet, tableZ, ...
#include "updates_batch.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix toMatrix(Rcpp::NumericVector x, int NR, int NC) {
  int nrow = NR ;
  int ncol = NC ;
  NumericMatrix Y(NR, NC) ;
  int iter = 0 ;
  for(int j = 0; j < NC; ++j){
    for(int i = 0; i < NR; ++i){
      Y(i, j) = x[iter] ;
      iter += 1 ;
    }
  }
  return Y ;
}

// Only one done by RS.
// [[Rcpp::export]]
Rcpp::NumericVector marginal_theta_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model_(xmod) ;
  Rcpp::S4 model = clone(model_) ;
  Rcpp::S4 params=model.slot("mcmc.params") ;
  Rcpp::S4 chains(model.slot("mcmc.chains")) ;  
  int S = params.slot("iter") ;
  List modes = model.slot("modes") ;
  NumericMatrix theta_ = as<NumericMatrix>(modes["theta"]) ;
  NumericMatrix thetastar=clone(theta_) ;
  int K = thetastar.ncol() ;
  NumericVector p_theta(S) ;
  NumericMatrix muc = chains.slot("mu") ;
  NumericMatrix tau2c = chains.slot("tau2") ;
  NumericMatrix sigma2 = chains.slot("sigma2") ;
  //NumericVector tauc = sqrt(tau2c) ;
  int B = thetastar.nrow() ;
  NumericVector tmp(1) ;

  IntegerMatrix Z = chains.slot("z") ;
  IntegerVector zz ;

  NumericVector tau2_tilde(K) ;
  NumericVector sigma2_tilde(K) ;

  // this should be updated for each iteration
  NumericMatrix data_mean(B,K) ;
  IntegerVector nn(K) ;
  double post_prec;
  double tau_n;
  double mu_n;
  double w1;
  double w2;
  double prod ;
  NumericVector tauc(K) ;
  NumericMatrix iSigma2(B,K) ;
  NumericVector invs2 ;
  NumericVector theta(1) ;

  for(int s=0; s < S; ++s){
    tauc = sqrt(tau2c(s, _)) ;
    zz = Z(s, _) ;
    model.slot("z") = zz ;
    nn = tableZ(K, zz) ;
    data_mean = compute_means_batch(model) ;
    tau2_tilde = 1/tau2c(s, _) ;
    invs2 = 1.0/sigma2(s, _) ;  // this is a vector of length B*K
    sigma2_tilde = as<Rcpp::NumericVector>(toMatrix(invs2, B, K));
    //tmp = dnorm(thetastar, muc[s], tauc[s]) ;
    double prod = 1.0 ;
    for(int k = 0; k < K; ++k) {
      for(int b = 0; b < B; ++b){
        post_prec = tau2_tilde[k] + sigma2_tilde(b, k) * nn[k] ;
        tau_n = sqrt(1/post_prec) ;
        w1 = tau2_tilde[k]/post_prec;
        w2 = nn[k] * sigma2_tilde(b, k)/post_prec;
        mu_n = w1*muc(s, k) + w2*data_mean(b, k);
        theta = thetastar(b, k) ;
        tmp = dnorm(theta, mu_n, tau_n) ;
        prod = prod * tmp[0] ;
      }
    }
    p_theta[s] = prod ;
  }
  return p_theta ;
}

