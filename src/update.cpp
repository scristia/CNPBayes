#include "update.h"
#include "miscfunctions.h"
#include "Rmath.h"

using namespace Rcpp;
// create accessor functions for s4 slots
// Update theta in marginal model
RcppExport SEXP update_theta(SEXP xmod, SEXP constraint) {
    // Rcpp::RNGScope scope;
    // initialize objects that are passed from R
    RNGScope scope;

    Rcpp::S4 model(xmod);
    NumericVector theta = model.slot("theta");
    NumericVector tau2 = model.slot("data.prec");
    NumericVector tau2_tilde = 1/tau2;
    NumericVector sigma2 = model.slot("sigma2");
    NumericVector data_mean = model.slot("data.mean");
    NumericVector sigma2_tilde = 1/sigma2;
    IntegerVector z = model.slot("z");
    int K = theta.size();
    //
    // Initialize nn, vector of component-wise sample size
    IntegerVector nn(K);
    for(int k = 0; k < K; k++) nn[k] = sum(z == k+1);

    NumericVector denom;
    for(int k = 0; k < K; k++)
        denom[k] = tau2_tilde[k] + sigma2_tilde[k] * nn[k];
    NumericVector tau2_n = 1/denom;
//    NumericVector denom = tau2_tilde + nn*sigma2_tilde;
    NumericVector w1 = tau2_tilde/denom;
    NumericVector w2 = as<NumericVector>(nn)*sigma2_tilde/denom;

    // Create accessor function
    NumericVector mu_prior = data_mean;
    NumericVector mu_n = w1*mu_prior + w2*data_mean;
    NumericVector thetas;
    for(int k = 0; k < K; k++)
        thetas[k] = as<double>(rnorm(1, mu_prior[k], sqrt(tau2_n[k])));
    return thetas;
/*
  tau2.n <- 1/tau2.n.tilde
  denom <- tau2.tilde + n.h*sigma2.tilde
  w1 <- tau2.tilde/denom
  w2 <- n.h*sigma2.tilde/denom
  mu.n <- w1*mu + w2*data.mean
  k <- length(tau2.n)
  thetas <- rnorm(k, mu.n, sqrt(tau2.n))
  */
    //Rcpp::Rcout << z[0] << " " << z[50];
    //return theta;
}
