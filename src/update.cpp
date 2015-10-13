#include "update.h"
#include "miscfunctions.h"
#include "Rmath.h"

using namespace Rcpp;
// create accessor functions for s4 slots
// Update theta in marginal model
Rcpp::NumericVector update_theta(Rcpp::S4 xmod) {
    RNGScope scope ;
    Rcpp::S4 model(xmod) ;
    NumericVector theta = model.slot("theta") ;
    double tau2 = model.slot("tau2") ;
    double tau2_tilde = 1/tau2 ;
    NumericVector sigma2 = model.slot("sigma2") ;
    NumericVector data_mean = model.slot("data.mean") ;
    NumericVector sigma2_tilde = 1.0/sigma2 ;
    IntegerVector z = model.slot("z");
    int K = getK(model.slot("hyperparams"));
    double mu_prior = model.slot("mu");
    //
    // Initialize nn, vector of component-wise sample size
    //
    IntegerVector nn = tableZ(K, z) ;
    double post_prec;
    double tau_n;
    double mu_n;
    double w1;
    double w2;
    NumericVector thetas(K);

    for(int k = 0; k < K; ++k) {
        post_prec = tau2_tilde + sigma2_tilde[k] * nn[k];
        if (post_prec == R_PosInf) {
            throw std::runtime_error("Bad simulation. Run again with different start.");
        }
        tau_n = sqrt(1/post_prec);
        w1 = tau2_tilde/post_prec;
        w2 = nn[k]*sigma2_tilde[k]/post_prec;
        mu_n = w1*mu_prior + w2*data_mean[k];
        thetas[k] = as<double>(rnorm(1, mu_n, tau_n));
    }
    return thetas;
}

// Accessors
Rcpp::IntegerVector getZ(Rcpp::S4 model) {
    IntegerVector z = model.slot("z");
    return z;
}

Rcpp::NumericVector getData(Rcpp::S4 model) {
    NumericVector y = model.slot("data");
    return y;
}

// FUNCTIONS FOR ACCESSING HYPERPARAMETERS
// [[Rcpp::export]]
int getK(Rcpp::S4 hyperparams) {
  int k = hyperparams.slot("k");
  return k;
}

Rcpp::NumericVector getMu(Rcpp::S4 hyperparams) {
  NumericVector mu = hyperparams.slot("mu");
  return mu;
}

Rcpp::NumericVector getTau2(Rcpp::S4 hyperparams) {
    NumericVector tau2 = hyperparams.slot("tau2");
    return tau2;
}

Rcpp::IntegerVector getAlpha(Rcpp::S4 hyperparams) {
    IntegerVector alpha = hyperparams.slot("alpha");
    return alpha;
}

Rcpp::LogicalVector nonZeroCopynumber(Rcpp::IntegerVector z) {
//nonZeroCopynumber <- function(object) as.integer(as.integer(z(object)) > 1)
 LogicalVector nz = z > 1;
 return nz;
}

