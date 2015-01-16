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
    int K = getK(model.slot("hyperparams"));
    NumericVector mu_prior = getMu(model.slot("hyperparams"));
    //
    // Initialize nn, vector of component-wise sample size
    IntegerVector nn(K);
    for(int k = 0; k < K; k++) nn[k] = sum(z == k+1);

    NumericVector denom = tau2_tilde + sigma2_tilde * as<NumericVector>(nn);

    NumericVector tau2_n = 1/denom;
//    NumericVector denom = tau2_tilde + nn*sigma2_tilde;
    NumericVector w1 = tau2_tilde/denom;
    NumericVector w2 = as<NumericVector>(nn)*sigma2_tilde/denom;

    NumericVector mu_n = w1*mu_prior + w2*data_mean;
    NumericVector thetas(K);
    for(int k = 0; k < K; k++) {
        thetas[k] = as<double>(rnorm(1, mu_n[k], sqrt(tau2_n[k])));
        //Rcpp::Rcout << mu_n[k] << " -- " << tau2_n[k] << " -- " << thetas[k] << "\n";
    }
    return thetas;
}

RcppExport SEXP update_sigma2(SEXP xmod) {
    // Rcpp::RNGScope scope;
    // initialize objects that are passed from R
    RNGScope scope;
    Rcpp::S4 model(xmod);
    int K = getK(model.slot("hyperparams"));
    IntegerVector z = getZ(model);
    IntegerVector zu = unique(z);
    return zu;
}

// Accessors
Rcpp::IntegerVector getZ(Rcpp::S4 model) {
    IntegerVector z = model.slot("z");
    return z;
}
// FUNCTIONS FOR ACCESSING HYPERPARAMETERS
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

