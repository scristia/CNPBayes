#ifndef _update_H
#define _update_H
#include <Rcpp.h>

Rcpp::NumericVector update_theta(Rcpp::S4 xmod);

// Access model values
Rcpp::IntegerVector getZ(Rcpp::S4 model);
Rcpp::NumericVector getData(Rcpp::S4 model);

// Access hyperparameters
int getK(Rcpp::S4 hyperparams);
Rcpp::NumericVector getMu(Rcpp::S4 hyperparams);
Rcpp::NumericVector getTau2(Rcpp::S4 hyperparams);
Rcpp::IntegerVector getAlpha(Rcpp::S4 hyperparams);

Rcpp::LogicalVector nonZeroCopynumber(Rcpp::IntegerVector z);
#endif
