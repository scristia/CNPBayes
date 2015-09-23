#ifndef _singlebatch_pooledvar_H
#define _singlebatch_pooledvar_H

Rcpp::NumericVector theta_pooled(Rcpp::S4 xmod);
Rcpp::NumericVector sigma2_pooled(Rcpp::S4 xmod);
Rcpp::NumericVector nu0_pooled(Rcpp::S4 xmod);
Rcpp::NumericVector sigma2_0_pooled(Rcpp::S4 xmod);
Rcpp::IntegerVector z_pooled(Rcpp::S4 xmod);

#endif
