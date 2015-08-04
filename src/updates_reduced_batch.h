#ifndef _update_red_batch_H
#define _update_red_batch_H

Rcpp::NumericMatrix toMatrix(Rcpp::NumericVector x, int NR, int NC);
Rcpp::NumericVector marginal_theta_batch(Rcpp::S4 xmod); 

#endif
