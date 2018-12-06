#ifndef _multibatch_reduced_H
#define _multibatch_reduced_H

double log_prob_theta(Rcpp::S4 xmod, Rcpp::NumericMatrix thetastar);
Rcpp::NumericVector marginal_theta_batch(Rcpp::S4 xmod);
double log_prob_sigma2(Rcpp::S4 model, Rcpp::NumericMatrix sigma2star);
Rcpp::NumericVector reduced_sigma_batch(Rcpp::S4 xmod) ;
double log_prob_pmix(Rcpp::S4 xmod, Rcpp::NumericVector pstar);
Rcpp::NumericVector reduced_pi_batch(Rcpp::S4 xmod) ;
double log_prob_pmix2(Rcpp::S4 xmod, Rcpp::NumericMatrix pstar);
Rcpp::NumericVector reduced_pi_batch2(Rcpp::S4 xmod) ;
double log_prob_mu(Rcpp::S4 xmod, Rcpp::NumericVector mustar) ;
Rcpp::NumericVector reduced_mu_batch(Rcpp::S4 xmod) ;
double log_prob_tau2(Rcpp::S4 xmod);
Rcpp::NumericVector reduced_tau_batch(Rcpp::S4 xmod) ;
double log_prob_nu0(Rcpp::S4 xmod, int nu0star);
Rcpp::NumericVector reduced_nu0_batch(Rcpp::S4 xmod);
double log_prob_s20(Rcpp::S4 xmod);
Rcpp::NumericVector reduced_s20_batch(Rcpp::S4 xmod) ;

#endif
