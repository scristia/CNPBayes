#ifndef _update_red_batch_H
#define _update_red_batch_H

Rcpp::NumericMatrix toMatrix(Rcpp::NumericVector x, int NR, int NC);
Rcpp::NumericVector marginal_theta_batch(Rcpp::S4 xmod); 
Rcpp::NumericVector p_theta_zpermuted_batch(Rcpp::S4 xmod);

Rcpp::NumericVector marginal_sigma2_batch(Rcpp::S4 xmod, Rcpp::S4 mcmcp);

Rcpp::S4 simulate_z_reduced1_batch(Rcpp::S4 object);
Rcpp::S4 simulate_z_reduced2_batch(Rcpp::S4 object);
Rcpp::S4 permutedz_reduced1_batch(Rcpp::S4 xmod);
Rcpp::S4 permutedz_reduced2_batch(Rcpp::S4 xmod);

Rcpp::NumericVector p_pmix_reduced_batch(Rcpp::S4 xmod);

Rcpp::S4 reduced_sigma_batch(Rcpp::S4 xmod);
Rcpp::NumericVector p_sigma_reduced_batch(Rcpp::S4 xmod);

Rcpp::S4 reduced_pi_batch(Rcpp::S4 xmod);

Rcpp::S4 reduced_mu_batch(Rcpp::S4 xmod);
Rcpp::NumericVector p_mu_reduced_batch(Rcpp::S4 xmod);

Rcpp::S4 reduced_tau_batch(Rcpp::S4 xmod);
Rcpp::NumericVector p_tau_reduced_batch(Rcpp::S4 xmod);

Rcpp::S4 reduced_nu0_batch(Rcpp::S4 xmod);
Rcpp::NumericVector p_nu0_reduced_batch(Rcpp::S4 xmod);

Rcpp::S4 reduced_s20_batch(Rcpp::S4 xmod);
Rcpp::NumericVector p_s20_reduced_batch(Rcpp::S4 xmod);

#endif
