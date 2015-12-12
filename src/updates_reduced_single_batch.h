#ifndef _update_red_marginal_H
#define _update_red_marginal_H

Rcpp::NumericVector marginal_theta(Rcpp::S4 xmod);
Rcpp::NumericVector p_theta_zpermuted(Rcpp::S4 xmod);

Rcpp::NumericVector marginal_sigma2(Rcpp::S4 xmod, Rcpp::S4 mcmcp);

Rcpp::S4 simulate_z_reduced1(Rcpp::S4 object);
Rcpp::S4 simulate_z_reduced2(Rcpp::S4 object);
Rcpp::S4 permutedz_reduced1(Rcpp::S4 xmod);
Rcpp::S4 permutedz_reduced2(Rcpp::S4 xmod);

Rcpp::NumericVector p_pmix_reduced(Rcpp::S4 xmod);

Rcpp::S4 reduced_sigma(Rcpp::S4 xmod);
Rcpp::NumericVector p_sigma_reduced(Rcpp::S4 xmod);

Rcpp::S4 reduced_pi(Rcpp::S4 xmod);

Rcpp::S4 reduced_mu(Rcpp::S4 xmod);
Rcpp::NumericVector p_mu_reduced(Rcpp::S4 xmod);

Rcpp::S4 reduced_tau(Rcpp::S4 xmod);
Rcpp::NumericVector p_tau_reduced(Rcpp::S4 xmod);

Rcpp::S4 reduced_nu0(Rcpp::S4 xmod);
Rcpp::NumericVector p_nu0_reduced(Rcpp::S4 xmod);

Rcpp::S4 reduced_s20(Rcpp::S4 xmod);
Rcpp::NumericVector p_s20_reduced(Rcpp::S4 xmod);

#endif
