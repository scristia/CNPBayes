#ifndef _update_batch_H
#define _update_batch_H

Rcpp::IntegerVector uniqueBatch(Rcpp::IntegerVector x);
Rcpp::NumericMatrix tableBatchZ(Rcpp::S4 xmod);
Rcpp::NumericVector compute_loglik_batch(Rcpp::S4 xmod);

Rcpp::NumericVector update_mu_batch(Rcpp::S4 xmod);
Rcpp::NumericVector update_tau2_batch(Rcpp::S4 xmod);
Rcpp::NumericVector update_sigma20_batch(Rcpp::S4 xmod);
Rcpp::NumericVector update_nu0_batch(Rcpp::S4 xmod);
Rcpp::NumericMatrix update_multinomialPr_batch(Rcpp::S4 xmod);
Rcpp::NumericVector update_p_batch(Rcpp::S4 xmod);
Rcpp::IntegerVector update_z_batch(Rcpp::S4 xmod);

Rcpp::NumericMatrix compute_means_batch(Rcpp::S4 xmod);
Rcpp::NumericMatrix compute_vars_batch(Rcpp::S4 xmod);
Rcpp::NumericMatrix compute_prec_batch(Rcpp::S4 xmod);
Rcpp::NumericVector compute_logprior_batch(Rcpp::S4 xmod);

Rcpp::NumericVector stageTwoLogLikBatch(Rcpp::S4 xmod);

Rcpp::NumericMatrix update_theta_batch(Rcpp::S4 xmod);
Rcpp::NumericMatrix update_sigma2_batch(Rcpp::S4 xmod);

Rcpp::IntegerVector order_(Rcpp::NumericVector x);

Rcpp::IntegerMatrix update_probz_batch(Rcpp::S4 xmod);

Rcpp::S4 mcmc_batch_burnin(Rcpp::S4 xmod, Rcpp::S4 mcmcp);
Rcpp::S4 mcmc_batch(Rcpp::S4 object, Rcpp::S4 mcmcp);

Rcpp::NumericVector p_theta_batch(Rcpp::S4 xmod);
Rcpp::NumericVector p_theta_zfixed_batch(Rcpp::S4 xmod);

Rcpp::S4 simulate_z_reduced1_batch(Rcpp::S4 object);

Rcpp::S4 reduced_z_theta_fixed(Rcpp::S4 object);
Rcpp::S4 simulate_z_reduced2_batch(Rcpp::S4 object);
Rcpp::NumericVector p_sigma2_batch(Rcpp::S4 xmod);

#endif
