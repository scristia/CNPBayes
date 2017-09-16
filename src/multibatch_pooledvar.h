#ifndef _multibatch_pooledvar_H
#define _multibatch_pooledvar_H

Rcpp::NumericVector loglik_multibatch_pvar(Rcpp::S4 xmod);
Rcpp::NumericVector sigma20_multibatch_pvar(Rcpp::S4 xmod);
Rcpp::NumericVector nu0_multibatch_pvar(Rcpp::S4 xmod);
Rcpp::NumericMatrix multinomialPr_multibatch_pvar(Rcpp::S4 xmod);
Rcpp::IntegerVector z_multibatch_pvar(Rcpp::S4 xmod);
Rcpp::NumericVector stagetwo_multibatch_pvar(Rcpp::S4 xmod);
Rcpp::NumericMatrix theta_multibatch_pvar(Rcpp::S4 xmod);
Rcpp::NumericVector sigma2_multibatch_pvar(Rcpp::S4 xmod);
Rcpp::S4 burnin_multibatch_pvar(Rcpp::S4 xmod, Rcpp::S4 mcmcp);
Rcpp::S4 mcmc_multibatch_pvar(Rcpp::S4 object, Rcpp::S4 mcmcp);

#endif
