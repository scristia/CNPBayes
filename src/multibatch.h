#ifndef _multibatch_H
#define _multibatch_H

//Rcpp::IntegerVector unique_batch(Rcpp::IntegerVector x);
//Rcpp::NumericMatrix tableBatchZ(Rcpp::S4 xmod);
Rcpp::NumericVector compute_loglik(Rcpp::S4 xmod);

Rcpp::NumericVector update_mu(Rcpp::S4 xmod);
Rcpp::NumericVector update_tau2(Rcpp::S4 xmod);
Rcpp::NumericVector update_sigma20(Rcpp::S4 xmod);
Rcpp::NumericVector update_nu0(Rcpp::S4 xmod);
Rcpp::NumericMatrix update_multinomialPr(Rcpp::S4 xmod);
Rcpp::NumericVector update_p(Rcpp::S4 xmod);
Rcpp::IntegerVector update_z(Rcpp::S4 xmod);

Rcpp::NumericMatrix compute_means(Rcpp::S4 xmod);
Rcpp::NumericMatrix compute_vars(Rcpp::S4 xmod);
Rcpp::NumericMatrix compute_prec(Rcpp::S4 xmod);
Rcpp::NumericVector compute_logprior(Rcpp::S4 xmod);

Rcpp::NumericVector stageTwoLogLikBatch(Rcpp::S4 xmod);

Rcpp::NumericMatrix update_theta(Rcpp::S4 xmod);
Rcpp::NumericMatrix update_sigma2(Rcpp::S4 xmod);

Rcpp::IntegerVector order_(Rcpp::NumericVector x);

Rcpp::IntegerMatrix update_probz(Rcpp::S4 xmod);

Rcpp::S4 cpp_burnin(Rcpp::S4 xmod, Rcpp::S4 mcmcp);
Rcpp::S4 cpp_mcmc(Rcpp::S4 object, Rcpp::S4 mcmcp);

#endif
