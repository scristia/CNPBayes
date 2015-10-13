#ifndef _update_marginal_H
#define _update_marginal_H

Rcpp::NumericVector loglik(Rcpp::S4 xmod); 
Rcpp::NumericVector log_ddirichlet_(Rcpp::NumericVector x_, 
                                    Rcpp::NumericVector alpha_);
Rcpp::NumericVector stageTwoLogLik(Rcpp::S4 xmod);

Rcpp::NumericVector update_sigma2(Rcpp::S4 xmod);
Rcpp::NumericVector update_mu(Rcpp::S4 xmod);
Rcpp::NumericVector update_tau2(Rcpp::S4 xmod);
Rcpp::NumericVector update_sigma2_0(Rcpp::S4 xmod);
Rcpp::NumericVector update_nu0(Rcpp::S4 xmod);
Rcpp::NumericVector update_p(Rcpp::S4 xmod);
Rcpp::NumericMatrix update_multinomialPr(Rcpp::S4 xmod);
Rcpp::IntegerVector update_z(Rcpp::S4 xmod);

Rcpp::NumericVector compute_means(Rcpp::S4 xmod);
Rcpp::NumericVector compute_vars(Rcpp::S4 xmod);
Rcpp::NumericVector compute_prec(Rcpp::S4 xmod);
Rcpp::NumericVector compute_logprior(Rcpp::S4 xmod);

Rcpp::NumericVector update_sigma2(Rcpp::S4 xmod);
Rcpp::IntegerVector ordertheta_(Rcpp::NumericVector x);
Rcpp::IntegerMatrix compute_probz(Rcpp::S4 xmod);

Rcpp::S4 mcmc_marginal_burnin(Rcpp::S4 xmod, Rcpp::S4 mcmcp);
Rcpp::S4 mcmc_marginal(Rcpp::S4 object, Rcpp::S4 mcmcp);

#endif
