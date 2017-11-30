#ifndef _miscfunctions_H
#define _miscfunctions_H
#include <Rcpp.h>


// Access model values
Rcpp::IntegerVector getZ(Rcpp::S4 model);
Rcpp::NumericVector getData(Rcpp::S4 model);

// Access hyperparameters
int getK(Rcpp::S4 hyperparams);
double getDf(Rcpp::S4 hyperparams);
Rcpp::NumericVector getMu(Rcpp::S4 hyperparams);
Rcpp::NumericVector getTau2(Rcpp::S4 hyperparams);
Rcpp::IntegerVector getAlpha(Rcpp::S4 hyperparams);

Rcpp::LogicalVector nonZeroCopynumber(Rcpp::IntegerVector z);

void rdirichlet(Rcpp::NumericVector, Rcpp::NumericVector);
double cons_normal(double, double, double, double);
double trunc_norm(double mean, double sd);
Rcpp::NumericVector dsn(Rcpp::NumericVector r, double xi,
        double omega, double alpha);
Rcpp::IntegerMatrix rMultinom(Rcpp::NumericMatrix probs, int m); 

// stuff from simulating from truncated normal
double norm_rs(double a, double b);
double half_norm_rs(double a, double b);
double unif_rs(double a, double b);
double exp_rs(double a, double b);
double rnorm_trunc (double mu, double sigma, double lower, double upper);

// Simulate from scale-location t-distribution
Rcpp::NumericVector dlocScale_t(Rcpp::NumericVector x, double df,
        double mu, double sigma);

// Rob
Rcpp::IntegerVector tableZ(int K, Rcpp::IntegerVector z) ;

Rcpp::NumericVector compute_u_sums(Rcpp::S4 xmod) ;
Rcpp::NumericVector compute_heavy_sums(Rcpp::S4 object) ;
Rcpp::NumericVector compute_heavy_sums(Rcpp::S4 object) ;
Rcpp::NumericVector compute_heavy_means(Rcpp::S4 xmod) ;
Rcpp::NumericVector dlocScale_t(Rcpp::NumericVector x, double df, double mu, double sigma);
#endif
