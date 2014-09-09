#ifndef _miscfunctions_H
#define _miscfunctions_H
#include <RcppArmadillo.h>

void rdirichlet(Rcpp::NumericVector, Rcpp::NumericVector);
double cons_normal(double, double, double, double);
double trunc_norm(double mean, double sd);
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
Rcpp::NumericVector dsn(Rcpp::NumericVector r, double xi,
        double omega, double alpha);

// stuff from simulating from truncated normal
double norm_rs(double a, double b);
double half_norm_rs(double a, double b);
double unif_rs(double a, double b);
double exp_rs(double a, double b);
double rnorm_trunc (double mu, double sigma, double lower, double upper);
#endif
