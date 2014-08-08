#ifndef _miscfunctions_H
#define _miscfunctions_H
#include <Rcpp.h>

void rdirichlet(Rcpp::NumericVector, Rcpp::NumericVector);
double cons_normal(double, double, double, double);

#endif
