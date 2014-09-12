#ifndef _MixtureModel_RCPP_GRID_H
#define _MixtureModel_RCPP_GRID_H
/* #include <Rcpp.h> */
#include <RcppArmadillo.h>

RcppExport SEXP rcpp_discrete(SEXP grid, SEXP meangrid, SEXP precgrid,
                              SEXP dtheta, SEXP dprec, SEXP y) ;
#endif
