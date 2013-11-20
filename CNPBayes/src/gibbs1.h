#ifndef _MixtureModel_RCPP_GIBBS_H
#define _MixtureModel_RCPP_GIBBS_H
#include <Rcpp.h>
RcppExport SEXP rcpp_gibbs(SEXP phi, SEXP nu0, SEXP mu0, SEXP tau20,
                           SEXP sigma20, SEXP ybar, SEXP s2, SEXP n) ;

#endif
