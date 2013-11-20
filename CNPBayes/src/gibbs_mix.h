#ifndef _MixtureModel_RCPP_GIBBS_MIX_H
#define _MixtureModel_RCPP_GIBBS_MIX_H
#include <Rcpp.h>
RcppExport SEXP gibbs_mix(SEXP r, SEXP means, SEXP precs, SEXP P, SEXP Z,
        SEXP nu0, SEXP mu0, SEXP kappa0, SEXP alpha,
        SEXP tau20, SEXP sigma20, SEXP rbar, SEXP s2,
        SEXP nn, SEXP delta, SEXP burnin) ;

#endif
