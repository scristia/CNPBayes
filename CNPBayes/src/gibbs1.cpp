#include "gibbs1.h"
#include "Rmath.h"

RcppExport SEXP rcpp_gibbs(SEXP phi, SEXP nu0, SEXP mu0, SEXP tau20,
                           SEXP sigma20, SEXP ybar, SEXP s2, SEXP n) {
    //Rcpp::RNGScope scope;
    using namespace Rcpp;
    RNGScope scope;
    Rcpp::NumericMatrix xphi(phi);
    Rcpp::NumericVector xnu0(nu0);
    Rcpp::NumericVector mun(1);
    Rcpp::NumericVector xmu0(mu0);
    Rcpp::NumericVector xtau20(tau20);
    Rcpp::NumericVector xsigma20(sigma20);
    Rcpp::NumericVector xybar(ybar);
    Rcpp::NumericVector xs2(s2);
    Rcpp::NumericVector xn(n);
    Rcpp::NumericVector tau2n(1);
    Rcpp::NumericVector nun(1);
    Rcpp::NumericVector s2n(1);
    Rcpp::NumericVector theta(1);
    Rcpp::NumericVector postprec(1);

    nun = xnu0 + xn;

    int n_mcmc = xphi.nrow();
    xphi(0,0) = xybar[0];
    xphi(0,1) = xs2[0];
    for(int i = 1; i < n_mcmc; i++) {
        mun = (xmu0/xtau20 + xn*xybar*1/xs2)/(1/xtau20 + xn*1/xs2);
        tau2n = 1/(1/xtau20 + xn*1/xs2);
        theta = rnorm(1, mun[0], sqrt(tau2n[0]));
        s2n = (xnu0*xsigma20 + (xn-1)*xs2 + xn*pow(xybar - theta, 2))/nun;
        postprec = rgamma(1, nun[0]/1, nun[0]/2*s2n[0]);
        xphi(i, 0) = theta[0];
        xphi(i, 1) = postprec[0];
    }
    return xphi;
}
