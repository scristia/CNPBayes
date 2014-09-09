#include "grid.h"
#include "Rmath.h"

RcppExport SEXP rcpp_discrete(SEXP grid, SEXP meangrid, SEXP precgrid,
                              SEXP dtheta, SEXP dprec, SEXP y) {
    //Rcpp::RNGScope scope;
    using namespace Rcpp;
    RNGScope scope;
    // convert SEXP inputs to C++ vectors
    NumericMatrix xgrid(grid);
    NumericVector xmeangrid(meangrid);
    NumericVector xprecgrid(precgrid);
    NumericVector xdtheta(dtheta);
    NumericVector xdprec(dprec);
    NumericVector xy(y);

    int G = xmeangrid.size();
    int H = xprecgrid.size();
    int N = xy.size();
    NumericVector tmp(N);
    double lik; // likelihood

    //dprec = dgamma(xprecgrid, xnu0[0]/2, xnu0[0]/2*xsigma20[0]);
    //dtheta = dnorm(xmeangrid, xmu0, sqrt(xtau20[0]));

    for(int g = 0; g < G; g++) {
        for(int h = 0; h < H; h++) {
            tmp = dnorm(xy, xmeangrid[g], 1/(sqrt(xprecgrid[h])));
            lik = 1;
            for(int n = 0; n < N; n++) {
                lik *= tmp[n];
            }
            xgrid(g, h) = xdtheta[g] * xdprec[h] * lik;
        }
    }
    return xgrid;
}
