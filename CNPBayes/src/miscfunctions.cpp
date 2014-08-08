#include "miscfunctions.h"

using namespace Rcpp;
// Function to simulate from dirichlet distribution
void rdirichlet(Rcpp::NumericVector a, Rcpp::NumericVector pr) {
    double sample[a.size()];
    double sample_sum = 0;
    for(int i=0; i<a.size(); i++) {
        sample[i] = as<double>(rgamma(1, a[i], 1));
        sample_sum += sample[i];
    }
    for(int i = 0; i<a.size(); i++) {
        pr[i] = sample[i] / sample_sum ;
    }
}


// Function for drawing from contrained normal distribution for theta
double cons_normal(double mean, double var, double a, double b) {
    double p = R::pnorm(a, mean, sqrt(var), 1, 0) + as<double>(runif(1)) *
        (R::pnorm(b, mean, sqrt(var), 1, 0) - R::pnorm(a, mean, sqrt(var), 1, 0));
    return R::qnorm(p, mean, sqrt(var), 1, 0);
}

