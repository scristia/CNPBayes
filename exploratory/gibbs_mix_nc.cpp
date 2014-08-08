/* Description
 *
 *
 * TO DO: - Add conditionals so don't take variance of single element,
 *          or mean of zero elements.
 *        - Add thinning parameter.
 *        - Check: priors and starting values.
 */

#include "gibbs_mix.h"
#include "Rmath.h"

using namespace Rcpp;

// Function to simulate from dirichlet distribution
void rdirichlet(NumericVector a, NumericVector pr) {
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

RcppExport SEXP gibbs_mix(SEXP r, SEXP means, SEXP precs, SEXP P, SEXP Z,
        SEXP nu0, SEXP mu0, SEXP kappa0, SEXP alpha,
        SEXP tau20, SEXP sigma20, SEXP rbar,
        SEXP s2, SEXP nn, SEXP delta, SEXP burnin) {
    // Rcpp::RNGScope scope;
    // initialize objects that are passed from R
    RNGScope scope;
    Rcpp::NumericVector xr(r);
    Rcpp::NumericMatrix xmeans(means);
    Rcpp::NumericMatrix xprecs(precs);
    Rcpp::NumericMatrix xP(P);
    Rcpp::NumericMatrix xZ(Z);
    Rcpp::NumericVector xnu0(nu0);
    Rcpp::NumericVector xmu0(mu0);
    Rcpp::NumericVector xkappa0(kappa0);
    Rcpp::NumericVector xalpha(alpha);
    Rcpp::NumericVector xtau20(tau20);
    Rcpp::NumericVector xsigma20(sigma20);
    Rcpp::NumericVector xrbar(rbar);
    Rcpp::NumericVector xs2(s2);
    Rcpp::IntegerVector xnn(nn);
    Rcpp::NumericVector xdelta(delta);
    Rcpp::IntegerVector xburnin(burnin);

    int n_mcmc = xmeans.nrow(); // length of mcmc chain
    int size = (int) std::accumulate(xnn.begin(), xnn.end(), 0.0);
    const int K = xmeans.ncol(); // number of components
    double a0 = min(xr);
    double b0 = max(xr);

    // initialize vectors that depend on size of K, and not passed from R
    Rcpp::NumericVector mun(K);
    Rcpp::NumericVector tau2n(K);
    Rcpp::NumericVector nun(K);
    Rcpp::NumericVector s2n(K);
    Rcpp::NumericVector pi(K);
    Rcpp::NumericVector nalpha(K); // alpha + nn
    Rcpp::NumericVector theta(K);
    Rcpp::NumericVector postprec(K);
    Rcpp::NumericMatrix d(size, K); // for finding latent variable assignments
    Rcpp::IntegerVector z(size);
    Rcpp::NumericVector u(size);

    for(int j=0; j<K; j++) {
        nun[j] = xnu0[0] + xnn[j];
        theta[j] = as<double>(rnorm(1, xmu0[j], sqrt(xtau20[0])));
        postprec[j] = 1/xsigma20[j];
    }
    // Reference row 1 of mean and prec matrix and store value
    NumericMatrix::Row meanrow = xmeans(0, _);
    meanrow = theta;
    NumericMatrix::Row precrow = xprecs(0, _);
    precrow = postprec;

    // check if homozygous deletion present
    bool homdel = 0;
    if(min(xr) < -1) homdel = 1;

    // GIBBS SAMPLER
    for(int s=1; s<n_mcmc; s++) {
        //should generate pi from dirichlet distribution
        // tau a single value
        for(int j=0; j<K; j++) {
            // update tau2n, nun, mun, s2n
            tau2n[j] = 1/(1/xtau20[0] + xnn[j]*postprec[j]);
            nun[j] = xnu0[0] + xnn[j];
            mun[j] = (1/xtau20[0])/(1/xtau20[0] + xnn[j]*postprec[j])*xmu0[j] +
                xnn[j]*postprec[j]/(1/xtau20[0] + xnn[j]*postprec[j])*xrbar[j];
            s2n[j] = 1/nun[j] * (xnu0[0]*xsigma20[j] + (xnn[j] - 1)*xs2[j] +
                    xkappa0[0]*xnn[j]/(xkappa0[0]*xnn[j]) *
                    pow(xrbar[j]-xmu0[j], 2));
            // theta[j] = as<double>(rnorm(1, mun[j], sqrt(tau2n[j])));
            // simulate from the precision's full conditional
            // if s2n is NA, rgamma returns NaN
            //
            // Metropolis update on varance
            postprec[j] = as<double>(rgamma(1, nun[j]/2, 1/(nun[j]/2*s2n[j])));
            // normal theta updates
            theta[j] = as<double>(rnorm(1, mun[j], sqrt(tau2n[j])));
        }

        // hierarchical update on precision
        // C = as<double>(rgamma(1, K*xnu0 - 1, 1/sum(1/postprec)));
        // b = as<double>(rnorm(1, 1/K * sum(theta, xtau20/K)));
        LogicalVector vec = s2n > -0.75;
        for(j=0; j<K; j++) {
            if(theta[j] < -0.75) {
            }
        // have different update rule
        //

        }

        // Update theta. Use constrained normal distribution
        std::vector<int> indices;
        indices.reserve(K);
        LogicalVector res = is_finite(theta);
        Rcpp::NumericVector endpoints(1, a0);
        for(int j=0; j<K; j++) {
            if(res[j]) {
                indices.push_back(j);
                endpoints.push_back(theta[j]);
            }
        }
        endpoints.push_back(b0);
        LogicalVector q = is_finite(mun);
        double a = 0;
        double b = 0;


        // store theta and precision in their respective matrices
        NumericMatrix::Row meanrow = xmeans(s, _);
        meanrow = theta;
        NumericMatrix::Row precrow = xprecs(s, _);
        precrow = postprec;

        // update alpha
        for(int j=0; j<K; j++) nalpha[j] = xalpha[j] + xnn[j];
        // simulate pi from its multinomial posterior
        rdirichlet(nalpha, pi);
        NumericMatrix::Row Prow = xP(s, _);
        Prow = pi;

        // Simulate latent variables
        // This can be made faster
        for(int j=0; j<K; j++) {
            NumericMatrix::Column dcol = d( _, j);
            dcol = pi[j]*dnorm(xr, theta[j], sqrt(1/postprec[j]));
        }
        // normalize d so rows sum to one
        // If d(i, j) not finite (as in the case null or singular components),
        // so need a check to prevent this.
        for(int i=0; i<size; i++) {
            double dsum = 0;
            // LogicalVector qs = is_finite( d( i, _));
            for(int j=0; j<K; j++) {
                dsum += d( i, j);
            }
            NumericMatrix::Row drow = d( i, _);
            drow = d( i, _) / dsum;
        }

        // generate multinomial samples
        std::fill(xnn.begin(), xnn.end(), 0);
        u = runif(size);
        double tmp;
        for(int i=0; i<size; i++) {
            tmp = d(i,0);
            if(u[i] < tmp) {
                z[i] = 0;
                xnn[0] = xnn[0] + 1; // recalculate nn
            }
            for(int j=1; j<K; j++) {
                if(tmp < u[i] && u[i] < tmp + d(i,j)) {
                    z[i] = j;
                    xnn[j] = xnn[j] + 1;
                }
                tmp += d(i,j);
            }
        }

        // Split data by components and take new sample means and variances.
        // Still needed: add checks for when component is empty or singular.
        std::fill(xrbar.begin(), xrbar.end(), 0);
        for(int i=0; i<size; i++) {
            for(int j=0; j<K; j++) {
                if(z[i] == j) xrbar[j] += xr[i] / xnn[j];
            }
        }
        // recalculate component sample variances
        std::fill(xsigma20.begin(), xsigma20.end(), 0);
        for(int i=0; i<size; i++) {
            for(int j=0; j<K; j++) {
                if(z[i] == j) xsigma20[j] += pow(xr[i]-xrbar[j], 2)/xnn[j];
            }
        }
        // Check for 1 and 0 member components, assign NAs accordingly
        LogicalVector is_one = xnn == 1;
        LogicalVector is_zero = xnn == 0;
        if( is_true(any(is_one)) ) {
            for(int j=0; j<K; j++) if(is_one[j] == 1) xsigma20[j] = NA_REAL;
        }
        if( is_true(any(is_zero)) ) {
            for(int j=0; j<K; j++) {
                if(is_zero[j] == 1) {
                    xsigma20[j] = NA_REAL;
                    xrbar[j] = NA_REAL;
                }
            }
        }

        // for each individual, save which copy number state
        // keep tally, find probability after loop terminates.
        if(s >= xburnin[0]) {
            for(int i=0; i<size; i++) {
                for(int j=0; j<K; j++) {
                    if(z[i] == j) xZ( i, j) += 1;
                }
            }
        }
        //  IntegerMatrix::Row Zrow = xZ(s, _);
        //  Zrow = z;
    } // end Gibbs


    /* Will return a list
     * Return: P, means, precs, Z
     * Do this using:
     * List ret; ret["x"] = xx; ret["y"] = yy; return(ret);
     *
     * Add thinning parameter? Will then have to wrap another for loop
     */

    List ret;
    ret["means"] = xmeans( Range( xburnin[0], n_mcmc-1), _);
    ret["precs"] = xprecs( Range( xburnin[0], n_mcmc-1), _);
    ret["P"] = xP( Range( xburnin[0], n_mcmc-1), _);
    ret["Z"] = xZ;
    ret["start"] = xmu0;
    return ret;
}



/* Function for finding marginal likelihood from mcmc chain
 * Use Bridge sampling as in Fruhwirth-Schnatter book
 */

