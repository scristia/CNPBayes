/* Description
 * REQUIRES RcppArmadillo
 *
 * TO DO: - Vectorize certain functions to reduce overhread from function calls
 *        - Return matrix of classifications
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gibbs_mix.h"
#include "miscfunctions.h"

using namespace Rcpp;
using namespace RcppArmadillo;

// Simulate from mutlivariate normal distribution.
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

/* Need: Hyperparams for alpha, mu, omega, pi (list?)
 *       Starting values for parameters -- decide on.
 *       Length of chain and burnin.
 *       Dimension of parameter space: K.
 *       Data: 1-d summary of log r ratios from probes.
 * */
// [[Rcpp::export]]
RcppExport SEXP skewnormal_mix(SEXP xr, SEXP xK, SEXP xS, SEXP xmu,
        SEXP xalpha, SEXP xomega2, SEXP xeta, SEXP sim) {
    RNGScope scope;
    NumericVector r(xr);
    IntegerVector S(xS); // initial assignments
    NumericVector mu(xmu); // initial mu
    NumericVector alpha(xalpha); // initial alpha (skew)
    NumericVector omega2(xomega2); // initial omega
    NumericVector eta(xeta); // initial eta (mixing probabilities)
    int K = as<int>(xK);
    int nsim = as<int>(sim);


                              /* PRIORS */
    /****************************************************************
     *                        MEAN AND SKEW                         *
     *                                                              *
     * MVN prior on mu and psi (skewness parameter)                 *
     * mean (beta0_mu, beta0_psi) with covariance matrix B0*sigma^2 *
     ****************************************************************/
    // MVN prior on mu and psi (skewness parameter)
    // mean (beta0_mu, beta0_psi) with covariance matrix B0*sigma^2
    double beta0_mu = mean(r);
    double beta0_psi = 0;
    arma::colvec beta0(2);
    beta0[0] = beta0_mu;
    beta0[1] = beta0_psi;

    double B0_mu = 0.1; //controls prior information on mu
    double B0_psi = 0.1; // likewise for psi
    arma::mat B0;
    B0.eye(2,2);
    B0(0,0) = B0_mu;
    B0(1,1) = B0_psi;

    /****************************************************************
     *                        VARIANCE                              *
     *                                                              *
     * Inverse gamma prior on sigma^2: sigma^2_k~InvGamma(c0, C0)   *
     ****************************************************************/
    double c0 = 2.5; // chosen to bound sigma^2 away from 0
    // we choose C0=phi*var(data), where phi influences prior expectation
    // of the amount of heterogeneity explained by differences in group means
    double phi = 0.5; // prior expection of 2/3 explained heterogeneity
    double C0 = phi * var(r);

    /* Can potentially put a hierarchical prior on C0, C0~Gamma(g0, G0)
    * which reduces sensitivity with respect to choosing prior component
    * specific scale paramters (see Stephens 1997) */


    /****************************************************************
     *                    MIXING PROBABILITIES                      *
     *                                                              *
     * eta~Dirichlet(a0, K)                                         *
     ****************************************************************/
    IntegerVector a0 = rep(5, K);
    NumericVector an; // updated parameter


                          /* STARTING VALUES */
    /* mu, psi, omega, eta, and assigments passed in from R.
     * Kmeans used to initialize means and assignments.
     * Omega initialized from mad of data, psi initalized to 0.
     * Transformations must be applied to parameters enable to ensure
     * full conditionals for gibbs sampling. */

    // Initialize nn, vector of component-wise sample size
    IntegerVector nn(K);
    for(int k = 0; k < K; k++) nn[k] = sum(S == k);

    int size = (int) std::accumulate(nn.begin(), nn.end(), 0.0);

    // TRANSFORMATIONS
    NumericVector omega = sqrt(omega2);
    NumericVector delta = alpha/sqrt(1 + pow(alpha, 2));
    NumericVector psi = omega * delta;
    NumericVector sigma2 = omega2*(1 - pow(delta, 2));
    NumericVector tau = 1/sigma2;

                          /* STORAGE */
    NumericMatrix MU(nsim, K);
    NumericMatrix ALPHA(nsim, K);
    NumericMatrix OMEGA(nsim, K);
    NumericMatrix ETA(nsim, K);

    arma::mat X, B, beta, mvdraw, rprd, Cn;
    NumericVector z, rr, u;
    LogicalVector indic;
    // matrix d used for simulating latent variables
    NumericMatrix d(size, K), p(size, K);
    double v, m, cc;
                         /* GIBBS SAMPLER */
    for(int s = 0; s < nsim; s++) {
        for(int k = 0; k < K; k++) {
            LogicalVector indic = (S==k);
            NumericVector rr = r[indic];
            //m = (double)v[k] + (double)tau[k]*(rr-(double)mu[k]);
            NumericVector z(rr.size());

            // may change to take vector size as argument
            for(int j = 0; j < rr.size(); j++) {
                v = 1/(1 + tau[k]*pow(psi[k], 2));
                m = v * tau[k] * psi[k] * (rr[j] - mu[k]);
                //z[j] = trunc_norm(m, sqrt(v[k]));
                z[j] = rnorm_trunc(m, sqrt(v), 0.0, INFINITY);
            }
            X = arma::join_horiz(arma::ones(z.size()), arma::colvec(z));

            B = (B0 + tau[k] * X.t() * X).i();
            beta = B * (B0 * beta0 + tau[k] * X.t() * arma::colvec(rr));

            mvdraw = mvrnormArma(1, beta, B);
            mu[k] = mvdraw(0,0);
            psi[k] = mvdraw(0,1);

            // Simulate tau from its posterior gamma distribution
            cc = c0 + nn[k]/2;

            rprd = arma::colvec(rr) - X * mvdraw.t();
            Cn = C0 + rprd.t() * rprd / 2;
            tau[k] = as<double>(rgamma(1, cc, 1/Cn(0,0)));

        }

        // TRANSFORMATIONS
        alpha = psi * sqrt(tau);
        omega = sqrt(1/tau + pow(psi, 2));


        // GENERATE MULTINOMIAL SAMPLES
        // Sample latent class variables
        // d <- matrix(NA, nrow = length(xx), ncol = K)
        for(int k = 0; k < K; k++) {
            NumericMatrix::Column dcol = d( _, k);
            dcol = eta[k] * dsn(r, mu[k], omega[k], alpha[k]);
        }
        // normalize d so rows sum to one
        // If d(i, j) not finite (as in the case null or singular components),
        // so need a check to prevent this.
        for(int i=0; i<size; i++) {
            double dsum = 0;
            for(int j=0; j<K; j++) {
                dsum += d( i, j);
            }
            NumericMatrix::Row drow = d( i, _);
            drow = d( i, _) / dsum;
        }

        std::fill(nn.begin(), nn.end(), 0); //reset nn
        // multinomial sampling procedure
        u = runif(size);
        double tmp;
        for(int i=0; i<size; i++) {
            tmp = d(i,0);
            if(u[i] < tmp) {
                S[i] = 0;
                nn[0] = nn[0] + 1; // recalculate nn
            }
            for(int j=1; j<K; j++) {
                if(tmp < u[i] && u[i] < tmp + d(i,j)) {
                    S[i] = j;
                    nn[j] = nn[j] + 1;
                }
                tmp += d(i,j);
            }
        }

        // UPDATE MIXING PROBABILITIES (eta)
        an = a0 + nn;
        // rdirichlet edits eta via reference
        rdirichlet(an, eta);


        // !!!STORE!!!
        NumericMatrix::Row mu_row = MU(s, _);
        mu_row = mu;
        NumericMatrix::Row alpha_row = ALPHA(s, _);
        alpha_row = alpha;
        NumericMatrix::Row omega_row = OMEGA(s, _);
        omega_row = omega;
        NumericMatrix::Row eta_row = ETA(s, _);
        eta_row = eta;

        Rprintf("%d\rs = ", s+1);
    }
    Rprintf("\n");

    List ret;
    ret["MU"] = MU;
    ret["ALPHA"] = ALPHA;
    ret["OMEGA"] = OMEGA;
    ret["ETA"] = ETA;
    return ret;
}
