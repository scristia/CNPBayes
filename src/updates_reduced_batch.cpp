#include "miscfunctions.h" // for rdirichlet, tableZ, ...
#include "updates_batch.h"
#include "updates_marginal.h" // for log_ddirichlet_
#include <Rcpp.h>
#include <Rmath.h>

Rcpp::NumericMatrix toMatrix(Rcpp::NumericVector x, int NR, int NC) {
    Rcpp::NumericMatrix Y(NR, NC);
    int iter = 0;
    for(int j = 0; j < NC; ++j) {
        for(int i = 0; i < NR; ++i) {
            Y(i, j) = x[iter];
            iter++;
        }
    }
    return Y;
}

// [[Rcpp::export]]
Rcpp::NumericVector marginal_theta_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params=model.slot("mcmc.params");
    Rcpp::S4 chains(model.slot("mcmc.chains")); 
    int S = params.slot("iter");
    Rcpp::List modes = model.slot("modes");
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix thetastar=clone(theta_);
    int K = thetastar.ncol();
    Rcpp::NumericVector p_theta(S);
    Rcpp::NumericMatrix muc = chains.slot("mu");
    Rcpp::NumericMatrix tau2c = chains.slot("tau2");
    Rcpp::NumericMatrix sigma2 = chains.slot("sigma2");
    int B = thetastar.nrow();
    Rcpp::NumericVector tmp(1);

    Rcpp::IntegerMatrix Z = chains.slot("z");
    Rcpp::IntegerVector zz;

    Rcpp::NumericVector tau2_tilde(K);
    Rcpp::NumericVector sigma2_tilde(K);

    // this should be updated for each iteration
    Rcpp::NumericMatrix data_mean(B, K);
    Rcpp::IntegerVector nn(K);
    double post_prec;
    double tau_n;
    double mu_n;
    double w1;
    double w2;
    Rcpp::NumericVector tauc(K);
    Rcpp::NumericMatrix iSigma2(B, K);
    Rcpp::NumericVector invs2;
    Rcpp::NumericVector theta(1);

    for (int s=0; s < S; ++s) {
        tauc = sqrt(tau2c(s, Rcpp::_));
        zz = Z(s, Rcpp::_);
        model.slot("z") = zz;
        nn = tableZ(K, zz);
        data_mean = compute_means_batch(model);
        tau2_tilde = 1.0 / tau2c(s, Rcpp::_);
        invs2 = 1.0 / sigma2(s, Rcpp::_);    // this is a vector of length B*K
        sigma2_tilde = Rcpp::as<Rcpp::NumericVector>(toMatrix(invs2, B, K));
        double prod = 1.0;

        for (int k = 0; k < K; ++k) {
            for (int b = 0; b < B; ++b) {
                post_prec = tau2_tilde[k] + sigma2_tilde(b, k) * nn[k];
                tau_n = sqrt(1/post_prec);
                w1 = tau2_tilde[k]/post_prec;
                w2 = nn[k] * sigma2_tilde(b, k)/post_prec;
                mu_n = w1*muc(s, k) + w2*data_mean(b, k);
                theta = thetastar(b, k);
                tmp = dnorm(theta, mu_n, tau_n);
                prod *= tmp[0];
            }
        }

        p_theta[s] = prod;
    }

    return p_theta;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_theta_zpermuted_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 mcmcp = model.slot("mcmc.params");
    int S = mcmcp.slot("iter");
    Rcpp::List modes = model.slot("modes");
    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix sigma2star=clone(sigma2_);
    Rcpp::NumericMatrix thetastar=clone(theta_);
    int K = thetastar.ncol();
    int B = thetastar.nrow();
    Rcpp::NumericVector logp_theta(S);
    Rcpp::S4 chains(model.slot("mcmc.chains"));
    Rcpp::NumericVector mu(K);
    Rcpp::NumericVector tau(K);
    Rcpp::NumericVector tmp(1);
    Rcpp::IntegerMatrix Z = chains.slot("z");
    Rcpp::NumericVector tau2(1);

    for (int s = 0; s < S; ++s) {
        // update parameters
        model.slot("z") = Z(s, Rcpp::_ );
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        model.slot("theta") = update_theta_batch(model);
        model.slot("sigma2") = update_sigma2_batch(model);
        model.slot("pi") = update_p_batch(model);
        model.slot("mu") = update_mu_batch(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);
        mu = model.slot("mu");
        tau2 = model.slot("tau2");
        tau = sqrt(tau2);
        
        // calculate probability
        double prod = 0.0;

        for (int k = 0; k < K; ++k) {
            for (int b = 0; b < B; ++b) {
                Rcpp::NumericVector theta = thetastar(b, k);
                tmp = dnorm(theta, mu[k], tau[0]);
                prod += log(tmp[0]);
            }
        }

        logp_theta[s] = prod;
    }

    return logp_theta;
}

Rcpp::NumericVector marginal_sigma2_batch(Rcpp::S4 xmod, Rcpp::S4 mcmcp) {
    Rcpp::RNGScope scope;
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_) ;    
    Rcpp::S4 params(mcmcp);
    int S = params.slot("iter");

    // Assume current values are the modes (in R, useModes(object) ensures this)
    Rcpp::List modes = model.slot("modes");
    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix sigma2star = clone(sigma2_);
    Rcpp::NumericMatrix thetastar = clone(theta_);

    int K = thetastar.ncol();
    int B = thetastar.nrow();
    Rcpp::NumericMatrix prec(B, K);

    for (int k = 0; k < K; ++k) {
        prec(Rcpp::_, k) = 1.0 / sigma2star(Rcpp::_, k);
    }  

    Rcpp::NumericVector logp_prec(S);

    //
    // Run reduced Gibbs    -- theta is fixed at modal ordinate
    //
    Rcpp::S4 chains(model.slot("mcmc.chains"));
    Rcpp::NumericVector tmp(K);
    Rcpp::NumericVector nu0 = chains.slot("nu.0");
    Rcpp::NumericVector s20 = chains.slot("sigma2.0");

    for (int s = 0; s < S; ++s) {
        // update parameters
        model.slot("z") = update_z_batch(model);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        // model.slot("theta") = update_theta(model) ;  Do not update theta!
        model.slot("sigma2") = update_sigma2_batch(model);
        model.slot("pi") = update_p_batch(model);
        model.slot("mu") = update_mu_batch(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);

        // calculate probability of sigma2
        nu0 = model.slot("nu.0");
        s20 = model.slot("sigma2.0");

        double total = 0.0 ;

        for(int b=0; b < B; ++b) {
            tmp = Rcpp::dgamma(prec(b, Rcpp::_), 
                               0.5 * nu0[0], 
                               2.0 / (nu0[0] * s20[0]));
            for(int k = 0; k < K; ++k){
                total += log(tmp[k]);
            }
        }

        logp_prec[s] = total;
    }

    return logp_prec;
}
    
// [[Rcpp::export]]
Rcpp::S4 simulate_z_reduced1_batch(Rcpp::S4 object) {
    Rcpp::RNGScope scope;

    // initialize S4 and list objects
    Rcpp::S4 model_(object);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params = model.slot("mcmc.params");
    Rcpp::S4 chains = model.slot("mcmc.chains");
    Rcpp::List modes = model.slot("modes");

    // fix theta at mode
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix thetastar = clone(theta_);

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    Rcpp::NumericVector nu0chain = chains.slot("nu.0");
    Rcpp::NumericVector s20chain = chains.slot("sigma2.0");
    model.slot("theta") = thetastar;
    int S = params.slot("iter");

    //
    // Run reduced Gibbs    -- theta is fixed at modal ordinate
    //  
    for (int s = 0; s < S; ++s) {
        // update parameters
        model.slot("z") = update_z_batch(model);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        //model.slot("theta") = update_theta(model) ; Do not update theta !
        model.slot("sigma2") = update_sigma2_batch(model);
        model.slot("pi") = update_p_batch(model);
        model.slot("mu") = update_mu_batch(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);

        nu0chain[s] = model.slot("nu.0");
        s20chain[s] = model.slot("sigma2.0");
        Z(s, Rcpp::_) = Rcpp::as<Rcpp::NumericVector>(model.slot("z"));
    }

    chains.slot("z") = Z;
    chains.slot("nu.0") = nu0chain;
    chains.slot("sigma2.0") = s20chain;
    model.slot("mcmc.chains") = chains;

    return model;
}

// [[Rcpp::export]]
Rcpp::S4 simulate_z_reduced2_batch(Rcpp::S4 object) {
    Rcpp::RNGScope scope;

    // initialize S4 and list objects
    Rcpp::S4 model_(object);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params=model.slot("mcmc.params");
    Rcpp::S4 chains=model.slot("mcmc.chains");
    Rcpp::List modes = model.slot("modes");

    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix sigma2star = clone(sigma2_);
    Rcpp::NumericMatrix thetastar = clone(theta_);

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    model.slot("theta") = thetastar;
    model.slot("sigma2") = sigma2star;  
    int S = params.slot("iter");

    //
    // Run reduced Gibbs:
    //   -- theta is fixed at modal ordinate
    //   -- sigma2 is fixed at modal ordinate
    for (int s = 0; s < S; ++s) {
        // update parameters
        model.slot("z") = update_z_batch(model);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        // model.slot("theta") = update_theta(model) ; Do not update theta !
        // model.slot("sigma2") = update_sigma2(model) ;
        model.slot("pi") = update_p_batch(model);
        model.slot("mu") = update_mu_batch(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);

        // store chain of Z's
        Z(s, Rcpp::_) = Rcpp::as<Rcpp::NumericVector>(model.slot("z"));
    }

    chains.slot("z") = Z;
    model.slot("mcmc.chains") = chains;

    return model;
}

// [[Rcpp::export]]
Rcpp::S4 permutedz_reduced1_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params = model.slot("mcmc.params");
    Rcpp::S4 chains = model.slot("mcmc.chains");
    Rcpp::List modes = model.slot("modes");

    Rcpp::NumericVector theta_ = Rcpp::as<Rcpp::NumericVector>(modes["theta"]);
    Rcpp::NumericVector thetastar = clone(theta_);

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    Rcpp::NumericVector nu0chain = chains.slot("nu.0");
    Rcpp::NumericVector s20chain = chains.slot("sigma2.0");
    model.slot("theta") = thetastar;
    int S = params.slot("iter");

    //
    // Run reduced Gibbs    -- theta is fixed at modal ordinate
    //  
    for (int s = 0; s < S; ++s) {
        //model.slot("z") = update_z(xmod) ;
        model.slot("z") = Z(s, Rcpp::_);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        //model.slot("theta") = update_theta(model) ; Do not update theta !
        model.slot("sigma2") = update_sigma2_batch(model);
        model.slot("pi") = update_p_batch(model);
        model.slot("mu") = update_mu_batch(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);

        nu0chain[s] = model.slot("nu.0");
        s20chain[s] = model.slot("sigma2.0");
    }

    chains.slot("nu.0") = nu0chain;
    chains.slot("sigma2.0") = s20chain;
    model.slot("mcmc.chains") = chains;

    return model ;
}


// [[Rcpp::export]]
Rcpp::S4 permutedz_reduced2_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params=model.slot("mcmc.params");
    Rcpp::S4 chains=model.slot("mcmc.chains");
    Rcpp::List modes = model.slot("modes");
    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix sigma2star = clone(sigma2_);
    Rcpp::NumericMatrix thetastar = clone(theta_);

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    Rcpp::NumericVector nu0chain = chains.slot("nu.0");
    Rcpp::NumericVector s20chain = chains.slot("sigma2.0");
    model.slot("theta") = thetastar;
    model.slot("sigma2") = sigma2star;
    int S = params.slot("iter");

    //
    // Run reduced Gibbs:
    //   -- theta is fixed at modal ordinate
    //   -- sigma2 is fixed at modal ordinate
    //  
    for (int s = 0; s < S; ++s) {
        //model.slot("z") = update_z(xmod) ;
        model.slot("z") = Z(s, Rcpp::_);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        // model.slot("theta") = update_theta(model) ; Do not update theta !
        // model.slot("sigma2") = update_sigma2(model) ;
        model.slot("pi") = update_p_batch(model);
        model.slot("mu") = update_mu_batch(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);
        nu0chain[s] = model.slot("nu.0");
        s20chain[s] = model.slot("sigma2.0");
    }

    chains.slot("nu.0") = nu0chain;
    chains.slot("sigma2.0") = s20chain;
    model.slot("mcmc.chains") = chains;

    return model;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_pmix_reduced_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    Rcpp::S4 model(xmod);
    Rcpp::S4 mcmcp = model.slot("mcmc.params");
    Rcpp::S4 chains = model.slot("mcmc.chains");
    Rcpp::S4 hypp = model.slot("hyperparams");
    Rcpp::List modes = model.slot("modes");

    Rcpp::NumericVector x = model.slot("data");      
    int K = hypp.slot("k");
    int S = mcmcp.slot("iter");    

    Rcpp::NumericVector p_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
    Rcpp::NumericVector pstar = clone(p_);
    Rcpp::NumericMatrix Z = chains.slot("z");
    Rcpp::IntegerVector alpha = hypp.slot("alpha");
    Rcpp::NumericVector pmix(S);

    //
    // Run reduced Gibbs    -- theta,sigma2 fixed at modal ordinates
    //
    Rcpp::NumericVector alpha_n(K);
    Rcpp::NumericVector tmp(1);

    for (int s = 0; s < S; ++s) {
        for (int k = 0 ; k < K; ++k) {
            alpha_n[k] = sum(Z(s, Rcpp::_) == k+1) + alpha[k];
        }

        tmp = log_ddirichlet_(pstar, alpha_n);
        pmix[s] = exp(tmp[0]);
    }

    return pmix;
}


// [[Rcpp::export]]
Rcpp::S4 reduced_sigma_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params=model.slot("mcmc.params");
    Rcpp::S4 chains=model.slot("mcmc.chains");
    Rcpp::List modes = model.slot("modes");

    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix thetastar = clone(theta_);
    model.slot("theta") = thetastar;

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    Rcpp::NumericVector nu0chain = chains.slot("nu.0");
    Rcpp::NumericVector s20chain = chains.slot("sigma2.0");
    Rcpp::NumericVector muchain = chains.slot("mu");
    Rcpp::NumericVector tauchain = chains.slot("tau2");
    Rcpp::NumericMatrix pichain = chains.slot("pi");
    Rcpp::NumericMatrix sigmachain = chains.slot("sigma2");
    int S = params.slot("iter");
    
    Rcpp::NumericVector sigma2 = model.slot("sigma2");
    Rcpp::NumericVector pi = model.slot("pi");
    Rcpp::NumericVector tau = model.slot("tau2");
    Rcpp::NumericVector mu = model.slot("mu");

    //
    // Run reduced Gibbs    -- theta is fixed at modal ordinate
    //  
    for(int s=0; s < S; ++s){
        model.slot("z") = update_z_batch(model);
        Z(s, Rcpp::_) = Rcpp::as<Rcpp::NumericVector>(model.slot("z"));
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        //model.slot("theta") = update_theta(model) ; Do not update theta !
        model.slot("sigma2") = update_sigma2_batch(model);
        model.slot("pi") = update_p_batch(model);
        model.slot("mu") = update_mu_batch(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);
        nu0chain[s] = model.slot("nu.0");
        s20chain[s] = model.slot("sigma2.0");
        // update the following chains for debugging small sigma2.0 values
        sigma2 = model.slot("sigma2");
        sigmachain(s, Rcpp::_) = sigma2;
        pi = model.slot("pi");
        pichain(s, Rcpp::_) = pi;
        tau = model.slot("tau2");
        tauchain[s] = tau[0];
        mu = model.slot("mu");
        muchain[s] = mu[0];
    }

    chains.slot("z") = Z;
    chains.slot("nu.0") = nu0chain;
    chains.slot("sigma2.0") = s20chain;

    // update the following chains for debugging
    chains.slot("pi") = pichain;
    chains.slot("sigma2") = sigmachain;
    chains.slot("tau2") = tauchain;
    chains.slot("mu") = muchain;

    model.slot("mcmc.chains") = chains;

    return model;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_sigma_reduced_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 chains = model.slot("mcmc.chains");
    Rcpp::S4 params = model.slot("mcmc.params");
    Rcpp::List modes = model.slot("modes");

    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix sigma2star = clone(sigma2_);
    Rcpp::NumericMatrix thetastar = clone(theta_);

    Rcpp::NumericVector x = model.slot("data");
    Rcpp::NumericMatrix tabz = tableBatchZ(xmod) ;
    Rcpp::IntegerVector batch = model.slot("batch") ;
    Rcpp::IntegerVector ub = uniqueBatch(batch) ;

    int n = x.size();
    int S = params.slot("iter");
    int K = thetastar.ncol();
    int B = thetastar.nrow();

    Rcpp::NumericMatrix prec(B, K);
    Rcpp::NumericVector p_prec(S);
    Rcpp::NumericVector tmp(K);
    Rcpp::NumericVector nu0(1);
    Rcpp::NumericVector s20(1);
    Rcpp::IntegerMatrix Z = chains.slot("z");
    Rcpp::IntegerVector zz;

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::NumericVector nu0chain = chains.slot("nu.0");
    Rcpp::NumericVector s20chain = chains.slot("sigma2.0");

    Rcpp::NumericVector nu_n(1);
    Rcpp::NumericVector sigma2_n(1);
    Rcpp::NumericVector shape(1);
    Rcpp::NumericVector rate(1);
    Rcpp::NumericVector sigma2_new(K) ;  

    //
    // Run reduced Gibbs    -- theta is fixed at modal ordinate
    //
    for (int k = 0; k < K; ++k) {
        prec(Rcpp::_, k) = 1.0 / sigma2star(Rcpp::_, k);
    }  

    for (int s = 0; s < S; ++s) {
        s20 = s20chain[s];
        nu0 = nu0chain[s];
        zz = Z(s, Rcpp::_);

        Rcpp::NumericMatrix ss(B, K);

        for (int i = 0; i < n; i++) {
            for (int b = 0; b < B; ++b) {
                if (batch[i] != ub[b]) {
                    continue;
                }

                for (int k = 0; k < K; ++k) {
                    if (zz[i] == k + 1) {
                        ss(b, k) += pow(x[i] - thetastar(b, k), 2);
                    }
                }
            }
        }

        Rcpp::NumericMatrix nn(B, K);
        double total = 1.0;

        for (int b = 0; b < B; ++b) {
            for (int k = 0; k < K; ++k) {
                // update nn
                nn(b, k) = sum((zz == (k + 1)) & (batch == ub[b]));

                // calculate nu_n and sigma2_n
                nu_n = nu0 + nn(b, k);
                sigma2_n = 1.0 / nu_n * (nu0 * s20 + ss(b, k));

                // calculate shape and rate
                shape = 0.5 * nu_n;
                rate = shape * sigma2_n;

                // calculate probability
                Rcpp::NumericVector prec_typed(1);
                prec_typed[0] = prec(b, k);

                tmp = Rcpp::dgamma(prec_typed, shape[0], 1.0 / rate[0]);

                total *= tmp[0];
            }
        }

        p_prec[s] = total;
    }

    return p_prec;
}


// [[Rcpp::export]]
Rcpp::S4 reduced_pi_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    // model objects and accessories
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::List modes = model.slot("modes");
    Rcpp::S4 params=model.slot("mcmc.params");
    Rcpp::S4 chains=model.slot("mcmc.chains");

    // modes
    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericMatrix sigma2star=clone(sigma2_);
    Rcpp::NumericMatrix thetastar=clone(theta_);

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    model.slot("theta") = thetastar;
    model.slot("sigma2") = sigma2star;
    int S = params.slot("iter");

    //
    // Run reduced Gibbs:
    //   -- theta is fixed at modal ordinate
    //   -- sigma2 is fixed at modal ordinate
    //  
    for (int s = 0; s < S; ++s) {
        // update parameters
        model.slot("z") = update_z_batch(model);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        // model.slot("theta") = update_theta(model) ; Do not update theta !
        // model.slot("sigma2") = update_sigma2(model) ;
        model.slot("pi") = update_p_batch(model);
        model.slot("mu") = update_mu_batch(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);

        // capture chain of Zs
        Z(s, Rcpp::_) = Rcpp::as<Rcpp::NumericVector>(model.slot("z"));
    }

    chains.slot("z") = Z;
    model.slot("mcmc.chains") = chains;

    return model;
}

// [[Rcpp::export]]
Rcpp::S4 reduced_mu_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    // model and accessories
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params=model.slot("mcmc.params");
    Rcpp::S4 chains=model.slot("mcmc.chains");
    Rcpp::List modes = model.slot("modes");

    // get modal ordinates
    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericVector pi_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
    Rcpp::NumericMatrix sigma2star=clone(sigma2_);
    Rcpp::NumericMatrix thetastar=clone(theta_);
    Rcpp::NumericVector pistar=clone(pi_);

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    model.slot("theta") = thetastar;
    model.slot("sigma2") = sigma2star;
    model.slot("pi") = pistar;
    int S = params.slot("iter");

    // keep chains for debugging
    Rcpp::NumericVector nu0chain = chains.slot("nu.0");
    Rcpp::NumericVector s20chain = chains.slot("sigma2.0");
    Rcpp::NumericVector muchain = chains.slot("mu");
    Rcpp::NumericVector tauchain = chains.slot("tau2");

    //
    // Run reduced Gibbs:
    //   -- theta is fixed at modal ordinate
    //   -- sigma2 is fixed at modal ordinate
    //  
    for (int s = 0; s < S; ++s) {
        // update parameters
        model.slot("z") = update_z_batch(model);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        // model.slot("theta") = update_theta(model) ; Do not update theta !
        // model.slot("sigma2") = update_sigma2(model) ;
        // model.slot("pi") = update_p(model) ;
        model.slot("mu") = update_mu_batch(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);

        // store chains
        Z(s, Rcpp::_) = Rcpp::as<Rcpp::NumericVector>(model.slot("z"));
    }

    // store chains
    chains.slot("z") = Z;
    model.slot("mcmc.chains") = chains;

    return model;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_mu_reduced_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;
    
    Rcpp::S4 model(xmod);
    Rcpp::S4 mcmcp = model.slot("mcmc.params");
    Rcpp::S4 chains = model.slot("mcmc.chains");
    Rcpp::S4 hypp = model.slot("hyperparams");
    Rcpp::List modes = model.slot("modes");

    Rcpp::NumericVector x = model.slot("data");
    Rcpp::IntegerVector batch = model.slot("batch") ;
    Rcpp::IntegerVector ub = uniqueBatch(batch) ;

    int S = mcmcp.slot("iter");
    int B = ub.size() ;

    // get hyperparameters
    int K = hypp.slot("k");
    double mu_0 = hypp.slot("mu.0");
    double tau2_0 = hypp.slot("tau2.0");
    double tau2_0_tilde = 1.0 / tau2_0;

    // get ordinal modes
    Rcpp::NumericVector p_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericVector mu_ = Rcpp::as<Rcpp::NumericVector>(modes["mu"]);
    Rcpp::NumericVector pstar = clone(p_);
    Rcpp::NumericVector mustar = clone(mu_);
    Rcpp::NumericMatrix thetastar = clone(theta_);

    // tau2
    Rcpp::NumericMatrix tau2chain = chains.slot("tau2");
    Rcpp::NumericMatrix tau2_tilde(S, K);
    Rcpp::NumericMatrix tau2_B_tilde(S, K);

    for (int k = 0; k < K; ++k) {
        tau2_tilde(Rcpp::_, k) = 1.0 / tau2chain(Rcpp::_, k);
        tau2_B_tilde(Rcpp::_, k) = tau2_0_tilde + B * tau2_tilde(Rcpp::_, k);
    }

    // initalize variables for calculation of p_mu
    Rcpp::NumericMatrix n_b = tableBatchZ(model);
    Rcpp::NumericVector p_mu(S);

    for (int s = 0; s < S; ++s) {
        // initialize variables
        Rcpp::NumericVector w1(K);
        Rcpp::NumericVector w2(K);
        Rcpp::NumericVector thetabar(K);
        Rcpp::NumericVector post_prec(K);
        Rcpp::NumericVector mu_k(K);
        Rcpp::NumericVector tau_k(K);
        double total = 1.0;

        for (int k = 0; k < K; ++k) {
            // calculate weights
            w1[k] = tau2_0_tilde / (tau2_0_tilde + B * tau2_tilde(s, k));
            w2[k] = B * tau2_tilde(s, k) / (tau2_0_tilde + B * tau2_tilde(s, k));

            // calculate thetabar
            double n_k = 0.0 ; // number of observations for component k
            double colsumtheta = 0.0;

            for (int b = 0; b < B; ++b) {
                colsumtheta += n_b(b, k) * thetastar(b, k);
                n_k += n_b(b, k);
            }

            thetabar[k] = colsumtheta / n_k;

            // calculate post_prec
            post_prec[k] = sqrt(1.0 / tau2_B_tilde(s, k));

            // calculate mu_k
            mu_k[k] = w1[k] * mu_0 + w2[k] * thetabar[k];

            // calculate p_mu[s]
            Rcpp::NumericVector mu(1);
            mu[0] = mustar[k];
            total *= Rcpp::dnorm(mu, mu_k[k], post_prec[k])[0];
        }

        p_mu[s] = total;
    }

    return p_mu ;
}

// [[Rcpp::export]]
Rcpp::S4 reduced_tau_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    // get model and accessories
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params=model.slot("mcmc.params");
    Rcpp::S4 chains=model.slot("mcmc.chains");
    Rcpp::List modes = model.slot("modes");

    // get modal ordinates
    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericVector pi_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
    Rcpp::NumericVector mu_ = Rcpp::as<Rcpp::NumericVector>(modes["mu"]);
    Rcpp::NumericMatrix sigma2star=clone(sigma2_);
    Rcpp::NumericVector thetastar=clone(theta_);
    Rcpp::NumericVector pistar=clone(pi_);
    Rcpp::NumericVector mustar=clone(mu_);

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    model.slot("theta") = thetastar;
    model.slot("sigma2") = sigma2star;
    model.slot("pi") = pistar;
    model.slot("mu") = mustar;

    int S = params.slot("iter");

    //
    // Run reduced Gibbs:
    //   -- theta is fixed at modal ordinate
    //   -- sigma2 is fixed at modal ordinate
    //  
    for (int s = 0; s < S; ++s) {
        // update parameters
        model.slot("z") = update_z_batch(model);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        // model.slot("theta") = update_theta(model) ; Do not update theta !
        // model.slot("sigma2") = update_sigma2(model);
        // model.slot("pi") = update_p(model);
        // model.slot("mu") = update_mu(model);
        model.slot("tau2") = update_tau2_batch(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);

        // store Z
        Z(s, Rcpp::_) = Rcpp::as<Rcpp::NumericVector>(model.slot("z"));
    }

    chains.slot("z") = Z;
    model.slot("mcmc.chains") = chains;

    return model;
}


// [[Rcpp::export]]
Rcpp::NumericVector p_tau_reduced_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    // get model and accessories
    Rcpp::S4 model(xmod);
    Rcpp::S4 mcmcp = model.slot("mcmc.params");
    Rcpp::S4 chains = model.slot("mcmc.chains");
    Rcpp::S4 hypp = model.slot("hyperparams");
    Rcpp::List modes = model.slot("modes");

    // get modal ordinates
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericVector mu_ = Rcpp::as<Rcpp::NumericVector>(modes["mu"]);
    Rcpp::NumericVector tau2_ = Rcpp::as<Rcpp::NumericVector>(modes["tau2"]);
    Rcpp::NumericVector mustar = clone(mu_);
    Rcpp::NumericVector tau2star = clone(tau2_);
    Rcpp::NumericMatrix thetastar = clone(theta_);

    // hyperparameters
    int K = hypp.slot("k");
    double m2_0 = hypp.slot("m2.0");
    double eta_0 = hypp.slot("eta.0");

    // batch and component parameters
    Rcpp::IntegerVector batch = model.slot("batch");
    Rcpp::IntegerVector ub = uniqueBatch(batch);
    int B = ub.size();

    // updated parameters
    double eta_B = eta_0 + B;
    Rcpp::NumericVector s2_k(K);
    Rcpp::NumericVector m2_k(K);

    double total = 1.0;
    Rcpp::NumericVector p_tau(1);

    for (int k = 0; k < K; ++k) {
        for (int b = 0; b < B; ++b) {
            s2_k[k] += pow(thetastar(b, k) - mustar[k], 2) ;
        }

        m2_k[k] = 1.0 / eta_B * (eta_0 * m2_0 + s2_k[k]);
        Rcpp::NumericVector tau2star_typed(1);
        tau2star_typed[0] = 1.0 / tau2star[k];

        total *= Rcpp::dgamma(tau2star_typed, 0.5 * eta_B,
                              1.0 / (0.5 * eta_B * m2_k[k]))[0];
    }

    p_tau[0] = total;

    return p_tau;
}

// [[Rcpp::export]]
Rcpp::S4 reduced_nu0_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    // get model and accessories
    Rcpp::S4 model_(xmod);
    Rcpp::S4 model = clone(model_);
    Rcpp::S4 params=model.slot("mcmc.params");
    Rcpp::S4 chains=model.slot("mcmc.chains");
    Rcpp::List modes = model.slot("modes");

    // get modal ordinates
    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix theta_ = Rcpp::as<Rcpp::NumericMatrix>(modes["theta"]);
    Rcpp::NumericVector pi_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
    Rcpp::NumericVector mu_ = Rcpp::as<Rcpp::NumericVector>(modes["mu"]);
    Rcpp::NumericVector tau2_ = Rcpp::as<Rcpp::NumericVector>(modes["tau2"]);
    Rcpp::NumericMatrix sigma2star=clone(sigma2_);
    Rcpp::NumericMatrix thetastar=clone(theta_);
    Rcpp::NumericVector pistar=clone(pi_);
    Rcpp::NumericVector mustar=clone(mu_);
    Rcpp::NumericVector tau2star=clone(tau2_);

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    int S = params.slot("iter");
    Rcpp::NumericVector s20chain(S) ;
    model.slot("theta") = thetastar;
    model.slot("sigma2") = sigma2star;
    model.slot("pi") = pistar;
    model.slot("mu") = mustar;
    model.slot("tau2") = tau2star;

    for (int s = 0; s < S; ++s) {
        model.slot("z") = update_z_batch(model);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        // model.slot("theta") = update_theta(model) ; Do not update theta !
        // model.slot("sigma2") = update_sigma2(model);
        // model.slot("pi") = update_p(model);
        // model.slot("mu") = update_mu(model);
        // model.slot("tau2") = update_tau2(model);
        model.slot("nu.0") = update_nu0_batch(model);
        model.slot("sigma2.0") = update_sigma20_batch(model);

        Z(s, Rcpp::_) = Rcpp::as<Rcpp::NumericVector>(model.slot("z"));
        s20chain[s] = model.slot("sigma2.0");
    }

    // update chains
    chains.slot("z") = Z;
    chains.slot("sigma2.0") = s20chain;
    model.slot("mcmc.chains") = chains;

    return model;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_nu0_reduced_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    // get model and accessories
    Rcpp::S4 model(xmod);
    Rcpp::S4 mcmcp = model.slot("mcmc.params");
    Rcpp::S4 chains = model.slot("mcmc.chains");
    Rcpp::S4 hypp = model.slot("hyperparams");
    Rcpp::List modes = model.slot("modes");

    // get modal ordinates
    Rcpp::IntegerVector nu0_ = Rcpp::as<Rcpp::IntegerVector>(modes["nu0"]);
    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericMatrix sigma2star = clone(sigma2_);
    int nu0star = clone(nu0_)[0];

    // hyperparameters
    int K = hypp.slot("k");
    double betas = hypp.slot("beta");
    int B = sigma2star.nrow() ;  

    int S = mcmcp.slot("iter") ;    
    Rcpp::NumericVector p_nu0(S);
    Rcpp::NumericVector s20chain = chains.slot("sigma2.0");

    //
    // compute p(nu0*, ) from *normalized* probabilities
    //
    Rcpp::NumericVector d(100) ;  // 100 is the maximum allowed value for nu_0
    Rcpp::NumericVector lpnu0(100);
    double prec = 0.0;
    double lprec = 0.0;
    d = Rcpp::seq_len(100);

    for (int b = 0; b < B; ++b) {
        for (int k = 0; k < K; ++k) {
            prec += 1.0 / sigma2star(b, k);
            lprec += log(1.0 / sigma2star(b, k));
        }
    }

    Rcpp::NumericVector y1(100);
    Rcpp::NumericVector y2(100);
    Rcpp::NumericVector y3(100);

    for (int s = 0; s < S; ++s) {
        y1 = B * K * (0.5 * d * log(s20chain[s] * 0.5 * d) - lgamma(d * 0.5));
        y2 = (0.5 * d - 1.0) * lprec;
        y3 = d * (betas + 0.5 * s20chain[s] * prec);

        lpnu0 = y1 + y2 - y3;

        Rcpp::NumericVector prob(100);
        prob = exp(lpnu0); // - maxprob);
        prob = prob / sum(prob);  // this is now normalized
        p_nu0[s] = prob[nu0star];
    }

    return p_nu0;
}

// [[Rcpp::export]]
Rcpp::S4 reduced_s20_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    // get model and accessories
    Rcpp::S4 model_(xmod) ;
    Rcpp::S4 model = clone(model_) ;
    Rcpp::S4 params = model.slot("mcmc.params") ;
    Rcpp::S4 chains = model.slot("mcmc.chains") ;
    Rcpp::List modes = model.slot("modes") ;

    // get modal ordinates
    Rcpp::NumericVector sigma2_ = Rcpp::as<Rcpp::NumericVector>(modes["sigma2"]);
    Rcpp::NumericVector theta_ = Rcpp::as<Rcpp::NumericVector>(modes["theta"]);
    Rcpp::NumericVector pi_ = Rcpp::as<Rcpp::NumericVector>(modes["mixprob"]);
    Rcpp::NumericVector mu_ = Rcpp::as<Rcpp::NumericVector>(modes["mu"]);
    Rcpp::NumericVector tau2_ = Rcpp::as<Rcpp::NumericVector>(modes["tau2"]);
    Rcpp::IntegerVector nu0_ = Rcpp::as<Rcpp::IntegerVector>(modes["nu0"]);
    Rcpp::NumericVector sigma2star=clone(sigma2_);
    Rcpp::NumericVector thetastar=clone(theta_);
    Rcpp::NumericVector pistar=clone(pi_);
    Rcpp::NumericVector mustar=clone(mu_);
    Rcpp::NumericVector tau2star=clone(tau2_);
    Rcpp::IntegerVector nu0star=clone(nu0_);

    //
    // We need to keep the Z|y,theta* chain
    //
    Rcpp::IntegerMatrix Z = chains.slot("z");
    model.slot("theta") = thetastar;
    model.slot("sigma2") = sigma2star;
    model.slot("pi") = pistar;
    model.slot("mu") = mustar;
    model.slot("tau2") = tau2star;
    model.slot("nu.0") = nu0star;

    int S = params.slot("iter");

    for (int s = 0; s < S; ++s) {
        // update parameters
        model.slot("z") = update_z_batch(model);
        model.slot("data.mean") = compute_means_batch(model);
        model.slot("data.prec") = compute_prec_batch(model);
        // model.slot("theta") = update_theta(model) ; Do not update theta !
        // model.slot("sigma2") = update_sigma2(model) ;
        // model.slot("pi") = update_p(model) ;
        // model.slot("mu") = update_mu(model) ;
        // model.slot("tau2") = update_tau2(model) ;
        // model.slot("nu.0") = update_nu0(model) ;
        model.slot("sigma2.0") = update_sigma20_batch(model);

        Z(s, Rcpp::_) = Rcpp::as<Rcpp::NumericVector>(model.slot("z"));
    }

    // return chains
    chains.slot("z") = Z;
    model.slot("mcmc.chains") = chains;

    return model;
}

// [[Rcpp::export]]
Rcpp::NumericVector p_s20_reduced_batch(Rcpp::S4 xmod) {
    Rcpp::RNGScope scope;

    // get model and accessories
    Rcpp::S4 model(xmod);
    Rcpp::S4 mcmcp = model.slot("mcmc.params");
    Rcpp::S4 chains = model.slot("mcmc.chains");
    Rcpp::S4 hypp = model.slot("hyperparams");
    Rcpp::List modes = model.slot("modes");

    Rcpp::IntegerVector batch = model.slot("batch");
    Rcpp::IntegerVector ub = uniqueBatch(batch);
    int B = ub.size();

    // get ordinal modes.
    Rcpp::IntegerVector nu0_ = Rcpp::as<Rcpp::IntegerVector>(modes["nu0"]);
    Rcpp::NumericMatrix sigma2_ = Rcpp::as<Rcpp::NumericMatrix>(modes["sigma2"]);
    Rcpp::NumericVector s20_ = Rcpp::as<Rcpp::NumericVector>(modes["sigma2.0"]);
    Rcpp::NumericMatrix sigma2star = clone(sigma2_);
    Rcpp::NumericVector s20star = clone(s20_);
    double nu0star = clone(nu0_)[0];
    
    // get hyperparameters
    int K = hypp.slot("k");
    double a = hypp.slot("a");
    double b = hypp.slot("b");  

    // calculate a_k, b_k
    double prec = 0.0;

    for (int b = 0; b < B; ++b) {
        for (int k = 0; k < K; ++k) {
            prec += 1.0 / sigma2star(b, k);
        }
    }

    double a_k = a + 0.5 * K * B * nu0star;
    double b_k = b + 0.5 * nu0star * prec;

    Rcpp::NumericVector p_s20 = dgamma(s20star, a_k, 1.0 / b_k);

    return p_s20;
}

