//function for skew normal gibbs sampler. Requires RcppArmadillo



// GIBBS
NumericVector group;
for(int s = 1; s < n_mcmc; s++) {
    // update mu, alpha, omega
    for(int k = 1; k < K; k++) {
        group = S[S == k];
        v[k] = 1/(1 + tau[k] * pow(psi[k], 2) );
        m[k] = v[k] * tau[k] * psi[k] * (r[S == k] - mu[k]);
        // simulate z from truncated normal distribution
        // z = constr_norm()?
        
    }



}


