#using Laplace-Metropolis estimator from gibbs output

#P(D|Mk) aprox= (2pi)^(d/2) |psi|^(1/2) Pr(D|theta_hat) Pr(Theta_hat)
# d = dimension of theta (number of parameters)
# theta_hat = posterior mode
# |psi| = minus the inverse hession of log(pr(D|theta)p(theta)) evaluated at theta = theta_hat

# BAYES FACTOR
mcmcbf <- function(theta, prec, p, data=NULL, method="likelihood") {
    ## Function adopted from Raftery 1997 in [______]
    niter <- length(theta[,1])
    d <- length(theta[1, ])*3 - 1
    k <- length(theta[1,])
    if(method == "likelihood") {
        h <- NULL
        for(t in 1:niter) h <- c(h, llik(theta[t,], prec[t,], p[t,] ,data, K=k)
                                 + lprior(theta[t,], prec[t,], p[t,] ))
        hmax <- max(h)
    }
    if(method == "median"){
        thetamed <- apply(theta,2,median)
        precmed <- apply(prec,2,median)
        pmed <- apply(p,2,median)
        hmax <- llik(thetamed, precmed, pmed, data, K=k)
        + lprior(thetamed, precmed, pmed)
    }
    if (d==1) logdetV <- 2*log(mad(theta[,1])) else
        logdetV <- sum(log(eigen(cov.mve(theta)$cov)$values ))
    return( hmax + 0.5 * d * log(2*pi) + 0.5 * logdetV )
}

## assume the parameters are all independent
lprior <- function(theta, prec, p) {
    v <- sum(log(dnorm(theta, mu0, tau20))) + sum(log(dgamma(prec, 1, 1)))
    + log(ddirichlet(p, alpha))
    return(v)
}

llik <- function(theta, prec, p, r=data, K) {
    lik.comp <- function(r, comp) {
        p[comp]*dnorm(r, mean=theta[comp], sd=1/sqrt(prec[comp]))
    }
    ## Create likelihood array with likelihoods of from each component
    liks <- sapply(1:K, lik.comp, r=r)

    d <- rowSums(liks, na.rm=TRUE)
    return(sum(log(d)))
}
