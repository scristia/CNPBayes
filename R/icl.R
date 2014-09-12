icl  <- function(r, post, K, bic) {
    ## find z_i,k (missing value variables)
    ## is.max: return vector of indicator variables, max(vec) = 1. If more than
    ## one max, choose one at random.
    is.max <- function(vec) {
        indic <- as.integer(vec == max(vec))
        #if(sum(indic) > 1)
        invisible(indic)
    }
    z.hat<- t(apply(post$Z, 1, is.max))

    ## find t_k(x_i|theta_hat). This is conditional probability that x_i in
    ## component k given the parameter estimates.

    p.hat <- colMeans(post$P)
    theta.hat <- colMeans(post$means)
    sigma.hat <- 1/sqrt(colMeans(post$precs))
    ## tk is n dimensional vector of conditional probabilities each sample
    ## x_i belongs to a particular component.
    lik.comp <- function(r, comp) {
        p.hat[comp]*dnorm(r, mean=theta.hat[comp], sd=sigma.hat[comp])
    }

    tk <- sapply(1:K, lik.comp, r=r)
    tk <- tk/rowSums(tk)

    if(K!= 1) {
        if(all(tk > 0 & is.finite(tk))) entropy <- sum(colSums(tk * log(tk)))
        else {
            ## tentative until something better
            v <- tk*log(tk)
            v[which(!is.finite(v))] <- 0
            entropy <- sum(colSums(v))
        }
    }
    else entropy <- 0

    invisible(bic - 2*entropy)
}
