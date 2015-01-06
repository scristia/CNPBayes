sntest <- function(r, K, nsim, burnin=500) {

    #mu <- rep(mean(xx), K)
    alpha0 <- rep(0, K) ## skewness parameter
    #alpha0 <- c(-3, 0)
    omega0 <- rep(mad(xx), K) ## scale parameter
    omega20 <- omega0^2

    pars <- kmeans(xx, centers=K, nstart=15)
    mu <- sort(pars$centers)
    S <- rep(NA, length(xx))
    for(i in 1:K) S[pars$cluster == order(pars$center)[i]] <- i
    S <- S-1L
    nn <- pars$size[order(pars$centers)]

    eta0 <- rep(1/K, K) ## intitial mixing params
    mat <- .Call("skewnormal_mix", r, K=K, S=S, centers=mu, alpha=alpha0,
                 omega2=omega20, eta=eta0, nsim)

    snmixture <- list("mu"=mat$MU[-burnin,], "omega"=mat$OMEGA[-burnin,],
                      "alpha"=mat$ALPHA[-burnin,], "P"=mat$ETA[-burnin,])
    loglik <- loglik.snmix(r, snmixture, K)
    bic <- -2*loglik + (4*K-1)*log(length(r))
    mat$BIC <- bic

    return(mat)
}
