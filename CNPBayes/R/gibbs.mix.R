gibbs.mix <- function(r, S=1000, k, delta=0.15, mu0, tau20,
                      nu0, sigma20, kappa0, burnin=100) {
    ## Check: if initial value vectors not of length k, STOP

    ## Check if initial values supplied (later)
    inits <- inits(r, k)
    mu0 <- inits$mu0
    sigma20 <- inits$sigma20
    nn <- inits$nn
    alpha <- rep(1, k)

    ## Initialize matrices for theta, prec, pi
    if(sum(nn) != length(r)) stop("nn must sum to length(r)")

    z <- rep(NA, length(r))
    s2n <- tau2n <- mun <- numeric(k)
    p <- matrix(NA, sum(nn), k)

    ## for sampling from contrained full conditionals
    a0 <- min(r)
    b0 <- max(r)

    ## For first iteration of mcmc, s2 is equal to sigma20
    s2 <- sigma20
    precs <- means <- PI <- matrix(NA, S, k)

    ## simulate from prior
    rbar <- means[1, ] <- rnorm(k, mu0, tau20)
    Z <- matrix(0, length(r), k)

    precs[1, ] <- 1/sigma20
    post <- .Call("gibbs_mix", r, means, precs, PI, Z, nu0, mu0, kappa0, alpha,
                  tau20, sigma20, rbar, s2, nn, delta, burnin)

    loglik <- loglik.normmix(r, post, K=k, burnin=burnin)
    bic <- -2*loglik + (3*k-1)*log(length(r))
    icl <- icl(r, post, K=k, bic)

    return(c(post, "loglik"=loglik, "bic"=bic, "icl"=icl, "K"=k))
}
