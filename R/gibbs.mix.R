gibbs.mix <- function(r, S=1000, k, delta=0.15, mu0, tau20,
                      nu0, sigma20, kappa0, burnin=100, outliers.rm=FALSE) {
    ## Check: if initial value vectors not of length k, STOP

  if(outliers.rm) {
        quant <- quantile(r, c(0.001, 0.999))
        rd <- r[r > quant[1] & r < quant[2] ]
        rd.low <- r[r <= quant[1]]
        rd.hi <- r[r >= quant[2]]
    }
    else rd <- r

    ## Check if initial values supplied (later)
    inits <- inits(rd, k)
    mu0 <- inits$mu0
#    sigma20 <- inits$sigma20
    nn <- inits$nn
    alpha <- rep(1, k)

    sigma20 <- rep(1, k)
    ## Initialize matrices for theta, prec, pi
    if(sum(nn) != length(rd)) stop("nn must sum to length(r)")

    z <- rep(NA, length(rd))
    s2n <- tau2n <- mun <- numeric(k)
    p <- matrix(NA, sum(nn), k)

    ## for sampling from constrained full conditionals
    a0 <- min(rd)
    b0 <- max(rd)

    ## For first iteration of mcmc, s2 is equal to sigma20
#    s2 <- sigma20
    s2 <- rep(mad(r), k)
    precs <- means <- PI <- matrix(NA, S, k)

    ## simulate from prior
#    rbar <- means[1, ] <- rnorm(k, mu0, tau20)
    rbar <- means[1, ] <- mu0
    Z <- matrix(0, length(rd), k, dimnames=list(names(rd), NULL))

    precs[1, ] <- 1/s2
    post <- .Call("gibbs_mix", r, means, precs, PI, Z, nu0, mu0, kappa0, alpha,
                  tau20, sigma20, rbar, s2, nn, delta, burnin)
#    post <- .Call("gibbs_mix_hier", r, means, precs, PI, Z, nu0, mu0, kappa0, alpha,
#                  tau20, sigma20, rbar, s2, nn, delta, burnin)

    if(outliers.rm) {
        post$Z <- rbind(post$Z, matrix(c(S-burnin, rep(0, k-1)),
                                       nrow=length(rd.low), ncol=k, byrow=TRUE,
                                       dimnames=list(names(rd.low), NULL)))
        post$Z <- rbind(post$Z, matrix(c(rep(0, k-1), S-burnin),
                                       nrow=length(rd.hi), ncol=k, byrow=TRUE,
                                       dimnames=list(names(rd.hi), NULL)))
    }

    MU <- colMeans(post$means)
    SD <- 1/sqrt(colMeans(post$precs))
    PI <- colMeans(post$P)
    CN <- getcn(MU)
    if(k > 1 & all(is.finite(SD) & is.finite(MU))) mixtures <- mixture(MU, SD, PI)
    else mixtures <- 1

    post <- c(post, list("CN"=CN))
    loglik <- loglik.normmix(rd, post, K=k, burnin=burnin)
    bic <- -2*loglik + (3*k-1)*log(length(rd))
    icl <- icl(rd, post, K=k, bic)

    return(c(post, "loglik"=loglik, "bic"=bic, "icl"=icl, "K"=k, "mix"=mixtures))
}

getcn <- function(mu) {
    if(any(is.na(mu))) return(NA)
    if(length(mu) == 1L) return(2)
    d <- diff(c(mu[1], mu[2]))
    if (mu < -0.7 || d > 0.6) {
        cn <- seq(0, length(mu)-1)
        return(cn)
    }
    else if (mu >= -0.7 && mu < -0.2 && d < 0.6) {
        cn <- seq(1, length(mu))
        return(cn)
    }
    else return(seq(2, length(mu)+1))
}
