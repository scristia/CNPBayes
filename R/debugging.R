debug.function <- function() {
    library(CNPBayes)

    eta.0 <- 1800
    m2.0 <- 1 / 60

    mp <- McmcParams(iter=5000, burnin=1000, thin=1)

    hypp <- Hyperparameters(type="batch", k=3,
                            eta.0=eta.0, m2.0=m2.0)

    k <- 3
    nbatch <- 3
    means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
    sds <- matrix(0.1, nbatch, k)
    N <- 400

    set.seed(42)
    sim.data <- simulateBatchData(N=N,
                                  batch=rep(letters[1:3], length.out=N),
                                  theta=means,
                                  sds=sds,
                                  p=c(1/5, 1/3, 1-1/3-1/5))

    model <- BatchModel(data=y(sim.data), k=3,
                        batch=batch(sim.data),
                        hypp=hypp,
                        mcmc.params=mp)

    set.seed(1)
    model <- posteriorSimulation(model)

    # Rcpp
    set.seed(1)
    .Call("update_tau2_batch", model)

    # MC sim
    set.seed(1)
    shape <- 0.5*eta.0
    rate <- 0.5*eta.0*m2.0
    1 / rgamma(k, shape, rate)

    # MCMC sim
    set.seed(1)
    theta <- theta(model)
    mu <- mu(model)

    s2.k <- numeric(3)
    for (i in 1:k) {
        s2.k[i] <- sum((theta[, i] - mu[i]) ^ 2)
    }

    eta.b <- eta.0 + length(unique(batch(sim.data)))
    m2.k <- 1 / eta.b * (eta.0 * m2.0 + s2.k)
    1 / rgamma(3, 0.5 * eta.b, (eta.b * m2.k) / 2.0)

    # posterior simulation with fixed true tau2
    set.seed(1337)
    model@tau2 <- rep(0.2 ^ 2, 3)
    u <- CNPBayes:::paramUpdates(model)
    u["tau2"] <- 0L
    CNPBayes:::paramUpdates(model) <- u
    model <- posteriorSimulation(model)

    # posterior simulation with fixed true tau2
    set.seed(1337)
    model <- BatchModel(data=y(sim.data), k=3,
                        batch=rep(letters[1:3], length.out=N),
                        hypp=hypp,
                        mcmc.params=mp)
    rownames(means) <- 1:3
    model@theta <- means
    u <- CNPBayes:::paramUpdates(model)
    u["theta"] <- 0L
    CNPBayes:::paramUpdates(model) <- u
    model <- posteriorSimulation(model)
}
