check_rounding_sim <- function() {
    library(dplyr)
    set.seed(1)
    N <- 7524

    X <- rnorm(N * 81, 0, 1)
    X_round <- as.integer(1000 * X) / 1000
    lrr <- data.frame(x=X_round, label=rep(1:N, each=81)) %>%
        group_by(label) %>%
        summarise(y=median(x)) %>%
        .$y

    mp <- McmcParams(iter=1000, burnin=1000, thin=5, nStarts=10)

    # model for data as seen
    model <- MarginalModel(data=lrr, mcmc.params=mp) 
    m.list <- posteriorSimulation(model, k=1:4)
    m.lik <- marginalLikelihood(m.list)

    # model for data w/ random noise
    model <- MarginalModel(data=lrr + runif(N, 0, 1) / 1000, 
                           mcmc.params=mp) 

    m.list <- posteriorSimulation(model, k=1:4)
    m.lik <- marginalLikelihood(m.list)
}

check_less_rounding_sim <- function() {
    library(dplyr)
    set.seed(1)
    N <- 7524

    X <- rnorm(N * 81, 0, 1)
    X_round <- as.integer(round(1000 * X)) / 1000
    lrr <- data.frame(x=X_round, label=rep(1:N, each=81)) %>%
        group_by(label) %>%
        summarise(y=median(x)) %>%
        .$y

    mp <- McmcParams(iter=1000, burnin=1000, thin=5, nStarts=10)

    # model for data as seen
    model <- MarginalModel(data=lrr, mcmc.params=mp) 
    m.list <- posteriorSimulation(model, k=1:4)
    m.lik <- marginalLikelihood(m.list)

    # model for data w/ random noise
    model <- MarginalModel(data=lrr + runif(N, 0, 1) / 100, 
                           mcmc.params=mp) 

    m.list <- posteriorSimulation(model, k=1:4)
    m.lik <- marginalLikelihood(m.list)
}
