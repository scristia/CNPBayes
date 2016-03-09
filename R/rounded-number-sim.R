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

check_rounding_sim_mean <- function() {
    library(dplyr)
    set.seed(1)
    N <- 7524

    X <- rnorm(N * 81, 0, 1)
    X_round <- as.integer(1000 * X) / 1000
    lrr <- data.frame(x=X_round, label=rep(1:N, each=81)) %>%
        group_by(label) %>%
        summarise(y=mean(x)) %>%
        .$y

    mp <- McmcParams(iter=1000, burnin=1000, thin=2, nStarts=1)

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

check_less_rounding_sim_mean <- function() {
    library(dplyr)
    set.seed(1)
    N <- 7524

    prec <- 1e5

    X <- rnorm(N * 81, 0, 1)
    X_round <- as.integer(round(prec * X)) / prec

    lrr.prec.5 <- data.frame(x=X_round, label=rep(1:N, each=81)) %>%
        group_by(label) %>%
        summarise(y=mean(x)) %>%
        .$y

    prec <- 1e3

    X <- rnorm(N * 81, 0, 1)
    X_round <- as.integer(round(prec * X)) / prec

    lrr.prec.3 <- data.frame(x=X_round, label=rep(1:N, each=81)) %>%
        group_by(label) %>%
        summarise(y=mean(x)) %>%
        .$y

    lrr.actual <- data.frame(x=X, label=rep(1:N, each=81)) %>%
        group_by(label) %>%
        summarise(y=mean(x)) %>%
        .$y

    par(mfrow=c(1, 4))
    hist(lrr.orig, breaks=N/10)
    hist(lrr.prec.3, breaks=N/10)
    hist(lrr.prec.5, breaks=N/10)
    hist(lrr.actual, breaks=N/10)

    mp <- McmcParams(iter=1000, burnin=1000, thin=5, nStarts=1)

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

check_no_rounding_sim_mean <- function() {
    library(CNPBayes)
    library(dplyr)
    set.seed(1)

    N <- 7524
    n <- 1e3
    X <- rnorm(N * n, 0, 1)

    lrr <- data.frame(x=X, label=rep(1:N, each=n)) %>%
        group_by(label) %>%
        summarise(y=median(x)) %>%
        .$y

    mp <- McmcParams(iter=1000, burnin=1000, thin=5, nStarts=1)

    # model for data as seen
    model <- MarginalModel(data=lrr, mcmc.params=mp) 
    m.list <- posteriorSimulation(model, k=1:4)
    m.lik <- marginalLikelihood(m.list)
}
