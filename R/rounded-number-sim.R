check_rounding_sim <- function() {
    set.seed(1)
    library(dplyr)
    X <- rnorm(7524 * 81, 0, 1)
    X_round <- as.integer(1000 * X) / 1000
    lrr <- data.frame(x=X_round, label=rep(1:7524, each=81)) %>%
        group_by(label) %>%
        summarise(y=median(x)) %>%
        .$y

    mp <- McmcParams(iter=1000, burnin=1000, thin=5, nStarts=10)
    model <- MarginalModel(data=lrr, mcmc.params=mp) 

    m.list <- posteriorSimulation(model, k=1:4)
}
