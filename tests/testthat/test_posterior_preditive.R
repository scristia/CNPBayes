context("Posterior predictive distribution")

test_that("posteriorPredictive", {
  model <- MarginalModelExample
  mp <- McmcParams(iter=500, burnin=50)
  mcmcParams(model) <- mp
  model <- posteriorSimulation(model)
  y <- posteriorPredictive(model)

  bmodel <- BatchModelExample
  mp <- McmcParams(iter=500, burnin=150, nStarts=20)
  mcmcParams(bmodel) <- mp
  bmodel <- posteriorSimulation(bmodel)
  batchy <- posteriorPredictive(bmodel)
  if(FALSE){
    ## overlay posterior predictive
    df.observed <- data.frame(y=observed(model))
    df.predict <- data.frame(y=y)
    ggplot(df.observed, aes(y)) +
      geom_density(adjust=1/10, fill="gray") +
      geom_density(data=df.predict, aes(y),
                   adjust=1/10,
                   fill="steelblue",
                   alpha=0.5, inherit.aes=FALSE)

    df.observed <- data.frame(y=observed(bmodel))
    df.predict <- data.frame(y=batchy)
    ggplot(df.observed, aes(y)) +
      geom_density(adjust=1/10, fill="beige") +
      geom_density(data=df.predict, aes(y), adjust=1/10,
                   fill="steelblue",
                   alpha=0.5, inherit.aes=FALSE)
  }
})
