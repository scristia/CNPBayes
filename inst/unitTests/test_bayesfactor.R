.test_bayesfactor <- function(){
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3), theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=1000, burnin=3500)
  p3 <- ModelParams("marginal", y=y(truth), k=3)
  model3 <- posteriorSimulation(p3, mcmcp)
  p4 <- ModelParams("marginal", y=y(truth), k=4)
  model4 <- posteriorSimulation(p4, mcmcp)
  pg3 <- partialGibbs(model3)
  stat3 <- partialGibbsSummary(model3, pg3)
  pg4 <- partialGibbs(model4)
  stat4 <- partialGibbsSummary(model4, pg4)
  bayesFactor(stat3, stat4)

  ##
  ## fit batch model to data without a batch effect
  ## - use Bayes factor as a means to compare models
  p3b <- ModelParams("batch", y=y(truth), k=3,
                     batch=rep(letters[1:3], length.out=2500))
  model3b <- posteriorSimulation(p3b, mcmcp)
  pg3b <- partialGibbs(model3b)
  trace(partialGibbsSummary, browser)
  stat3b <- partialGibbsSummary(model3b, pg3b)




  ## p(theta | y ) = p(y|theta)*p(theta)/m.y
  ## log(m.y) = log(p(y.theta)) + log(p(theta)) - log(p(theta|y))
  ##
  ## If model has not converged, the estimate of the posterior
  ## probability in the denominator can be very small.  Small
  ## p(theta|y) causes the estimate of the marginal density to blow
  ## up.
  ##


  ## "moderately" difficult
  set.seed(100)
  truth <- simulateData(N=2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.05, 0.1, 0.8))
  if(FALSE) plot(truth, use.current=TRUE)
  mcmcp <- McmcParams(iter=1000, burnin=3500, nStarts=20, nStartIter=200)
  p3 <- ModelParams("marginal", y=y(truth), k=3)
  model3 <- posteriorSimulation(p3, mcmcp)
  p4 <- ModelParams("marginal", y=y(truth), k=4)
  model4 <- posteriorSimulation(p4, mcmcp)

  ## neither model has converged
  pg3 <- partialGibbs(model3)
  stat3 <- partialGibbsSummary(model3, pg3)
  pg4 <- partialGibbs(model4)
  stat4 <- partialGibbsSummary(model4, pg4)
  bayesFactor(stat3, stat4)

  ##
  ## p(theta1, theta2 | y) = p(theta2 | y, theta1) * p(theta1 | y)
  ##
  ## p(theta1* | y) = Int p(theta1*, theta2, z | y) d(theta2, z)
  ##                = Int p(theta1*| theta2, z, y) p(theta2, z | y) d(theta2, z)
  ##                ~  mean full conditional for theta1 with complete Gibbs
  ##
  ## p(theta2* | y, theta1*) = Int p(theta2*  | y, theta1*, z) p(z | y, theta1*) dz
  ##                         ~ mean of full conditional for theta2 with reduced Gibbs (theta1 is fixed at theta1*)
  ##
  ## Question: In the Gibbs' sampler we are sampling theta from the full conditional.
  ##           Here, are we computing
  ##            dnorm(theta*, mu(g), tau(g)), for each iteration g
  ##           where mu and tau are sampled from the full conditional at iteration g.  Is this correct?
  ##
  ## Notation
  ##
  ## Let psi = [theta, sigma2, mu, tau2, nu.0, sigma2.0]
  ##
  ## Let psi[-theta] = [sigma2, mu, tau2, nu.0, sigma2.0] (i.e., all
  ## parameters except theta)
  ##
  ## Let psi* denote the value of theta giving rise to the mode of
  ## p(psi|y)p(psi)
  ##
  ## p(theta* | y) is computed as the ergodic average of the full
  ## conditional density with the posterior draws of psi[-theta]
  ## leading to the estimate. Because the z-values are not saved, we
  ## need to run T additional iterations.
  ##




  if(FALSE){
    op <- par(mfrow=c(1,2),las=1)
    plot(truth, use.current=T)
    plot(model)
    par(op)
  }
}
