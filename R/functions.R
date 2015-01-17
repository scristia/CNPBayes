posteriorPrecisionConjugateNormal <- function(prior.precision, data.precision) {
  prior.precision+data.precision
}

posteriorMeanConjugateNormal <- function(prior.precision, data.precision, posterior.precision, prior.mean, data.mean){
  (prior.precision/posterior.precision)*prior.mean + (data.precision/posterior.precision) * data.mean
}

##sumOfSquares <- function(ylist, means){
##  mapply(function(y, theta) sum((y-theta)^2), y=ylist, theta=means)
##}

stopif <- function(x) stopifnot(!x)


precision <- function(x) 1/var(x, na.rm=TRUE)


gammaShapeRate <- function(mn, sd){
  ##shape/rate   = mn (1)
  ##shape/rate^2 = sd (2)
  ##shape = mn * rate
  ##mn * rate / rate^2 = sd
  ##mn/rate = sd
  ##mn/sd = rate
  rate <- mn/sd
  shape <- mn*rate
  setNames(c(shape, rate), c("shape", "rate"))
}
