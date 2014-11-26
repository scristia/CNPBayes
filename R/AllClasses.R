#' @include AllGenerics.R

setClassUnion("numericOrMatrix", c("numeric", "matrix"))

## Better to have a field for each parameter
setClass("Hyperparameters", representation(k="integer",
                                           mu.0="numeric",
                                           tau2.0="numeric",
                                           eta.0="numeric",
                                           m2.0="numeric",
                                           alpha="numeric",
                                           beta="numeric",
                                           a="numeric",
                                           b="numeric"))

setClass("McmcChains", representation(theta="matrix",
                                      sigma2="matrix",
                                      pi="matrix",
                                      mu="numericOrMatrix",
                                      tau2="numericOrMatrix",
                                      nu.0="numeric",
                                      sigma2.0="numeric",
                                      logpotential="numeric"))

setClass("MixtureModel", representation("VIRTUAL",
                                        hyperparams="Hyperparameters",
                                        theta="numericOrMatrix",
                                        sigma2="numericOrMatrix",
                                        mu="numericOrMatrix",
                                        tau2="numericOrMatrix",
                                        nu.0="numeric",
                                        sigma2.0="numeric",
                                        pi="numeric",
                                        data="numeric",
                                        data.mean="numericOrMatrix",
                                        data.prec="numericOrMatrix",
                                        z="factor",
                                        logpotential="numeric",
                                        mcmc.chains="McmcChains",
                                        batch="vector"))

setClass("BatchModel", contains="MixtureModel")

setClass("MarginalModel", contains="MixtureModel")

setClass("UnivariateBatchModel", contains="BatchModel")

setClass("UnivariateMarginalModel", contains="MarginalModel")

## Better to have a field for each variable
setClass("Posterior", contains="MixtureModel")

setClass("McmcParams", representation(thin="numeric",
                                      iter="numeric",
                                      burnin="numeric"))
