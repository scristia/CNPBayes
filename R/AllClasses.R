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

setClass("HyperparametersMarginal", representation(mu="numeric", tau2="numeric"),
         contains="Hyperparameters")

setClass("HyperparametersBatch", representation(mu.0="numeric", tau2.0="numeric"),
         contains="Hyperparameters")

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
                                        ##mu="numericOrMatrix",
                                        ##tau2="numericOrMatrix",
                                        nu.0="numeric",
                                        sigma2.0="numeric",
                                        pi="numeric",
                                        data="numeric",
                                        data.mean="numericOrMatrix",
                                        data.prec="numericOrMatrix",
                                        z="factor",
                                        ## the posterior probability for each copy number state (sample specific)
                                        probz="matrix",
                                        logpotential="numeric",
                                        mcmc.chains="McmcChains",
                                        batch="vector",
                                        hwe="numeric"))


setClass("BatchModel", contains="MixtureModel", representation(mu="numericOrMatrix", tau2="numericOrMatrix"))

setClass("MarginalModel", contains="MixtureModel")

setClass("UnivariateBatchModel", contains="BatchModel")

setClass("UnivariateMarginalModel", contains="MarginalModel")

## Better to have a field for each variable
setClass("Posterior", contains="MixtureModel")

setClass("McmcParams", representation(thin="numeric",
                                      iter="numeric",
                                      burnin="numeric",
                                      ## whether to constrain theta after the burnin
                                      constrainTheta="logical"))

setClass("ModelParams", representation(type="character",
                                       k="numeric",
                                       data="numeric",
                                       batch="character",
                                       mcmc.params="McmcParams"))
