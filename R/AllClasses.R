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
setClass("HyperparametersMarginal", contains="Hyperparameters")

setClass("HyperparametersBatch",  contains="Hyperparameters")

setClass("McmcChains", representation(theta="matrix",
                                      sigma2="matrix",
                                      pi="matrix",
                                      mu="numericOrMatrix",
                                      tau2="numericOrMatrix",
                                      nu.0="numeric",
                                      sigma2.0="numeric",
                                      logprior="numeric",
                                      loglik="numeric",
                                      zfreq="matrix",
                                      z="matrix"))

setClass("McmcParams", representation(thin="numeric",
                                      iter="numeric",
                                      burnin="numeric",
                                      nstarts="numeric",
                                      nstart_iter="numeric",
                                      check_labels="logical",
                                      param_updates="integer"))

setClass("MixtureModel", representation("VIRTUAL",
                                        k = "integer",
                                        hyperparams="Hyperparameters",
                                        theta="numericOrMatrix",
                                        sigma2="numericOrMatrix",
                                        nu.0="numeric",
                                        sigma2.0="numeric",
                                        pi="numeric",
                                        mu="numericOrMatrix",
                                        tau2="numericOrMatrix",
                                        data="numeric",
                                        data.mean="numericOrMatrix",
                                        data.prec="numericOrMatrix",
                                        z="integer",
                                        zfreq="integer",
                                        probz="matrix",
                                        ##logpotential="numeric",
                                        logprior="numeric",
                                        loglik="numeric",
                                        mcmc.chains="McmcChains",
                                        batch="integer",
                                        batchElements="integer",
                                        hwe="numeric",
                                        modes="list",
                                        theta_order="numeric",
                                        m.y="numeric",
                                        mcmc.params="McmcParams"))



setClass("BatchModel", contains="MixtureModel")


##
setClass("BatchModelPlusHom", contains="BatchModel")
setClass("BatchModelNoHom", contains="BatchModel")

setClass("MarginalModel", contains="MixtureModel")

setClass("UnivariateBatchModel", contains="BatchModel")

setClass("UnivariateMarginalModel", contains="MarginalModel")

## Better to have a field for each variable
setClass("Posterior", contains="MixtureModel")

setClass("ModelParams", representation(type="character",
                                       k="numeric",
                                       data="numeric",
                                       batch="factor",
                                       mcmc.params="McmcParams"))

setClass("PosteriorFiles",
         representation(
             isMarginalModel="logical",
             model="character",
             post1="character"))


setClass("ModelList", representation(model_list="list", names="character",
                                     data="numeric"))
##                                     ##elementType="character",
##                                     names="character")
##                                     ##mode_index="list",
##                                     ##mode_list="list",
##                                     ##maxperm="numeric"))
setClass("MarginalModelList", contains="ModelList")
setClass("BatchModelList", contains="ModelList")

setClass("PosteriorSummary", representation(p_theta="matrix", chib="numeric", berkhof="numeric",
                                            marginal="numeric", delta_marginal="numeric"))

setClass("PosteriorSummaryList", representation(data="list", names="character"))
