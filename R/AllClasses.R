#' @include AllGenerics.R

setClassUnion("numericOrMatrix", c("numeric", "matrix"))

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

setClass("MarginalModel", contains="MixtureModel")
setClass("UnivariateBatchModel", contains="BatchModel")
setClass("UnivariateMarginalModel", contains="MarginalModel")

setClass("ModelParams", representation(type="character",
                                       k="numeric",
                                       data="numeric",
                                       batch="factor",
                                       mcmc.params="McmcParams"))


setClass("ModelList", representation(model_list="list", names="character",
                                     data="numeric"))

setClass("MarginalModelList", contains="ModelList")
setClass("BatchModelList", contains="ModelList")

setClass("Posterior", contains="MixtureModel")

setClass("PosteriorSummary", representation(p_theta="matrix", chib="numeric", berkhof="numeric",
                                            marginal="numeric", delta_marginal="numeric"))

setClass("PosteriorSummaryList", representation(data="list", names="character"))
