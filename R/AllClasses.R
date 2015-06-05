#' @include AllGenerics.R

setClassUnion("numericOrMatrix", c("numeric", "matrix"))

#' An object to specify the hyperparameters of a model.
#'
#' @slot k Number of components
#' @slot mu.0 a priori mean?
#' @slot tau2.0 a priori precision?
#' @slot eta.0 ?
#' @slot m2.0 ?
#' @slot alpha gamma parameter?
#' @slot beta gamma parameter?
#' @slot a ?
#' @slot b ?
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

#' An object to specify MCMC options for a later simulation
#'
#' @slot thin A one length numeric to specify thinning. A value of n indicates that every nth sample should be saved. Thinning helps to reduce autocorrelation.
#' @slot iter A one length numeric to specify how many MCMC iterations should be sampled.
#' @slot burnin A one length numeric to specify burnin. The first $n$ samples will be discarded.
#' @slot nstarts A one length numeric to specify the number of chains in a simulation.
#' @slot nstart_iter Not sure about this one. RS: Is this used?
#' @slot check_labels Not sure about this one either. RS: Is this used?
#' @slot param_updates Still looking.
#' @examples 
#' McmcParams()
#' McmcParams(k=3)
#' mp <- McmcParams()
#' iter(mp)
#' @aliases iter,burnin,nStarts,McmcParams-method
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


#' The 'MarginalModel' class 
#'
#' @slot k An integer value specifying the number of latent classes.
#' @slot hyperparams An object of class `Hyperparameters` used to specify the hyperparameters of....
#' @slot theta A
#' @slot sigma2 A
#' @slot nu.0 A
#' @slot sigma2.0 A
#' @slot pi A
#' @slot mu A
#' @slot tau2 A
#' @slot data A
#' @slot data.mean A
#' @slot data.prec A
#' @slot z A
#' @slot zfreq A
#' @slot probz A
#' @slot logprior A
#' @slot loglik A
#' @slot mcmc.chains A
#' @slot batch A
#' @slot batchElements A
#' @slot hwe A
#' @slot modes A
#' @slot theta_order A
#' @slot m.y A
#' @slot mcmc.params An object of class 'McmcParams'
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
