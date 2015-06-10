#' @include AllGenerics.R

setClassUnion("numericOrMatrix", c("numeric", "matrix"))

#' An object to specify the hyperparameters of a model.
#'
#' @slot k Number of components
#' @slot mu.0 Prior mean for mu.
#' @slot tau2.0 prior variance on mu
#' @slot eta.0 rate paramater for tau2
#' @slot m2.0 shape parameter for tau2
#' @slot alpha mixture probabilities
#' @slot beta parameter for nu.0 distribution
#' @slot a shape for sigma2.0
#' @slot b rate for sigma2.0
#' @aliases k,Hyperparameters-method
setClass("Hyperparameters", representation(k="integer",
                                           mu.0="numeric",
                                           tau2.0="numeric",
                                           eta.0="numeric",
                                           m2.0="numeric",
                                           alpha="numeric",
                                           beta="numeric",
                                           a="numeric",
                                           b="numeric"))

#' An object to specify the hyperparameters of a marginal model.
#'
#' This class inherits from the Hyperparameters class. This class is for hyperparameters which are marginal over the batches.
#' @slot k Number of components
#' @slot mu.0 Prior mean for mu.
#' @slot tau2.0 prior variance on mu
#' @slot eta.0 rate paramater for tau2
#' @slot m2.0 shape parameter for tau2
#' @slot alpha mixture probabilities
#' @slot beta parameter for nu.0 distribution
#' @slot a shape for sigma2.0
#' @slot b rate for sigma2.0
setClass("HyperparametersMarginal", contains="Hyperparameters")

#' An object to specify the hyperparameters of a batch effect model.
#'
#' This class inherits from the Hyperparameters class. This class is for hyperparameters which are hierachical over the batches.
#' @slot k Number of components
#' @slot mu.0 Prior mean for mu.
#' @slot tau2.0 prior variance on mu
#' @slot eta.0 rate paramater for tau2
#' @slot m2.0 shape parameter for tau2
#' @slot alpha mixture probabilities
#' @slot beta parameter for nu.0 distribution
#' @slot a shape for sigma2.0
#' @slot b rate for sigma2.0
setClass("HyperparametersBatch",  contains="Hyperparameters")

#' An object to hold estimated paraeters.
#'
#' An object of this class holds estimates of each parameter at each iteration of the MCMC simulation.
#' @slot theta means of each batch and component
#' @slot sigma2 variances of each batch and component
#' @slot pi mixture probabilities
#' @slot mu overall mean in a marginal. In batch model, averaged across batches
#' @slot tau2 overall variance in a marginal model. In a batch model, weighted average by precision across batches.
#' @slot nu.0 shape parameter for sigma.2 distribution
#' @slot sigma2.0 rate parameter for sigma.2 distribution
#' @slot logprior log likelihood of prior.
#' @slot loglik log likelihood.
#' @slot zfreq table of z.
#' @slot z latent variables
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
#' @slot param_updates Still looking.
#' @examples 
#' McmcParams()
#' McmcParams(iter=1000)
#' mp <- McmcParams()
#' iter(mp)
#' @aliases iter,burnin,nStarts,McmcParams-method
setClass("McmcParams", representation(thin="numeric",
                                      iter="numeric",
                                      burnin="numeric",
                                      nstarts="numeric",
                                      param_updates="integer"))

#' An object for running MCMC simulations.
#'
#' BatchModel and MarginalModel both inherit from this class.
#' @slot k An integer value specifying the number of latent classes.
#' @slot hyperparams An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @slot theta the means of each component and batch
#' @slot sigma2 the variances of each component and batch
#' @slot nu.0 the shape parameter for sigma2
#' @slot sigma2.0 the rate parameter for sigma2
#' @slot pi mixture probabilities which are assumed to be the same for all batches
#' @slot mu overall mean
#' @slot tau2 overall variance
#' @slot data the data for the simulation.
#' @slot data.mean the empirical means of the components
#' @slot data.prec the empirical precisions
#' @slot z latent variables
#' @slot zfreq table of latent variables
#' @slot probz n x k matrix of probabilities
#' @slot logprior log likelihood of prior: log(p(sigma2.0)p(nu.0)p(mu))
#' @slot loglik log likelihood: $\sum p_k \Phi(\theta_k, \sigma_k)
#' @slot mcmc.chains an object of class 'McmcChains' to store MCMC samples
#' @slot batch a vector of the different batch numbers
#' @slot batchElements a vector labeling from which batch each observation came from
#' @slot hwe Hardy Weinberg Equilibrium
#' @slot modes the values of parameters from the iteration which maximizes log likelihood and log prior
#' @slot theta_order to be removed
#' @slot mcmc.params An object of class 'McmcParams'
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
#' @slot mcmc.params An object of class 'McmcParams'
setClass("MarginalModel", contains="MixtureModel")
setClass("UnivariateBatchModel", contains="BatchModel")
# setClass("UnivariateMarginalModel", contains="MarginalModel")

#' An object to specify the parameters for MCMC simulation.
#'
#' This class has slots for parameters and data needed for simulation.
#' @slot type Indicates whether the model should be marginal over the batches, or hierarchical over the batches.
#' @slot k Number of components.
#' @slot data A numeric vector containing the data
#' @slot batch A vector indicating from which batch each observation came.
#' @slot mcmc.params An object of class McmcParams indicating number of iterations, etc.
#' @examples 
#' ModelParams()
#' ModelParams(k=3)
setClass("ModelParams", representation(type="character",
                                       k="numeric",
                                       data="numeric",
                                       batch="factor",
                                       mcmc.params="McmcParams"))


setClass("ModelList", representation(model_list="list", names="character",
                                     data="numeric"))

setClass("MarginalModelList", contains="ModelList")
setClass("BatchModelList", contains="ModelList")

setClass("PosteriorSummary", representation(p_theta="matrix", chib="numeric", berkhof="numeric",
                                            marginal="numeric", delta_marginal="numeric"))
