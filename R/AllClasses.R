#' @include AllGenerics.R
NULL

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
#' @slot dfr positive number for t-distribution degrees of freedom
#' @aliases k,Hyperparameters-method
setClass("Hyperparameters", representation(k="integer",
                                           mu.0="numeric",
                                           tau2.0="numeric",
                                           eta.0="numeric",
                                           m2.0="numeric",
                                           alpha="numeric",
                                           beta="numeric",
                                           a="numeric",
                                           b="numeric",
                                           dfr="numeric"))


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
#' @slot dfr positive number for t-distribution degrees of freedom
setClass("HyperparametersSingleBatch", contains="Hyperparameters")

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
#' @slot dfr positive number for t-distribution degrees of freedom
setClass("HyperparametersMultiBatch",  contains="Hyperparameters")

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
setClass("McmcChains", representation(theta="matrix",
                                      sigma2="matrix",
                                      pi="matrix",
                                      mu="numericOrMatrix",
                                      tau2="numericOrMatrix",
                                      nu.0="numeric",
                                      sigma2.0="numeric",
                                      logprior="numeric",
                                      loglik="numeric",
                                      zfreq="matrix"))

#' An object to specify MCMC options for a later simulation
#'
#' @slot thin A one length numeric to specify thinning. A value of n indicates that every nth sample should be saved. Thinning helps to reduce autocorrelation.
#' @slot iter A one length numeric to specify how many MCMC iterations should be sampled.
#' @slot burnin A one length numeric to specify burnin. The first $n$ samples will be discarded.
#' @slot nstarts A one length numeric to specify the number of chains in a simulation.
#' @slot param_updates Indicates whether each parameter should be updated (1) or fixed (0).
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
#' @slot u chi-square draws for controlling t-distribution
#' @slot logprior log likelihood of prior: log(p(sigma2.0)p(nu.0)p(mu))
#' @slot loglik log likelihood: \eqn{\sum p_k \Phi(\theta_k, \sigma_k)}
#' @slot mcmc.chains an object of class 'McmcChains' to store MCMC samples
#' @slot batch an integer-vector numbering the different batches. Must the same length as \code{data}.
#' @slot batchElements a vector labeling from which batch each observation came from
#' @slot modes the values of parameters from the iteration which maximizes log likelihood and log prior
#' @slot mcmc.params An object of class 'McmcParams'
#' @slot label_switch length-one logical indicating problems with label switching
#' @slot .internal.constraint Constraint on parameters. For internal use only.
#' @slot .internal.counter For internal use only.
#' @slot marginal_lik scalar for marginal likelihood
#' @export
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
                                        u="numeric",
                                        logprior="numeric",
                                        loglik="numeric",
                                        mcmc.chains="McmcChains",
                                        batch="integer",
                                        batchElements="integer",
                                        modes="list",
                                        mcmc.params="McmcParams",
                                        label_switch="logical",
                                        marginal_lik="numeric",
                                        .internal.constraint="numeric",
                                        .internal.counter="integer"))


#' An object for running MCMC simulations.
#'
#' Run hierarchical MCMC for batch model.
#' @slot k An integer value specifying the number of latent classes.
#' @slot hyperparams An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @slot theta the means of each component and batch
#' @slot sigma2 the variances of each component and batch
#' @slot nu.0 the shape parameter for sigma2
#' @slot sigma2.0 the rate parameter for sigma2
#' @slot pi mixture probabilities which are assumed to be the same for all batches
#' @slot mu means from batches, averaged across batches
#' @slot tau2 variances from batches,  weighted by precisions
#' @slot data the data for the simulation.
#' @slot data.mean the empirical means of the components
#' @slot data.prec the empirical precisions
#' @slot z latent variables
#' @slot zfreq table of latent variables
#' @slot probz n x k matrix of probabilities
#' @slot logprior log likelihood of prior: log(p(sigma2.0)p(nu.0)p(mu))
#' @slot loglik log likelihood: \eqn{\sum p_k \Phi(\theta_k, \sigma_k)}
#' @slot mcmc.chains an object of class 'McmcChains' to store MCMC samples
#' @slot batch a vector of the different batch numbers
#' @slot batchElements a vector labeling from which batch each observation came from
#' @slot modes the values of parameters from the iteration which maximizes log likelihood and log prior
#' @slot label_switch length-one logical vector indicating whether label-switching occurs (possibly an overfit model)
#' @slot mcmc.params An object of class 'McmcParams'
#' @slot .internal.constraint Constraint on parameters. For internal use only.
#' @export
setClass("MultiBatchModel", contains="MixtureModel")


#' The 'SingleBatchModel' class
#'
#' Run marginal MCMC simulation
#' @slot k An integer value specifying the number of latent classes.
#' @slot hyperparams An object of class `Hyperparameters` used to specify the hyperparameters of the model.
#' @slot theta the means of each component and batch
#' @slot sigma2 the variances of each component and batch
#' @slot nu.0 the shape parameter for sigma2
#' @slot sigma2.0 the rate parameter for sigma2
#' @slot pi mixture probabilities which are assumed to be the same for all batches
#' @slot mu overall mean
#' @slot tau2 overall variance
#' @slot u latent chi square variable for t-distribution
#' @slot data the data for the simulation.
#' @slot data.mean the empirical means of the components
#' @slot data.prec the empirical precisions
#' @slot z latent variables
#' @slot zfreq table of latent variables
#' @slot probz n x k matrix of probabilities
#' @slot logprior log likelihood of prior: log(p(sigma2.0)p(nu.0)p(mu))
#' @slot loglik log likelihood: \eqn{\sum p_k \Phi(\theta_k, \sigma_k)}
#' @slot mcmc.chains an object of class 'McmcChains' to store MCMC samples
#' @slot batch a vector of the different batch numbers
#' @slot batchElements a vector labeling from which batch each observation came from
#' @slot modes the values of parameters from the iteration which maximizes log likelihood and log prior
#' @slot mcmc.params An object of class 'McmcParams'
#' @slot label_switch length-one logical vector indicating whether label-switching occurs (possibly an overfit model)
#' @slot .internal.constraint Constraint on parameters. For internal use only.
#' @export
setClass("SingleBatchModel", contains="MixtureModel")


setClass("SingleBatchPooled", contains="SingleBatchModel")

## SB ## single batch
## SBP ## single batch pooled
## SBt ## single batch, t-
setClass("SBPt", contains="SingleBatchPooled")


setClass("MultiBatchPooled", contains="MultiBatchModel")

setClass("UnivariateBatchModel", contains="MultiBatchModel")

#' Mixture model container where mixture components have been genotyped
#'
#' The components in a mixture model need not correspond to distinct copy number states. For example, when batch or copy number does not explain skewness or heavy-tails.
#'
#' @slot mapping character string vector indicating the copy number states.  Typically '0', '1', '2', '3', or '4'.
#'
#' @details Suppose a mixture model with four components is selected, where the 3rd and 4th components both correspond to the diploid state.  The mapping slot will be the vector "0", "1", "2", and "2".
#' @rdname CopyNumber-classes
#' @export
setClass("SingleBatchCopyNumber", contains="SingleBatchModel",
         representation(mapping="character"))

#' @rdname CopyNumber-classes
#' @export
setClass("MultiBatchCopyNumber", contains="MultiBatchModel",
         representation(mapping="character"))

#' @rdname CopyNumber-classes
#' @export
setClass("MultiBatchCopyNumberPooled", contains="MultiBatchModel",
         representation(mapping="character"))
