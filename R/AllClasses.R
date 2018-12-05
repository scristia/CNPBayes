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

#' An object to specify the hyperparameters of a model with additional parameters for trios
#' 
#' This class inherits from the hyperparameters class. 
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
setClass("HyperparametersTrios", contains="Hyperparameters")

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
#' @slot predictive posterior predictive distribution
#' @slot zstar needed for plotting posterior predictive distribution
#' @slot k integer specifying number of components
#' @slot iter integer specifying number of MCMC simulations
#' @slot B integer specifying number of batches
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
                                      predictive="matrix",
                                      zstar="matrix",
                                      k="integer",
                                      iter="integer",
                                      B="integer"))

setClass("McmcChainsTrios", contains="McmcChains",
         slots=c(pi_parents="matrix",
                 zfreq_parents="matrix",
                 is_mendelian="integer"))
##family_member="character"))


#' An object to specify MCMC options for a later simulation
#'
#' @slot thin A one length numeric to specify thinning. A value of n indicates that every nth sample should be saved. Thinning helps to reduce autocorrelation.
#' @slot iter A one length numeric to specify how many MCMC iterations should be sampled.
#' @slot burnin A one length numeric to specify burnin. The first $n$ samples will be discarded.
#' @slot nstarts A one length numeric to specify the number of chains in a simulation.
#' @slot param_updates Indicates whether each parameter should be updated (1) or fixed (0).
#' @slot min_GR minimum value of multivariate Gelman Rubin statistic for diagnosing convergence. Default is 1.2.
#' @slot min_effsize  the minimum mean effective size of the chains. Default is 1/3 * iter.
#' @slot max_burnin The maximum number of burnin iterations before we give up and return the existing model.
#' @slot min_chains minimum number of independence MCMC chains used for assessing convergence. Default is 3.
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
                                      param_updates="integer",
                                      min_GR="numeric",
                                      min_effsize="numeric",
                                      max_burnin="numeric",
                                      min_chains="numeric"))

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
#' @slot marginal_lik the marginal likelihood of the model
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
                                        pi="matrix",
                                        mu="numericOrMatrix",
                                        tau2="numericOrMatrix",
                                        predictive="numeric",
                                        zstar="numeric",
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
#' @slot is_mendelian integer vector equal in length to the number of trios
#' @slot .internal.constraint Constraint on parameters. For internal use only.
#' @export
setClass("MultiBatchModel", contains="MixtureModel")

setClass("TrioBatchModel", contains="MultiBatchModel",
         slots=c(triodata="tbl_df",
                 mprob="matrix",
                 father="numeric",
                 mother="numeric",
                 maplabel="numeric",
                 pi_parents="numeric",
                 zfreq_parents="integer",
                 probz_par="matrix",
                 is_mendelian="integer"))

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
#' @rdname SingleBatchModel-class
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
#' @rdname SingleBatchModel-class
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


#' An object for running MCMC simulations.
#'
#' BatchModel and MarginalModel both inherit from this class.
#' @slot data a tibble with one-dimensional summaries (oned), id, and batch
#' @slot parameters list of parameters
#' @slot chains object of class McmcChains
#' @slot current_values current value of each chain
#' @slot summaries list of empirical data and model summaries
#' @slot flags list of model flags
#' @export
setClass("MultiBatch", representation(data="tbl_df",
                                      specs="tbl_df",
                                      parameters="list",
                                      chains="McmcChains",
                                      current_values="list",
                                      summaries="list",
                                      flags="list"))

setClass("CnList", contains="MultiBatch")

setClass("MultiBatchP", contains="MultiBatch")

setClass("Trios", contains="MultiBatch")
