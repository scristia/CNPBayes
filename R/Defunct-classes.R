#' @include Defunct-classes.R
NULL


#' Defunct classes in CNPBayes
#'
#' The following classes in CNPBayes are deprecated and are provided only for compatability.
#'
#'
#' @slot component The component densities.
#' @slot overall The overall (marginal across batches and components) estimate of the density.
#' @slot modes A numeric vector providing the estimated modes in the
#' overall density.  The modes are defined by a crude estimate of the
#' first derivative of the overall density (see \code{findModes}).
#' @slot data A numeric vector containing the data
#' @slot clusters A vector providing the k-means clustering of the
#' component means using the modes as centers.  If an object of class
#' \code{DensityModel} is instantiated with \code{merge=FALSE}, this
#' slot takes values 1, ..., K, where K is the number of components.
#' @aliases DensityBatchModel-class
#' @export
#' @rdname Defunct-classes
setClass("DensityModel", representation(component="list",
                                        overall="numeric",
                                        modes="numeric",
                                        clusters="numeric",
                                        data="numeric",
                                        quantiles="numeric"))

#' @rdname Defunct-classes
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
#' @slot loglik log likelihood: \eqn{\sum p_k \Phi(\theta_k, \sigma_k)}
#' @slot mcmc.chains an object of class 'McmcChains' to store MCMC samples
#' @slot batch a vector of the different batch numbers
#' @slot batchElements a vector labeling from which batch each observation came from
#' @slot modes the values of parameters from the iteration which maximizes log likelihood and log prior
#' @slot mcmc.params An object of class 'McmcParams'
#' @slot label_switch length-one logical vector indicating whether label-switching occurs (possibly an overfit model)
#' @slot .internal.constraint Constraint on parameters. For internal use only.
setClass("MarginalModel", contains="MixtureModel")

#' @rdname Defunct-classes
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
setClass("BatchModel", contains="MixtureModel")

#' @rdname Defunct-classes
#' @export
setClass("DensityBatchModel", representation(batch="list"), contains="DensityModel")


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
#' @rdname Defunct-classes
setClass("BatchModel", contains="MixtureModel")

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
#' @rdname Defunct-classes
setClass("MarginalModel", contains="MixtureModel")


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
