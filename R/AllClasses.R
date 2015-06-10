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

#' An object to store estimated mixture model densities
#'
#' Instances of DensityModel store the estimated densities for each
#' component and the overall (marginal) estimate of the density. The
#' derived class DensityBatchModel additionally stores the density for
#' each batch / component combination (i.e., if there are 3 components
#' and 10 batches, there are 30 estimated densities).  The intended
#' use-case of the DensityModel class is to faciliate visualization of
#' the estimated densities (see examples) as well as to provide an
#' estimate of the number of modes in the overall density. If the
#' number of estimated modes is smaller than the number of components
#' of the best-fitting mixture model, post-hoc merging of components
#' may be useful.
#' @slot component The component densities.
#' @slot overall The overall (marginal across batches and components) estimate of the density.
#' @slot modes A numeric vector providing the estimated modes in the
#' overall density.  The modes are defined by a crude estimate of the
#' first derivative of the overall density (see \code{findModes}).
#' @slot clusters A vector providing the k-means clustering of the
#' component means using the modes as centers.  If an object of class
#' \code{DensityModel} is instantiated with \code{merge=FALSE}, this
#' slot takes values 1, ..., K, where K is the number of components.
#' @examples
#' ## marginal model
#' truth <- simulateData(N=2500, p=rep(1/3, 3),
#'                       theta=c(-1, 0, 1),
#'                       sds=rep(0.1, 3))
#' dm <- DensityModel(truth)
#' print(dm)
#' dm.merged <- DensityModel(truth, merge=TRUE)
#' print(dm.merged)
#' ## here, because there are 3 distinct modes, specifying merge=TRUE
#' ## does not change the resulting clusters
#' identical(clusters(dm), clusters(dm.merged))
#' ## These objects can be plotted
#' plot(dm, oned(truth))
#' ## Note that calling plot on a MixtureModel-derived object returns
#' ## a density object as a side-effect of the plotting
#' dm2 <- CNPBayes::plot(truth)
#' identical(dm, dm2)
#' ## batch model
#' k <- 3
#' nbatch <- 3
#' means <- matrix(c(-1.2, -1.0, -0.8,
#'                  -0.2, 0, 0.2,
#'                   0.8, 1, 1.2), nbatch, k, byrow=FALSE)
#' sds <- matrix(0.1, nbatch, k)
#' N <- 1500
#' truth <- simulateBatchData(N=N,
#'                            batch=rep(letters[1:3], length.out=N),
#'                            theta=means,
#'                            sds=sds,
#'                            p=c(1/5, 1/3, 1-1/3-1/5))
#' dm <- DensityModel(truth)
#' dm.merged <- DensityModel(truth, merge=TRUE)
#' print(dm)
#' dm2 <- CNPBayes::plot(truth)
#' identical(dm, dm2)
#' ## suppress plotting of the batch-specific densities
#' CNPBayes::plot(dm2, oned(truth), show.batch=FALSE)
#' @aliases DensityBatchModel-class
#' @seealso \code{\link{DensityModel}}
#' @export
setClass("DensityModel", representation(component="list",
                                        overall="numeric",
                                        modes="numeric",
                                        clusters="numeric",
                                        quantiles="numeric"))

#' @export
setClass("DensityBatchModel", representation(batch="list"), contains="DensityModel")
