#' @include help.R
NULL

#' Number of components.
#'
#' This function retrieves the number of a priori components.
#' @param object see \code{showMethods(k)}
#' @export
#' @docType methods
#' @rdname k-method
setGeneric("k", function(object) standardGeneric("k"))

setGeneric("k<-", function(object, value) standardGeneric("k<-"))

setGeneric("z<-", function(object, value) standardGeneric("z<-"))

setGeneric("theta<-", function(object, value) standardGeneric("theta<-"))
setGeneric("sigma2<-", function(object, value) standardGeneric("sigma2<-"))
setGeneric("p<-", function(object, value) standardGeneric("p<-"))

#' Retrieve overall mean
#'
#' @param object see \code{showMethods(mu)}
#' @export
#' @docType methods
#' @rdname mu-method
setGeneric("mu", function(object) standardGeneric("mu"))
setGeneric("mu<-", function(object, value) standardGeneric("mu<-"))
setGeneric("tau2", function(object) standardGeneric("tau2"))
setGeneric("tau2<-", function(object, value) standardGeneric("tau2<-"))
setGeneric("nu.0<-", function(object, value) standardGeneric("nu.0<-"))
setGeneric("sigma2.0<-", function(object, value) standardGeneric("sigma2.0<-"))
setGeneric("logpotential<-", function(object, value) standardGeneric("logpotential<-"))
setGeneric("dataMean<-", function(object, value) standardGeneric("dataMean<-"))
setGeneric("dataPrec<-", function(object, value) standardGeneric("dataPrec<-"))

#' Replace model's simulated MCMC chain.
#' @param object see \code{showMethods(mcmcChains)}
#' @param value object with which to replace the mcmcChains
#' @export
#' @docType methods
#' @rdname mcmcChains-method
setGeneric("mcmcChains<-", function(object, value) standardGeneric("mcmcChains<-"))

#' Retrieve simulated chains from model object.
#' 
#' @param object \code{showMethods(mcmcChains)}
#' @export
#' @docType methods
#' @rdname mcmcChains-method
setGeneric("mcmcChains", function(object) standardGeneric("mcmcChains"))

#' Retrieve simulated chains from model object.
#' 
#' @param object \code{showMethods(chains)}
#' @export
#' @docType methods
#' @rdname chains-method
setGeneric("chains", function(object) standardGeneric("chains"))

#' Retried specified hyper parameters from model.
#'
#' @param object see \code{showMethods(hyperParams)}
#' @export
#' @docType methods
#' @rdname hyperParams-method
setGeneric("hyperParams", function(object) standardGeneric("hyperParams"))

#' Replace hyper paramaters of model.
#'
#' @param object see \code{showMethods(hyperParams)}
#' @param value an object of class 'Hyperparameters'
#' @export
#' @docType methods
#' @rdname hyperParams-method
setGeneric("hyperParams<-", function(object,value) standardGeneric("hyperParams<-"))

#' Initialize empty chain for model.
#'
#' @param object see \code{showMethods(McmcChains)}
#' @export
#' @docType methods
#' @rdname McmcChains-method
setGeneric("McmcChains", function(object) standardGeneric("McmcChains"))

setGeneric("hist")

## RS: might need some help documenting this one
#' @export
setGeneric("plot")

#' Retrieve batches from object.
#'
#' The batches are represented as a vector of integers.
#' @param object see \code{showMethods(batch)}
#' @export
#' @docType methods
#' @rdname batch-method
setGeneric("batch", function(object) standardGeneric("batch"))
setGeneric("batch<-", function(object,value) standardGeneric("batch<-"))

setGeneric("startingValues", function(object, params, zz) standardGeneric("startingValues"))

setGeneric("computeMeans", function(object) standardGeneric("computeMeans"))
setGeneric("computeVars", function(object) standardGeneric("computeVars"))

setGeneric("updateZ", function(object) standardGeneric("updateZ"))

setGeneric("computePotential", function(object) standardGeneric("computePotential"))

setGeneric("dat", function(object) standardGeneric("dat"))
setGeneric("dat<-", function(object,value) standardGeneric("dat<-"))

setGeneric("updateMu", function(object) standardGeneric("updateMu"))

setGeneric("initializeSigma2.0", function(object) standardGeneric("initializeSigma2.0"))

setGeneric("initializeSigma2", function(object) standardGeneric("initializeSigma2"))

setGeneric("initializeMu", function(object) standardGeneric("initializeMu"))

setGeneric("initializeTheta", function(object) standardGeneric("initializeTheta"))

setGeneric("initializeTau2", function(object) standardGeneric("initializeTau2"))

setGeneric("posteriorMultinomial", function(object) standardGeneric("posteriorMultinomial"))

setGeneric("simulateY", function(object, N) standardGeneric("simulateY"))
setGeneric("batchCorrect", function(object) standardGeneric("batchCorrect"))

setGeneric("moveChain", function(object, s) standardGeneric("moveChain"))

setGeneric("updateThetaCpp", function(object, constrain) standardGeneric("updateThetaCpp"))
setGeneric("updateSigma2Cpp", function(object) standardGeneric("updateSigma2Cpp"))
setGeneric("updateTheta", function(object, constrain) standardGeneric("updateTheta"))
setGeneric("updateSigma2", function(object) standardGeneric("updateSigma2"))
setGeneric("updateSigma2.0", function(object) standardGeneric("updateSigma2.0"))
setGeneric("updateNu.0", function(object) standardGeneric("updateNu.0"))
setGeneric("updateTau2", function(object) standardGeneric("updateTau2"))

setGeneric("alpha<-", function(object, value) standardGeneric("alpha<-"))

setGeneric("updateWithPosteriorMeans", function(object) standardGeneric("updateWithPosteriorMeans"))

#' Calculate BIC of a model
#' 
#' @param object see \code{showMethods(bic)}
#' @docType methods
#' @rdname bic-method
#' @export
setGeneric("bic", function(object) standardGeneric("bic"))

#' Retrieve theta.
#'
#' This function retrieves theta from an object. Theta is a vector or matrix representing the means of each component/batch.
#' @param object see \code{showMethods(theta)}
#' @export
#' @docType methods
#' @rdname theta-method
setGeneric("theta", function(object) standardGeneric("theta"))

#' Retrieve the variances of each component and batch distribution
#'
#' For a MarginalModel, this function returns a vector of variances. For a BatchModel, returns a matrix of size number of batches by number of components.
#' @param object see \code{showMethods(sigma2)}
#' @export
#' @docType methods
#' @rdname sigma2-method
setGeneric("sigma2", function(object) standardGeneric("sigma2"))

setGeneric("reorderComponents", function(object, new_levels) standardGeneric("reorderComponents"))

setGeneric("hwe", function(object) standardGeneric("hwe"))

#' @export
setGeneric("probz", function(object) standardGeneric("probz"))

setGeneric("probz<-", function(object, value) standardGeneric("probz<-"))

setGeneric("fitMixtureModels", function(object, mcmcp, K=1:5, batch) standardGeneric("fitMixtureModels"))

#' Retrieve the shape parameter for the sigma.2 distribution.
#'
#' @param object see \code{showMethods(nu.0)}
#' @export
#' @docType methods
#' @rdname nu.0-method
setGeneric("nu.0", function(object) standardGeneric("nu.0"))

#' Retrieve the rate parameter for the sigma.2 distribution.
#'
#' @param object see \code{showMethods(sigma2.0)}
#' @export
#' @docType methods
#' @rdname sigma2.0-method
setGeneric("sigma2.0", function(object) standardGeneric("sigma2.0"))

#' Retrieve data.
#'
#' @param object see \code{showMethods(y)}
#' @export
#' @docType methods
#' @rdname y-method
setGeneric("y", function(object) standardGeneric("y"))


## TODO: oned is for 1-dimensional summary. Perhaps more informative
## than 'y'.

#' Retrieve data.
#'
#' @param object see \code{showMethods(oned)}
#' @export
#' @docType methods
#' @rdname oned-method
#' @export
setGeneric("oned", function(object) standardGeneric("oned"))

#' @export
setGeneric("z", function(object) standardGeneric("z"))

setGeneric("modalParameters", function(object) standardGeneric("modalParameters"))

#' Calculates the modes for each of the parameters estimated by a model.
#'
#' @param object see \code{showMethods(computeModes)}
#' @export
#' @docType methods
#' @rdname computeModes-method
setGeneric("computeModes", function(object) standardGeneric("computeModes"))

setGeneric("switchLabels", function(object) standardGeneric("switchLabels"))

setGeneric("computeDistance", function(object) standardGeneric("computeDistance"))

#' @export
setGeneric("modes", function(object) standardGeneric("modes"))

#' @export
setGeneric("modes<-", function(object,value) standardGeneric("modes<-"))


setGeneric("mu.0", function(object) standardGeneric("mu.0"))
setGeneric("tau2.0", function(object) standardGeneric("tau2.0"))
setGeneric("eta.0", function(object) standardGeneric("eta.0"))
setGeneric("eta.0<-", function(object,value) standardGeneric("eta.0<-"))
setGeneric("m2.0", function(object) standardGeneric("m2.0"))
setGeneric("m2.0<-", function(object,value) standardGeneric("m2.0<-"))

setGeneric("showMeans", function(object) standardGeneric("showMeans"))
setGeneric("showSigmas", function(object) standardGeneric("showSigmas"))

#' @export
setGeneric("collapseBatch", function(object, plate, THR=0.1) standardGeneric("collapseBatch"))

#' @export
setGeneric("thetac", function(object) standardGeneric("thetac"))

setGeneric("thetaMean", function(object) standardGeneric("thetaMean"))

#' Calculate mean variance of each component/batch.
#'
#' Mean is taken across simulated MCMC chain.
#' @param object see \code{showMethods(sigmaMean)}
#' @export
#' @docType methods
#' @rdname sigmaMean-method
setGeneric("sigmaMean", function(object) standardGeneric("sigmaMean"))

#' Calculate mean mixture proportion of components.
#'
#' @param object see showMethods(pMean)
#' @export
#' @docType methods
#' @rdname pMean-method
setGeneric("pMean", function(object) standardGeneric("pMean"))

#' @export
setGeneric("tracePlot", function(object, name, ...) standardGeneric("tracePlot"))

setGeneric("tablez", function(object) standardGeneric("tablez"))

#' Number of MCMC chains.
#'
#' This function retrieves the number of chains used for an MCMC simulation.
#' @param object see showMethods(nStarts)
#' @export
#' @docType methods
#' @rdname nStarts-method
setGeneric("nStarts", function(object) standardGeneric("nStarts"))

#' Reset number of MCMC chains in simulation.
#'
#' This function changes the number of chains used for an MCMC simulation.
#' @param object see showMethods(nStarts)
#' @param value new number of chains
#' @export
#' @docType methods
#' @rdname nStarts-method
setGeneric("nStarts<-", function(object, value) standardGeneric("nStarts<-"))

setGeneric("updateAlpha", function(object) standardGeneric("updateAlpha"))
setGeneric("alpha", function(object) standardGeneric("alpha"))

setGeneric("orderTheta<-", function(object, value) standardGeneric("orderTheta<-"))
setGeneric("orderTheta", function(object) standardGeneric("orderTheta"))

#' Retrieve log likelihood.
#'
#' @param object see showMethods(logLik)
#' @export
#' @docType methods
#' @rdname logLik-method
setGeneric("logLik", function(object) standardGeneric("logLik"))

setGeneric("logLik<-", function(object,value) standardGeneric("logLik<-"))

#' Compute log likelihood for a BatchModel or MarginalModel
#'
#' @param object see \code{showMethods(computeLoglik)}
#' @export
#' @docType methods
#' @rdname computeLoglik-method
setGeneric("computeLoglik", function(object) standardGeneric("computeLoglik"))

#' Number of burnin iterations.
#'
#' This function retrieves the number of burnin simulations to be discarded.
#' @param object see \code{showMethods(burnin)}
#' @export
#' @docType methods
#' @rdname burnin-method
setGeneric("burnin", function(object) standardGeneric("burnin"))

#' Reset number of burnin iterations.
#'
#' This function changes the number of burnin simulations to be discarded.
#' @param object see \code{showMethods("burnin<-")}
#' @param value new number of burnin iterations
#' @export
#' @docType methods
#' @rdname burnin-method
setGeneric("burnin<-", function(object, value) standardGeneric("burnin<-"))

#' Reset number of iterations.
#'
#' This function changes the number of simulations.
#' @param object see \code{showMethods("iter<-")}
#' @param force Allow changing of the size of the elements?
#' @param value new number of iterations
#' @export
#' @docType methods
#' @rdname iter-method
setGeneric("iter<-", function(object, force=FALSE, value) standardGeneric("iter<-"))

#' Number of MCMC iterations.
#'
#' This function retrieves the number of iterations of an MCMC simulation.
#' @param object see \code{showMethods(iter)}
#' @export
#' @docType methods
#' @rdname iter-method
setGeneric("iter", function(object) standardGeneric("iter"))

#' Number of thinning intervals.
#'
#' This function retrieves the number of thinning intervals used for an MCMC simulation.
#' @param object see showMethods(thin)
#' @export
#' @docType methods
#' @rdname thin-method
setGeneric("thin", function(object) standardGeneric("thin"))


#' Run the MCMC simulation.
#'
#' nStarts chains are run. b burnin iterations are run and then discarded. Next, s iterations are run in each train. The mode of the MCMC simulation is also calculated.
#' @param object see showMethods(posteriorSimulation)
#' @export
#' @docType methods
#' @rdname posteriorSimulation-method
setGeneric("posteriorSimulation", function(object) standardGeneric("posteriorSimulation"))

setGeneric("initializeModel", function(params, hypp) standardGeneric("initializeModel"))

setGeneric("updateLabels", function(object) standardGeneric("updateLabels"))

setGeneric("updateMixProbs", function(object) standardGeneric("updateMixProbs"))

#' @export
setGeneric("priorProbs", function(object) standardGeneric("priorProbs"))

setGeneric("posteriorTheta", function(object, mcmcp) standardGeneric("posteriorTheta"))
setGeneric("posteriorSigma2", function(object, mcmcp) standardGeneric("posteriorSigma2"))

#' @export
setGeneric("model", function(object) standardGeneric("model"))
setGeneric("post1", function(object) standardGeneric("post1"))
setGeneric("post2", function(object) standardGeneric("post2"))
setGeneric("post3", function(object) standardGeneric("post3"))
setGeneric("isMarginalModel", function(object) standardGeneric("isMarginalModel"))

setGeneric("computePrior", function(object) standardGeneric("computePrior"))
setGeneric("computeLogLikxPrior", function(object) standardGeneric("computeLogLikxPrior"))

setGeneric("postFiles", function(object) standardGeneric("postFiles"))

#' @export
setGeneric("relabel", function(object, zindex) standardGeneric("relabel"))

setGeneric("HyperParameterList", function(hypp, K) standardGeneric("HyperParameterList"))

#' Creates a list of ModelParams for differing number of components
#'
#' Simplifies creating multiple objects for model parameters for multiple component sizes. Used when number of components is not known. Can be used with Marginal and Batch models.
#' @param hypp An object of class 'Hyperparameters'
#' @param K An integer vector containing the different component sizes for which you wish to create ModelParams.
#' @param data A vector containing the data.
#' @param mcmcp An object of class 'ModelParams'. Can have varying numbers of iterations, etc for each component size.
#' @param batch A vector indicating the batches from which each observation came from.
#' @export
#' @docType methods
#' @rdname ModelParamList-method
setGeneric("ModelParamList", function(hypp, K, data, mcmcp, batch) standardGeneric("ModelParamList"))

setGeneric("runBurnin", function(object, mcmcp) standardGeneric("runBurnin"))
setGeneric("runMcmc", function(object, mcmcp) standardGeneric("runMcmc"))

#' Change MCMC parameter status.
#'
#' This function changes the status of MCMC parameters for a model.
#' @param object see showMethods("paramUpdates<-")
#' @param value new MCMC parameter status.
setGeneric("paramUpdates<-", function(x, value) standardGeneric("paramUpdates<-"))

#' @export
setGeneric("zFreq", function(object) standardGeneric("zFreq"))
setGeneric("zFreq<-", function(object,value) standardGeneric("zFreq<-"))

#' @export
setGeneric("mcmcParams", function(object) standardGeneric("mcmcParams"))
#' @export
setGeneric("mcmcParams<-", function(object, force=FALSE, value) standardGeneric("mcmcParams<-"))

#' Calculate log likelihood of prior for model
#' 
#' @param object see \code{showMethods(logPrior)}
#' @export
#' @docType methods
#' @rdname logPrior-method
setGeneric("logPrior", function(object) standardGeneric("logPrior"))

setGeneric("logPrior<-", function(object,value) standardGeneric("logPrior<-"))

setGeneric("fullGibbs", function(object, mcmcp) standardGeneric("fullGibbs"))

#' MCMC parameter status.
#'
#' This function retrieves the status of MCMC parameters for a model.
#' @param object see showMethods(paramUpdates)
setGeneric("paramUpdates", function(x) standardGeneric("paramUpdates"))

setGeneric("computePrec", function(object) standardGeneric("computePrec"))

#' @export
setGeneric("marginal", function(object, batch, mcmc.params, K=1:4, maxperm=5, ...)
  standardGeneric("marginal"))

#' @export
setGeneric("rowMarginal", function(object, batch, mcmc.params,
                                   model.files, K=1:4, maxperm=5, ...)
  standardGeneric("rowMarginal"))

setGeneric("modelList", function(object) standardGeneric("modelList"))

setGeneric("modelList<-", function(object, value) standardGeneric("modelList<-"))


setGeneric("computeMarginalEachK", function(object, K=1:4, hypp, mcmcp=McmcParams(), maxperm=5)
  standardGeneric("computeMarginalEachK"))

setGeneric("computeMarginalEachK2",
           function(object, batch, maxperm=3,
                    K=1:4, mcmcp=McmcParams(),
                    hypp) standardGeneric("computeMarginalEachK2"))

setGeneric("best", function(object) standardGeneric("best"))

setGeneric("posteriorP", function(object) standardGeneric("posteriorP"))

setGeneric("zChain", function(object) standardGeneric("zChain"))
setGeneric("zChain<-", function(object,value) standardGeneric("zChain<-"))

setGeneric("permuteModes", function(object, ix) standardGeneric("permuteModes"))

setGeneric("pThetaStar", function(kmod, maxperm=5, T2) standardGeneric("pThetaStar"))

setGeneric("pSigma2", function(object) standardGeneric("pSigma2"))

setGeneric("pMixProb", function(object) standardGeneric("pMixProb"))

setGeneric("reducedGibbsThetaFixed", function(object)
  standardGeneric("reducedGibbsThetaFixed"))

setGeneric("reducedGibbsThetaSigmaFixed", function(object)
  standardGeneric("reducedGibbsThetaSigmaFixed"))

setGeneric("reducedGibbsZThetaFixed", function(object) standardGeneric("reducedGibbsZThetaFixed"))

setGeneric("pTheta_Zfixed", function(object) standardGeneric("pTheta_Zfixed"))

setGeneric("pTheta", function(object) standardGeneric("pTheta"))

setGeneric("updateMultinomialProb", function(object) standardGeneric("updateMultinomialProb"))

setGeneric("component", function(object) standardGeneric("component"))
setGeneric("overall", function(object) standardGeneric("overall"))

setGeneric("densities", function(object) standardGeneric("densities"))
setGeneric("densitiesCluster", function(object) standardGeneric("densitiesCluster"))

#' Accessor for extracting the kmeans clusters from a DensityModel
#' instance
#'
#' @seealso \code{\link{DensityModel-class}}
#' @examples
#' truth <- simulateData(N=2500, p=rep(1/3, 3),
#'                       theta=c(-1, 0, 1),
#'                       sds=rep(0.1, 3))
#' dm <- DensityModel(truth)
#' clusters(dm)
#' @export
setGeneric("clusters", function(object) standardGeneric("clusters"))
setGeneric("quantiles", function(object) standardGeneric("quantiles"))
