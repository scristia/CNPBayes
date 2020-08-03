#' @include help.R
NULL

#' Number of observations
setGeneric("numberObs", function(model) standardGeneric("numberObs"))

#' Number of components.
setGeneric("k", function(object) standardGeneric("k"))

#' Change the number of components.
setGeneric("k<-", function(object, value) standardGeneric("k<-"))

setGeneric("z<-", function(object, value) standardGeneric("z<-"))

setGeneric("theta<-", function(object, value) standardGeneric("theta<-"))

setGeneric("sigma2<-", function(object, value) standardGeneric("sigma2<-"))
setGeneric("p<-", function(object, value) standardGeneric("p<-"))
setGeneric("pp<-", function(object, value) standardGeneric("pp<-"))

setGeneric("sigma<-", function(object, value) standardGeneric("sigma<-"))

setGeneric("mu", function(object) standardGeneric("mu"))
setGeneric("mu<-", function(object, value) standardGeneric("mu<-"))

setGeneric("tau2", function(object) standardGeneric("tau2"))

setGeneric("tau2<-", function(object, value) standardGeneric("tau2<-"))
setGeneric("nu.0<-", function(object, value) standardGeneric("nu.0<-"))
setGeneric("sigma2.0<-", function(object, value) standardGeneric("sigma2.0<-"))

setGeneric("dataMean", function(object) standardGeneric("dataMean"))
setGeneric("dataPrec", function(object) standardGeneric("dataPrec"))
setGeneric("dataMean<-", function(object, value) standardGeneric("dataMean<-"))
setGeneric("dataPrec<-", function(object, value) standardGeneric("dataPrec<-"))

setGeneric("chains<-", function(object, value) standardGeneric("chains<-"))

#' Retrieve simulated chains from model object.
#'
#' The method \code{chains} applied to a \code{MixtureModel}-derived
#' class will return an object of class \code{McmcChains} that
#' contains the chains for all simulated parameters. Typically, chains
#' is called in conjunction with an accessor for one of these
#' parameters.
#' @examples
#'     theta.chain <- theta(chains(MultiBatchModelExample))
#'     dim(theta.chain)
#'     plot.ts(theta.chain, plot.type="single",
#'             col=seq_len(k(SingleBatchModelExample)))
#' @param object \code{showMethods(chains)}
#' @return The simulated chains.
#' @export
#' @docType methods
#' @rdname chains-method
setGeneric("chains", function(object) standardGeneric("chains"))

#' Accessor for Hyperparameters object for a MixtureModel-derived object
#'
#' @examples
#'     hyperParams(SingleBatchModelExample)
#' @param object see \code{showMethods(hyperParams)}
#' @return The Hyperparameters of a MixtureModel
#' @export
#' @docType methods
#' @rdname hyperParams-method
setGeneric("hyperParams", function(object) standardGeneric("hyperParams"))

#' Replace the hyperparameters for a \code{MixtureModel}-derived object
#'
#' @examples
#' hypp <- Hyperparameters(type="marginal",
#'                         k=k(SingleBatchModelExample),
#'                         alpha=c(9, 9, 10))
#' hyperParams(SingleBatchModelExample) <- hypp
#'
#' @param value an object of class 'Hyperparameters'
#' @export
#' @docType methods
#' @rdname hyperParams-method
setGeneric("hyperParams<-", function(object,value) standardGeneric("hyperParams<-"))

setGeneric("McmcChains", function(object) standardGeneric("McmcChains"))

setGeneric("McmcChainsTrios", function(object) standardGeneric("McmcChainsTrios"))


setGeneric("hist")

#' Retrieve batches from object.
#'
#' The batches are represented as a vector of integers.
#' @examples
#'      batch(MultiBatchModelExample)
#' @param object see \code{showMethods(batch)}
#' @return The batch of each data element.
#' @export
#' @docType methods
#' @rdname batch-method
setGeneric("batch", function(object) standardGeneric("batch"))
setGeneric("batch<-", function(object,value) standardGeneric("batch<-"))

setGeneric("startingValues", function(object, params, zz) standardGeneric("startingValues"))

setGeneric("computeMeans", function(object) standardGeneric("computeMeans"))
setGeneric("computeVars", function(object) standardGeneric("computeVars"))

setGeneric("computePotential", function(object) standardGeneric("computePotential"))

setGeneric("initializeSigma2", function(object) standardGeneric("initializeSigma2"))

setGeneric("initializeTheta", function(object) standardGeneric("initializeTheta"))

setGeneric("alpha<-", function(object, value) standardGeneric("alpha<-"))

#' Calculate BIC of a model
#'
#' @examples
#' mb <- MultiBatchModelExample
#' mcmcParams(mb) <- McmcParams(iter=100, burnin=50)
#' mb <- posteriorSimulation(mb)
#' bic(mb)
#' @param object see \code{showMethods(bic)}
#' @return The BIC of the model.
#' @docType methods
#' @rdname bic-method
#' @export
setGeneric("bic", function(object) standardGeneric("bic"))

#' Accessor for the theta parameter in the hierarchical mixture model
#'
#' The interpretation of \code{theta} depends on whether \code{object}
#' is a \code{MarginalModel} or a \code{BatchModel}.  For
#' \code{BatchModel}, \code{theta} is a matrix of size B x K, where B is
#' the number of batches and K is the number of components.
#' Each column of the \code{theta} matrix can be interpreted as the
#' batch means for a particular component. For objects of class
#' \code{MarginalModel} (assumes no batch effect), \code{theta} is a
#' vector of length K. Each element of \code{theta} can be interpreted
#' as the mean for a component. See the following examples for accessing
#' the current value of \code{theta} from a \code{MixtureModel}-derived
#' object, and for plotting the chain of \code{theta} values.
#' @examples
#' ## MarginalModel
#' \dontrun{
#' k(SingleBatchModelExample)
#' theta(SingleBatchModelExample)
#' plot.ts(theta(chains(SingleBatchModelExample)))
#' ## BatchModel
#' k(MultiBatchModelExample)
#' length(unique(batch(MultiBatchModelExample)))
#' theta(MultiBatchModelExample)
#' ## Plot means for batches in one component
#' plot.ts(theta(chains(MultiBatchModelExample))[, 1:3])
#' }
#' @param object see \code{showMethods(theta)}
#' @return A vector of length number of components or a matrix of size
#' number of batches x number of components
#' @export
#' @docType methods
#' @rdname theta-method
setGeneric("theta", function(object) standardGeneric("theta"))

#' @rdname theta-method
#' @export
setGeneric("theta<-", function(object, value) standardGeneric("theta<-"))

#' Retrieve the variances of each component and batch distribution
#'
#' For a MarginalModel, this function returns a vector of variances. For a BatchModel, returns a matrix of size number of batches by number of components.
#' @examples
#'      sigma2(SingleBatchModelExample)
#' @param object see \code{showMethods(sigma2)}
#' @return A vector of length number of components or a matrix of size 
#' number of batches x number of components
#' @export
#' @docType methods
#' @rdname sigma2-method
setGeneric("sigma2", function(object) standardGeneric("sigma2"))

setGeneric("sigma_", function(object) standardGeneric("sigma_"))

#' Retrieve the probability of latent variable membership by observation.
#'
#' @examples
#'      head(probz(SingleBatchModelExample))
#' @param object see \code{showMethods(probz)}
#' @return A matrix of size number of observations x number of components
#' @export
#' @docType methods
#' @rdname probz-method
setGeneric("probz", function(object) standardGeneric("probz"))

setGeneric("probz<-", function(object, value) standardGeneric("probz<-"))

setGeneric("nu.0", function(object) standardGeneric("nu.0"))

setGeneric("sigma2.0", function(object) standardGeneric("sigma2.0"))

setGeneric("y", function(object) standardGeneric("y"))

setGeneric("y<-", function(object, value) standardGeneric("y<-"))

#' Retrieve data.
#'
#' @param object see \code{showMethods(oned)}
#' @return A vector the length of the data
#' @export
#' @docType methods
#' @rdname oned-method
#' @export
setGeneric("oned", function(object) standardGeneric("oned"))

setGeneric("oned<-", function(object, value) standardGeneric("oned<-"))

#' Retrieve latent variable assignments.
#'
#' Retrieves the simulated latent variable assignments of each observation at each MCMC simulation.
#' @examples
#'      head(z(SingleBatchModelExample))
#' @param object see \code{showMethods(z)}
#' @return A vector the length of the data
#' @export
#' @docType methods
#' @rdname z-method
setGeneric("z", function(object) standardGeneric("z"))

setGeneric("modalParameters", function(object) standardGeneric("modalParameters"))

setGeneric("computeModes", function(object) standardGeneric("computeModes"))

setGeneric("modes", function(object) standardGeneric("modes"))

setGeneric("modes<-", function(object,value) standardGeneric("modes<-"))


setGeneric("mu.0", function(object) standardGeneric("mu.0"))
setGeneric("tau2.0", function(object) standardGeneric("tau2.0"))

setGeneric("eta.0", function(object) standardGeneric("eta.0"))
setGeneric("eta.0<-", function(object,value) standardGeneric("eta.0<-"))

setGeneric("m2.0", function(object) standardGeneric("m2.0"))

setGeneric("m2.0<-", function(object,value) standardGeneric("m2.0<-"))

setGeneric("showMeans", function(object) standardGeneric("showMeans"))
setGeneric("showSigmas", function(object) standardGeneric("showSigmas"))

setGeneric("collapseBatch", function(object, provisional_batch, THR=0.1, nchar=8) standardGeneric("collapseBatch"))



setGeneric("thetac", function(object) standardGeneric("thetac"))

setGeneric("thetaMean", function(object) standardGeneric("thetaMean"))

setGeneric("sigmaMean", function(object) standardGeneric("sigmaMean"))

setGeneric("pMean", function(object) standardGeneric("pMean"))

setGeneric("tablez", function(object) standardGeneric("tablez"))

setGeneric("nStarts", function(object) standardGeneric("nStarts"))

setGeneric("nStarts<-", function(object, value) standardGeneric("nStarts<-"))

setGeneric("alpha", function(object) standardGeneric("alpha"))


#' Retrieve log likelihood.
#'
#' @examples
#' log_lik(SingleBatchModelExample)
#' @param object see showMethods(log_lik)
#' @return The log likelihood
#' @export
#' @docType methods
#' @rdname log_lik-method
setGeneric("log_lik", function(object) standardGeneric("log_lik"))

setGeneric("log_lik<-", function(object,value) standardGeneric("log_lik<-"))

setGeneric("computeLoglik", function(object) standardGeneric("computeLoglik"))


#' Number of burnin iterations.
#'
#' This function retrieves the number of burnin simulations to be discarded.
#' @examples
#' burnin(SingleBatchModelExample)
#' mp <- mcmcParams(SingleBatchModelExample)
#' burnin(mp)
#' @param object see \code{showMethods(burnin)}
#' @return The number of burnin simulations.
#' @export
#' @docType methods
setGeneric("burnin", function(object) standardGeneric("burnin"))


setGeneric("min_GR", function(object) standardGeneric("min_GR"))

setGeneric("min_effsize", function(object)
  standardGeneric("min_effsize"))

setGeneric("max_burnin", function(object)
  standardGeneric("max_burnin"))

setGeneric("min_chains", function(object)
  standardGeneric("min_chains"))

#' Reset number of burnin iterations.
#'
#' This function changes the number of burnin simulations to be discarded.
#' @param value new number of burnin iterations
#' @export
#' @docType methods
#' @rdname burnin-method
setGeneric("burnin<-", function(object, value) standardGeneric("burnin<-"))

#' Reset number of iterations.
#'
#' This function changes the number of simulations.
#' @param force Allow changing of the size of the elements?
#' @param value new number of iterations
#' @export
#' @docType methods
#' @rdname iter-method
setGeneric("iter<-", function(object, value) standardGeneric("iter<-"))

#' Number of MCMC iterations.
#'
#' This function retrieves the number of iterations of an MCMC simulation.
#' @examples
#'      iter(SingleBatchModelExample)
#' @param object see \code{showMethods(iter)}
#' @return The number of MCMC iterations
#' @export
#' @docType methods
#' @rdname iter-method
setGeneric("iter", function(object) standardGeneric("iter"))

#' Get or set the number of thinning intervals.
#'
#' This function gets or sets the number of thinning intervals used for an MCMC
#' simulation.
#' @examples
#'      thin(SingleBatchModelExample)
#' @param object see showMethods(thin)
#' @return An integer of the number of thinning intervals
#' @export
#' @docType methods
#' @rdname thin-method
setGeneric("thin", function(object) standardGeneric("thin"))

#' @examples
#'      thin(SingleBatchModelExample) <- 10L
#' @export
#' @docType methods
#' @rdname thin-method
setGeneric("thin<-", function(object, value) standardGeneric("thin<-"))

#' Run MCMC simulation.
#'
#' nStarts chains are run. b burnin iterations are run and then discarded.
#' Next, s iterations are run in each crain. The user can also specify
#' an alternative number of components.
#' The mode of the MCMC simulation is also calculated.
#' @examples
#' # Fit model with pre-specified number of components (k=3)
#' set.seed(123)
#' ## specify small number of iterations so that the example runs quickly
#' mp <- McmcParams(iter=2, burnin=0, nStarts=3)
#' sb <- SingleBatchModelExample
#' mcmcParams(sb) <- mp
#' posteriorSimulation(sb)
#'
#' # Run additional iterations, but set nStart = 0 so that the last value of the chain is the first value of the next chain
#' mcmcParams(sb) <- McmcParams(iter=5, nStarts=0, burnin=0)
#' posteriorSimulation(sb)
#' @param object see showMethods(posteriorSimulation)
#' @param k The number of a priori components. This is optional and if not
#' specified, the stored k model components are used. This parameters is
#' useful for running multiple models of varying components.
#' @return An object of class 'MarginalModel' or 'BatchModel'
#' @seealso  See \code{\link{ggMixture}} for plotting the model-based densities.
#' @export
#' @docType methods
#' @rdname posteriorSimulation-method
setGeneric("posteriorSimulation", function(object, k) standardGeneric("posteriorSimulation"))

setGeneric("isSB", function(object) standardGeneric("isSB"))

setGeneric("computePrior", function(object) standardGeneric("computePrior"))

setGeneric("relabel", function(object, zindex) standardGeneric("relabel"))

setGeneric("runBurnin", function(object, mcmcp) standardGeneric("runBurnin"))
setGeneric("runMcmc", function(object, mcmcp) standardGeneric("runMcmc"))

setGeneric("paramUpdates<-", function(x, value) standardGeneric("paramUpdates<-"))

setGeneric("zFreq", function(object) standardGeneric("zFreq"))
setGeneric("zFreq<-", function(object,value) standardGeneric("zFreq<-"))

setGeneric("zFreqPar", function(object) standardGeneric("zFreqPar"))
setGeneric("zFreqPar<-", function(object,value) standardGeneric("zFreqPar<-"))

setGeneric("triodata_lrr", function(object) standardGeneric("triodata_lrr"))

#' Retrieve MCMC parameters from model.
#'
#' View number of iterations, burnin, etc.
#' @examples
#'      mcmcParams(SingleBatchModelExample)
#' @param object see \code{showMethods(mcmcParams)}
#' @return An object of class 'McmcParams'
#' @export
#' @docType methods
#' @rdname mcmcParams-method
setGeneric("mcmcParams", function(object) standardGeneric("mcmcParams"))

#' Replace MCMC parameters of model.
#'
#' Replace number of iterations, burnin, etc. Any update of the MCMC parameters will trigger an update of the chains. However, if iter (the number of MCMC iterations) is set to a nonpositive value, the chains will not be updated and kept as is.
#' @param force logical value. If false (default) the update will not proceed.
#' @param value an object of class 'McmcParams' containing the new number of iterations, etc.
#' @export
#' @docType methods
#' @rdname mcmcParams-method
setGeneric("mcmcParams<-", function(object, value) standardGeneric("mcmcParams<-"))

setGeneric("logPrior", function(object) standardGeneric("logPrior"))

setGeneric("logPrior<-", function(object,value) standardGeneric("logPrior<-"))

setGeneric("paramUpdates", function(x) standardGeneric("paramUpdates"))

setGeneric("computePrec", function(object) standardGeneric("computePrec"))

setGeneric("marginal", function(object, batch, mcmc.params, K=1:4, maxperm=5, ...)
  standardGeneric("marginal"))

setGeneric("zChain", function(object) standardGeneric("zChain"))

setGeneric("updateMultinomialProb", function(object) standardGeneric("updateMultinomialProb"))

setGeneric("component", function(object) standardGeneric("component"))
setGeneric("overall", function(object) standardGeneric("overall"))

setGeneric("densities", function(object) standardGeneric("densities"))
setGeneric("densitiesCluster", function(object) standardGeneric("densitiesCluster"))

setGeneric("quantiles", function(object) standardGeneric("quantiles"))

#' Compute the marginal likelihood of a converged model.
#'
#' The recommended function for fitting mixture models and evaluating convergence is through the `gibbs` function. This function will return a list of models ordered by the marginal likelihood. The marginal likelihood is computed using the Chib's estimator (JASA, Volume 90 (435), 1995).
#'
#'
#' @examples
#' ## In practice, run a much longer burnin and increase the number of
#' ## iterations to save after burnin
#'    mm <- SingleBatchModelExample
#'    mcmcParams(mm) <- McmcParams(iter=50, burnin=0, nStarts=0)
#'    mm <- posteriorSimulation(mm)
#'    marginalLikelihood(mm)
#' @param model An object of class \code{MarginalModel}, or a list of
#'        \code{MarginalModel}'s. Can also be an object of \code{BatchModel} or
#'        a list of such models.
#' @param params A list containing parameters for marginalLikelihood computation. See \code{mlParams} for details.
#'
#' @seealso See \code{\link{mlParams}} for parameters related to computing the log marginal likelihood via Chib's estimator. See \code{\link{gibbs}} for fitting multiple mixture models and returning a list sorted by the marginal likelihood.  See \code{\link{marginal_lik}} for the accessor.
#'
#' Note: currently thinning of the reduced MCMC chains is not allowed.
#'
#' @return A vector of the marginal likelihood of the model(s)
#' @export
#' @docType methods
#' @rdname marginalLikelihood-method
setGeneric("marginalLikelihood", function(model, params=mlParams()) standardGeneric("marginalLikelihood"))

setGeneric("chromosome", function(object, ...) standardGeneric("chromosome"))

setMethod("chromosome", "GenomicRanges", function(object, ...){
  as.character(seqnames(object))
})

setGeneric("label_switch", function(object) standardGeneric("label_switch"))


setGeneric("label_switch<-", function(object, value) standardGeneric("label_switch<-"))

#' Accessor for the log marginal likelihood of a SB, SBP, MB, or MBP model
#'
#' The marginal likelihood is computed by Chib's estimator (JASA, Volume 90 (435), 1995).
#'
#' @seealso See \code{\link{marginalLikelihood}} for computing the marginal likelihood of a mixture model.
#' @param object a SB, SBP, MB, or MBP model
#' @export
#' @examples
#' sb <- SingleBatchModelExample
#' marginal_lik(sb)
setGeneric("marginal_lik", function(object) standardGeneric("marginal_lik"))

#' @export
#' @rdname marginal_lik
setGeneric("marginal_lik<-", function(object, value) standardGeneric("marginal_lik<-"))

setGeneric("updateZ", function(object) standardGeneric("updateZ"))

setGeneric("gatherChains", function(object) standardGeneric("gatherChains"))

setGeneric("sortComponentLabels", function(model) standardGeneric("sortComponentLabels"))

setGeneric("isOrdered", function(object) standardGeneric("isOrdered"))

setGeneric("ggChains", function(model) standardGeneric("ggChains"))

#' Plotting posterior predictive densities from mixture model
#' 
#' @param bins a length-one numeric vector indicating the number of bins -- passed to \code{geom_hist}
#' @param mixtheme a ggplot theme
#' @export
#' @rdname ggplot-functions
setGeneric("ggMixture", function(model, bins=100, mixtheme, shift_homozygous) standardGeneric("ggMixture"))


setGeneric("CopyNumberModel", function(model, params=mapParams()) standardGeneric("CopyNumberModel"))

#' Map mixture components to copy number states
#'
#' @return a character vector of length k (number of components) indicating the inferred copy number of each component.  The mapping is not neccessarily one-to-one as multiple mixture components can be mapped to a single copy number state.
#' @export
#' @param object a `MultiBatch` model
#' @rdname mapping
setGeneric("mapping", function(object) standardGeneric("mapping"))

#' @param value a k-length numeric vector with values in {1, 2, ..., k}, where k is the number of mixture components
#' @export
#' @rdname mapping
setGeneric("mapping<-", function(object, value) standardGeneric("mapping<-"))

setGeneric("numberStates", function(model) standardGeneric("numberStates"))

setGeneric("probCopyNumber", function(model) standardGeneric("probCopyNumber"))

setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))

setGeneric("mergeComponents", function(model, j) standardGeneric("mergeComponents"))

### t-distribution stuff

setGeneric("dfr", function(object) standardGeneric("dfr"))

setGeneric("dfr<-", function(object, value) standardGeneric("dfr<-"))

setGeneric("u", function(object) standardGeneric("u"))


setGeneric("probzpar", function(object) standardGeneric("probzpar"))

setGeneric("probzpar<-", function(object, value) standardGeneric("probzpar<-"))

setGeneric("p", function(object) standardGeneric("p"))

setGeneric("sigma_<-", function(object, value) standardGeneric("sigma_<-"))

setGeneric("flags", function(object) standardGeneric("flags"))
setGeneric("flags<-", function(object, value) standardGeneric("flags<-"))
setGeneric("current_values<-", function(object, value) standardGeneric("current_values<-"))

setGeneric("current_values", function(object, value) standardGeneric("current_values"))

setGeneric("current_values2<-", function(object, value) standardGeneric("current_values2<-"))

setGeneric("current_values2", function(object, value) standardGeneric("current_values2"))

setGeneric("u<-", function(object, value) standardGeneric("u<-"))

setGeneric("parameters", function(object) standardGeneric("parameters"))
setGeneric("parameters<-", function(object, value) standardGeneric("parameters<-"))

#' @export
setGeneric("summaries", function(object) standardGeneric("summaries"))

setGeneric("summaries<-", function(object, value) standardGeneric("summaries<-"))

setGeneric("summaries2", function(object) standardGeneric("summaries2"))
setGeneric("summaries2<-", function(object, value) standardGeneric("summaries2<-"))

setGeneric("dataSd", function(object) standardGeneric("dataSd"))

setGeneric("findSurrogates", function(object, THR=0.1, min_oned=-1) standardGeneric("findSurrogates"))

setGeneric("isSampled", function(object, THR=0.1) standardGeneric("isSampled"))

setGeneric("down_sample", function(object) standardGeneric("down_sample"))

setGeneric("down_sample<-", function(object, value) standardGeneric("down_sample<-"))

setGeneric("downSampledData", function(object) standardGeneric("downSampledData"))

setGeneric("downSampledData<-", function(x, value) standardGeneric("downSampledData<-"))

setGeneric("probability_z", function(object) standardGeneric("probability_z"))

setGeneric("mcmc2", function(object, guide) standardGeneric("mcmc2"))

setGeneric("setModes", function(object) standardGeneric("setModes"))

setGeneric("upsample_z", function(object) standardGeneric("upsample_z"))

setGeneric("specs", function(object) standardGeneric("specs"))

setGeneric("specs<-", function(object, value) standardGeneric("specs<-"))

setGeneric("downSampleModel", function(object, N, i) standardGeneric("downSampleModel") )

setGeneric("upSampleModel", function(downsampled.model, data.full) standardGeneric("upSampleModel"))

setGeneric("listChains", function(object) standardGeneric("listChains"))

setGeneric("predictive", function(object) standardGeneric("predictive"))

setGeneric("zstar", function(object) standardGeneric("zstar"))

setGeneric("predictive<-", function(object, value) standardGeneric("predictive<-"))

setGeneric("revertBack", function(object, mbm) standardGeneric("revertBack"))

setGeneric("numBatch", function(object) standardGeneric("numBatch"))

setGeneric("numBatch<-", function(object, value) standardGeneric("numBatch<-"))

setGeneric("max_burnin<-", function(object, value) standardGeneric("max_burnin<-"))

setGeneric("useModes", function(object) standardGeneric("useModes"))

setGeneric("compute_marginal_lik", function(object, params) standardGeneric("compute_marginal_lik"))

#' Accessor for name of MultiBatch model
#'
#' MultiBatch objects are named according to whether the mixture components are stratified by batch (MB) or single batch (SB), the variance across all mixture components within a batch is pooled (P), and the number of mixture components.  For example, MB3 is a multi-batch 3-component mixture model without pooled variance, while SBP4 is a 4-component single-batch mixture model using a pooled estimate of the variance.  
#' @param object a `MultiBatch` object
#' @export
setGeneric("modelName", function(object) standardGeneric("modelName"))

setGeneric("singleBatchGuided", function(x, guide) standardGeneric("singleBatchGuided"))

setGeneric("convergence", function(object) standardGeneric("convergence"))

setGeneric("isMendelian", function(object) standardGeneric("isMendelian"))

setGeneric("augmentData2", function(object) standardGeneric("augmentData2"))

setGeneric("baf_loglik", function(object, snpdat) standardGeneric("baf_loglik"))

setGeneric("toSingleBatch", function(object) standardGeneric("toSingleBatch"))

#' @export
setGeneric("isSimulated", function(object)
  standardGeneric("isSimulated"))

setGeneric(".compute_loglik", function(object) standardGeneric(".compute_loglik"))
